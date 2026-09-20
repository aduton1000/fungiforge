// FungiForge Layer-1 subworkflow — route each isolate through the 15 stages, then
// aggregate every per-stage result.json into one per-isolate report + master row.
include { DB_CHECK }    from '../modules/stage00_dbcheck.nf'
include { BASECALL }    from '../modules/stage00_basecall.nf'
include { READ_QC }     from '../modules/stage01_readqc.nf'
include { READ_TRIAGE } from '../modules/stage00b_triage.nf'
include { KMER_PROFILE } from '../modules/stage01b_kmer.nf'
include { ASSEMBLE }    from '../modules/stage02_assemble.nf'
include { SR_ASSEMBLE } from '../modules/stage02b_srassemble.nf'
include { PURGE_DUPS } from '../modules/stage02c_purge.nf'
include { MEDAKA }      from '../modules/stage03_polish.nf'
include { SRPOLISH }    from '../modules/stage03b_srpolish.nf'
include { DECONTAM }    from '../modules/stage04_decontam.nf'
include { ORGANELLE }   from '../modules/stage04b_organelle.nf'
include { ASSEMBLY_QC } from '../modules/stage05_assembly_qc.nf'
include { GATE as GATE_READS; GATE as GATE_ASSEMBLY } from '../modules/stage05b_gate.nf'
include { REPEATMASK }  from '../modules/stage06_repeatmask.nf'
include { PREDICT }     from '../modules/stage07_predict.nf'
include { EGGNOG }      from '../modules/stage07b_eggnog.nf'
include { INTERPROSCAN } from '../modules/stage07c_interproscan.nf'
include { ANNOTATE }    from '../modules/stage07_annotate.nf'
include { IDENTIFY }    from '../modules/stage08_identify.nf'
include { BUSCO_LINEAGE } from '../modules/stage08b_busco_lineage.nf'
include { RESISTANCE }  from '../modules/stage09_resistance.nf'
include { MOBILE }      from '../modules/stage10_mobile.nf'
include { BGC }         from '../modules/stage11_bgc.nf'
include { NOVELTY }     from '../modules/stage12_novelty.nf'
include { EXTRAS }      from '../modules/stage13_extras.nf'
include { REPORT }      from '../modules/stage14_report.nf'
include { PROVENANCE }  from '../modules/stage15_provenance.nf'
include { COHORT }      from '../modules/stage16_cohort.nf'
include { SUMMARY }     from '../modules/stage17_summary.nf'

workflow FUNGIFORGE {
  take:
    samples          // tuple(meta, filesMap[ont,r1,r2])
    run_info         // path: run_info.json (pipeline/nextflow/params/container map, from main.nf)

  main:
    // 00. database check — one task per run, BEFORE any compute: every isolate's first stage
    //     waits for it, so a missing database fails the run in seconds, not hours.
    DB_CHECK()
    reads0 = samples.combine(DB_CHECK.out.manifest).map { meta, f, _manifest -> tuple(meta, f.ont, f.r1, f.r2) }
    // W6.2: the RNA-seq pair travels separately to stage 07a. It is not read before prediction,
    // so it stays out of the read channel that the QC, assembly and polishing stages consume.
    rna_ch = samples.map { meta, f -> tuple(meta, f.rna1, f.rna2) }

    // 0. optional Dorado basecalling (ont column = pod5 dir when --basecall).
    //    Only long-read isolates carry ONT/pod5; short-read-only rows pass through.
    if (params.basecall) {
      lr0 = reads0.filter { row -> row[0].assembly_mode == 'longread' }
      sr0 = reads0.filter { row -> row[0].assembly_mode == 'shortread' }
      reads_bc = BASECALL(lr0).reads.mix(sr0)
      basecall_json = BASECALL.out.json
    } else {
      reads_bc = reads0
      basecall_json = channel.empty()
    }

    // 1. read QC + filtering (ONT chopper/NanoPlot; Illumina fastp — mode-aware). Emits the
    //    filtered ONT reads and the TRIMMED Illumina reads: everything downstream uses these.
    READ_QC(reads_bc)

    // 00b. read-level triage (W2.2): Kraken2 on a subsample -> early verdict; an isolate that
    //      is non_fungal/human at read level is stopped here (reason read_triage:<verdict>)
    //      and never assembled. --force_all disables it; --read_triage false skips the stage.
    if (params.read_triage) {
      READ_TRIAGE(READ_QC.out.reads)
      triage_json = READ_TRIAGE.out.json
      read_gate = READ_TRIAGE.out.json.map { meta, tj ->
        def v = new groovy.json.JsonSlurper().parseText(tj.text)?.verdict ?: 'not_run'
        tuple(meta, (!params.force_all && v in ['non_fungal', 'human']) ? "read_triage:${v}" : null)
      }
      GATE_READS(read_gate.filter { _meta, reason -> reason != null })
      gate_reads_json = GATE_READS.out.json
      reads_ok = READ_QC.out.reads.join(read_gate.filter { _meta, reason -> reason == null }.map { meta, _reason -> tuple(meta) })
    } else {
      triage_json = channel.empty(); gate_reads_json = channel.empty()
      reads_ok = READ_QC.out.reads
    }

    // 01b. k-mer profile (W2.2): genome size, heterozygosity, ploidy hint, coverage
    if (params.kmer_profile) { KMER_PROFILE(reads_ok); kmer_json = KMER_PROFILE.out.json }
    else                     { kmer_json = channel.empty() }

    // Split isolates by assembly mode: long-read (ONT ± hybrid) vs short-read (Illumina only).
    reads_lr = reads_ok.filter { row -> row[0].assembly_mode == 'longread' }
    reads_sr = reads_ok.filter { row -> row[0].assembly_mode == 'shortread' }

    // 2. assembly — long-read (Flye --nano-hq + purge_dups) OR short-read (SPAdes)
    ASSEMBLE(reads_lr)
    SR_ASSEMBLE(reads_sr)

    // 3. polishing: Medaka (ONT long-read only) then Illumina Polypolish hybrid when
    //    short reads are present. Short-read-only assemblies skip Medaka (no ONT) and
    //    flow straight into Stage 03b, which records mode=illumina_only (high conf).
    ont_ch  = reads_lr.map { meta, ont, _r1, _r2 -> tuple(meta, ont) }
    ilmn_ch = reads_ok.map { meta, _ont, r1, r2 -> tuple(meta, r1, r2) }
    // 2c. haplotig purging in the base image (L30): long-read assemblies only
    PURGE_DUPS(ASSEMBLE.out.assembly.join(ont_ch))
    MEDAKA(PURGE_DUPS.out.assembly.join(ont_ch))
    pre_srpolish = MEDAKA.out.assembly.mix(SR_ASSEMBLE.out.assembly)
    SRPOLISH(pre_srpolish.join(ilmn_ch))

    // 4. decontamination + organelle split -> nuclear / mito
    DECONTAM(SRPOLISH.out.assembly)

    // 4b. mitochondrial genome: core genes, introns, circularity, copy number, heteroplasmy (W2.7)
    ORGANELLE(DECONTAM.out.mito.join(DECONTAM.out.nuclear).join(reads_ok))

    // 5. assembly QC + completeness (compleasm/BUSCO + contiguity)
    ASSEMBLY_QC(DECONTAM.out.nuclear)

    // 5b. the gate (W2.1): an isolate that is not a fungus (Kraken2 verdict non_fungal or human)
    //     or whose assembly failed QC does not go through the fungal stages 06-13. It gets a
    //     `gate` JSON (status skipped + reason) and still reaches REPORT, so its master row
    //     exists with the reason. --force_all sends every isolate through regardless.
    gate = DECONTAM.out.json.join(ASSEMBLY_QC.out.json).map { meta, dj, qj ->
      def d = new groovy.json.JsonSlurper().parseText(dj.text)
      def q = new groovy.json.JsonSlurper().parseText(qj.text)
      def verdict = d?.verdict ?: 'not_run'
      def reason = null
      if (!params.force_all) {
        if (verdict in ['non_fungal', 'human']) reason = "verdict:${verdict}"
        else if (q?.qc_pass == false)          reason = "qc_pass:false"
      }
      tuple(meta, reason)
    }
    pass_meta = gate.filter { _meta, reason -> reason == null }.map { meta, _reason -> tuple(meta) }
    GATE_ASSEMBLY(gate.filter { _meta, reason -> reason != null })
    nuclear_ok = DECONTAM.out.nuclear.join(pass_meta)     // (meta, nuclear) for isolates that passed the gate

    // 6. repeat modeling + soft-masking (RepeatModeler2 / RepeatMasker)
    REPEATMASK(nuclear_ok)

    // 8. identification (ITS + CaM/BenA/TEF1/RPB2/LSU concordance, MLST, optional genome ANI).
    //    Runs before prediction so the species call can seed the annotation training (W2.5).
    IDENTIFY(nuclear_ok)

    // 7a. gene prediction (Funannotate predict; GeneMark-ES when licensed; species-aware training)
    PREDICT(REPEATMASK.out.masked.join(IDENTIFY.out.species).join(rna_ch))

    // 7b/7c. eggNOG-mapper and InterProScan on the predicted proteins (separate images);
    //        off switches hand ANNOTATE an empty placeholder so funannotate runs without them.
    if (!params.skip_eggnog) { EGGNOG(PREDICT.out.proteins); eggnog_ch = EGGNOG.out.annotations; eggnog_json = EGGNOG.out.json }
    else { eggnog_ch = PREDICT.out.proteins.map { m, _p -> tuple(m, file("${projectDir}/assets/NO_FILE")) }; eggnog_json = channel.empty() }
    if (params.run_interproscan) { INTERPROSCAN(PREDICT.out.proteins); ips_ch = INTERPROSCAN.out.xml; ips_json = INTERPROSCAN.out.json }
    else { ips_ch = PREDICT.out.proteins.map { m, _p -> tuple(m, file("${projectDir}/assets/NO_FILE")) }; ips_json = channel.empty() }

    // 7. functional annotation (Funannotate annotate + eggNOG/InterProScan results)
    ANNOTATE(PREDICT.out.results.join(eggnog_ch).join(ips_ch))

    // 8b. species-aware BUSCO with the lineage chosen from the species call (W2.3)
    BUSCO_LINEAGE(nuclear_ok.join(IDENTIFY.out.lineage))

    // 9. antifungal resistance (bespoke panel + cyp51A TR34/TR46 module).
    //    needs proteins (substitutions) + species + nuclear & GBK (promoter TR locus) + reads (W2.4)
    RESISTANCE(ANNOTATE.out.proteins
                 .join(IDENTIFY.out.species)
                 .join(DECONTAM.out.nuclear)
                 .join(ANNOTATE.out.gbk)
                 .join(SRPOLISH.out.json)
                 .join(reads_ok))                      // W2.4: read-level genotyping of the hotspots

    // 10-13: skippable stages — a skip flag truly skips the process (empty channel),
    // rather than running it emptily (which would also drag in its container).
    mobile_ch = channel.empty()
    if (!params.skip_mge)     { MOBILE(REPEATMASK.out.telib.join(REPEATMASK.out.tbl).join(DECONTAM.out.nuclear).join(DECONTAM.out.mito)); mobile_ch = MOBILE.out.json }   // W2.8

    bgc_ch = channel.empty()
    if (!params.skip_bgc)     { BGC(ANNOTATE.out.gbk); bgc_ch = BGC.out.json }

    novelty_ch = channel.empty()
    if (!params.skip_novelty) { NOVELTY(DECONTAM.out.nuclear.join(IDENTIFY.out.markers).join(IDENTIFY.out.json)); novelty_ch = NOVELTY.out.json }

    extras_ch = channel.empty()
    if (!params.skip_extras)  { EXTRAS(ANNOTATE.out.proteins.join(ANNOTATE.out.gbk).join(RESISTANCE.out.bam)); extras_ch = EXTRAS.out.json }   // W2.6: GenBank for MAT synteny, stage-09 BAM for ploidy

    // 14. aggregate every per-stage result.json per isolate -> report + master row
    all_json = READ_QC.out.json
      .mix(basecall_json, triage_json, gate_reads_json, kmer_json, ASSEMBLE.out.json, SR_ASSEMBLE.out.json, MEDAKA.out.json,
           PURGE_DUPS.out.json, SRPOLISH.out.json, DECONTAM.out.json, ORGANELLE.out.json, ASSEMBLY_QC.out.json, GATE_ASSEMBLY.out.json, REPEATMASK.out.json,
           PREDICT.out.json, eggnog_json, ips_json, ANNOTATE.out.json, IDENTIFY.out.json, BUSCO_LINEAGE.out.json, RESISTANCE.out.json,
           mobile_ch, bgc_ch, novelty_ch, extras_ch)
      .map { meta, j -> tuple(meta.id, meta, j) }
      .groupTuple(by: 0)
      .map { _id, metas, jsons -> tuple(metas[0], jsons) }
    REPORT(all_json)

    // 16. run-level cohort phylogenomics + clonality over the isolates that passed the gates (W3.1)
    cohort_json = channel.empty()
    if (!params.skip_cohort) {
      regions_ch = params.skip_bgc ? channel.empty() : BGC.out.regions.map { _m, f -> f }
      COHORT(nuclear_ok.map { _m, f -> f }.collect(),
             ASSEMBLY_QC.out.busco_sc.join(pass_meta).map { _m, f -> f }.collect(),
             IDENTIFY.out.species.map { _m, f -> f }.collect(),
             ASSEMBLY_QC.out.json.join(pass_meta).map { _m, f -> f }.collect(),
             regions_ch.collect().ifEmpty([]),
             bgc_ch.map { _m, f -> f }.collect().ifEmpty([]))
      cohort_json = COHORT.out.json
    }

    // 15. run-level provenance: every stage JSON of every isolate + DB manifest + run metadata.
    //     Waits for all reports (masters.collect()) so it is the last task of the run.
    PROVENANCE(run_info,
               DB_CHECK.out.manifest,
               all_json.map { _meta, jsons -> jsons }.flatten().collect(),
               REPORT.out.master.collect())

    // 17. run-level cohort summary: merged master table (schema-validated), cohort report,
    //     MultiQC custom content (W3.3). Waits for every isolate's master row.
    SUMMARY(REPORT.out.master.collect(),
            cohort_json.ifEmpty { file("${projectDir}/assets/NO_COHORT") },
            PROVENANCE.out.provenance.ifEmpty { file("${projectDir}/assets/NO_PROVENANCE") })

  emit:
    cohort     = cohort_json
    summary    = SUMMARY.out.summary
    master     = REPORT.out.master
    reports    = REPORT.out.report
    provenance = PROVENANCE.out.provenance
}
