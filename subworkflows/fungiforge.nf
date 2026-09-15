// FungiForge Layer-1 subworkflow — route each isolate through the 15 stages, then
// aggregate every per-stage result.json into one per-isolate report + master row.
include { DB_CHECK }    from '../modules/stage00_dbcheck.nf'
include { BASECALL }    from '../modules/stage00_basecall.nf'
include { READ_QC }     from '../modules/stage01_readqc.nf'
include { ASSEMBLE }    from '../modules/stage02_assemble.nf'
include { SR_ASSEMBLE } from '../modules/stage02b_srassemble.nf'
include { MEDAKA }      from '../modules/stage03_polish.nf'
include { SRPOLISH }    from '../modules/stage03b_srpolish.nf'
include { DECONTAM }    from '../modules/stage04_decontam.nf'
include { ASSEMBLY_QC } from '../modules/stage05_assembly_qc.nf'
include { REPEATMASK }  from '../modules/stage06_repeatmask.nf'
include { ANNOTATE }    from '../modules/stage07_annotate.nf'
include { IDENTIFY }    from '../modules/stage08_identify.nf'
include { RESISTANCE }  from '../modules/stage09_resistance.nf'
include { MOBILE }      from '../modules/stage10_mobile.nf'
include { BGC }         from '../modules/stage11_bgc.nf'
include { NOVELTY }     from '../modules/stage12_novelty.nf'
include { EXTRAS }      from '../modules/stage13_extras.nf'
include { REPORT }      from '../modules/stage14_report.nf'
include { PROVENANCE }  from '../modules/stage15_provenance.nf'

workflow FUNGIFORGE {
  take:
    samples          // tuple(meta, filesMap[ont,r1,r2])
    run_info         // path: run_info.json (pipeline/nextflow/params/container map, from main.nf)

  main:
    // 00. database check — one task per run, BEFORE any compute: every isolate's first stage
    //     waits for it, so a missing database fails the run in seconds, not hours.
    DB_CHECK()
    reads0 = samples.combine(DB_CHECK.out.manifest).map { meta, f, _manifest -> tuple(meta, f.ont, f.r1, f.r2) }

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

    // 1. read QC + filtering (ONT chopper/NanoPlot; Illumina fastp — mode-aware)
    READ_QC(reads_bc)

    // Split isolates by assembly mode: long-read (ONT ± hybrid) vs short-read (Illumina only).
    reads_lr = READ_QC.out.reads.filter { row -> row[0].assembly_mode == 'longread' }
    reads_sr = READ_QC.out.reads.filter { row -> row[0].assembly_mode == 'shortread' }

    // 2. assembly — long-read (Flye --nano-hq + purge_dups) OR short-read (SPAdes)
    ASSEMBLE(reads_lr)
    SR_ASSEMBLE(reads_sr)

    // 3. polishing: Medaka (ONT long-read only) then Illumina Polypolish hybrid when
    //    short reads are present. Short-read-only assemblies skip Medaka (no ONT) and
    //    flow straight into Stage 03b, which records mode=illumina_only (high conf).
    ont_ch  = reads_lr.map { meta, ont, _r1, _r2 -> tuple(meta, ont) }
    ilmn_ch = READ_QC.out.reads.map { meta, _ont, r1, r2 -> tuple(meta, r1, r2) }
    MEDAKA(ASSEMBLE.out.assembly.join(ont_ch))
    pre_srpolish = MEDAKA.out.assembly.mix(SR_ASSEMBLE.out.assembly)
    SRPOLISH(pre_srpolish.join(ilmn_ch))

    // 4. decontamination + organelle split -> nuclear / mito
    DECONTAM(SRPOLISH.out.assembly)

    // 5. assembly QC + completeness (QUAST + compleasm/BUSCO)
    ASSEMBLY_QC(DECONTAM.out.nuclear)

    // 6. repeat modeling + soft-masking (RepeatModeler2 / RepeatMasker)
    REPEATMASK(DECONTAM.out.nuclear)

    // 7. eukaryotic gene prediction + functional annotation (Funannotate)
    ANNOTATE(REPEATMASK.out.masked)

    // 8. identification (ITS/LSU + genome ANI + MLST, GCPSR multi-locus)
    IDENTIFY(DECONTAM.out.nuclear)

    // 9. antifungal resistance (bespoke panel + cyp51A TR34/TR46 module).
    //    needs proteins (substitutions) + species + nuclear & GBK (promoter TR locus)
    RESISTANCE(ANNOTATE.out.proteins
                 .join(IDENTIFY.out.species)
                 .join(DECONTAM.out.nuclear)
                 .join(ANNOTATE.out.gbk)
                 .join(SRPOLISH.out.json))

    // 10-13: skippable stages — a skip flag truly skips the process (empty channel),
    // rather than running it emptily (which would also drag in its container).
    mobile_ch = channel.empty()
    if (!params.skip_mge)     { MOBILE(REPEATMASK.out.telib.join(DECONTAM.out.mito)); mobile_ch = MOBILE.out.json }

    bgc_ch = channel.empty()
    if (!params.skip_bgc)     { BGC(ANNOTATE.out.gbk); bgc_ch = BGC.out.json }

    novelty_ch = channel.empty()
    if (!params.skip_novelty) { NOVELTY(DECONTAM.out.nuclear.join(IDENTIFY.out.markers).join(IDENTIFY.out.json)); novelty_ch = NOVELTY.out.json }

    extras_ch = channel.empty()
    if (!params.skip_extras)  { EXTRAS(ANNOTATE.out.proteins.join(READ_QC.out.reads)); extras_ch = EXTRAS.out.json }

    // 14. aggregate every per-stage result.json per isolate -> report + master row
    all_json = READ_QC.out.json
      .mix(basecall_json, ASSEMBLE.out.json, SR_ASSEMBLE.out.json, MEDAKA.out.json,
           SRPOLISH.out.json, DECONTAM.out.json, ASSEMBLY_QC.out.json, REPEATMASK.out.json,
           ANNOTATE.out.json, IDENTIFY.out.json, RESISTANCE.out.json,
           mobile_ch, bgc_ch, novelty_ch, extras_ch)
      .map { meta, j -> tuple(meta.id, meta, j) }
      .groupTuple(by: 0)
      .map { _id, metas, jsons -> tuple(metas[0], jsons) }
    REPORT(all_json)

    // 15. run-level provenance: every stage JSON of every isolate + DB manifest + run metadata.
    //     Waits for all reports (masters.collect()) so it is the last task of the run.
    PROVENANCE(run_info,
               DB_CHECK.out.manifest,
               all_json.map { _meta, jsons -> jsons }.flatten().collect(),
               REPORT.out.master.collect())

  emit:
    master     = REPORT.out.master
    reports    = REPORT.out.report
    provenance = PROVENANCE.out.provenance
}
