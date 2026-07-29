// FungiForge Layer-1 subworkflow — route each isolate through the 15 stages, then
// aggregate every per-stage result.json into one per-isolate report + master row.
include { BASECALL }    from '../modules/stage00_basecall.nf'
include { READ_QC }     from '../modules/stage01_readqc.nf'
include { ASSEMBLE }    from '../modules/stage02_assemble.nf'
include { POLISH }      from '../modules/stage03_polish.nf'
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

workflow FUNGIFORGE {
  take:
    samples          // tuple(meta, filesMap[ont,r1,r2])

  main:
    reads0 = samples.map { meta, f -> tuple(meta, f.ont, f.r1, f.r2) }

    // 0. optional Dorado basecalling (ont column = pod5 dir when --basecall)
    reads_bc = params.basecall ? BASECALL(reads0).reads : reads0

    // 1. read QC + filtering (ONT chopper/NanoPlot; Illumina fastp)
    READ_QC(reads_bc)

    // 2. assembly (Flye --nano-hq by default) + purge_dups
    ASSEMBLE(READ_QC.out.reads)

    // 3. polishing: Medaka (+ optional Illumina Polypolish/POLCA hybrid branch)
    POLISH(ASSEMBLE.out.assembly.join(READ_QC.out.reads))

    // 4. decontamination + organelle split -> nuclear / mito
    DECONTAM(POLISH.out.assembly)

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
                 .join(ANNOTATE.out.gbk))

    // 10-13: skippable stages — a skip flag truly skips the process (empty channel),
    // rather than running it emptily (which would also drag in its container).
    mobile_ch = Channel.empty()
    if (!params.skip_mge)     { MOBILE(REPEATMASK.out.telib.join(DECONTAM.out.mito)); mobile_ch = MOBILE.out.json }

    bgc_ch = Channel.empty()
    if (!params.skip_bgc)     { BGC(ANNOTATE.out.gbk); bgc_ch = BGC.out.json }

    novelty_ch = Channel.empty()
    if (!params.skip_novelty) { NOVELTY(DECONTAM.out.nuclear.join(IDENTIFY.out.markers)); novelty_ch = NOVELTY.out.json }

    extras_ch = Channel.empty()
    if (!params.skip_extras)  { EXTRAS(ANNOTATE.out.proteins.join(READ_QC.out.reads)); extras_ch = EXTRAS.out.json }

    // 14. aggregate every per-stage result.json per isolate -> report + master row
    all_json = READ_QC.out.json
      .mix(POLISH.out.json, DECONTAM.out.json, ASSEMBLY_QC.out.json, REPEATMASK.out.json,
           ANNOTATE.out.json, IDENTIFY.out.json, RESISTANCE.out.json,
           mobile_ch, bgc_ch, novelty_ch, extras_ch)
      .map { meta, j -> tuple(meta.id, meta, j) }
      .groupTuple(by: 0)
      .map { id, metas, jsons -> tuple(metas[0], jsons) }
    REPORT(all_json)

  emit:
    master  = REPORT.out.master
    reports = REPORT.out.report
}
