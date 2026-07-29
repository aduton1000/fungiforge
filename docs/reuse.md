# fungiforge — reused vs bespoke

The engineering stance (as in callforge): stand on validated tools; add the bespoke
layers fungal One Health work actually needs. Per stage:

| Stage | Reused / adapted | Built bespoke |
|---|---|---|
| 00 basecall | Dorado (ONT) | — |
| 01 read QC | NanoPlot, chopper, fastp, Kraken2, GenomeScope2 | QC JSON contract |
| 02 assembly | Flye / Canu / Raven, purge_dups, minimap2 | assembler switch + size-vs-expectation check |
| 03 polish | Medaka, Polypolish, POLCA | **hybrid branch + polishing-confidence flag feeding resistance** |
| 04 decontam | Kraken2/tiara, BlobTools, GetOrganelle/oatk | nuclear/mito split contract |
| 05 QC | QUAST, compleasm/BUSCO, Merqury | MIMAG-style QC-pass gate |
| 06 repeats | RepeatModeler2, RepeatMasker (dfam/tetools) | TE-class + RIP summary for Stage 10 |
| 07 annotation | **Funannotate** (GeneMark/Augustus/SNAP/EVM), eggNOG, dbCAN, InterProScan | evidence-mode wrapper, InterProScan toggle |
| 08 identify | ITSx, barrnap, sourmash, skani, MLST, UNITE | **GCPSR multi-locus concordance call** (`id_classify.py`) |
| 09 resistance | BLAST/DIAMOND/HMMER, FungAMR/MARDy references | **the whole caller**: curated panel + ortholog/substitution/LoF/GOF logic + **cyp51A TR34/TR46 structural detector** (`af_resistance.py`, `cyp51a_TR.py`) |
| 10 mobile | RepeatModeler classes, geNomad, RVDB | fungal-MGE merge (TE + mycovirus + mito introns) |
| 11 BGC | fungiSMASH (antiSMASH 8), BiG-SCAPE, MIBiG | per-isolate BGC summary + cross-isolate novelty layer |
| 12 novelty | skani, IQ-TREE, BUSCO SCOs | ANI+ITS+phylo novelty verdict with GCPSR caveat |
| 13 extras | nQuire, Smudgeplot, dbCAN, EffectorP, SignalP | mating-type + secretome/virulence tally |
| 14 report | jinja2 | **master-table schema** (the Layer-2 handoff) |
| Layer 2 | vegan, lme4, ggplot2 (as bacterial r-analysis) | 8 fungal objectives + compartment framing |

There is **no ResFinder/GTDB for fungi** — Stages 08/09/12 are necessarily bespoke, which is where fungiforge adds the most value.
