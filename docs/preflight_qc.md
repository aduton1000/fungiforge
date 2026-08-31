# Preflight QC (`bin/preflight_qc.sh`)

A fast, **read-only** triage you run on freshly arrived Oxford Nanopore (ONT)
data **before** committing to a long fungiforge run. On this Apple-Silicon Mac
the pipeline is emulated (linux/amd64) and a real run is multi-hour, so the goal
is to answer three questions per sample in a few minutes:

1. **Is there enough usable data to assemble a genome?** — GO / MARGINAL / NO-GO.
2. **Is this actually whole-genome shotgun, or an ITS amplicon run** that
   fungiforge is the wrong tool for?
3. **What is each isolate?** — a read-based species ID that works **even when
   coverage is too low for assembly to succeed** (the salvage step).

Nothing here launches fungiforge or modifies the reads/DB; it only reads data,
merges per-barcode chunks, and writes a summary + a samplesheet stub.

---

## Prerequisites

- **Docker** running, with the base image `aduton1000/fungiforge:0.1.0`
  (`docker images | grep fungiforge`). All tools run inside this image, so you
  do **not** need seqkit/vsearch on the host.
- The **fungiforge reference DB**, pointed to by `--db` or `$FUNGIFORGE_DB`,
  containing `unite/`, `funannotate/`, `antismash/`, `busco/lineages/`,
  `fungamr/reference_proteins.faa`. Build it with `bin/fetch_references.sh`.
- Enough free disk for the merged reads (roughly the size of the raw fastqs).

**Machine-specific facts baked into the script** (don't remove them): every
`docker run` passes `--platform linux/amd64` and `-e HOME=/tmp -e USER=fungiforge`
(the arm64/amd64 emulation + host-UID gotcha); binaries are called by full path
under `/opt/conda/bin/`; there is **no kraken2** (ID uses vsearch vs UNITE) and
**no bc** (arithmetic is done in awk); globs are expanded inside the container.

---

## Quick start

```bash
export FUNGIFORGE_DB=/path/to/fungiforge_db

bin/preflight_qc.sh \
    --input "/path/to/Passed samples" \
    --merged-out ./preflight_merged \
    --outdir ./preflight_qc \
    --genome-size 30000000        # molds ~30 Mb; yeasts ~15 Mb
```

- The **`--input`** path may contain spaces (ONT often writes a `Passed samples`
  dir). It should contain one subdir per barcode (`barcode74/…*.fastq.gz`).
- **`--merged-out` and `--outdir` must be space-free** — Flye/funannotate/
  antiSMASH break on spaces, so the merged reads fungiforge will consume must
  live on a clean path.
- Already have merged per-sample fastqs? Add **`--skip-merge`** and point
  `--input` at the dir of `*.fastq.gz`.
- Re-running is safe; existing merged files are skipped unless you pass
  **`--force`**.

See `bin/preflight_qc.sh --help` for all options (`--sample-prefix`,
`--amplicon-subn`, `--id-max-reads`, `--image`, …) and `head -80
bin/preflight_qc.sh` for the inline reference.

---

## The 8 checks and how to read them

1. **Environment readiness.** Confirms the docker image exists, the DB root has
   every required subdir, resolves the versioned UNITE fasta
   (`sh_general_release_dynamic_*.fasta`, ignoring the `_dev` copy), reports free
   disk, and warns **loudly** if the DB or merged-reads path contains a space. It
   also counts macOS AppleDouble `._*` files in the DB (they break
   antiSMASH/funannotate) — clean them with `dot_clean -m "$FUNGIFORGE_DB"`.
   *Read it as a go/no-go gate: any `[MISSING]` or space warning must be fixed
   before a real run.*

2. **Merge per-barcode chunks.** ONT writes many `*.fastq.gz` per barcode;
   gzip members concatenate with plain `cat`, so each barcode dir becomes one
   `<SAMPLE>.ont.fastq.gz`. Default naming: `barcode74 -> BC74` (change the `BC`
   with `--sample-prefix`). *Read it as a record of what got merged into what.*

3. **Read stats** (`seqkit stats -a`). num_seqs, sum_len, avg_len, **read N50**,
   Q20%, Q30%, GC%. *The headline "how much / how good" table.* Long read N50 and
   high Q20/Q30 are good; a tiny N50 or low Q means basecalling/library trouble.

4. **Read-length buckets** (`<300 / 300-500 / 500-800 / 800-1200 / 1200-3000 /
   >3000` bp). *Read the **shape**:* a single narrow spike ⇒ amplicon-like; a
   broad smear reaching long reads ⇒ shotgun; a big pile under 300 bp ⇒
   short/junk reads that won't assemble and inflate counts.

5. **Amplicon-vs-WGS discriminator.** Subsamples ~2000 reads (≥200 bp) and maps
   them to UNITE at `--id 0.75`, reporting the **fraction that hit ITS** and how
   many distinct species those hits span.
   - **High hit fraction (>50%) concentrated on one species ⇒ ITS AMPLICON run**
     — fungiforge is the **wrong tool**; use an ITS/metabarcoding workflow.
   - **Very low fraction (<~5%), scattered across many taxa ⇒ WGS shotgun** —
     correct input for fungiforge.
   - Anything in between is flagged **AMBIGUOUS**; cross-read checks 4 and 7.

6. **Assembly-usable yield & rough coverage.** Reports reads/bp at ≥1000 and
   ≥3000 bp and fold-coverage for `--genome-size`. The verdict uses the **≥1 kb**
   coverage (short reads don't help ONT assembly). **Heuristic** thresholds:
   `≥30x ⇒ GO`, `10–30x ⇒ MARGINAL` (expect a fragmented assembly), `<10x ⇒
   NO-GO`. These are rules of thumb, not guarantees; a yeast genome is ~15 Mb, so
   for yeasts halve `--genome-size` (or read the coverage as roughly double).

7. **Read-based species ID (the salvage step).** Reads ≥300 bp are mapped to
   UNITE at `--id 0.85`; the top 3 target species by hit count are printed. This
   gives an identity **even for NO-GO samples**, because fungal rDNA/ITS is
   tandemly repeated ~50–200×, so ITS-bearing reads are enriched and present long
   before genome coverage is assembly-worthy. *Read the top species + its hit
   count; a clear #1 with many hits is a confident ID, a near-tie between genera
   means treat it as tentative.* (Species labels come from the first `|`-field of
   the UNITE header, e.g. `Aspergillus_fumigatus`.)

8. **Artifacts.** Writes `preflight_summary.tsv` (one row per sample),
   `samplesheet.stub.csv` (fungiforge columns, `ont_fastq` pre-filled to the
   merged files, metadata as `NA` for you to fill), and a **saveable QC report**
   in two formats (see below).

---

## Outputs (what lands in `--outdir`)

| file | what it is |
|---|---|
| `preflight_qc_report.md` | Human-readable report: run overview, environment readiness (PASS/WARN flags), per-sample QC table, figures, amplicon-vs-WGS calls, and a prominent **recommendations** section. Good for pasting into a repo / lab notes. |
| `preflight_qc_report.html` | The **same report, self-contained** — every figure is embedded as a base64 data-URI, so the single file is portable (no external assets) and prints/emails cleanly. Light/dark aware. |
| `figs/*.png` | The rendered figures, also saved as standalone PNGs (referenced by the Markdown report). |
| `preflight_summary.tsv` | Machine-readable one-row-per-sample summary. |
| `samplesheet.stub.csv` | Ready-to-edit fungiforge samplesheet (fill in `compartment/facility/season`, add Illumina columns for the hybrid branch). |

**Figures** (via matplotlib inside the base image; `NanoPlot` is also present in
the image if you want richer per-read plots):

- per-sample **read-length distribution** (the six buckets, as a bar chart);
- a cross-sample **coverage** bar chart (all reads vs ≥1 kb, with the GO
  threshold drawn as a reference line);
- a cross-sample **read-N50** bar chart (reference line at the ~5 kb "good ONT"
  mark);
- a cross-sample **%ITS-hit** bar chart (the amplicon-vs-WGS visual).

If matplotlib is unavailable for any reason, figure rendering is skipped and the
report **degrades gracefully to compact ASCII bar charts** — the run never fails
because plotting is missing.

The **recommendations** section states, per sample, a GO / MARGINAL / NO-GO
verdict *with its reason and the exact thresholds behind it*, then an overall
"how to proceed" narrative: which samples to assemble, which to attempt only
with relaxed settings (e.g. a lower `--ont_min_len`), which to treat as ID-only,
and which to flag for re-sequencing (short read-N50 or low yield) or route away
from fungiforge entirely (amplicon libraries).

---

## Worked example — reading the summary table

```
sample  num_reads  sum_bp     read_N50  GC     pct_ITS_hits  cov_all  cov_ge1kb  assembly_verdict  top_species
BC74    412033     2.98e9     11840     49.1   0.8           99.3     94.1       GO                Aspergillus_fumigatus(3187)
BC75    38150      210000000  6420      47.6   1.9           7.0      5.4        NO-GO             Aspergillus_flavus(214)
BC76    950120     640000000  520       52.0   71.4          21.3     2.1        NO-GO             Penicillium_citrinum(6620)
```

- **BC74 — GO.** ~94× usable (≥1 kb) coverage, long N50, ITS fraction <1%
  (shotgun). This is a strong assembly candidate; run fungiforge and trust the
  `Aspergillus fumigatus` call to be confirmed by the genome-based ID.
- **BC75 — NO-GO but identified.** Only ~5× usable coverage, so de novo assembly
  will likely fail or fragment badly — **don't spend the multi-hour run hoping
  for a genome.** But check 7 still IDs it as *Aspergillus flavus* from 214 ITS
  reads: you get an identity without assembly. Re-sequence for more depth if a
  genome is needed.
- **BC76 — NO-GO, and the wrong assay.** `pct_ITS_hits` is 71% with a read-length
  profile dominated by a ~520 bp N50 (short, amplicon-shaped) concentrated on one
  species: this is an **ITS amplicon library**, not shotgun. fungiforge is the
  wrong tool here — route it to an ITS/metabarcoding workflow. The high "coverage"
  is an artifact of amplicon read counts, not genome breadth.

**Rule of thumb:** run fungiforge on the **GO** (and case-by-case **MARGINAL**)
shotgun samples; for **NO-GO** samples, take the check-7 identity and decide
whether to re-sequence; for anything flagged **amplicon**, don't use fungiforge
at all.
