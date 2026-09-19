#!/usr/bin/env python3
"""Functional-annotation coverage from funannotate's annotations table (W2.5, stage 07).

  annotate_stats.py --sample ID --annotations fun/annotate_results/ID.annotations.txt
                    [--proteins ID.proteins.faa] [--training training.json] [--eggnog yes|no|skipped]
                    [--interproscan yes|no|skipped] [--genemark yes|no] --out ID.annotate.json

Counts proteins and the fraction with a PFAM domain, an InterPro entry, GO terms, an eggNOG
orthologue, a product name other than "hypothetical protein", plus secreted / CAZyme / protease /
BUSCO marks, so annotation quality is comparable between isolates and between runs."""
from __future__ import annotations
import argparse, csv, json, os


# Record types in funannotate's annotations table that are not protein-coding transcripts.
NONCODING_FEATURES = {"trna", "rrna", "ncrna", "snrna", "snorna", "tmrna", "misc_rna", "repeat_region"}


def count_proteins(faa):
    if not faa or not os.path.exists(faa):
        return None
    return sum(1 for l in open(faa) if l.startswith(">"))


def annotation_stats(path):
    """Coverage from funannotate's annotations table, with the reason when there is none.

    Returning a bare {} left every pct_* column as NA in the master row with nothing to say
    whether the table was missing, empty, or filtered away — indistinguishable from a stage that
    never ran. Every return now carries `annotations_note`, and the header actually seen.
    """
    if not path or not os.path.exists(path):
        return {"annotations_note": "annotations table not found: %s" % (path or "<none>")}
    with open(path, encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = list(reader.fieldnames or [])
        all_rows = list(reader)
    # The "Feature" column names the record type, and which value marks a protein-coding row
    # depends on the funannotate build: 1.8.17 writes "mRNA" (the table is per transcript, hence
    # the TranscriptID column), older notes say "CDS". Keeping only "CDS" discarded every row of a
    # real 9,844-gene table and left all four coverage columns NA. Exclude the non-coding types
    # instead, which stays correct whichever label the build uses.
    rows = [r for r in all_rows
            if (r.get("Feature") or "").strip().lower() not in NONCODING_FEATURES]
    n = len(rows)
    if not n:
        return {"n_annotated": 0, "annotations_header": header,
                "annotations_note": "no usable rows in %s (%d line(s) read, header: %s)"
                                    % (path, len(all_rows), ", ".join(header) or "<empty>")}

    def has(col):
        return sum(1 for r in rows if (r.get(col) or "").strip())

    def pct(k):
        return round(100.0 * k / n, 1)
    named = sum(1 for r in rows if (r.get("Product") or "").strip().lower() not in ("", "hypothetical protein"))
    return {"n_annotated": n, "annotations_header": header, "pct_pfam": pct(has("PFAM")), "pct_interpro": pct(has("InterPro")), "pct_go": pct(has("GO Terms")),
            "pct_eggnog": pct(has("EggNog")), "pct_named_product": pct(named), "n_secreted": has("Secreted"),
            "n_cazyme": has("CAZyme"), "n_protease": has("Protease"), "n_busco": has("BUSCO"), "n_ec": has("EC_number")}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True); ap.add_argument("--annotations", default=None); ap.add_argument("--proteins", default=None)
    ap.add_argument("--training", default=None); ap.add_argument("--eggnog", default="no"); ap.add_argument("--interproscan", default="no")
    ap.add_argument("--genemark", default="no"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    doc = {"sample": a.sample, "stage": "annotate", "proteins": f"{a.sample}.proteins.faa",
           "n_proteins": count_proteins(a.proteins), "eggnog": a.eggnog, "interproscan": a.interproscan, "genemark": a.genemark}
    doc.update(annotation_stats(a.annotations))
    if a.training and os.path.exists(a.training):
        try:
            doc["training"] = json.load(open(a.training))
        except Exception:  # noqa: BLE001
            pass
    json.dump(doc, open(a.out, "w"), indent=2)
    print(f"[annotate_stats] {a.sample}: {doc['n_proteins']} proteins; PFAM {doc.get('pct_pfam', 'NA')}%, GO {doc.get('pct_go', 'NA')}%, "
          f"eggNOG {doc.get('pct_eggnog', 'NA')}%, InterPro {doc.get('pct_interpro', 'NA')}%")


if __name__ == "__main__":
    main()
