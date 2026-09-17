#!/usr/bin/env python3
"""Stage 12 — novelty verdict (W2.9). Genome ANI (skani against the staged reference genome set)
is the strong signal; the ITS identity and the secondary-locus concordance from stage 08 back it.

  novelty_call.py --sample ID --skani skani.tsv --identify-json ID.identify.json [--manifest manifest.tsv]
                  [--min-af 15] --out ID.novelty.json

Rules (nearest genome = highest ANI with aligned fraction >= --min-af on either genome):
  ANI >= 95 %                   known_species (the nearest genome's species; a flag when stage 08
                                named a different species)
  90 <= ANI < 95                candidate_novel_species (sister of the nearest species)
  ANI < 90, or AF too low       candidate_novel_or_unrepresented (nothing close in the set: a novel
                                lineage, or a genus the set lacks — the manifest tells which)
  no reference set              ITS fallback: < 98.5 % => candidate novel, else known; undetermined without ITS
Species-level agreement of ITS + a secondary locus (stage 08 confidence high) with the nearest
genome's species raises the verdict to `known_species` even at ANI 93-95 % (intraspecific ANI
in some fungi is broad) and is stated in `basis`."""
from __future__ import annotations
import argparse, csv, json, os

KNOWN, NOVEL, FAR = 95.0, 90.0, 90.0


def read_skani(path, min_af=15.0):
    """skani dist TSV -> rows sorted by ANI (desc), with the aligned-fraction filter applied."""
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            try:
                ani = float(p[ci.get("ANI", 2)])
                af_r = float(p[ci["Align_fraction_ref"]]) if "Align_fraction_ref" in ci else None
                af_q = float(p[ci["Align_fraction_query"]]) if "Align_fraction_query" in ci else None
            except (ValueError, IndexError):
                continue
            ref_file = p[ci.get("Ref_file", 0)]
            ref_name = p[ci["Ref_name"]] if "Ref_name" in ci and len(p) > ci["Ref_name"] else ref_file
            rows.append({"ani": ani, "af_ref": af_r, "af_query": af_q, "ref_file": ref_file, "ref_name": ref_name,
                         "af_ok": max(af_r or 0, af_q or 0) >= min_af})
    rows.sort(key=lambda r: -r["ani"])
    return rows


def top_skani(path):
    """Backwards-compatible: (ANI, ref_name) of the best hit, or None."""
    rows = read_skani(path, min_af=0)
    return (rows[0]["ani"], rows[0]["ref_name"]) if rows else None


def read_manifest(path):
    m = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                m[r.get("accession", "")] = r
                if r.get("file"):
                    m[r["file"]] = r
    return m


def describe_ref(row, manifest):
    base = os.path.basename(row["ref_file"])
    acc = base.replace(".fna", "").replace(".fa", "").replace(".fasta", "")
    info = manifest.get(acc) or manifest.get(base) or {}
    return {"accession": acc, "organism": info.get("organism") or row["ref_name"], "species": info.get("species") or None, "category": info.get("category")}


def verdict(rows, ident, manifest):
    its_pid = ident.get("its_identity")
    species = ident.get("species")
    conf = ident.get("confidence")
    usable = [r for r in rows if r["af_ok"]]
    out = {"species": species, "its_identity": its_pid, "id_confidence": conf, "genome_ani": None, "nearest": None, "flags": []}
    if rows:
        best = usable[0] if usable else rows[0]
        near = describe_ref(best, manifest)
        out.update(genome_ani=round(best["ani"], 2), aligned_fraction=round(max(best["af_ref"] or 0, best["af_query"] or 0), 1), nearest=near,
                   n_references_hit=len(rows))
        same = bool(species and near["species"] and species.lower() == near["species"].lower())
        if not usable:
            out.update(verdict="candidate_novel_or_unrepresented", basis="genome_ani(skani): aligned fraction too low for every reference",
                       note=f"best ANI {best['ani']:.1f} % but aligned fraction below the floor: no reference genome resembles this isolate (novel lineage, or a genus the reference set lacks)")
        elif best["ani"] >= KNOWN:
            out.update(verdict="known_species", basis="genome_ani(skani)", note=f"ANI {best['ani']:.1f} % to {near['organism']} ({near['accession']})")
            if species and near["species"] and not same:
                out["flags"].append(f"identification_disagrees:{species}!={near['species']}")
        elif best["ani"] >= NOVEL:
            if same and conf == "high":
                out.update(verdict="known_species", basis="genome_ani(skani)+ITS+secondary_loci",
                           note=f"ANI {best['ani']:.1f} % to {near['organism']} is below 95 % but ITS and a secondary locus agree on {species}: broad intraspecific ANI, not a new species")
            else:
                out.update(verdict="candidate_novel_species", basis="genome_ani(skani)",
                           note=f"ANI {best['ani']:.1f} % to the nearest reference {near['organism']} ({near['accession']}): sister lineage of that species; confirm by multi-locus GCPSR + polyphasic description")
        else:
            out.update(verdict="candidate_novel_or_unrepresented", basis="genome_ani(skani)",
                       note=f"ANI {best['ani']:.1f} % to the nearest reference {near['organism']}: no close genome in the set (novel lineage, or an unrepresented genus)")
        if its_pid is not None and its_pid < 98.5 and out["verdict"] == "known_species":
            out["flags"].append(f"its_below_species_threshold:{its_pid}")
    elif its_pid is not None:
        novel = its_pid < 98.5
        out.update(verdict="candidate_novel_species" if novel else "known_species", basis="ITS_distance(UNITE)",
                   note=(f"best ITS identity {its_pid}% (< 98.5% species threshold) — candidate novel; confirm with genome ANI + multi-locus GCPSR" if novel
                         else f"ITS identity {its_pid}% is within the species range"))
    else:
        out.update(verdict="undetermined", basis="none", note="no ITS hit and no genome-ANI reference set staged")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True); ap.add_argument("--skani"); ap.add_argument("--identify-json"); ap.add_argument("--manifest", default=None)
    ap.add_argument("--min-af", type=float, default=15.0); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    ident = {}
    if a.identify_json and os.path.exists(a.identify_json):
        try:
            ident = json.load(open(a.identify_json))
        except Exception:  # noqa: BLE001
            ident = {}
    rows = read_skani(a.skani, a.min_af)
    out = {"sample": a.sample, "stage": "novelty", **verdict(rows, ident, read_manifest(a.manifest))}
    out["top_references"] = [dict(describe_ref(r, read_manifest(a.manifest)), ani=r["ani"], af_ref=r["af_ref"], af_query=r["af_query"]) for r in rows[:5]]
    json.dump(out, open(a.out, "w"), indent=2)
    print(f"[novelty_call] {a.sample}: {out['verdict']} ({out['basis']})" + (f" nearest {out['nearest']['organism']} ANI {out['genome_ani']}" if out.get('nearest') else ""))


if __name__ == "__main__":
    main()
