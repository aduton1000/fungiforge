#!/usr/bin/env python3
"""Stage 08 — resolve a species call by multi-locus concordance (W2.3).

Evidence lines, each parsed into a per-locus call {species, genus, pident, ...}:
  ITS         vsearch blast6 against UNITE (the primary barcode; species-level >= 98.5 %)
  secondary   blastn (outfmt 6 + slen stitle) of each extracted locus (CaM, BenA, TEF1, RPB2,
              LSU D1/D2; bin/extract_markers.py) against the type-material reference set of that
              locus (<data_dir>/markers/<locus>.fasta, NCBI "sequence from type" records)
  genome      sourmash gather top match (optional, --genome_id)
  MLST        `mlst` scheme/ST/alleles (reported; scheme genus must agree with the call)

Concordance rule (GCPSR-style, the rule the plan calls "true concordance"):
  high    ITS gives a species AND at least one secondary line (locus or genome) gives the same
          species, and no secondary line gives a different species of the same genus
  medium  ITS-only species (no secondary line could be evaluated), or secondary-only species
          (no ITS), or ITS 94-98.5 %, or every secondary line is a tie that includes the ITS species
  low     genus only (ITS below species threshold and no secondary species), or the secondary
          lines contradict ITS (reported as "Genus sp." with a `discordant` flag)
  none    nothing classified
A tie is when reference species within --tie-margin % identity of the best hit differ; a tied
locus never confirms a species on its own (A. flavus / A. oryzae on CaM and BenA is the
canonical case) but is reported with its candidates.

Also chooses the species-aware BUSCO lineage (fungiforge/resources/busco_lineages.tsv) for
stage 08b and keeps the top ITS identity for novelty (stage 12)."""
from __future__ import annotations
import argparse, csv, json, os, re

ITS_SPECIES, ITS_GENUS = 98.5, 94.0
# species-level identity per locus; LSU D1/D2 is genus-level evidence only (it does not separate
# sibling species in Aspergillus, Penicillium or Fusarium: CEA10's D1/D2 is 99.9 % to several section
# Fumigati type strains), so it never carries a species call
LOCUS_SPECIES = {"CaM": 99.0, "BenA": 99.0, "TEF1": 99.0, "RPB2": 98.5, "LSU": None}
LOCUS_GENUS = 95.0
NON_SPECIES = re.compile(r"(^|\s)(sp|spp|cf|aff|nom)\.?$|Incertae|unidentified|uncultured|unclassified", re.I)
MLST_SCHEME_GENUS = {"afumigatus": "Aspergillus", "calbicans": "Candida", "cglabrata": "Candida", "ctropicalis": "Candida",
                     "ckrusei": "Candida", "cauris": "Candida"}


def parse_unite_b6(path):
    """Best UNITE hit from a vsearch --blast6out file -> (species, genus, pident, sh)."""
    if not path or not os.path.exists(path):
        return None
    best = None
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        try:
            pid = float(f[2])
        except ValueError:
            continue
        target = f[1]
        if best is None or pid > best[0]:
            best = (pid, target)
    if not best:
        return None
    pid, target = best
    sp = re.search(r"s__([A-Za-z0-9_.\-]+)", target)
    ge = re.search(r"g__([A-Za-z0-9_.\-]+)", target)
    sh = re.search(r"(SH\d+\.\d+FU)", target)
    species = sp.group(1).replace("_", " ") if sp else None
    # "Genus sp" / "Genus_sp." / Incertae sedis / unidentified are not species-level names
    if species and NON_SPECIES.search(species):
        species = None
    return {"pident": round(pid, 2),
            "species": species,
            "genus": ge.group(1) if ge else None,
            "sh": sh.group(1) if sh else None}


def top_gather(path):
    if not path or not os.path.exists(path):
        return None
    try:
        rows = list(csv.DictReader(open(path)))
    except Exception:
        return None
    if not rows:
        return None
    name = rows[0].get("name") or rows[0].get("match_name") or ""
    m = re.search(r"([A-Z][a-z]+ [a-z]+)", name)
    return {"species": m.group(1) if m else name[:60], "raw": name[:120]}


def species_from_title(title):
    """'Aspergillus flavus strain CBS 100927 calmodulin (CaM) gene' -> ('Aspergillus flavus', 'Aspergillus');
    'Aspergillus sp. XYZ ...' -> (None, 'Aspergillus'). Brackets ([Candida] auris) are stripped."""
    t = re.sub(r"[\[\]]", "", title or "").strip()
    t = re.sub(r"^(UNVERIFIED(_[A-Z]+)?|TPA(_[a-z]+)?|PREDICTED|MAG):\s*", "", t)   # record-status prefixes
    t = re.sub(r"^[A-Za-z_]*\d[\w.]*\s+", "", t)   # leading accession (stitle in BLAST output starts with it)
    m = re.match(r"([A-Z][a-z]+)\s+([a-z][a-z\-]+)", t)
    if not m:
        return None, None
    genus, epithet = m.group(1), m.group(2)
    if NON_SPECIES.search(epithet) or epithet in ("sp", "spp", "cf", "aff"):
        return None, genus
    return f"{genus} {epithet}", genus


def parse_locus_b6(path, locus, tie_margin=0.3, min_scov=0.5):
    """Best hits of one extracted locus vs its type-material set.
    outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen stitle".
    Returns None when the file is absent/empty; otherwise a locus call with tie handling."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    hits = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            try:
                pid, alen, bits, slen = float(f[2]), int(f[3]), float(f[9]), int(f[11])
            except ValueError:
                continue
            scov = alen / slen if slen else 0.0
            sp, ge = species_from_title(f[12])
            hits.append({"accession": f[1], "pident": pid, "alen": alen, "bitscore": bits, "scov": round(scov, 3),
                         "species": sp, "genus": ge, "title": f[12][:100]})
    hits = [h for h in hits if h["scov"] >= min_scov] or hits
    if not hits:
        return None
    hits.sort(key=lambda h: (-h["bitscore"], -h["pident"]))
    best = hits[0]
    thr = LOCUS_SPECIES.get(locus, 99.0)   # None: genus-level marker
    # species named within tie_margin of the best identity (species-level hits only)
    top_pid = max(h["pident"] for h in hits)
    cands = []
    for h in hits:
        if h["species"] and h["pident"] >= top_pid - tie_margin and h["species"] not in cands:
            cands.append(h["species"])
    call = {"locus": locus, "pident": round(top_pid, 2), "best_accession": best["accession"], "best_title": best["title"],
            "scov": best["scov"], "n_hits": len(hits), "candidates": cands, "species": None, "genus": best["genus"], "level": "none"}
    if thr is not None and top_pid >= thr and cands:
        if len(cands) == 1:
            call["species"], call["level"] = cands[0], "species"
        else:
            call["level"] = "tie"
    elif top_pid >= LOCUS_GENUS and best["genus"]:
        call["level"] = "genus"
    return call


def parse_mlst(path):
    """`mlst` TSV line: file, scheme, ST, locus(allele)...; '-' when no scheme matched."""
    if not path or not os.path.exists(path):
        return None
    with open(path) as fh:
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    if not rows:
        return None
    f = rows[0]
    if len(f) < 3:
        return None
    scheme, st = f[1], f[2]
    alleles = {}
    for a in f[3:]:
        m = re.match(r"([^()]+)\(([^)]*)\)", a)
        if m:
            alleles[m.group(1)] = m.group(2)
    return {"scheme": None if scheme in ("-", "") else scheme, "st": None if st in ("-", "") else st, "alleles": alleles,
            "n_loci": len(alleles), "n_alleles_found": sum(1 for v in alleles.values() if v and v != "-")}


def choose_lineage(species, genus, map_path, default="fungi_odb10"):
    """Most specific BUSCO lineage for the call: species row, then genus row, else default."""
    if not map_path or not os.path.exists(map_path):
        return default, "no map"
    rows = {}
    with open(map_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] != "taxon":
                rows[f[0]] = f[1]
    if species and species in rows:
        return rows[species], f"species:{species}"
    if genus and genus in rows:
        return rows[genus], f"genus:{genus}"
    return default, "default"


def concordance(unite, loci, gather):
    """Apply the rule in the module docstring. Returns species, confidence, method, flags, evidence."""
    flags, lines = [], []          # lines: (name, species, genus) for every secondary line with a species or genus
    its_species = unite["species"] if unite and unite["pident"] >= ITS_GENUS else None
    its_genus = unite["genus"] if unite else None
    its_strong = bool(its_species and unite["pident"] >= ITS_SPECIES)
    for name, call in loci.items():
        if not call:
            continue
        if call["level"] == "species":
            lines.append((name, call["species"], call["genus"]))
        elif call["level"] == "tie":
            lines.append((name, None, call["genus"]))
            flags.append(f"tie:{name}={'/'.join(call['candidates'])}")
        elif call["level"] == "genus":
            lines.append((name, None, call["genus"]))
    if gather and gather.get("species"):
        gs = gather["species"]
        lines.append(("genome", gs, gs.split()[0]))
    sec_species = [(n, s) for n, s, _ in lines if s]
    agree = [n for n, s in sec_species if its_species and s.lower() == its_species.lower()]
    disagree = [(n, s) for n, s in sec_species if its_species and s.lower() != its_species.lower()]
    tie_with_its = [n for n, c in loci.items() if c and c["level"] == "tie" and its_species
                    and any(x.lower() == its_species.lower() for x in c["candidates"])]
    evidence = {"its_species": its_species, "its_identity": unite["pident"] if unite else None,
                "secondary_agree": agree, "secondary_disagree": [f"{n}={s}" for n, s in disagree], "ties_including_its": tie_with_its}
    # 1. ITS species with concordant secondary evidence
    if its_species and agree and not disagree:
        conf = "high" if its_strong or len(agree) >= 2 else "medium"
        return its_species, conf, f"ITS+{'+'.join(agree)}(concordant)", flags, evidence
    # 2. secondary lines contradict ITS -> genus with a flag
    if its_species and disagree:
        genus = its_genus or its_species.split()[0]
        flags.append("discordant:" + ",".join(f"{n}={s}" for n, s in disagree))
        return f"{genus} sp.", "low", f"ITS={its_species};" + ";".join(f"{n}={s}" for n, s in disagree) + "(discordant)", flags, evidence
    # 3. ITS species alone (no secondary species-level line)
    if its_species:
        if tie_with_its:
            flags.append("secondary_ties_include_its_species")
            method = f"ITS+{'+'.join(tie_with_its)}(tie)"
        elif lines:
            flags.append("no_secondary_species_call")
            method = "ITS/UNITE(secondary loci unresolved)"
        else:
            flags.append("secondary_loci_unavailable")
            method = "ITS/UNITE(no secondary loci)"
        if not its_strong:
            method += "(below-species-threshold)"
        return its_species, "medium", method, flags, evidence
    # 4. no ITS species: secondary lines alone
    if sec_species:
        names = {}
        for n, s in sec_species:
            names.setdefault(s.lower(), []).append((n, s))
        if len(names) == 1:
            (n_list,) = names.values()
            sp = n_list[0][1]
            conf = "medium" if len(n_list) >= 2 else "low" if not its_genus else "medium"
            flags.append("no_its_species")
            return sp, conf, "+".join(n for n, _ in n_list) + "(secondary-only)", flags, evidence
        flags.append("discordant:" + ",".join(f"{n}={s}" for n, s in sec_species))
        genus = its_genus or sec_species[0][1].split()[0]
        return f"{genus} sp.", "low", ";".join(f"{n}={s}" for n, s in sec_species) + "(discordant)", flags, evidence
    # 5. genus only
    genera = [g for _, _, g in lines if g] + ([its_genus] if its_genus else [])
    if genera:
        g = its_genus or genera[0]
        src = "ITS/UNITE" if its_genus else "+".join(n for n, _, gg in lines if gg)
        return f"{g} sp.", "low", f"{src}(genus-only)", flags, evidence
    return "unknown", "none", "no_reference", flags, evidence


def locus_doc(name, call, extracted):
    if call:
        return dict(call, extracted=extracted.get(name, {}).get("found"))
    return {"locus": name, "extracted": extracted.get(name, {}).get("found", False), "level": "none", "species": None}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--its")                  # ITS fasta (presence flag)
    ap.add_argument("--unite-b6")             # vsearch blast6 vs UNITE
    ap.add_argument("--gather")               # optional sourmash gather csv
    ap.add_argument("--markers-json")         # bin/extract_markers.py summary
    ap.add_argument("--locus-b6", action="append", default=[], metavar="LOCUS=FILE",
                    help="blastn outfmt-6(+qlen slen stitle) of one extracted locus vs its reference set (repeatable)")
    ap.add_argument("--mlst")                 # mlst TSV output
    ap.add_argument("--lineage-map")          # fungiforge/resources/busco_lineages.tsv
    ap.add_argument("--tie-margin", type=float, default=0.3)
    ap.add_argument("--out-species", required=True)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-lineage", default=None)
    a = ap.parse_args()

    its_present = bool(a.its and os.path.exists(a.its) and os.path.getsize(a.its) > 0)
    unite = parse_unite_b6(a.unite_b6)
    gather = top_gather(a.gather)
    extracted = {}
    if a.markers_json and os.path.exists(a.markers_json):
        try:
            extracted = json.load(open(a.markers_json)).get("loci", {})
        except Exception:  # noqa: BLE001
            extracted = {}
    loci = {}
    for kv in a.locus_b6:
        if "=" not in kv:
            continue
        name, path = kv.split("=", 1)
        loci[name] = parse_locus_b6(path, name, a.tie_margin)
    mlst = parse_mlst(a.mlst)

    species, conf, method, flags, evidence = concordance(unite, loci, gather)
    genus = None if species in (None, "unknown") else species.split()[0]
    if mlst and mlst.get("scheme"):
        exp = MLST_SCHEME_GENUS.get(mlst["scheme"])
        if exp and genus and exp != genus:
            flags.append(f"mlst_scheme_genus_mismatch:{mlst['scheme']}")
    lineage, lineage_basis = choose_lineage(None if species.endswith(" sp.") else species, genus, a.lineage_map)

    open(a.out_species, "w").write(f"{species}\t{conf}\n")
    if a.out_lineage:
        open(a.out_lineage, "w").write(lineage + "\n")
    doc = {"sample": a.sample, "stage": "identify", "species": species, "confidence": conf, "method": method,
           "flags": flags, "its_present": its_present, "its_identity": unite["pident"] if unite else None,
           "unite_hit": unite, "gather_top": gather,
           "loci": {n: locus_doc(n, loci.get(n), extracted) for n in sorted(set(loci) | set(extracted))},
           "concordance": evidence, "mlst": mlst,
           "busco_lineage": lineage, "busco_lineage_basis": lineage_basis}
    json.dump(doc, open(a.out_json, "w"), indent=2)
    loci_txt = ", ".join("%s:%s" % (n, c.get("level")) for n, c in doc["loci"].items()) or "none"
    mlst_txt = "%s ST%s" % (mlst["scheme"], mlst["st"]) if mlst and mlst.get("scheme") else "none"
    print("[id_classify] %s: %s (%s, %s, ITS%%=%s; loci %s; MLST %s; lineage %s)"
          % (a.sample, species, conf, method, doc["its_identity"], loci_txt, mlst_txt, lineage))


if __name__ == "__main__":
    main()
