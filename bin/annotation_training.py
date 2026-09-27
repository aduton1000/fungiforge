#!/usr/bin/env python3
"""Choose funannotate's training start for an isolate from its species call (W2.5).

  annotation_training.py --species ID.species.txt --map annotation_training.tsv
                         [--augustus-config $AUGUSTUS_CONFIG_PATH] [--funannotate-db $FUNANNOTATE_DB]
                         [--default-species anidulans] [--default-busco dikarya] [--json out.json]

Prints two shell assignments (AUGUSTUS_SPECIES=…, BUSCO_DB=…). A mapped Augustus species is used
only when its parameter directory exists under the Augustus config; a mapped BUSCO set only when
it is staged under the funannotate database; otherwise the defaults, with the reason recorded."""
from __future__ import annotations
import argparse, json, os


def read_map(path):
    rows = {}
    if not path or not os.path.exists(path):
        return rows
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[0] != "taxon":
                rows[f[0]] = (f[1], f[2])
    return rows


def read_species(path):
    try:
        first = open(path).readline().strip()
    except Exception:  # noqa: BLE001
        return "unknown"
    return first.split("\t")[0].strip() or "unknown"


def choose(species, rows, augustus_config=None, funannotate_db=None, default_species="anidulans", default_busco="dikarya"):
    genus = species.split()[0] if species and species != "unknown" else ""
    basis = "default"
    aug, bus = default_species, default_busco
    if species in rows:
        aug, bus = rows[species]; basis = f"species:{species}"
    elif genus in rows:
        aug, bus = rows[genus]; basis = f"genus:{genus}"
    notes = []
    if augustus_config and not os.path.isdir(os.path.join(augustus_config, "species", aug)):
        notes.append(f"augustus species {aug} not in {augustus_config}/species; using {default_species}")
        aug = default_species
    if funannotate_db and not os.path.isdir(os.path.join(funannotate_db, bus)):
        notes.append(f"BUSCO set {bus} not staged under {funannotate_db}; using {default_busco}")
        bus = default_busco
    return {"species": species, "augustus_species": aug, "busco_db": bus, "basis": basis, "notes": notes}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True); ap.add_argument("--map", required=True)
    ap.add_argument("--augustus-config", default=None); ap.add_argument("--funannotate-db", default=None)
    ap.add_argument("--default-species", default="anidulans"); ap.add_argument("--default-busco", default="dikarya")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    r = choose(read_species(a.species), read_map(a.map), a.augustus_config, a.funannotate_db, a.default_species, a.default_busco)
    if a.json:
        json.dump(r, open(a.json, "w"), indent=2)
    print(f"AUGUSTUS_SPECIES={r['augustus_species']}\nBUSCO_DB={r['busco_db']}")


if __name__ == "__main__":
    main()
