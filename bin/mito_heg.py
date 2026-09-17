#!/usr/bin/env python3
"""Mitochondrial homing-endonuclease / maturase ORFs (W2.8, stage 10).

  mito_heg.py --sample ID --mito ID.mito.fasta [--pfam-hmm Pfam-A.hmm | --domtbl heg.domtbl] --out heg.json
              [--min-orf 100] [--threads 4] [--workdir heg_work]

Six-frame ORFs (>= --min-orf aa, genetic code 4: UGA = Trp) of the mitochondrial contigs are
searched with the Pfam profiles LAGLIDADG_1/2/3 (PF00961, PF03161, PF14528), GIY-YIG (PF01541) and
HNH (PF01844) fetched from Pfam-A.hmm; ORFs with a hit (E <= 1e-5) are the mobile intron-encoded
elements the mitogenome carries (the plan's "mito mobile introns / homing endonucleases")."""
from __future__ import annotations
import argparse, collections, json, os, re, shutil, subprocess

PFAM = {"LAGLIDADG_1": "PF00961", "LAGLIDADG_2": "PF03161", "LAGLIDADG_3": "PF14528", "GIY-YIG": "PF01541", "HNH": "PF01844"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_B = "TCAG"; _AA = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"   # code 4 (mold mitochondrial): TGA = W
CODON = {_B[i // 16] + _B[(i // 4) % 4] + _B[i % 4]: a for i, a in enumerate(_AA)}


def read_fasta(path):
    seqs, name, buf = {}, None, []
    if not path or not os.path.exists(path):
        return seqs
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def translate(s):
    return "".join(CODON.get(s[i:i + 3].upper(), "X") for i in range(0, len(s) - 2, 3))


def six_frame_orfs(seqs, min_len=100):
    """ORFs between stops (Met not required: mitochondrial intron ORFs are often in-frame with the exon)."""
    orfs = []
    for name, seq in seqs.items():
        for strand, s in (("+", seq), ("-", seq.translate(COMP)[::-1])):
            for frame in range(3):
                prot = translate(s[frame:])
                pos = 0
                for seg in prot.split("*"):
                    if len(seg) >= min_len:
                        nt_start = frame + pos * 3
                        nt_end = nt_start + len(seg) * 3
                        if strand == "+":
                            gs, ge = nt_start + 1, nt_end
                        else:
                            gs, ge = len(seq) - nt_end + 1, len(seq) - nt_start
                        orfs.append({"id": f"{name}_{strand}{frame}_{gs}_{ge}", "contig": name, "start": gs, "end": ge, "strand": strand, "aa": len(seg), "seq": seg})
                    pos += len(seg) + 1
    return orfs


def parse_domtbl(path):
    hits = collections.defaultdict(dict)
    if not path or not os.path.exists(path):
        return hits
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 13:
                continue
            try:
                ev = float(f[12])
            except ValueError:
                continue
            if ev <= 1e-5:
                fam = f[3]
                if fam not in hits[f[0]] or ev < hits[f[0]][fam]:
                    hits[f[0]][fam] = ev
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True); ap.add_argument("--mito", required=True); ap.add_argument("--pfam-hmm", default=None)
    ap.add_argument("--domtbl", default=None); ap.add_argument("--min-orf", type=int, default=100); ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--workdir", default="heg_work"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    seqs = read_fasta(a.mito)
    orfs = six_frame_orfs(seqs, a.min_orf)
    doc = {"sample": a.sample, "n_mito_contigs": len(seqs), "n_orfs": len(orfs), "min_orf_aa": a.min_orf}
    domtbl = a.domtbl
    if not domtbl and orfs and a.pfam_hmm and os.path.exists(a.pfam_hmm) and shutil.which("hmmfetch") and shutil.which("hmmsearch"):
        faa = os.path.join(a.workdir, "orfs.faa")
        with open(faa, "w") as fh:
            for o in orfs:
                fh.write(f">{o['id']}\n{o['seq']}\n")
        keys, sub = os.path.join(a.workdir, "heg.keys"), os.path.join(a.workdir, "heg.hmm")
        ok = False
        for names in (list(PFAM.values()), list(PFAM.keys())):
            open(keys, "w").write("\n".join(names) + "\n")
            p = subprocess.run(["hmmfetch", "-f", "-o", sub, a.pfam_hmm, keys], capture_output=True, text=True)
            if p.returncode == 0 and os.path.getsize(sub) > 0:
                ok = True; break
        if ok:
            domtbl = os.path.join(a.workdir, "heg.domtbl")
            p = subprocess.run(["hmmsearch", "--cpu", str(a.threads), "--domtblout", domtbl, "-E", "1e-5", "--noali", sub, faa], capture_output=True, text=True)
            if p.returncode != 0:
                domtbl = None; doc["note"] = "hmmsearch failed"
        else:
            doc["note"] = "Pfam profiles not found in Pfam-A.hmm"
    if domtbl:
        hits = parse_domtbl(domtbl)
        heg = [dict({k: v for k, v in o.items() if k != "seq"}, families=sorted(hits[o["id"]])) for o in orfs if o["id"] in hits]
        doc["heg_orfs"] = heg; doc["n_heg_orfs"] = len(heg)
        doc["by_family"] = dict(collections.Counter(f for h in heg for f in h["families"]))
        doc["status"] = "ok"
    else:
        doc["heg_orfs"] = None; doc["n_heg_orfs"] = None
        doc.setdefault("note", "no Pfam-A.hmm (funannotate database) or hmmer: homing-endonuclease ORFs not scanned")
        doc["status"] = "skipped"
    json.dump(doc, open(a.out, "w"), indent=2)
    print(f"[mito_heg] {a.sample}: {len(orfs)} ORFs >= {a.min_orf} aa; HEG ORFs: {doc['n_heg_orfs']}")


if __name__ == "__main__":
    main()
