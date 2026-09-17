#!/usr/bin/env python3
"""Stage 13 — eukaryote biology extras (W2.6). Runs what is available and records what is not:

  secretome   SignalP 6 (site-built image; `signalp6` on PATH) -> signal-peptide proteins;
              EffectorP 3 (--effectorp-dir, Java + bundled WEKA) on the secretome -> effectors
  CAZymes     hmmsearch against the dbCAN family HMMs (--dbcan-hmm), run_dbcan's HMMER filter
              (E <= 1e-15, HMM coverage >= 0.35) -> families per protein
  virulence   diamond blastp against PHI-base (--phibase), best hit per protein with the PHI-base
              phenotype (from the header) -> counts by phenotype
  mating type MATalpha_HMGbox (PF04769) = MAT1-1-1 / MTLalpha1; an HMG_box (PF00505) protein in
              the MAT locus context (next to APN2 Exo_endo_phos PF03372 / SLA2 I_LWEQ PF01608 in
              the GenBank, or the only HMG-box protein without other domains) = MAT1-2-1
  ploidy      nQuire on the read BAM (--bam): lrdmodel free/di/tri/tetraploid log-likelihoods
              -> best model; GenomeScope's k-mer hint (stage 01b) stays the second opinion

Every tool is optional: missing tool or database -> that block is `null` with a reason, the
stage records `partial`, never fails. Tests drive the parsers on canned outputs.
"""
from __future__ import annotations
import argparse, collections, json, os, re, shutil, subprocess, sys, tempfile

PFAM = {"MATalpha_HMGbox": "PF04769", "HMG_box": "PF00505", "Exo_endo_phos": "PF03372", "I_LWEQ": "PF01608", "Homeodomain": "PF00046"}
CAZY_CLASS = re.compile(r"^(GH|GT|PL|CE|AA|CBM)\d+")


# ---------- generic ----------
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


def run(cmd, log=None, **kw):
    p = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if log:
        with open(log, "a") as fh:
            fh.write(f"$ {' '.join(map(str, cmd))}\n{p.stderr}\n")
    return p


def parse_domtbl(path):
    """hmmsearch --domtblout -> [(target, query_name, query_acc, evalue, tlen, qlen, hmm_from, hmm_to, ali_from, ali_to)]."""
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 22:
                continue
            try:
                rows.append({"target": f[0], "tlen": int(f[2]), "query": f[3], "acc": f[4], "qlen": int(f[5]), "evalue": float(f[12]),
                             "hmm_from": int(f[15]), "hmm_to": int(f[16]), "ali_from": int(f[17]), "ali_to": int(f[18])})
            except ValueError:
                continue
    return rows


# ---------- CAZymes ----------
def cazyme_families(domtbl_rows, max_evalue=1e-15, min_cov=0.35):
    """run_dbcan's HMMER filter: per protein, the families with E <= 1e-15 and HMM coverage >= 0.35."""
    per = collections.defaultdict(set)
    for r in domtbl_rows:
        cov = (r["hmm_to"] - r["hmm_from"] + 1) / r["qlen"] if r["qlen"] else 0
        if r["evalue"] <= max_evalue and cov >= min_cov:
            fam = r["query"].replace(".hmm", "")
            per[r["target"]].add(fam)
    return {k: sorted(v) for k, v in per.items()}


def cazyme_summary(fams):
    by_class = collections.Counter()
    for f in fams.values():
        for fam in f:
            m = CAZY_CLASS.match(fam)
            by_class[m.group(1) if m else "other"] += 1
    return {"n_cazymes": len(fams), "by_class": dict(sorted(by_class.items())),
            "families": dict(sorted(collections.Counter(x for f in fams.values() for x in f).items()))}


# ---------- secretome ----------
def parse_signalp6(path):
    """prediction_results.txt of signalp6: ID<tab>Prediction<tab>OTHER<tab>SP(Sec/SPI)...<tab>CS Position."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            pred = f[1].strip()
            cs = None
            m = re.search(r"CS pos: (\d+)-(\d+)", line)
            if m:
                cs = int(m.group(1))
            out[f[0].split()[0]] = {"prediction": pred, "signal_peptide": pred.upper().startswith("SP"), "cs": cs}
    return out


def parse_effectorp(path):
    """EffectorP 3 -o table: identifier, cytoplasmic, apoplastic, non-effector, prediction (tab or spaces)."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip() or line.startswith("-"):
                continue
            f = [x.strip() for x in re.split(r"\t|\s{2,}", line.rstrip("\n")) if x.strip()]
            if len(f) < 5:
                continue
            pred = f[-1]
            out[f[0].split()[0]] = {"prediction": pred, "effector": "effector" in pred.lower() and "non-effector" not in pred.lower()}
    return out


# ---------- virulence ----------
def parse_phibase_hits(path, min_pident=40.0, min_cov=0.5):
    """diamond outfmt 6 (qseqid sseqid pident length qlen slen evalue bitscore stitle) vs PHI-base;
    best hit per query; phenotype parsed from the PHI-base header (…#phenotype)."""
    best = {}
    if not path or not os.path.exists(path):
        return best
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            try:
                pid, alen, qlen, slen, bits = float(f[2]), int(f[3]), int(f[4]), int(f[5]), float(f[7])
            except ValueError:
                continue
            if pid < min_pident or alen < min_cov * min(qlen, slen):
                continue
            if f[0] not in best or bits > best[f[0]]["bitscore"]:
                title = f[8]
                parts = title.split("#")
                pheno = parts[-1].strip() if len(parts) > 1 else "unknown"
                # PHI-base header: accession#PHI:id#gene#taxid#organism#phenotype
                best[f[0]] = {"hit": f[1], "pident": pid, "bitscore": bits, "phenotype": pheno.lower().replace(" ", "_"),
                              "gene": parts[2].strip() if len(parts) > 3 else None, "organism": parts[4].strip() if len(parts) > 5 else None}
    return best


def virulence_summary(hits):
    return {"n_phibase_hits": len(hits), "by_phenotype": dict(sorted(collections.Counter(h["phenotype"] for h in hits.values()).items()))}


# ---------- mating type ----------
def gbk_gene_order(gbk_path):
    """{contig: [(protein_id_or_locus, start)]} in genomic order, plus {name: contig}."""
    order, where = {}, {}
    try:
        from Bio import SeqIO
    except Exception:  # noqa: BLE001
        return order, where
    for rec in SeqIO.parse(gbk_path, "genbank"):
        genes = []
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            q = feat.qualifiers
            names = [x.replace("ncbi:", "") for x in q.get("protein_id", [])] + q.get("locus_tag", [])
            if not names:
                continue
            genes.append((names, int(feat.location.start)))
        genes.sort(key=lambda g: g[1])
        order[rec.id] = genes
        for i, (names, _) in enumerate(genes):
            for n in names:
                where[n] = (rec.id, i)
    return order, where


def _name_variants(name):
    base = re.sub(r"-T\d+$", "", name)
    return {name, base, f"{base}-T1"}


def mating_type(domtbl_rows, gbk_path=None, all_domains=None, window=6):
    """MAT1-1 (alpha box) / MAT1-2 (HMG box in the MAT context) / both / undetermined."""
    by_acc = collections.defaultdict(set)
    for r in domtbl_rows:
        if r["evalue"] <= 1e-5:
            by_acc[r["acc"].split(".")[0]].add(r["target"])
    alpha = sorted(by_acc.get(PFAM["MATalpha_HMGbox"], set()))
    hmg = sorted(by_acc.get(PFAM["HMG_box"], set()))
    apn2 = by_acc.get(PFAM["Exo_endo_phos"], set())
    sla2 = by_acc.get(PFAM["I_LWEQ"], set())
    hmg_in_context, context = [], {}
    if hmg and gbk_path and os.path.exists(gbk_path):
        order, where = gbk_gene_order(gbk_path)
        flank = {n for n in apn2 | sla2}
        for h in hmg:
            loc = next((where[v] for v in _name_variants(h) if v in where), None)
            if not loc:
                continue
            contig, i = loc
            genes = order[contig]
            neighbours = {n for names, _ in genes[max(0, i - window):i + window + 1] for n in names}
            near = [f for f in flank if _name_variants(f) & neighbours]
            if near:
                hmg_in_context.append(h); context[h] = sorted(near)
    # fallback without synteny: an HMG-box protein that carries no other Pfam domain and is short
    if not hmg_in_context and hmg and all_domains:
        for h in hmg:
            doms = all_domains.get(h, set())
            if doms <= {PFAM["HMG_box"]}:
                hmg_in_context.append(h); context[h] = ["single-domain HMG-box (no synteny evidence)"]
    if alpha and hmg_in_context:
        call, note = "MAT1-1/MAT1-2", "both idiomorphs present (homothallic, or a heterozygous diploid MTLa/alpha)"
    elif alpha:
        call, note = "MAT1-1", "alpha-box protein (PF04769) present"
    elif hmg_in_context:
        call, note = "MAT1-2", "HMG-box protein in the MAT locus context (" + "; ".join(f"{h}: {', '.join(context[h])}" for h in hmg_in_context) + ")"
    else:
        call, note = "undetermined", "no alpha-box protein and no HMG-box protein in a MAT context"
    return {"mating_type": call, "alpha_box_proteins": alpha, "mat1_2_candidates": hmg_in_context, "n_hmg_box_proteins": len(hmg), "note": note}


# ---------- ploidy ----------
def parse_nquire_lrdmodel(text):
    """nQuire lrdmodel output: header + one line per BAM with free / dip / tri / tet log-likelihoods
    (and their deltas). Best model = the fixed-ploidy model with the highest log-likelihood."""
    lines = [l for l in text.splitlines() if l.strip() and not l.startswith("file")]
    if not lines:
        return None
    hdr = [h.strip() for h in [l for l in text.splitlines() if l.startswith("file")][0].split("\t")] if any(l.startswith("file") for l in text.splitlines()) else None
    f = lines[-1].split("\t")
    try:
        if hdr and "dip" in hdr:
            vals = {k: float(v) for k, v in zip(hdr[1:], f[1:]) if k in ("free", "dip", "tri", "tet", "d_dip", "d_tri", "d_tet")}
        else:
            vals = {"free": float(f[1]), "dip": float(f[2]), "tri": float(f[3]), "tet": float(f[4])}
    except (ValueError, IndexError):
        return None
    models = {k: vals[k] for k in ("dip", "tri", "tet") if k in vals}
    if not models:
        return None
    best = max(models, key=models.get)
    ploidy = {"dip": 2, "tri": 3, "tet": 4}[best]
    return {"loglik": vals, "best_model": best, "ploidy_call": ploidy,
            "delta_to_free": round(vals.get("free", 0) - models[best], 2) if "free" in vals else None}


def ploidy_from_nquire(bam, prefix, min_sites=1000, log=None):
    """nQuire create (denoised) -> lrdmodel; a haploid genome gives few heterozygous sites, which is
    reported as such (nQuire cannot model ploidy 1: no biallelic sites)."""
    if not bam or not os.path.exists(bam) or not shutil.which("nQuire"):
        return None
    p = run(["nQuire", "create", "-b", bam, "-o", prefix, "-q", "20", "-c", "10"], log)
    if p.returncode != 0 or not os.path.exists(prefix + ".bin"):
        return {"error": "nQuire create failed"}
    d = run(["nQuire", "denoise", prefix + ".bin", "-o", prefix + ".denoised"], log)
    binfile = prefix + ".denoised.bin" if d.returncode == 0 and os.path.exists(prefix + ".denoised.bin") else prefix + ".bin"
    v = run(["nQuire", "view", binfile], log)
    n_sites = sum(1 for l in v.stdout.splitlines() if l.strip()) if v.returncode == 0 else None
    res = {"n_sites": n_sites, "denoised": binfile.endswith("denoised.bin")}
    if n_sites is not None and n_sites < min_sites:
        res.update(ploidy_call=1, best_model="haploid_like", note=f"only {n_sites} biallelic sites after denoising: no heterozygosity signal (haploid or homozygous)")
        return res
    m = run(["nQuire", "lrdmodel", binfile], log)
    parsed = parse_nquire_lrdmodel(m.stdout) if m.returncode == 0 else None
    if not parsed:
        res["error"] = "nQuire lrdmodel failed"
        return res
    res.update(parsed)
    return res


# ---------- driver ----------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True); ap.add_argument("--proteins", required=True); ap.add_argument("--gbk", default=None)
    ap.add_argument("--bam", default=None); ap.add_argument("--data-dir", default="")
    ap.add_argument("--dbcan-hmm", default=None); ap.add_argument("--phibase", default=None); ap.add_argument("--effectorp-dir", default=None)
    ap.add_argument("--pfam-hmm", default=None, help="Pfam-A.hmm (default <data_dir>/funannotate/Pfam-A.hmm)")
    ap.add_argument("--threads", type=int, default=4); ap.add_argument("--workdir", default="extras_work")
    ap.add_argument("--signalp-results", default=None, help="pre-computed signalp6 prediction_results.txt (tests)")
    ap.add_argument("--effectorp-results", default=None); ap.add_argument("--dbcan-domtbl", default=None); ap.add_argument("--pfam-domtbl", default=None)
    ap.add_argument("--phibase-hits", default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    log = os.path.join(a.workdir, "extras.log")
    prots = read_fasta(a.proteins)
    doc = {"sample": a.sample, "stage": "extras", "n_proteins": len(prots), "tools": {}}
    dd = a.data_dir or ""
    dbcan = a.dbcan_hmm or (os.path.join(dd, "dbcan", "dbCAN-HMMdb.txt") if dd else None)
    phibase = a.phibase or (os.path.join(dd, "phibase", "phi-base_current.fas") if dd else None)
    effdir = a.effectorp_dir or (os.path.join(dd, "effectorp", "EffectorP-3.0") if dd else None)
    pfam = a.pfam_hmm or (os.path.join(dd, "funannotate", "Pfam-A.hmm") if dd else None)
    clean = os.path.join(a.workdir, "proteins.faa")
    with open(clean, "w") as fh:
        for n, s in prots.items():
            fh.write(f">{n}\n{s.rstrip('*')}\n")

    # CAZymes
    domtbl = a.dbcan_domtbl
    if not domtbl and dbcan and os.path.exists(dbcan) and shutil.which("hmmsearch") and prots:
        domtbl = os.path.join(a.workdir, "dbcan.domtbl")
        p = run(["hmmsearch", "--cpu", str(a.threads), "--domtblout", domtbl, "-E", "1e-5", "--noali", dbcan, clean], log)
        if p.returncode != 0:
            domtbl = None; doc["tools"]["dbcan"] = "hmmsearch failed"
    if domtbl:
        fams = cazyme_families(parse_domtbl(domtbl))
        doc["cazymes"] = cazyme_summary(fams); doc["cazymes"]["per_protein"] = fams; doc["tools"]["dbcan"] = "ok"
    else:
        doc["cazymes"] = None; doc["tools"].setdefault("dbcan", "dbCAN HMMs not staged (fetch_references.sh dbcan)")

    # secretome + effectors
    sp = None
    if a.signalp_results:
        sp = parse_signalp6(a.signalp_results)
    elif shutil.which("signalp6") and prots:
        outd = os.path.join(a.workdir, "signalp6")
        p = run(["signalp6", "--fastafile", clean, "--organism", "eukarya", "--output_dir", outd, "--format", "txt", "--mode", "fast"], log)
        res = os.path.join(outd, "prediction_results.txt")
        sp = parse_signalp6(res) if p.returncode == 0 and os.path.exists(res) else None
        if sp is None:
            doc["tools"]["signalp6"] = "signalp6 failed"
    if sp is not None:
        secreted = sorted(n for n, v in sp.items() if v["signal_peptide"])
        doc["secretome"] = {"n_signal_peptide": len(secreted), "pct_signal_peptide": round(100.0 * len(secreted) / len(prots), 1) if prots else None, "proteins": secreted}
        doc["tools"]["signalp6"] = "ok"
        eff = None
        if a.effectorp_results:
            eff = parse_effectorp(a.effectorp_results)
        elif secreted and effdir and os.path.exists(os.path.join(effdir, "EffectorP.py")) and shutil.which("java"):
            sec_fa = os.path.join(a.workdir, "secreted.faa")
            with open(sec_fa, "w") as fh:
                for n in secreted:
                    fh.write(f">{n}\n{prots[n].rstrip('*')}\n")
            res = os.path.join(a.workdir, "effectorp.txt")
            p = run([sys.executable, os.path.join(effdir, "EffectorP.py"), "-i", os.path.abspath(sec_fa), "-o", os.path.abspath(res)], log, cwd=effdir)
            eff = parse_effectorp(res) if p.returncode == 0 and os.path.exists(res) else None
            if eff is None:
                doc["tools"]["effectorp"] = "EffectorP failed"
        if eff is not None:
            effectors = sorted(n for n, v in eff.items() if v["effector"])
            doc["effectors"] = {"n_effectors": len(effectors), "by_prediction": dict(sorted(collections.Counter(v["prediction"] for v in eff.values()).items())), "proteins": effectors}
            doc["tools"]["effectorp"] = "ok"
        else:
            doc["effectors"] = None; doc["tools"].setdefault("effectorp", "EffectorP 3 not staged (fetch_references.sh effectorp) or no secreted proteins")
    else:
        doc["secretome"] = None; doc["effectors"] = None
        doc["tools"].setdefault("signalp6", "signalp6 not in this image (site-built SignalP 6 image: bin/hpc_install.sh --signalp)")

    # virulence (PHI-base)
    hits_path = a.phibase_hits
    if not hits_path and phibase and os.path.exists(phibase) and shutil.which("diamond") and prots:
        dbp = os.path.join(a.workdir, "phibase")
        p1 = run(["diamond", "makedb", "--in", phibase, "-d", dbp, "--quiet"], log)
        hits_path = os.path.join(a.workdir, "phibase.tsv")
        p2 = run(["diamond", "blastp", "-q", clean, "-d", dbp, "-o", hits_path, "--quiet", "-p", str(a.threads), "-e", "1e-10", "-k", "5",
                  "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore", "stitle"], log)
        if p1.returncode != 0 or p2.returncode != 0:
            hits_path = None; doc["tools"]["phibase"] = "diamond failed"
    if hits_path:
        hits = parse_phibase_hits(hits_path)
        doc["virulence"] = virulence_summary(hits); doc["virulence"]["per_protein"] = hits; doc["tools"]["phibase"] = "ok"
    else:
        doc["virulence"] = None; doc["tools"].setdefault("phibase", "PHI-base not staged (fetch_references.sh phibase)")

    # mating type (Pfam profiles fetched from Pfam-A.hmm of the funannotate database)
    pf_domtbl = a.pfam_domtbl
    if not pf_domtbl and pfam and os.path.exists(pfam) and shutil.which("hmmfetch") and prots:
        keys = os.path.join(a.workdir, "mat.keys"); sub = os.path.join(a.workdir, "mat.hmm")
        open(keys, "w").write("\n".join(PFAM.values()) + "\n")
        names_ok = run(["hmmfetch", "-f", "-o", sub, pfam, keys], log).returncode == 0 and os.path.getsize(sub) > 0
        if not names_ok:   # Pfam-A.hmm keys by NAME in some builds
            open(keys, "w").write("\n".join(PFAM.keys()) + "\n")
            names_ok = run(["hmmfetch", "-f", "-o", sub, pfam, keys], log).returncode == 0 and os.path.getsize(sub) > 0
        if names_ok:
            pf_domtbl = os.path.join(a.workdir, "mat.domtbl")
            if run(["hmmsearch", "--cpu", str(a.threads), "--domtblout", pf_domtbl, "-E", "1e-5", "--noali", sub, clean], log).returncode != 0:
                pf_domtbl = None
    if pf_domtbl:
        rows = parse_domtbl(pf_domtbl)
        # accession may be '-' when profiles were fetched by name: map names back to accessions
        for r in rows:
            if r["acc"] in ("-", "") and r["query"] in PFAM:
                r["acc"] = PFAM[r["query"]]
        all_dom = collections.defaultdict(set)
        for r in rows:
            if r["evalue"] <= 1e-5:
                all_dom[r["target"]].add(r["acc"].split(".")[0])
        doc["mating"] = mating_type(rows, a.gbk, all_dom); doc["tools"]["mating_hmm"] = "ok"
    else:
        doc["mating"] = {"mating_type": "undetermined", "note": "Pfam-A.hmm not available (funannotate database) or hmmer missing"}
        doc["tools"].setdefault("mating_hmm", "Pfam profiles unavailable")

    # ploidy
    if a.bam and os.path.exists(a.bam) and os.path.getsize(a.bam) > 0:
        doc["ploidy_reads"] = ploidy_from_nquire(a.bam, os.path.join(a.workdir, "nquire"), log=log)
        doc["tools"]["nquire"] = "ok" if doc["ploidy_reads"] and "error" not in doc["ploidy_reads"] else (doc["ploidy_reads"] or {}).get("error", "nQuire not in this image")
    else:
        doc["ploidy_reads"] = None; doc["tools"]["nquire"] = "no read BAM (stage 09 read genotyping off or failed)"

    # flat summary for the master table
    doc["mating_type"] = doc["mating"]["mating_type"]
    doc["ploidy"] = (doc["ploidy_reads"] or {}).get("ploidy_call", "NA") if doc["ploidy_reads"] else "NA"
    doc["n_cazymes"] = doc["cazymes"]["n_cazymes"] if doc["cazymes"] else None
    doc["n_secreted"] = doc["secretome"]["n_signal_peptide"] if doc["secretome"] else None
    doc["n_effectors"] = doc["effectors"]["n_effectors"] if doc["effectors"] else None
    doc["n_phibase_hits"] = doc["virulence"]["n_phibase_hits"] if doc["virulence"] else None
    doc["skipped"] = sorted(k for k, v in doc["tools"].items() if v != "ok")
    json.dump(doc, open(a.out, "w"), indent=2)
    print(f"[extras] {a.sample}: mating={doc['mating_type']} ploidy={doc['ploidy']} cazymes={doc['n_cazymes']} "
          f"secreted={doc['n_secreted']} effectors={doc['n_effectors']} phibase={doc['n_phibase_hits']}; skipped: {', '.join(doc['skipped']) or 'none'}")


if __name__ == "__main__":
    main()
