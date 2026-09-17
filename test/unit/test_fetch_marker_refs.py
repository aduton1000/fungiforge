"""Unit tests for bin/fetch_marker_refs.py and bin/fetch_mlst_schemes.py — the NCBI/PubMLST
download helpers, with the HTTP layer replaced by canned responses."""
import json, os
import fetch_marker_refs as fm  # noqa: E402
import fetch_mlst_schemes as fs  # noqa: E402


def test_fetch_locus_batches_and_counts(tmp_path, monkeypatch):
    calls = []

    def fake_get(url, params, **kw):
        calls.append((url.rsplit("/", 1)[-1], dict(params)))
        if url.endswith("esearch.fcgi"):
            return json.dumps({"esearchresult": {"count": "7", "webenv": "WE1", "querykey": "1"}})
        start = int(params["retstart"]); n = min(3, 7 - start)
        return "".join(f">ACC{start + i}.1 Aspergillus flavus strain T calmodulin gene\nACGT\n" for i in range(n))
    monkeypatch.setattr(fm, "_get", fake_get)
    monkeypatch.setattr(fm.time, "sleep", lambda *_: None)
    out = tmp_path / "CaM.fasta"
    count, n = fm.fetch_locus("CaM", fm.QUERIES["CaM"], str(out), batch=3)
    assert (count, n) == (7, 7) and open(out).read().count(">") == 7 and not (tmp_path / "CaM.fasta.part").exists()
    assert calls[0][0] == "esearch.fcgi" and calls[0][1]["usehistory"] == "y"
    assert [c[1]["retstart"] for c in calls[1:]] == [0, 3, 6] and all(c[1]["WebEnv"] == "WE1" for c in calls[1:])


def test_get_adds_api_key_and_tool(monkeypatch):
    seen = {}

    class R:
        def __init__(self, url): seen["url"] = url
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self): return b"ok"
    monkeypatch.setattr(fm.urllib.request, "urlopen", lambda url, timeout=0: R(url))
    monkeypatch.setenv("NCBI_API_KEY", "K123"); monkeypatch.setenv("NCBI_EMAIL", "x@y.z")
    assert fm._get("https://h/x.fcgi", {"db": "nuccore"}) == "ok"
    assert "api_key=K123" in seen["url"] and "tool=fungiforge" in seen["url"] and "email=x%40y.z" in seen["url"]


def test_queries_are_type_material_only():
    for locus, q in fm.QUERIES.items():
        assert '"sequence from type"[Filter]' in q and "fungi[Organism]" in q and "[Sequence Length]" in q and "NOT unverified[Title]" in q, locus
    assert set(fm.QUERIES) == {"CaM", "BenA", "TEF1", "RPB2", "LSU"}


def test_fetch_scheme_writes_mlst_layout(tmp_path, monkeypatch):
    base = "https://rest.pubmlst.org/db/pubmlst_afumigatus_seqdef"

    def fake_get(url, **kw):
        if url == f"{base}/schemes/1":
            return json.dumps({"description": "MLST", "loci": [f"{base}/loci/ANX4", f"{base}/loci/MAT1_2"], "profiles_csv": f"{base}/schemes/1/profiles_csv"})
        if url.endswith("/loci/ANX4/alleles_fasta"):
            return ">ANX4_1\nACGT\n>ANX4_2\nACGA"
        if url.endswith("/loci/MAT1_2/alleles_fasta"):
            return ">MAT1_2_1\nTTTT\n"
        if url.endswith("/profiles_csv"):
            return "ST\tANX4\tMAT1_2\n1\t1\t1\n2\t2\t1\n"
        raise AssertionError(url)
    monkeypatch.setattr(fs, "_get", fake_get)
    info = fs.fetch_scheme("afumigatus", str(tmp_path), log=lambda *_: None)
    d = tmp_path / "pubmlst" / "afumigatus"
    assert info == {"loci": ["ANX4", "MAT1_2"], "profiles": 2, "description": "MLST"}
    assert open(d / "ANX4.tfa").read() == ">ANX4_1\nACGT\n>ANX4_2\nACGA\n" and (d / "MAT1_2.tfa").exists()
    assert open(d / "afumigatus.txt").read().splitlines()[0] == "ST\tANX4\tMAT1_2"
    assert fs.DEFAULT.split(",") == ["afumigatus", "calbicans", "cglabrata", "ctropicalis", "ckrusei"]
