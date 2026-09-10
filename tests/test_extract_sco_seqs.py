"""
Stage 2 test: 02_extract_sco_seqs.py against a synthetic 3-species fixture.

Fixture (tests/fixtures/stage2/) exercises every path:
    gene1  all IDs exact                       -> written, no flag
    gene2  fx3 protein ID ends ".p", CDS lacks it -> strip:.p match, written
    gene3  fx2 has no CDS record for it        -> whole gene dropped
    gene4  fx1's CDS is the wrong (short) one  -> written but length-flagged
Expected: 3 genes written, 1 dropped, 1 length flag, fx3 uses strip:.p once.
"""
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIXTURE = ROOT / "tests" / "fixtures" / "stage2"

spec = importlib.util.spec_from_file_location(
    "extract_sco_seqs", ROOT / "codes" / "core" / "02_extract_sco_seqs.py")
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


def test_extract(tmp_path):
    rc = mod.main([
        str(FIXTURE / "common_scos.tsv"),
        "--prt-dir", str(FIXTURE / "prt"),
        "--cds-dir", str(FIXTURE / "cds"),
        "-o", str(tmp_path),
    ])
    assert rc == 0

    prt = sorted(p.name for p in (tmp_path / "prt").glob("*.faa"))
    cds = sorted(p.name for p in (tmp_path / "cds").glob("*.fna"))
    assert prt == ["gene1.faa", "gene2.faa", "gene4.faa"], prt
    assert cds == ["gene1.fna", "gene2.fna", "gene4.fna"], cds

    # each per-gene file has one record per species, header ">{sp} {busco_id}"
    g1 = (tmp_path / "prt" / "gene1.faa").read_text().splitlines()
    assert [ln for ln in g1 if ln.startswith(">")] == [
        ">fx1 gene1", ">fx2 gene1", ">fx3 gene1"]

    dropped = (tmp_path / "stage2_dropped.tsv").read_text().splitlines()
    assert dropped[0] == "busco_id\treason"
    assert len(dropped) == 2
    assert dropped[1].startswith("gene3\tfx2:no_cds")

    man = (tmp_path / "stage2_manifest.tsv").read_text().splitlines()
    header, *rows = man
    cols = header.split("\t")
    fi, ri = cols.index("length_flag"), cols.index("match_rule")
    flags = [r.split("\t")[fi] for r in rows]
    assert flags.count("1") == 1, flags          # only gene4/fx1
    rules = [r.split("\t")[ri] for r in rows if r.split("\t")[1] == "fx3"]
    assert rules.count("strip:.p") == 1, rules   # gene2/fx3
    assert rules.count("exact") == 2, rules      # gene1, gene4

    print("stage 2 OK: 3 written, 1 dropped, 1 length flag, fx3 strip:.p x1")


if __name__ == "__main__":
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        test_extract(Path(td))
    print("ALL TESTS PASSED")
