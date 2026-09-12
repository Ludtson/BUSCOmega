"""
Stage 1 test: 01_common_scos.py against a synthetic 3-species fixture.

Fixture (tests/fixtures/stage1_busco/) exercises every branch:
    geneA  Complete    Complete    Complete      -> common set
    geneB  Complete    Complete    Complete      -> common set
    geneC  Complete    Complete    Duplicated    -> sp3 sole-blocks (dup)
    geneD  Complete    Missing     Complete      -> sp2 sole-blocks (missing)
    geneE  Fragmented  Complete    Complete      -> sp1 sole-blocks (frag)
    geneF  Complete    Complete    Complete      -> common set
    geneG  Missing     Fragmented  Complete      -> shared loss, no sole blocker
Expected common set: {geneA, geneB, geneF}
Expected recoverable: sp1=1 (geneE), sp2=1 (geneD), sp3=1 (geneC)
"""
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIXTURE = ROOT / "tests" / "fixtures" / "stage1_busco"

spec = importlib.util.spec_from_file_location(
    "common_scos",
    ROOT / "codes" / "core" / "01_common_scos.py",
)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


def test_common_set_and_summary(tmp_path):
    rc = mod.main([str(FIXTURE), "-o", str(tmp_path)])
    assert rc == 0

    common = (tmp_path / "common_scos.tsv").read_text().splitlines()
    header, *rows = common
    assert header.split("\t") == ["busco_id", "sp1", "sp2", "sp3"]
    ids = [r.split("\t")[0] for r in rows]
    assert ids == ["geneA", "geneB", "geneF"], ids

    row_by_id = {r.split("\t")[0]: r.split("\t") for r in rows}
    assert row_by_id["geneA"] == ["geneA", "sp1_protA", "sp2_protA", "sp3_protA"]
    assert row_by_id["geneF"][3] == "sp3_protF"

    summary = (tmp_path / "busco_status_summary.tsv").read_text().splitlines()
    s_header, *s_rows = summary
    assert s_header.split("\t") == [
        "species", "complete", "duplicated", "fragmented", "missing", "total",
        "in_common_set", "sole_block_dup", "sole_block_frag",
        "sole_block_missing", "recoverable",
    ]
    srow = {r.split("\t")[0]: r.split("\t")[1:] for r in s_rows}
    #                     C  D  F  M  tot ic  sbD sbF sbM rec
    assert srow["sp1"] == ["5", "0", "1", "1", "7", "3", "0", "1", "0", "1"]
    assert srow["sp2"] == ["5", "0", "1", "1", "7", "3", "0", "0", "1", "1"]
    assert srow["sp3"] == ["6", "2", "0", "0", "8", "3", "1", "0", "0", "1"]

    print("stage 1 OK: common set =", ids,
          "| recoverable:", {s: srow[s][-1] for s in ("sp1", "sp2", "sp3")})


def test_exclude_drops_species_from_consideration(tmp_path):
    """--exclude sp3 removes it from both the species set and the
    intersection -- geneC, which only sp3 sole-blocked (Duplicated), should
    now join the common set since it's Complete in the remaining sp1/sp2."""
    rc = mod.main([str(FIXTURE), "-o", str(tmp_path), "--exclude", "sp3"])
    assert rc == 0
    common = (tmp_path / "common_scos.tsv").read_text().splitlines()
    header, *rows = common
    assert header.split("\t") == ["busco_id", "sp1", "sp2"]
    ids = [r.split("\t")[0] for r in rows]
    assert ids == ["geneA", "geneB", "geneC", "geneF"], ids


if __name__ == "__main__":
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        test_common_set_and_summary(Path(td))
    print("ALL TESTS PASSED")
