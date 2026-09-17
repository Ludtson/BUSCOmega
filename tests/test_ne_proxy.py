"""
Stage 7 test: 07_ne_proxy.py — the pooling arithmetic, the dS-ceiling and
ds_floor filtering, bootstrap determinism, and an end-to-end run on a tiny
synthetic records dir.
"""
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location(
    "ne_proxy", ROOT / "codes" / "core" / "07_ne_proxy.py")
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)

HEADER = ("gene_id\tdataset_index\ttip\tbranch\tt\tN\tS\tomega\tdN\tdS\t"
          "lnL\tkappa\tis_terminal\tqc_flag")


def _row(gene, tip, N, S, dN, dS, qc="ok"):
    om = dN / dS if dS else 999.0
    return (f"{gene}\t1\t{tip}\tb\t0.2\t{N}\t{S}\t{om}\t{dN}\t{dS}\t"
            f"-100\t2\tTrue\t{qc}")


def _write(path, rows):
    path.write_text(HEADER + "\n" + "\n".join(rows) + "\n", encoding="utf-8")


def test_pooled_omega():
    rows = [mod.Row(dict(zip(HEADER.split("\t"), r.split("\t"))))
            for r in [_row("g1", "spA", 300, 100, 0.01, 0.10),
                      _row("g2", "spA", 300, 100, 0.02, 0.10)]]
    # dN_pool = (300*0.01 + 300*0.02)/600 = 0.015 ; dS_pool = 20/200 = 0.10
    # omega = 0.015 / 0.10 = 0.15
    assert abs(mod.pooled_omega(rows) - 0.15) < 1e-9
    # unequal gene sizes: the big gene dominates the pool
    r2 = [mod.Row(dict(zip(HEADER.split("\t"), x.split("\t")))) for x in [
        _row("g1", "spA", 900, 300, 0.02, 0.10),      # big, omega 0.2
        _row("g2", "spA", 90, 30, 0.05, 0.10)]]       # small, omega 0.5
    # dN = (900*0.02 + 90*0.05)/990 = (18+4.5)/990 = 0.022727
    # dS = (300*0.1 + 30*0.1)/330 = 33/330 = 0.1
    assert abs(mod.pooled_omega(r2) - 0.227273) < 1e-5


def test_filter_focal():
    rows = [mod.Row(dict(zip(HEADER.split("\t"), r.split("\t")))) for r in [
        _row("g1", "spA", 300, 100, 0.01, 0.10),                 # keep
        _row("g2", "spA", 300, 100, 0.00, 0.10, qc="omega_boundary"),  # keep
        _row("g3", "spA", 300, 100, 0.02, 0.0005, qc="ds_floor"),      # drop
        _row("g4", "spA", 300, 100, 0.03, 2.5),                  # drop: dS>1.5
        _row("g5", "spB", 300, 100, 0.01, 0.10),                 # other species
    ]]
    used, annot, n_floor, n_ceil = mod.filter_focal(rows, "spA", 1.5)
    assert [r.gene_id for r in used] == ["g1", "g2"]
    assert n_floor == 1 and n_ceil == 1
    assert len(annot) == 4                       # spB excluded entirely


def test_bootstrap_deterministic():
    rows = [mod.Row(dict(zip(HEADER.split("\t"), r.split("\t")))) for r in [
        _row(f"g{i}", "spA", 300, 100, 0.01 + 0.001 * i, 0.10)
        for i in range(20)]]
    a = mod.bootstrap_ci(rows, 200, seed=42)
    b = mod.bootstrap_ci(rows, 200, seed=42)
    c = mod.bootstrap_ci(rows, 200, seed=7)
    assert a == b and a != c
    assert a[0] <= mod.pooled_omega(rows) <= a[1]


def test_end_to_end(tmp_path):
    rd = tmp_path / "06_parsed"
    rd.mkdir()
    _write(rd / "m0_records.tsv",
           [_row("g1", "spA", 300, 100, 0.005, 0.10),
            _row("g1", "spB", 300, 100, 0.005, 0.10),
            _row("g2", "spA", 300, 100, 0.006, 0.12),
            _row("g2", "spB", 300, 100, 0.006, 0.12)])
    _write(rd / "two_ratio.spA_records.tsv",
           [_row("g1", "spA", 300, 100, 0.004, 0.10),
            _row("g2", "spA", 300, 100, 0.006, 0.11),
            _row("g3", "spA", 300, 100, 0.02, 0.0004, qc="ds_floor")])
    _write(rd / "two_ratio.spB_records.tsv",
           [_row("g1", "spB", 300, 100, 0.010, 0.10),
            _row("g2", "spB", 300, 100, 0.012, 0.11)])

    out = tmp_path / "07"
    rc = mod.main(["--records-dir", str(rd), "-o", str(out),
                   "--bootstrap", "50"])
    assert rc == 0

    proxy = (out / "ne_proxy.tsv").read_text().splitlines()
    hdr = proxy[0].split("\t")
    recs = {r.split("\t")[0]: dict(zip(hdr, r.split("\t"))) for r in proxy[1:]}
    assert set(recs) == {"spA", "spB"}
    assert recs["spA"]["n_genes"] == "2"           # g3 dropped
    assert recs["spA"]["n_excl_ds_floor"] == "1"
    # spB more relaxed than spA -> higher pooled omega
    assert float(recs["spB"]["omega_pooled"]) > float(recs["spA"]["omega_pooled"])
    assert recs["spA"]["omega_M0"] == recs["spB"]["omega_M0"]  # one clade value

    pg = (out / "per_gene_omega.tsv").read_text().splitlines()
    assert any("\tg3\t" in ln and "\t0\tds_floor" in ln for ln in pg)
    for name in ("ne_proxy_forest.svg", "omega_per_gene_dist.svg",
                 "omega_vs_divergence.svg"):
        assert (out / "plots" / name).exists()
    print("stage 7 OK: per-species pooled omega, filtering, plots")


def test_free_ratio_column_optional(tmp_path):
    """free_ratio_records.tsv is optional: absent -> unchanged schema
    (no free_ratio columns at all); present -> one shared analysis pooled
    per species the same way two_ratio is, added as extra columns."""
    rd = tmp_path / "06_parsed"
    rd.mkdir()
    _write(rd / "m0_records.tsv",
           [_row("g1", "spA", 300, 100, 0.005, 0.10),
            _row("g1", "spB", 300, 100, 0.005, 0.10)])
    _write(rd / "two_ratio.spA_records.tsv",
           [_row("g1", "spA", 300, 100, 0.004, 0.10)])
    _write(rd / "two_ratio.spB_records.tsv",
           [_row("g1", "spB", 300, 100, 0.010, 0.10)])

    out_no_fr = tmp_path / "07_no_fr"
    mod.main(["--records-dir", str(rd), "-o", str(out_no_fr), "--bootstrap", "20"])
    hdr_no_fr = (out_no_fr / "ne_proxy.tsv").read_text().splitlines()[0]
    assert "free_ratio" not in hdr_no_fr

    # one shared free_ratio analysis: both species' branches in one table
    _write(rd / "free_ratio_records.tsv",
           [_row("g1", "spA", 300, 100, 0.003, 0.10),
            _row("g1", "spB", 300, 100, 0.009, 0.10)])
    out = tmp_path / "07"
    rc = mod.main(["--records-dir", str(rd), "-o", str(out), "--bootstrap", "20"])
    assert rc == 0

    proxy = (out / "ne_proxy.tsv").read_text().splitlines()
    hdr = proxy[0].split("\t")
    assert {"omega_free_ratio", "omega_free_ratio_ci_lo",
           "omega_free_ratio_ci_hi", "n_genes_free_ratio"} <= set(hdr)
    recs = {r.split("\t")[0]: dict(zip(hdr, r.split("\t"))) for r in proxy[1:]}
    assert abs(float(recs["spA"]["omega_free_ratio"]) - 0.03) < 1e-6   # 0.003/0.10
    assert abs(float(recs["spB"]["omega_free_ratio"]) - 0.09) < 1e-6   # 0.009/0.10
    assert recs["spA"]["n_genes_free_ratio"] == "1"
    print("stage 7 OK: free_ratio column optional, pooled per species like two_ratio")


if __name__ == "__main__":
    import tempfile
    test_pooled_omega()
    test_filter_focal()
    test_bootstrap_deterministic()
    with tempfile.TemporaryDirectory() as td:
        test_end_to_end(Path(td))
    print("ALL TESTS PASSED")
