"""
Stage 5 test: 05_run_codeml.py orchestration.

codeml is external, so `_one_codeml` (which sets up a batch dir and runs
it) is replaced with a fake that reports a dataset count from the gene
list it is handed:
    - a batch containing "genebad" comes up one short  -> triggers retry
    - "genebad" alone produces nothing                 -> logged as failed
    - every other gene / batch succeeds
The rest is the real code path: plan reading, batching, the retry loop,
the manifest and the failed log. `concat_pml`, `fill_template` and
`n_datasets` are covered as unit tests.
"""
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location(
    "run_codeml", ROOT / "codes" / "core" / "05_run_codeml.py")
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


def _fake_one_codeml(workdir, template, tree_src, genes, aln_dir, out_name):
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / out_name).write_text("fake mlc\n")
    if "genebad" in genes:
        return len(genes) - 1          # 0 when alone, N-1 in a batch
    return len(genes)


def _mk_stage4(d: Path):
    (d / "ctl").mkdir(parents=True)
    (d / "trees").mkdir()
    (d / "ctl" / "m0.ctl").write_text(
        "seqfile = __SEQFILE__\ntreefile = __TREEFILE__\n"
        "outfile = __OUTFILE__\nndata = __NDATA__\nmodel = 0\n")
    (d / "trees" / "m0.nwk").write_text("3 1\n(a,(b,c));\n")
    (d / "stage4_analyses.tsv").write_text(
        "analysis\ttype\tmodel\tnssites\ttreefile\tctl_template\n"
        "m0\tm0\t0\t0\ttrees/m0.nwk\tctl/m0.ctl\n")


def _mk_alns(d: Path, genes):
    d.mkdir(parents=True)
    for g in genes:
        (d / f"{g}.pml").write_text(f"  3 6\na\nATGAAA\nb\nATGAAA\nc\nATGAAA\n")


def test_unit_helpers(tmp_path):
    assert mod.n_datasets("x\nDataset 1\ny\nDataset 2\n") == 2
    assert mod.n_datasets("CODONML (in paml version 4.10)\n") == 1
    assert mod.n_datasets("nothing here") == 0
    t = "a=__SEQFILE__ b=__TREEFILE__ c=__OUTFILE__ d=__NDATA__"
    assert mod.fill_template(t, seqfile="S", treefile="T", outfile="O",
                             ndata=5) == "a=S b=T c=O d=5"
    ad = tmp_path / "aln"
    _mk_alns(ad, ["g1", "g2"])
    dest = tmp_path / "batch.pml"
    mod.concat_pml(["g1", "g2"], ad, dest)
    assert dest.read_text().count("  3 6") == 2


def test_all_ok(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, "_one_codeml", _fake_one_codeml)
    s4 = tmp_path / "s4"; _mk_stage4(s4)
    aln = tmp_path / "aln"; _mk_alns(aln, ["g1", "g2", "g3", "g4"])
    out = tmp_path / "out"
    rc = mod.main(["-o", str(out), str(s4), "--aln-dir", str(aln),
                   "--batch-size", "2", "--skip-tool-check"])
    assert rc == 0
    man = (out / "stage5_manifest.tsv").read_text().splitlines()
    assert man[0].split("\t") == ["analysis", "batch", "n_genes", "n_ok",
                                  "n_datasets", "status", "wall_s", "mlc"]
    rows = [ln.split("\t") for ln in man[1:]]
    assert [r[1] for r in rows] == ["batch_0000", "batch_0001"]
    assert all(r[5] == "ok" for r in rows), rows
    assert (out / "stage5_failed.tsv").read_text().strip() == \
        "analysis\tbusco_id\treason"
    assert (out / "stage5.log").read_text().count("stage5") == 2


def test_retry_and_fail(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, "_one_codeml", _fake_one_codeml)
    s4 = tmp_path / "s4"; _mk_stage4(s4)
    aln = tmp_path / "aln"; _mk_alns(aln, ["g1", "genebad", "g3"])
    out = tmp_path / "out"
    rc = mod.main(["-o", str(out), str(s4), "--aln-dir", str(aln),
                   "--batch-size", "3", "--skip-tool-check"])
    assert rc == 0
    rows = [ln.split("\t") for ln in
            (out / "stage5_manifest.tsv").read_text().splitlines()[1:]]
    assert len(rows) == 1
    r = rows[0]
    assert r[2] == "3" and r[3] == "2"        # 3 genes, 2 recovered
    assert r[5] == "partial"
    failed = (out / "stage5_failed.tsv").read_text().splitlines()[1:]
    assert failed == ["m0\tgenebad\tcodeml produced no result block"]


if __name__ == "__main__":
    import tempfile

    class _MP:
        def setattr(self, o, n, v): setattr(o, n, v)

    with tempfile.TemporaryDirectory() as td:
        test_unit_helpers(Path(td))
    with tempfile.TemporaryDirectory() as td:
        test_all_ok(Path(td), _MP())
    with tempfile.TemporaryDirectory() as td:
        test_retry_and_fail(Path(td), _MP())
    print("ALL TESTS PASSED")
