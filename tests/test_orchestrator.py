"""
Orchestrator test: the audit join and run summary in run_buscomega.py,
built from a synthetic run directory (no subprocess, no binaries).
"""
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location(
    "run_buscomega", ROOT / "run_buscomega.py")
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


def test_dry_run_wires_analyses_and_free_ratio_concat(tmp_path, capsys):
    """--analyses passes through to Stage 4; --free-ratio-concat inserts
    Stage 8 between Stage 6 and Stage 7 and feeds Stage 7 --concat-tsv."""
    busco = tmp_path / "busco_out"
    for sp in ("sp1", "sp2", "sp3"):
        d = busco / sp / "run_x"
        d.mkdir(parents=True)
        (d / "full_table.tsv").write_text("g1\tComplete\tp1\n", encoding="utf-8")
    tree = tmp_path / "tree.nwk"
    tree.write_text("(sp1,sp2,sp3);\n")
    out = tmp_path / "run"

    rc = mod.main(["--busco-dir", str(busco), "--prt-dir", str(tmp_path / "prt"),
                  "--cds-dir", str(tmp_path / "cds"), "--tree", str(tree),
                  "-o", str(out), "--analyses", "m0,two_ratio,free_ratio",
                  "--free-ratio-concat", "--dry-run"])
    assert rc == 0
    log = capsys.readouterr().out

    assert "--analyses m0,two_ratio,free_ratio" in log
    stage8_pos = log.index("stage 8:")
    stage6_pos = log.index("stage 6:")
    stage7_pos = log.index("stage 7:")
    assert stage6_pos < stage8_pos < stage7_pos      # 8 runs between 6 and 7
    assert "08_free_ratio_concat.py" in log
    assert "--concat-tsv" in log and "free_ratio_concat.tsv" in log
    print("orchestrator OK: --analyses and --free-ratio-concat wired correctly")


def _mkrun(root: Path):
    (root / "01_common_scos").mkdir(parents=True)
    (root / "02_sequences").mkdir()
    (root / "03_alignments").mkdir()
    (root / "05_codeml_out").mkdir()
    (root / "07_ne_proxy").mkdir()

    (root / "01_common_scos" / "common_scos.tsv").write_text(
        "busco_id\tspA\tspB\n"
        "g1\tp1\tp2\ng2\tp3\tp4\ng3\tp5\tp6\ng4\tp7\tp8\n")
    (root / "02_sequences" / "stage2_manifest.tsv").write_text(
        "busco_id\tspecies\ng1\tspA\ng2\tspA\ng3\tspA\ng4\tspA\n")
    (root / "02_sequences" / "stage2_dropped.tsv").write_text("busco_id\treason\n")
    (root / "03_alignments" / "stage3_manifest.tsv").write_text(
        "busco_id\tn_taxa\ng1\t2\ng2\t2\ng4\t2\n")
    (root / "03_alignments" / "stage3_dropped.tsv").write_text(
        "busco_id\tstep\treason\ng3\tpal2nal\tmismatch\n")
    (root / "05_codeml_out" / "stage5_failed.tsv").write_text(
        "analysis\tbusco_id\treason\ntwo_ratio.spB\tg4\tno result block\n")
    (root / "07_ne_proxy" / "per_gene_omega.tsv").write_text(
        "species\tgene_id\tN\tS\tdN\tdS\tomega\tt\tqc_flag\tused\texcl_reason\n"
        "spA\tg1\t300\t100\t0.01\t0.1\t0.1\t0.2\tok\t1\t\n"
        "spA\tg2\t300\t100\t0.01\t0.1\t0.1\t0.2\tok\t1\t\n"
        "spB\tg1\t300\t100\t0.02\t0.1\t0.2\t0.2\tok\t1\t\n")
    for n, name in mod.STAGE_DIRS.items():
        (root / name).mkdir(exist_ok=True)
        (root / name / f"stage{n}.log").write_text(
            f"[t] stage{n} start\n[t] stage{n} done  elapsed_s=1\n")


def test_audit_and_summary(tmp_path):
    run = tmp_path / "run"
    _mkrun(run)
    mod.build_summary(run)
    mod.build_audit(run)

    gt = (run / "gene_tracking.tsv").read_text().splitlines()
    assert gt[0].split("\t") == ["busco_id", "stage1", "stage2", "stage3",
                                 "stage5_failed_for", "stage7_used_in"]
    rows = {r.split("\t")[0]: r.split("\t") for r in gt[1:]}
    assert set(rows) == {"g1", "g2", "g3", "g4"}
    assert rows["g3"][3].startswith("dropped:pal2nal")        # stage3 drop
    assert rows["g1"][5] == "spA,spB"                         # used in both
    assert rows["g2"][5] == "spA"
    assert rows["g4"][4] == "two_ratio.spB"                   # stage5 failure
    assert rows["g3"][5] == "-"                               # never reached 7

    summ = (run / "run_summary.tsv").read_text().splitlines()
    assert summ[0] == "stage\tlog_line"
    assert any(ln.startswith("7\t") for ln in summ)
    print("orchestrator OK: audit join + summary")


if __name__ == "__main__":
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        test_audit_and_summary(Path(td))
    print("ALL TESTS PASSED")
