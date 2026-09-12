"""
Test the optional prep_optional/run_busco.py orchestrator: species
discovery, --resume skipping, per-species drop-and-continue on a BUSCO
failure, the manifest, and complete_pct parsing -- with `busco` itself
replaced by a fake (the binary is not assumed present in CI).
"""
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location(
    "run_busco", ROOT / "codes" / "prep_optional" / "run_busco.py")
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)

FULL_TABLE = """# BUSCO full table
# Busco id\tStatus\tSequence
g1\tComplete\tp1
g2\tComplete\tp2
g3\tMissing
g4\tFragmented\tp4
"""


def _fake_run_one(species, fasta, out_dir, lineage, mode, threads, log_dir):
    """Stand-in for run_one(): 'bad_species' fails, everyone else succeeds
    and gets a real full_table.tsv written where complete_pct() expects it."""
    import time
    if species == "bad_species":
        (log_dir / f"{species}.log").write_text("boom", encoding="utf-8")
        return (species, "failed", 0.1, None)
    run_dir = out_dir / species / "run_lineage"
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "full_table.tsv").write_text(FULL_TABLE, encoding="utf-8")
    return (species, "ok", 0.2, mod.complete_pct(out_dir, species))


def test_find_fasta(tmp_path):
    for name in ["a_halleri.faa", "a_thaliana.fa", "c_grandiflora.fna",
                "notes.txt"]:
        (tmp_path / name).write_text(">x\nMK\n")
    found = mod.find_fasta(tmp_path)
    assert [s for s, _ in found] == ["a_halleri", "a_thaliana", "c_grandiflora"]


def test_complete_pct(tmp_path):
    run_dir = tmp_path / "sp1" / "run_lineage"
    run_dir.mkdir(parents=True)
    (run_dir / "full_table.tsv").write_text(FULL_TABLE, encoding="utf-8")
    assert mod.already_done(tmp_path, "sp1")
    assert not mod.already_done(tmp_path, "sp2")
    pct = mod.complete_pct(tmp_path, "sp1")
    assert abs(pct - 50.0) < 1e-6          # 2 Complete / 4 total


def test_complete_pct_counts_duplicated_markers_once(tmp_path):
    """A Duplicated marker gets one full_table.tsv row per retained copy
    (e.g. a polyploid's homeologs) -- complete_pct() must count it as one
    marker, not one row per copy, or a real polyploid's single-copy
    fraction is understated the more copies it retains."""
    polyploid = """# Busco id\tStatus\tSequence
g1\tComplete\tp1
g2\tDuplicated\tp2a
g2\tDuplicated\tp2b
g2\tDuplicated\tp2c
g3\tMissing
"""
    run_dir = tmp_path / "poly" / "run_lineage"
    run_dir.mkdir(parents=True)
    (run_dir / "full_table.tsv").write_text(polyploid, encoding="utf-8")
    pct = mod.complete_pct(tmp_path, "poly")
    assert abs(pct - (100 / 3)) < 1e-6     # 1 Complete / 3 markers (g1,g2,g3)


def test_end_to_end(tmp_path, monkeypatch):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    for sp in ("sp1", "sp2", "bad_species"):
        (fasta_dir / f"{sp}.faa").write_text(">x\nMK\n")
    lineage = tmp_path / "lineage"
    lineage.mkdir()

    monkeypatch.setattr(mod, "check_tools", lambda: None)
    monkeypatch.setattr(mod, "run_one", _fake_run_one)

    out = tmp_path / "out"
    rc = mod.main([str(fasta_dir), "-o", str(out), "--lineage", str(lineage),
                  "--mode", "protein"])
    assert rc == 0                          # 2 ok >= 2 usable

    man = (out / "busco_run_manifest.tsv").read_text().splitlines()
    rows = {ln.split("\t")[0]: ln.split("\t")[1] for ln in man[1:]}
    assert rows == {"bad_species": "failed", "sp1": "ok", "sp2": "ok"}
    assert (out / "busco_run.log").read_text().count("busco") == 2
    assert (out / "logs" / "bad_species.log").exists()

    print("run_busco OK: 2 ok, 1 failed, manifest + logs written")


def test_resume_skips(tmp_path, monkeypatch):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    (fasta_dir / "sp1.faa").write_text(">x\nMK\n")
    lineage = tmp_path / "lineage"
    lineage.mkdir()
    out = tmp_path / "out"
    run_dir = out / "sp1" / "run_lineage"
    run_dir.mkdir(parents=True)
    (run_dir / "full_table.tsv").write_text(FULL_TABLE, encoding="utf-8")

    calls = []
    monkeypatch.setattr(mod, "check_tools", lambda: None)
    monkeypatch.setattr(mod, "run_one",
                        lambda *a, **k: calls.append(1) or (a[0], "ok", 0.1, 100.0))

    rc = mod.main([str(fasta_dir), "-o", str(out), "--lineage", str(lineage),
                  "--mode", "protein", "--resume"])
    assert rc == 1                          # only 1 usable species (< 2)
    assert calls == []                      # run_one never called -- skipped


if __name__ == "__main__":
    import tempfile

    class _MP:
        def setattr(self, obj, name, val): setattr(obj, name, val)

    test_find_fasta(Path(tempfile.mkdtemp()))
    test_complete_pct(Path(tempfile.mkdtemp()))
    test_end_to_end(Path(tempfile.mkdtemp()), _MP())
    test_resume_skips(Path(tempfile.mkdtemp()), _MP())
    print("ALL TESTS PASSED")
