"""
Stage 3 test: 03_codon_align.py orchestration, against a 4-gene fixture.

MAFFT and PAL2NAL are external binaries, so the two functions that call them
(`run_mafft`, `run_pal2nal`) are replaced with fakes that copy pre-made
fixture alignments or return a non-zero code. Everything else -- gene
discovery, per-gene dispatch, the drop rules, the manifest, the dropped
log, gap-fraction parsing -- is the real code path.

Fixture (tests/fixtures/stage3/):
    geneA  clean 3-taxon gene                 -> aligned, in manifest
    geneB  PAL2NAL reports a translation mismatch -> dropped at pal2nal
    geneC  protein file present, no CDS file  -> dropped at input
    geneD  MAFFT exits non-zero               -> dropped at mafft
Expected: 1 aligned, 3 dropped (one per step).
"""
import importlib.util
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIX = ROOT / "tests" / "fixtures" / "stage3"

spec = importlib.util.spec_from_file_location(
    "codon_align", ROOT / "codes" / "core" / "03_codon_align.py")
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


def fake_mafft(faa_in, aln_out, threads):
    gid = Path(aln_out).stem
    if gid == "geneD":
        return 1, "mafft: internal error, exiting"
    shutil.copy(FIX / "aln" / f"{gid}.aln", aln_out)
    return 0, ""


def fake_pal2nal(prot_aln, cds, pml_out):
    gid = Path(pml_out).stem
    if gid == "geneB":
        return 1, ("Error: inconsistency between the following pep and nuc "
                   "seqs: sp2")
    shutil.copy(FIX / "pml" / f"{gid}.pml", pml_out)
    return 0, ""


def test_codon_align(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, "run_mafft", fake_mafft)
    monkeypatch.setattr(mod, "run_pal2nal", fake_pal2nal)

    rc = mod.main([str(FIX), "-o", str(tmp_path), "--skip-tool-check"])
    assert rc == 0

    aln = sorted(p.name for p in (tmp_path / "prot_aln").glob("*.aln"))
    assert aln == ["geneA.aln"], aln
    pml = sorted(p.name for p in (tmp_path / "codon_aln").glob("*.pml"))
    assert pml == ["geneA.pml"], pml

    man = (tmp_path / "stage3_manifest.tsv").read_text().splitlines()
    assert man[0].split("\t") == ["busco_id", "n_taxa", "prot_aln_aa",
                                  "codon_aln_bp", "gap_frac", "mafft_rc",
                                  "pal2nal_rc"]
    assert len(man) == 2
    a = dict(zip(man[0].split("\t"), man[1].split("\t")))
    assert a["busco_id"] == "geneA"
    assert a["n_taxa"] == "3"
    assert a["prot_aln_aa"] == "13"
    assert a["codon_aln_bp"] == "39"
    # 9 gap chars / 117 cells
    assert abs(float(a["gap_frac"]) - 9 / 117) < 1e-4, a["gap_frac"]

    drops = (tmp_path / "stage3_dropped.tsv").read_text().splitlines()
    assert drops[0] == "busco_id\tstep\treason"
    by_gene = {ln.split("\t")[0]: ln.split("\t")[1] for ln in drops[1:]}
    assert by_gene == {"geneB": "pal2nal", "geneC": "input", "geneD": "mafft"}, \
        by_gene

    print("stage 3 OK: 1 aligned, 3 dropped (input/mafft/pal2nal)")


def test_paml_stats():
    n_taxa, n_sites, gap = mod.paml_stats(
        (FIX / "pml" / "geneA.pml").read_text())
    assert (n_taxa, n_sites) == (3, 39)
    assert abs(gap - 9 / 117) < 1e-4


def test_strip_stop():
    assert mod.strip_stop("MKTAYI*") == "MKTAYI"
    assert mod.strip_stop("MKTAYI.") == "MKTAYI"
    assert mod.strip_stop("MKTAYI") == "MKTAYI"


if __name__ == "__main__":
    import tempfile

    class _MP:
        def setattr(self, obj, name, val): setattr(obj, name, val)

    test_strip_stop()
    test_paml_stats()
    with tempfile.TemporaryDirectory() as td:
        test_codon_align(Path(td), _MP())
    print("ALL TESTS PASSED")
