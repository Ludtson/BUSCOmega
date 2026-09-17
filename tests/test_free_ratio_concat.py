"""
Stage 8 test: 08_free_ratio_concat.py -- concatenation with skip reporting,
and end-to-end parsing against the real 3-species pilot free-ratio mlc
(examples/3species_pilot/pilot_free_model.mlc), with codeml itself faked
(the .pml content doesn't matter once codeml is bypassed; only the tree's
tip set does, since that drives which species main() looks up afterward).
"""
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location(
    "free_ratio_concat", ROOT / "codes" / "core" / "08_free_ratio_concat.py")
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)

REAL_MLC = ROOT / "examples" / "3species_pilot" / "pilot_free_model.mlc"
TIPS = ["a_halleri", "a_thaliana", "c_grandiflora"]


def _write_pml(path: Path, taxa_seqs: dict[str, str]) -> None:
    width = len(next(iter(taxa_seqs.values())))
    lines = [f" {len(taxa_seqs)} {width}"]
    for t, s in taxa_seqs.items():
        lines.append(t)
        lines.append(s)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_concatenate_skips_mismatched_genes(tmp_path):
    aln = tmp_path / "codon_aln"
    aln.mkdir()
    _write_pml(aln / "g1.pml", {t: "ATG" * 4 for t in TIPS})
    _write_pml(aln / "g2.pml", {t: "ATG" * 6 for t in TIPS})
    # missing a tip -> skipped
    _write_pml(aln / "g3_bad_tips.pml", {t: "ATG" * 4 for t in TIPS[:2]})
    # unequal lengths within the gene -> skipped
    _write_pml(aln / "g4_bad_len.pml", {"a_halleri": "ATG" * 4,
                                        "a_thaliana": "ATG" * 5,
                                        "c_grandiflora": "ATG" * 4})

    dest = tmp_path / "concat.pml"
    taxa_order, n_used, n_sites, skipped = mod.concatenate(aln, dest, set(TIPS))
    assert n_used == 2
    assert n_sites == 12 + 18            # g1 + g2 widths
    assert {g for g, _ in skipped} == {"g3_bad_tips", "g4_bad_len"}
    assert dest.exists()


def test_end_to_end_against_real_free_ratio_mlc(tmp_path, monkeypatch):
    aln = tmp_path / "codon_aln"
    aln.mkdir()
    _write_pml(aln / "g1.pml", {t: "ATG" * 4 for t in TIPS})
    _write_pml(aln / "g2.pml", {t: "ATG" * 6 for t in TIPS})

    tree = tmp_path / "tree.nwk"
    tree.write_text("((a_halleri,a_thaliana),c_grandiflora);\n")

    # isolate just gene 1's block from the real pilot free-ratio output --
    # codeml itself is faked; only the mlc content downstream parsing sees
    real_text = REAL_MLC.read_text(errors="replace")
    gene1_block = next(mod.PARSE.iter_gene_blocks(real_text))

    def _fake_run_codeml(workdir, seqfile, treefile):
        workdir.mkdir(parents=True, exist_ok=True)
        mlc = workdir / "concat.mlc"
        mlc.write_text(gene1_block, encoding="utf-8")
        return mlc

    monkeypatch.setattr(mod, "run_codeml", _fake_run_codeml)

    out = tmp_path / "08"
    rc = mod.main(["--aln-dir", str(aln), "--tree", str(tree), "-o", str(out)])
    assert rc == 0

    rows = (out / "free_ratio_concat.tsv").read_text().splitlines()
    hdr = rows[0].split("\t")
    recs = {r.split("\t")[0]: dict(zip(hdr, r.split("\t"))) for r in rows[1:]}
    assert set(recs) == set(TIPS)
    # known values from the real fixture (see test_parse_codeml_output.py)
    assert abs(float(recs["c_grandiflora"]["omega"]) - 0.0001) < 1e-6
    assert abs(float(recs["a_halleri"]["omega"]) - 0.05296) < 1e-4
    assert abs(float(recs["a_thaliana"]["omega"]) - 0.09462) < 1e-4

    skipped = (out / "skipped_genes.tsv").read_text().splitlines()
    assert skipped == ["busco_id\treason"]           # both genes used, none skipped
    print("stage 8 OK: concatenate + real free-ratio mlc -> per-species omega, no pooling")


if __name__ == "__main__":
    import tempfile

    class _MP:
        def setattr(self, obj, name, val): setattr(obj, name, val)

    test_concatenate_skips_mismatched_genes(Path(tempfile.mkdtemp()))
    test_end_to_end_against_real_free_ratio_mlc(Path(tempfile.mkdtemp()), _MP())
    print("ALL TESTS PASSED")
