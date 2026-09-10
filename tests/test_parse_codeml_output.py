"""
Validation test for 06_parse_codeml_output.py against the real 3-species
pilot free-ratio output (examples/3species_pilot/pilot_free_model.mlc).

Gene #1 in that file (368-gene batched run) is known by manual inspection
to report:
    w (dN/dS) for branches:  0.00010 999.00000 0.05296 0.09462
    w ratios as node labels: c_grandiflora #0.0001, a_halleri #0.0529584,
                              a_thaliana #0.0946202
So this test checks the parser recovers exactly that, with the two boundary
values correctly flagged.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "codes" / "core"))
import importlib.util

spec = importlib.util.spec_from_file_location(
    "parse_codeml_output",
    Path(__file__).resolve().parent.parent / "codes" / "core" / "06_parse_codeml_output.py",
)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)

MLC = Path(__file__).resolve().parent.parent / "examples" / "3species_pilot" / "pilot_free_model.mlc"


def test_gene1_branches():
    text = MLC.read_text(errors="replace")
    blocks = list(mod.iter_gene_blocks(text))
    assert len(blocks) == 368, f"expected 368 gene blocks, got {len(blocks)}"

    records = mod.parse_free_ratio_block(blocks[0], gene_id="gene1", dataset_index=1)
    assert len(records) == 4, f"expected 4 branches (3 tips + 1 internal), got {len(records)}"

    by_tip = {r.tip: r for r in records if r.tip}
    assert set(by_tip) == {"c_grandiflora", "a_halleri", "a_thaliana"}, by_tip

    cg = by_tip["c_grandiflora"]
    assert abs(cg.omega - 0.0001) < 1e-6
    assert cg.qc_flag == "omega_boundary", cg

    ah = by_tip["a_halleri"]
    assert abs(ah.omega - 0.05296) < 1e-4
    assert ah.qc_flag == "ok", ah

    at = by_tip["a_thaliana"]
    assert abs(at.omega - 0.09462) < 1e-4
    assert at.qc_flag == "ok", at

    internal = [r for r in records if r.tip is None]
    assert len(internal) == 1
    assert internal[0].omega >= 999.0
    assert internal[0].qc_flag == "omega_boundary"

    print("gene 1 OK:", {t: (r.omega, r.qc_flag) for t, r in by_tip.items()})


def test_full_file_parses_all_genes():
    records = mod.parse_mlc_file(MLC, gene_id=None)
    # 368 genes x 4 branches = 1472, allow for any gene whose block failed to match
    n_genes = len({r.gene_id for r in records})
    assert n_genes >= 360, f"only parsed {n_genes}/368 genes cleanly"
    n_tips_resolved = sum(1 for r in records if r.tip is not None)
    n_terminal_expected = n_genes * 3  # 3 species per gene
    print(f"parsed {n_genes} genes, {len(records)} branch rows, "
          f"{n_tips_resolved} tip-resolved (~{n_terminal_expected} expected)")
    # Deterministic (sequence-order) tip resolution should hit every
    # terminal branch in every cleanly-parsed gene block, not just ~95%
    # the way value-matching on omega did (ties between branches).
    assert n_tips_resolved == n_terminal_expected, \
        f"expected exact tip resolution: {n_tips_resolved} != {n_terminal_expected}"


if __name__ == "__main__":
    test_gene1_branches()
    test_full_file_parses_all_genes()
    print("ALL TESTS PASSED")
