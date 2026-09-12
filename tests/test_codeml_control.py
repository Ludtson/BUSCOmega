"""
Stage 4 test: 04_codeml_control.py against the 3-species pilot tree.

Checks the model configuration is right (this is where the pilot's
"M0 and free-ratio weren't comparable" bug lived -- cleandata), the
foreground labelling is exact, and the plan table lists the right analyses.
"""
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

spec = importlib.util.spec_from_file_location(
    "codeml_control", ROOT / "codes" / "core" / "04_codeml_control.py")
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)

TREE = "(c_grandiflora,(a_halleri,a_thaliana));"


def test_newick_helpers():
    assert mod.strip_to_topology("3 1\n(c,(a,b));") == "(c,(a,b));"
    assert mod.strip_to_topology("(c:0.1,(a:0.2,b:0.3):0.4);") == "(c,(a,b));"
    assert sorted(mod.tip_labels(TREE)) == ["a_halleri", "a_thaliana",
                                            "c_grandiflora"]
    # named internal node labels (e.g. RAxML/IQ-TREE "node_labels": N0, N11,
    # ...) must be stripped like numeric support -- a real bug found on a
    # 22-taxon tree: these were being picked up as bogus extra tips.
    labelled_internal = ("(Tcacao:0.19,((A:0.12,(B:0.06,C:0.03)N4:0.01)N2:"
                         "0.09,D:0.14)N1:0.19)N0;")
    topo = mod.strip_to_topology(labelled_internal)
    assert topo == "(Tcacao,((A,(B,C)),D));", topo
    assert mod.tip_labels(topo) == ["Tcacao", "A", "B", "C", "D"]
    lab = mod.label_foreground(TREE, "a_halleri")
    assert lab == "(c_grandiflora,(a_halleri #1,a_thaliana));", lab
    assert lab.count("#1") == 1
    # a tip whose name is a substring of another must still match exactly once
    t2 = "(sp1,(sp10,sp1x));"
    assert mod.label_foreground(t2, "sp1").count("#1") == 1
    try:
        mod.label_foreground(TREE, "not_a_tip")
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for missing tip")


def test_write_control_families():
    m0 = mod.write_control(seqfile="S", treefile="T", outfile="O", ndata="N",
                           model=0, nssites="0")
    assert "model = 0" in m0
    assert "cleandata = 1" in m0
    tr = mod.write_control(seqfile="S", treefile="T", outfile="O", ndata="N",
                           model=2, nssites="0")
    assert "model = 2" in tr
    assert "cleandata = 1" in tr
    # cleandata override is honoured
    c0 = mod.write_control(seqfile="S", treefile="T", outfile="O", ndata="N",
                           model=0, nssites="0", cleandata=0)
    assert "cleandata = 0" in c0


def test_main_pilot(tmp_path):
    tree = tmp_path / "sp.nwk"
    tree.write_text(TREE + "\n")
    out = tmp_path / "04_codeml"
    rc = mod.main(["--tree", str(tree), "-o", str(out)])
    assert rc == 0

    m0 = (out / "ctl" / "m0.ctl").read_text()
    assert "model = 0" in m0 and "cleandata = 1" in m0
    for ph in (mod.PH_SEQ, mod.PH_TREE, mod.PH_OUT, mod.PH_NDATA):
        assert ph in m0, ph
    tr = (out / "ctl" / "two_ratio.ctl").read_text()
    assert "model = 2" in tr

    m0_tree = (out / "trees" / "m0.nwk").read_text().splitlines()
    assert m0_tree[0] == "3 1", m0_tree          # <ntax> <ntrees=1> for codeml
    assert m0_tree[1] == TREE
    for sp in ("a_halleri", "a_thaliana", "c_grandiflora"):
        t = (out / "trees" / "two_ratio" / f"{sp}.nwk").read_text().splitlines()
        assert t[0] == "3 1"
        assert t[1].count("#1") == 1
        assert f"{sp} #1" in t[1]

    plan = (out / "stage4_analyses.tsv").read_text().splitlines()
    assert plan[0].split("\t") == ["analysis", "type", "model", "nssites",
                                   "treefile", "ctl_template"]
    names = {ln.split("\t")[0] for ln in plan[1:]}
    assert names == {"m0", "two_ratio.a_halleri", "two_ratio.a_thaliana",
                     "two_ratio.c_grandiflora"}, names

    params = dict(ln.split("\t") for ln in
                  (out / "stage4_params.tsv").read_text().splitlines()[1:])
    assert params["cleandata"] == "1"
    assert (out / "stage4.log").read_text().count("stage4") == 2

    print("stage 4 OK: 1 M0 + 3 two_ratio, cleandata=1, foreground exact")


def test_main_focal_subset(tmp_path):
    tree = tmp_path / "sp.nwk"
    tree.write_text(TREE + "\n")
    out = tmp_path / "o"
    rc = mod.main(["--tree", str(tree), "-o", str(out),
                   "--focal", "a_halleri", "--analyses", "two_ratio"])
    assert rc == 0
    plan = (out / "stage4_analyses.tsv").read_text().splitlines()[1:]
    assert [ln.split("\t")[0] for ln in plan] == ["two_ratio.a_halleri"]
    assert not (out / "ctl" / "m0.ctl").exists()


def test_tree_data_mismatch(tmp_path, capsys=None):
    tree = tmp_path / "sp.nwk"
    tree.write_text("(c_grandiflora,(a_halleri,a_lyrata));\n")
    s3 = tmp_path / "s3"
    (s3 / "codon_aln").mkdir(parents=True)
    (s3 / "stage3_manifest.tsv").write_text("busco_id\tn_taxa\n")
    (s3 / "codon_aln" / "g1.pml").write_text(
        "  3 6\na_halleri\nATGAAA\na_thaliana\nATGAAA\nc_grandiflora\nATGAAA\n")
    try:
        mod.main(["--tree", str(tree), "-o", str(tmp_path / "o"),
                  "--stage3-dir", str(s3)])
    except SystemExit as e:
        assert e.code != 0
    else:
        raise AssertionError("expected the tree/data mismatch to abort")


if __name__ == "__main__":
    import tempfile
    test_newick_helpers()
    test_write_control_families()
    with tempfile.TemporaryDirectory() as td:
        test_main_pilot(Path(td))
    with tempfile.TemporaryDirectory() as td:
        test_main_focal_subset(Path(td))
    with tempfile.TemporaryDirectory() as td:
        test_tree_data_mismatch(Path(td))
    print("ALL TESTS PASSED")
