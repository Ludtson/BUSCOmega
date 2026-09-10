"""
Test the optional species_tree.py helper: the stdlib parts only
(tree cleaning, rename, tip/data check, supermatrix concat). The IQ-TREE
call is not exercised.
"""
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location(
    "species_tree", ROOT / "codes" / "prep_optional" / "species_tree.py")
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


def test_clean_newick():
    # branch lengths + support + a NEXUS-ish wrapper + quoting
    raw = ("#NEXUS\nbegin trees;\n tree t = "
           "(('A. thaliana':0.1,B_sp:0.2)95:0.05,C_sp:0.3);\nend;")
    assert mod.clean_newick(raw) == "((A._thaliana,B_sp),C_sp);"
    # [&&NHX] comments
    assert mod.clean_newick("(A[&&NHX:x=1]:0.4,(B:0.1,C:0.1):0.2);") \
        == "(A,(B,C));"
    assert sorted(mod.tip_labels("((a_hal,b_sp),c_sp);")) == \
        ["a_hal", "b_sp", "c_sp"]


def test_rename_and_match():
    topo = "((Arabidopsis_thaliana,Arabidopsis_halleri),Capsella_grandiflora);"
    mp = {"Arabidopsis_thaliana": "a_thaliana",
          "Arabidopsis_halleri": "a_halleri",
          "Capsella_grandiflora": "c_grandiflora"}
    renamed = mod.apply_rename(topo, mp)
    assert set(mod.tip_labels(renamed)) == {"a_thaliana", "a_halleri",
                                            "c_grandiflora"}


def test_prepare_cli(tmp_path):
    intree = tmp_path / "pub.nwk"
    intree.write_text("((At:0.1,Ah:0.1)80:0.2,Cg:0.3);\n")
    ren = tmp_path / "map.tsv"
    ren.write_text("At\ta_thaliana\nAh\ta_halleri\nCg\tc_grandiflora\n")
    scos = tmp_path / "common_scos.tsv"
    scos.write_text("busco_id\ta_halleri\ta_thaliana\tc_grandiflora\n"
                    "g1\tx\ty\tz\n")
    out = tmp_path / "clean.nwk"
    rc = mod.main(["prepare", str(intree), "-o", str(out),
                   "--rename", str(ren), "--match-to", str(scos)])
    assert rc == 0
    txt = out.read_text().strip()
    assert txt == "((a_thaliana,a_halleri),c_grandiflora);"
    assert not txt.startswith("3 ")     # no <ntax> header -- Stage 4 adds it

    # a mismatch returns nonzero
    bad = tmp_path / "bad.nwk"
    bad.write_text("((a_thaliana,a_halleri),a_lyrata);\n")
    rc2 = mod.main(["prepare", str(bad), "-o", str(tmp_path / "o.nwk"),
                    "--match-to", str(scos)])
    assert rc2 == 1


def test_build_supermatrix(tmp_path):
    d = tmp_path / "prot_aln"
    d.mkdir()
    (d / "g1.aln").write_text(">a_hal\nMKT-\n>b_sp\nMKTA\n>c_sp\nMK-A\n")
    (d / "g2.aln").write_text(">a_hal\nGGG\n>b_sp\nGGA\n>c_sp\nGGA\n")
    phy, parts, ntax, nchar = mod.build_supermatrix(d, tmp_path / "super")
    assert (ntax, nchar) == (3, 7)
    body = phy.read_text().splitlines()
    assert body[0].strip() == "3 7"
    rows = dict(ln.split() for ln in body[1:])
    assert rows["a_hal"] == "MKT-GGG"
    assert rows["b_sp"] == "MKTAGGA"
    ptxt = parts.read_text().splitlines()
    assert ptxt[0] == "AA, g1 = 1-4"
    assert ptxt[1] == "AA, g2 = 5-7"


if __name__ == "__main__":
    import tempfile
    test_clean_newick()
    test_rename_and_match()
    with tempfile.TemporaryDirectory() as td:
        test_prepare_cli(Path(td))
    with tempfile.TemporaryDirectory() as td:
        test_build_supermatrix(Path(td))
    print("ALL TESTS PASSED")
