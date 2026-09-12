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
    # named internal node labels (RAxML/IQ-TREE "node_labels" style: N0, N11,
    # ...), not just numeric support -- a real bug: these must be stripped
    # the same as `)95`, or they get picked up as bogus extra tips.
    assert mod.clean_newick(
        "(Tcacao:0.19,((A:0.12,(B:0.06,C:0.03)N4:0.01)N2:0.09,D:0.14)N1:0.19)N0;"
    ) == "(Tcacao,((A,(B,C)),D));"
    assert mod.tip_labels(mod.clean_newick(
        "(Tcacao:0.19,((A:0.12,(B:0.06,C:0.03)N4:0.01)N2:0.09,D:0.14)N1:0.19)N0;"
    )) == ["Tcacao", "A", "B", "C", "D"]
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


def test_looks_like_orthofinder(tmp_path):
    of_tree = ("(A_longest_isoforms:0.1,((B_longest_isoforms:0.05,"
              "C_longest_isoforms:0.05)N2:0.01,D_longest_isoforms:0.06)"
              "N1:0.02)N0;")
    assert mod.looks_like_orthofinder(Path("SpeciesTree_rooted.txt"), of_tree)
    assert mod.looks_like_orthofinder(Path("anything.nwk"), of_tree)  # by content
    plain = "((A:0.1,B:0.1)80:0.05,C:0.2);"
    assert not mod.looks_like_orthofinder(Path("tree.nwk"), plain)


def test_guess_matches():
    targets = {"SpeciesA", "SpeciesB", "Foo", "FooBar"}
    tips = ["SpeciesA", "speciesb", "SpeciesA_longest_isoforms",
           "Foo_v1", "Xyz"]
    g = mod.guess_matches(tips, targets)
    assert g["SpeciesA"] == ("SpeciesA", "exact", [])
    assert g["speciesb"] == ("SpeciesB", "case-insensitive", [])
    assert g["SpeciesA_longest_isoforms"] == ("SpeciesA", "suffix-stripped", [])
    # "Foo_v1" is a prefix match for "Foo", and NOT for "FooBar" (FooBar
    # does not start with Foo_v1, nor is Foo_v1 a prefix of it) -> unique
    assert g["Foo_v1"] == ("Foo", "unique-prefix", [])
    assert g["Xyz"] == (None, "no-match", [])
    # a genuinely ambiguous case: two targets both prefix-compatible
    amb = mod.guess_matches(["Sp"], {"Species1", "Species2"})
    assert amb["Sp"][1] == "ambiguous"
    assert sorted(amb["Sp"][2]) == ["Species1", "Species2"]


def test_cmd_match_never_auto_applies_uncertain(tmp_path):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    (fasta_dir / "Foo.faa").write_text(">x\nMK\n")
    (fasta_dir / "Bar.faa").write_text(">x\nMK\n")
    tree = tmp_path / "SpeciesTree_rooted_node_labels.txt"
    tree.write_text("(Foo_v1:0.1,(Bar_v1:0.05,Unrelated:0.05)N1:0.02)N0;\n")
    out = tmp_path / "map.tsv"

    rc = mod.main(["match", str(tree), "--to", str(fasta_dir), "-o", str(out)])
    assert rc == 1                              # nothing confident -> review needed

    text = out.read_text()
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        raise AssertionError(f"an uncertain match was left active: {line!r}")
    assert "# VERIFY" in text and "Foo_v1" in text
    assert "# NO MATCH" in text and "Unrelated" in text

    # read_rename must ignore every commented line -- applying this map
    # as-is renames nothing
    mapping = mod.read_rename(out)
    assert mapping == {}


def test_cmd_match_confident_map_is_usable(tmp_path):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    for sp in ("SpeciesA", "SpeciesB", "SpeciesC"):
        (fasta_dir / f"{sp}.faa").write_text(">x\nMK\n")
    tree = tmp_path / "tree.nwk"
    tree.write_text("(SpeciesA_longest_isoforms:0.1,"
                    "(SpeciesB_longest_isoforms:0.1,"
                    "SpeciesC_longest_isoforms:0.1)N1:0.1);\n")
    map_out = tmp_path / "map.tsv"
    rc = mod.main(["match", str(tree), "--to", str(fasta_dir),
                  "-o", str(map_out)])
    assert rc == 0

    clean_out = tmp_path / "clean.nwk"
    rc2 = mod.main(["prepare", str(tree), "-o", str(clean_out),
                    "--rename", str(map_out), "--match-to", str(fasta_dir)])
    assert rc2 == 0
    assert clean_out.read_text().strip() == "(SpeciesA,(SpeciesB,SpeciesC));"


if __name__ == "__main__":
    import tempfile
    test_clean_newick()
    test_rename_and_match()
    with tempfile.TemporaryDirectory() as td:
        test_prepare_cli(Path(td))
    with tempfile.TemporaryDirectory() as td:
        test_build_supermatrix(Path(td))
    with tempfile.TemporaryDirectory() as td:
        test_looks_like_orthofinder(Path(td))
    test_guess_matches()
    with tempfile.TemporaryDirectory() as td:
        test_cmd_match_never_auto_applies_uncertain(Path(td))
    with tempfile.TemporaryDirectory() as td:
        test_cmd_match_confident_map_is_usable(Path(td))
    print("ALL TESTS PASSED")
