import pytest

from audit.prototypes import busco_lineage as bl

WHITELIST = bl.load_whitelist([
    "fungi_odb12", "basidiomycota_odb12", "agaricomycetes_odb12", "agaricales_odb12",
    "ascomycota_odb12", "eurotiomycetes_odb12", "eurotiales_odb12", "# comment", "",
])


def test_family_absent_falls_to_order():
    ranks = {"family": "Hymenogastraceae", "order": "Agaricales",
             "class": "Agaricomycetes", "phylum": "Basidiomycota"}
    lin, rank = bl.select_lineage(ranks, WHITELIST, "fungi_odb12")
    assert lin == "agaricales_odb12" and rank == "order"


def test_walk_skips_blank_ranks():
    ranks = {"family": "", "order": None, "class": "Agaricomycetes", "phylum": "Basidiomycota"}
    assert bl.select_lineage(ranks, WHITELIST, "fungi_odb12") == ("agaricomycetes_odb12", "class")


def test_fallback_when_nothing_matches():
    ranks = {"family": "Nope", "order": "Nope", "class": "Nope", "phylum": "Chytridiomycota"}
    assert bl.select_lineage(ranks, WHITELIST, "fungi_odb12") == ("fungi_odb12", None)


def test_case_and_space_insensitive():
    ranks = {"order": "  EUROTIALES "}
    assert bl.select_lineage(ranks, WHITELIST, "fungi_odb12")[0] == "eurotiales_odb12"


def test_pair_for_egap_tsv():
    ranks = {"order": "Agaricales", "phylum": "Basidiomycota"}
    spec, gen = bl.select_lineage_pair(ranks, WHITELIST, "basidiomycota_odb12")
    assert (spec, gen) == ("agaricales_odb12", "basidiomycota_odb12")
    assert bl.strip_odb(spec) == "agaricales"


def test_pair_never_returns_generic_as_specific_when_alternative_exists():
    ranks = {"phylum": "Basidiomycota", "class": "Agaricomycetes"}
    spec, gen = bl.select_lineage_pair(ranks, WHITELIST, "basidiomycota_odb12")
    assert spec == "agaricomycetes_odb12"


def test_pair_collapses_to_generic_when_nothing_else():
    ranks = {"phylum": "Basidiomycota"}
    spec, gen = bl.select_lineage_pair(ranks, WHITELIST, "basidiomycota_odb12")
    assert spec == gen == "basidiomycota_odb12"


def test_whitelist_file_roundtrip(tmp_path):
    p = tmp_path / "lineages.txt"
    p.write_text("fungi_odb12\n# c\nAgaricales_odb12\n\n")
    assert bl.load_whitelist(p) == {"fungi_odb12", "agaricales_odb12"}


def test_parse_busco_listing():
    text = """################################################
    Datasets available to be used with BUSCO v5:
    eukaryota_odb12
        - fungi_odb12
            - ascomycota_odb12
                - eurotiomycetes_odb12
    """
    got = bl.whitelist_from_busco_listing(text)
    assert got == {"eukaryota_odb12", "fungi_odb12", "ascomycota_odb12", "eurotiomycetes_odb12"}


def test_matches_fsp_reference_samplesheet_rows():
    """The two rows in fsp's own assets/samplesheet.csv, run against a
    whitelist copied from assets/fungi_busco_lineages.txt (subset)."""
    wl = bl.load_whitelist(["fungi_odb12", "basidiomycota_odb12", "agaricomycetes_odb12",
                            "cantharellales_odb12", "ascomycota_odb12",
                            "eurotiomycetes_odb12", "eurotiales_odb12"])
    s1 = {"family": "Hydnaceae", "order": "Cantharellales", "class": "Agaricomycetes",
          "phylum": "Basidiomycota"}
    s2 = {"family": "Trichocomaceae", "order": "Eurotiales", "class": "Eurotiomycetes",
          "phylum": "Ascomycota"}
    assert bl.select_lineage(s1, wl, "fungi_odb12")[0] == "cantharellales_odb12"
    assert bl.select_lineage(s2, wl, "fungi_odb12")[0] == "eurotiales_odb12"
