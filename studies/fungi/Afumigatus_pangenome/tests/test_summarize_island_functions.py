import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from summarize_island_functions import parse_domtblout, load_eligible_background, domain_enrichment


def _domtblout_line(target, query, i_evalue):
    # domtblout columns 1,4,13 are target name, query name, domain i-Evalue
    # (1-indexed); fill the rest with plausible placeholders.
    cols = ["-"] * 22
    cols[0] = target
    cols[1] = "PFxxxxx"
    cols[2] = "100"
    cols[3] = query
    cols[12] = str(i_evalue)
    return " ".join(cols)


def test_parse_domtblout_filters_by_independent_evalue(tmp_path):
    path = tmp_path / "hits.domtblout"
    path.write_text(
        _domtblout_line("DomainA", "famA", 1e-10) + "\n"
        + _domtblout_line("DomainB", "famA", 0.5) + "\n"  # too weak, filtered
        + "# comment\n"
    )
    hits = parse_domtblout([str(path)], max_ievalue=1e-3)
    assert hits == {"famA": {"DomainA"}}


def test_parse_domtblout_merges_multiple_files(tmp_path):
    path1 = tmp_path / "a.domtblout"
    path1.write_text(_domtblout_line("DomainA", "famA", 1e-10) + "\n")
    path2 = tmp_path / "b.domtblout"
    path2.write_text(_domtblout_line("DomainB", "famB", 1e-10) + "\n")
    hits = parse_domtblout([str(path1), str(path2)])
    assert hits == {"famA": {"DomainA"}, "famB": {"DomainB"}}


def test_load_eligible_background_excludes_core_and_singleton(tmp_path):
    path = tmp_path / "freq.tsv"
    path.write_text(
        "family\tfrequency\tstrain_count\tbin\n"
        "famA\t1.0\t10\tcore\n"
        "famB\t0.5\t5\tshell\n"
        "famC\t0.1\t1\tcloud\n"
        "famD\t0.01\t1\tsingleton\n"
    )
    background = load_eligible_background(str(path))
    assert background == {"famB", "famC"}


def test_domain_enrichment_finds_domain_overrepresented_in_islands():
    # DomainX is in ALL island-member families but only half the wider
    # background -- should come back strongly enriched.
    background = {"famA", "famB", "famC", "famD"}
    island_members = {"famA", "famB"}
    family_domains = {
        "famA": {"DomainX"}, "famB": {"DomainX"},
        "famC": {"DomainX"}, "famD": set(),
    }
    rows = domain_enrichment(island_members, background, family_domains)
    assert len(rows) == 1
    assert rows[0]["domain"] == "DomainX"
    assert rows[0]["n_with_domain_in_islands"] == 2
    assert rows[0]["n_with_domain_in_background"] == 3
    assert rows[0]["fisher_p"] < 1.0


def test_domain_enrichment_gives_p_one_when_uniformly_distributed():
    # DomainY appears in exactly the same PROPORTION inside and outside
    # the island-member set -- should not look enriched (high p-value).
    background = {"famA", "famB", "famC", "famD"}
    island_members = {"famA", "famB"}
    family_domains = {"famA": {"DomainY"}, "famC": {"DomainY"}}
    rows = domain_enrichment(island_members, background, family_domains)
    assert rows[0]["fisher_p"] > 0.4


def test_domain_enrichment_excludes_families_outside_background_not_crashes(capsys):
    # famZ (e.g. a singleton family an island happened to include) isn't
    # in the eligible shell+cloud background at all -- must be excluded
    # from the test, not crash scipy with an inconsistent 2x2 table.
    background = {"famA", "famB"}
    island_members = {"famA", "famZ"}
    family_domains = {"famA": {"DomainX"}, "famB": {"DomainX"}, "famZ": {"DomainX"}}
    rows = domain_enrichment(island_members, background, family_domains)
    assert rows[0]["n_island_families"] == 1  # famZ excluded, only famA remains
    err = capsys.readouterr().err
    assert "1 island-member families are outside" in err
