import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

from pangenome_matrix import PresenceMatrix, PRESENT, GENOME_ONLY, ABSENT
from extract_rescue_positions import (
    parse_tblastn_best_hit_positions, merge_best_hits, _parse_one_file, main,
)


def _outfmt6_line(qseqid, sseqid, pident, sstart, send, bitscore, qcovs):
    # std outfmt6 columns: qseqid sseqid pident length mismatch gapopen
    # qstart qend sstart send evalue bitscore + trailing qcovs
    return "\t".join([
        qseqid, sseqid, str(pident), "100", "0", "0", "1", "100",
        str(sstart), str(send), "1e-50", str(bitscore), str(qcovs),
    ])


def test_parse_tblastn_best_hit_positions_normalizes_minus_strand_coordinates():
    # sstart > send indicates a minus-strand hit -- start must normalize to
    # the smaller value, matching build_gene_positions.py's GFF3 convention.
    lines = [_outfmt6_line("famA", "s1|contig1", 100.0, 500, 300, 200.0, 100.0)]
    hits = parse_tblastn_best_hit_positions(lines)
    assert hits[("famA", "s1")] == ("contig1", 300, 200.0)


def test_parse_tblastn_best_hit_positions_keeps_highest_bitscore_per_pair():
    lines = [
        _outfmt6_line("famA", "s1|contig1", 100.0, 100, 200, 150.0, 100.0),
        _outfmt6_line("famA", "s1|contig2", 99.0, 400, 500, 300.0, 100.0),  # higher bitscore
    ]
    hits = parse_tblastn_best_hit_positions(lines)
    assert hits[("famA", "s1")] == ("contig2", 400, 300.0)


def test_parse_tblastn_best_hit_positions_filters_below_threshold():
    lines = [_outfmt6_line("famA", "s1|contig1", 50.0, 100, 200, 150.0, 100.0)]
    hits = parse_tblastn_best_hit_positions(lines, min_pident=90.0, min_qcov=80.0)
    assert hits == {}


def test_merge_best_hits_keeps_highest_bitscore_across_files():
    # Same key from two different "files" -- the reduction must pick the
    # higher bitscore regardless of which dict it came from or what order
    # they're merged in (parallel-safe: workers can finish in any order).
    per_file = [
        {("famA", "s1"): ("contig1", 100, 150.0)},
        {("famA", "s1"): ("contig2", 400, 300.0)},  # higher bitscore, different file
        {("famB", "s1"): ("contig1", 50, 90.0)},
    ]
    merged = merge_best_hits(per_file)
    assert merged[("famA", "s1")] == ("contig2", 400, 300.0)
    assert merged[("famB", "s1")] == ("contig1", 50, 90.0)


def test_merge_best_hits_order_independent():
    per_file = [
        {("famA", "s1"): ("contig1", 100, 150.0)},
        {("famA", "s1"): ("contig2", 400, 300.0)},
    ]
    forward = merge_best_hits(per_file)
    backward = merge_best_hits(list(reversed(per_file)))
    assert forward == backward


def test_parse_one_file_matches_direct_parse(tmp_path):
    # _parse_one_file (the Pool worker unit) must produce the same result
    # as calling parse_tblastn_best_hit_positions directly on the same
    # content -- it's just a picklable wrapper around file I/O + parsing.
    path = tmp_path / "hits.tsv"
    path.write_text(_outfmt6_line("famA", "s1|contig1", 100.0, 100, 200, 150.0, 100.0) + "\n")
    direct = parse_tblastn_best_hit_positions(path.read_text().splitlines())
    via_worker = _parse_one_file((str(path), 90.0, 80.0))
    assert direct == via_worker


def test_main_parallel_processes_matches_sequential_result(tmp_path):
    # End-to-end: --processes 2 must produce byte-identical output to the
    # default sequential (--processes 1) path on the same real inputs --
    # the whole point of the parallelization is that results don't depend
    # on how the work was split.
    strains = ["s1", "s2"]
    pm = PresenceMatrix(families=["famA", "famB"], strains=strains)
    pm.set_call("famA", "s1", GENOME_ONLY)
    pm.set_call("famB", "s2", GENOME_ONLY)
    matrix_path = tmp_path / "matrix.tsv"
    pm.to_tsv(matrix_path)

    file1 = tmp_path / "hits1.tsv"
    file1.write_text(_outfmt6_line("famA", "s1|contig1", 100.0, 100, 200, 150.0, 100.0) + "\n")
    file2 = tmp_path / "hits2.tsv"
    file2.write_text(_outfmt6_line("famB", "s2|contig2", 100.0, 300, 400, 150.0, 100.0) + "\n")

    def run(processes):
        output_path = tmp_path / f"out_{processes}.tsv"
        old_argv = sys.argv
        sys.argv = [
            "extract_rescue_positions.py",
            "--matrix", str(matrix_path),
            "--tblastn_tsv", str(file1),
            "--tblastn_tsv", str(file2),
            "--processes", str(processes),
            "--output", str(output_path),
        ]
        try:
            main()
        finally:
            sys.argv = old_argv
        return output_path.read_text()

    assert run(1) == run(2)


def test_main_only_writes_positions_for_genome_only_calls(tmp_path):
    # famA in s1 is GENOME_ONLY -> should get a position row.
    # famA in s2 is already PRESENT (protein-level) -> a hit there is real
    # but must be SKIPPED, since s2 already has a GFF3-based position and
    # doesn't need a rescue one (also guards against silently overwriting
    # a real position with a possibly-worse genomic-hit-derived one).
    strains = ["s1", "s2"]
    pm = PresenceMatrix(families=["famA"], strains=strains)
    pm.set_call("famA", "s1", GENOME_ONLY)
    pm.set_call("famA", "s2", PRESENT)
    matrix_path = tmp_path / "matrix.tsv"
    pm.to_tsv(matrix_path)

    tblastn_path = tmp_path / "hits.tsv"
    tblastn_path.write_text(
        _outfmt6_line("famA", "s1|contig1", 100.0, 100, 200, 150.0, 100.0) + "\n"
        + _outfmt6_line("famA", "s2|contig1", 100.0, 300, 400, 150.0, 100.0) + "\n"
    )

    output_path = tmp_path / "rescue_positions.tsv"
    old_argv = sys.argv
    sys.argv = [
        "extract_rescue_positions.py",
        "--matrix", str(matrix_path),
        "--tblastn_tsv", str(tblastn_path),
        "--output", str(output_path),
    ]
    try:
        main()
    finally:
        sys.argv = old_argv

    rows = output_path.read_text().splitlines()
    assert rows[0] == "Short\tfamily\tcontig\tstart"
    assert rows[1:] == ["s1\tfamA\tcontig1\t100"]
