"""Tests for bin/hmm_domtbl_to_annotations.py (HMMER Pfam domtblout -> annotation TSV)."""
import csv
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
import hmm_domtbl_to_annotations as h2a  # noqa: E402

HEAD = (
    "#                                                                            --- full sequence --- "
    "-------------- this domain -------------   hmm coord   ali coord   env coord\n"
    "# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  "
    "of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target\n"
    "#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- "
    "--- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------\n"
)


def trailer(program, cut_ga=True):
    opt = "--cut_ga " if cut_ga else ""
    if program == "hmmsearch":
        q, t = "/db/2026-01-27-Pfam38.2/Pfam-A.hmm", "in.fa"
        cmd = f"hmmsearch -o /dev/null --domtblout x.domtbl --noali {opt}--cpu 8 {q} {t}"
    else:
        q, t = "in.fa", "/db/2026-01-27-Pfam38.2/Pfam-A.hmm"
        cmd = f"hmmscan -o /dev/null --domtblout x.domtbl --noali {opt}--cpu 8 {t} {q}"
    return ("#\n"
            f"# Program:         {program}\n"
            "# Version:         3.4 (Aug 2023)\n"
            f"# Pipeline mode:   {'SEARCH' if program == 'hmmsearch' else 'SCAN'}\n"
            f"# Query file:      {q}\n"
            f"# Target file:     {t}\n"
            f"# Option settings: {cmd} \n"
            "# Current dir:     /scratch/x\n"
            "# Date:            Tue Oct 13 16:25:23 2026\n"
            "# [ok]\n")


# (protein, pfam name, pfam acc with version, i-Evalue, domain score, env from)
HITS = [
    ("P1.1", "Kinase", "PF00069.32", "1e-50", "180.5", 40),
    ("P1.1", "SH2", "PF00017.30", "2e-10", "35.0", 5),       # earlier start: listed first
    ("P1.1", "Kinase", "PF00069.32", "3e-20", "70.2", 300),  # second Kinase domain
    ("P2.1", "14-3-3", "PF00244.27", "6e-100", "333.4", 1),
]


def search_line(prot, name, acc, iev, score, start, desc="kinase protein [Org x]"):
    return (f"{prot:<20} -            400 {name:<20} {acc:<10}   223   1e-100  366.1   5.7   1   2  "
            f"7.6e-114 {iev:>9} {score:>6}   5.7     1   223  {start:>4}   223  {start:>4}   250 1.00 {desc}\n")


def scan_line(prot, name, acc, iev, score, start):
    return (f"{name:<20} {acc:<10}   223 {prot:<20} -            400   1e-100  366.1   5.7   1   2  "
            f"7.6e-114 {iev:>9} {score:>6}   5.7     1   223  {start:>4}   223  {start:>4}   250 1.00 "
            "Some Pfam family description\n")


def write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


@pytest.fixture
def fasta(tmp_path):
    return write(tmp_path / "p.fa", ">P1.1 kinase protein [Org x]\nMKV\n>P0.1 hypothetical protein\nMKK\n"
                                    ">P2.1 14-3-3 protein [Org x]\nMAA\n")


def read_out(path: Path):
    with gzip.open(path, "rt") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        return rdr.fieldnames, list(rdr)


def test_hmmsearch_layout_multidomain(tmp_path, fasta):
    tbl = write(tmp_path / "s.domtbl", HEAD + "".join(search_line(*h) for h in HITS) + trailer("hmmsearch"))
    out = tmp_path / "ann" / "GCA_1.tsv.gz"
    s = h2a.convert([tbl], out, fasta=fasta, taxon_id="42", argv=["--domtbl", str(tbl)])
    cols, rows = read_out(out)
    assert cols[:11] == h2a.UNIPROT_COLUMNS and cols[11:] == h2a.EXTRA_COLUMNS
    assert [r["accession"] for r in rows] == ["P1.1", "P0.1", "P2.1"]      # FASTA order, all proteins
    p1 = rows[0]
    assert p1["pfam_ids"] == "PF00017|PF00069"                            # version stripped, start order
    assert p1["pfam_names"] == "SH2|Kinase"
    assert p1["pfam_best_ievalue"] == "2e-10|1e-50"                        # best of the two Kinase domains
    assert p1["pfam_best_domain_score"] == "35.0|180.5"
    assert p1["description"] == "kinase protein" and p1["taxon_id"] == "42"
    assert rows[1]["pfam_ids"] == "" and rows[1]["description"] == "hypothetical protein"
    assert all(r["go_ids"] == "" and r["interpro_ids"] == "" for r in rows)
    assert s == {"n_rows": 3, "n_proteins_with_pfam": 2, "n_domain_rows": 4, "n_pfam_families": 3,
                 "n_fasta_proteins": 3}
    prov = yaml.safe_load(Path(str(out) + ".provenance.yaml").read_text())
    src = prov["source_domtbl"][0]
    assert src["layout"] == "hmmsearch" and src["hmmer_version"] == "3.4 (Aug 2023)"
    assert src["cut_ga_in_command"] is True and src["sha256"] == h2a.sha256_of(tbl)
    assert prov["pfam_hmm"] == "/db/2026-01-27-Pfam38.2/Pfam-A.hmm"
    assert prov["pfam_release"].startswith("Pfam 38.2")
    assert prov["vocabularies_not_available"] == ["go", "interpro"]
    assert "not available" in prov["statement"]
    assert prov["checksum_sha256"] == h2a.sha256_of(out)
    assert prov["protein_fasta"]["sha256"] == h2a.sha256_of(fasta)
    assert "hmm_domtbl_to_annotations.py" in prov["derived_by"]


def test_hmmscan_layout_gives_same_table(tmp_path, fasta):
    t1 = write(tmp_path / "s.domtbl", HEAD + "".join(search_line(*h) for h in HITS) + trailer("hmmsearch"))
    t2 = write(tmp_path / "c.domtbl", HEAD + "".join(scan_line(*h) for h in reversed(HITS)) + trailer("hmmscan"))
    o1, o2 = tmp_path / "a.tsv.gz", tmp_path / "b.tsv.gz"
    h2a.convert([t1], o1, fasta=fasta, argv=[])
    h2a.convert([t2], o2, fasta=fasta, argv=[])
    assert o1.read_bytes() == o2.read_bytes()                              # deterministic, layout-independent
    prov = yaml.safe_load(Path(str(o2) + ".provenance.yaml").read_text())
    assert prov["source_domtbl"][0]["layout"] == "hmmscan"
    assert prov["pfam_hmm"] == "/db/2026-01-27-Pfam38.2/Pfam-A.hmm"        # hmmscan: target file


def test_layout_from_columns_without_trailer(tmp_path):
    t = write(tmp_path / "c.domtbl", HEAD + scan_line(*HITS[3]))
    lines = h2a.read_text(t)
    assert h2a.detect_layout(lines, h2a.parse_trailer(lines), "auto") == "hmmscan"
    t = write(tmp_path / "s.domtbl", HEAD + search_line(*HITS[3]))
    lines = h2a.read_text(t)
    assert h2a.detect_layout(lines, h2a.parse_trailer(lines), "auto") == "hmmsearch"


def test_tool_conflict_stops(tmp_path):
    t = write(tmp_path / "s.domtbl", HEAD + search_line(*HITS[3]) + trailer("hmmsearch"))
    with pytest.raises(SystemExit, match="trailer"):
        h2a.convert([t], tmp_path / "o.tsv.gz", tool="hmmscan", argv=[])


def test_without_fasta_only_hit_proteins(tmp_path):
    t = write(tmp_path / "s.domtbl", HEAD + "".join(search_line(*h) for h in HITS) + trailer("hmmsearch"))
    out = tmp_path / "o.tsv.gz"
    h2a.convert([t], out, argv=[])
    _, rows = read_out(out)
    assert [r["accession"] for r in rows] == ["P1.1", "P2.1"]
    assert rows[0]["description"] == "kinase protein"                      # from hmmsearch col 23


def test_protein_not_in_fasta_stops(tmp_path):
    fa = write(tmp_path / "p.fa", ">P1.1 x\nM\n")
    t = write(tmp_path / "s.domtbl", HEAD + "".join(search_line(*h) for h in HITS) + trailer("hmmsearch"))
    with pytest.raises(SystemExit, match="not in the FASTA"):
        h2a.convert([t], tmp_path / "o.tsv.gz", fasta=fa, argv=[])


@pytest.mark.skipif(shutil.which("zstd") is None, reason="zstd not in PATH")
def test_zstd_and_gz_input_and_chunks(tmp_path, fasta):
    # Two chunk tables (one zstd, one gzip) read as one proteome.
    t1 = write(tmp_path / "c0.domtbl", HEAD + "".join(search_line(*h) for h in HITS[:3]) + trailer("hmmsearch"))
    t2 = tmp_path / "c1.domtbl.gz"
    with gzip.open(t2, "wt") as fh:
        fh.write(HEAD + search_line(*HITS[3]) + trailer("hmmsearch"))
    subprocess.run(["zstd", "-q", str(t1), "-o", str(tmp_path / "c0.domtbl.zst")], check=True)
    out = tmp_path / "o.tsv.gz"
    s = h2a.convert([tmp_path / "c0.domtbl.zst", t2], out, fasta=fasta, argv=[])
    assert s["n_proteins_with_pfam"] == 2 and s["n_domain_rows"] == 4
    ref = tmp_path / "ref.tsv.gz"
    h2a.convert([write(tmp_path / "all.domtbl", HEAD + "".join(search_line(*h) for h in HITS)
                       + trailer("hmmsearch"))], ref, fasta=fasta, argv=[])
    assert out.read_bytes() == ref.read_bytes()
    prov = yaml.safe_load(Path(str(out) + ".provenance.yaml").read_text())
    assert len(prov["source_domtbl"]) == 2


def test_empty_result(tmp_path, fasta):
    t = write(tmp_path / "e.domtbl", HEAD + trailer("hmmsearch"))
    out = tmp_path / "o.tsv.gz"
    s = h2a.convert([t], out, fasta=fasta, argv=[])
    _, rows = read_out(out)
    assert len(rows) == 3 and all(r["pfam_ids"] == "" for r in rows)
    assert s["n_proteins_with_pfam"] == 0 and s["n_domain_rows"] == 0
    out2 = tmp_path / "o2.tsv.gz"
    h2a.convert([t], out2, argv=[])                                        # no FASTA: header only
    cols, rows = read_out(out2)
    assert rows == [] and cols == h2a.COLUMNS


def test_enrichment_loader_reads_output(tmp_path, fasta):
    """The reader in enrich_novelties_v2.py keys on 'accession' and reads pfam_ids/go_ids/interpro_ids
    by name with csv.DictReader; check the same calls work on our output."""
    t = write(tmp_path / "s.domtbl", HEAD + "".join(search_line(*h) for h in HITS) + trailer("hmmsearch"))
    out = tmp_path / "o.tsv.gz"
    h2a.convert([t], out, fasta=fasta, argv=[])
    with gzip.open(out, "rt") as fh:
        ann = {r["accession"]: (r.get("go_ids", ""), r.get("pfam_ids", ""), r.get("interpro_ids", ""))
               for r in csv.DictReader(fh, delimiter="\t")}
    assert ann["P2.1"] == ("", "PF00244", "")


def test_cli_main(tmp_path, fasta):
    t = write(tmp_path / "s.domtbl", HEAD + search_line(*HITS[3]) + trailer("hmmsearch"))
    out = tmp_path / "o.tsv.gz"
    assert h2a.main(["--domtbl", str(t), "--out", str(out), "--fasta", str(fasta)]) == 0
    assert out.is_file() and Path(str(out) + ".provenance.yaml").is_file()
