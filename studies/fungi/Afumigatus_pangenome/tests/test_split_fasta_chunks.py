import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from split_fasta_chunks import split_fasta_chunks


def test_split_fasta_chunks_round_robin_distributes_evenly(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text(
        ">seq1\nAAAA\n"
        ">seq2\nCCCC\n"
        ">seq3\nGGGG\n"
        ">seq4\nTTTT\n"
    )
    chunks = split_fasta_chunks(str(fasta), n_chunks=2)
    assert len(chunks) == 2
    chunk0_text = "".join(chunks[0])
    chunk1_text = "".join(chunks[1])
    assert ">seq1" in chunk0_text and ">seq3" in chunk0_text
    assert ">seq2" in chunk1_text and ">seq4" in chunk1_text


def test_split_fasta_chunks_preserves_all_sequences(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text("".join(f">seq{i}\nAAAA\n" for i in range(7)))
    chunks = split_fasta_chunks(str(fasta), n_chunks=3)
    total_seqs = sum("".join(c).count(">") for c in chunks)
    assert total_seqs == 7


def test_split_fasta_chunks_handles_multiline_sequences(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text(">seq1\nAAAA\nCCCC\n>seq2\nGGGG\n")
    chunks = split_fasta_chunks(str(fasta), n_chunks=2)
    assert "".join(chunks[0]) == ">seq1\nAAAA\nCCCC\n"
    assert "".join(chunks[1]) == ">seq2\nGGGG\n"


def test_split_fasta_chunks_empty_input(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text("")
    chunks = split_fasta_chunks(str(fasta), n_chunks=3)
    assert chunks == [[], [], []]
