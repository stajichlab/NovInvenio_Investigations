"""Tests for bin/pangenome_gain_loss.py on tiny synthetic pangenome runs."""
import csv
import gzip
import hashlib
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "bin" / "pangenome_gain_loss.py"
spec = importlib.util.spec_from_file_location("pangenome_gain_loss", SCRIPT)
gl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gl)


def _write_run(d: Path, n_out: int, fams: dict, pairs=None, nonrep_out=0):
    """fams: name -> (bin, ingroup rep carriers 0..4, outgroup carriers 0..n_out).
    Ingroup: 4 reps (i1..i4, clades A,A,B,B) + 1 non-rep (i5). Outgroup: n_out strains,
    the last `nonrep_out` of them non-representative."""
    pdir = d / "cfg" / "pangenome"
    pdir.mkdir(parents=True)
    ins = ["i1", "i2", "i3", "i4", "i5"]
    outs = [f"o{k}" for k in range(1, n_out + 1)]
    strains = outs + ins
    with open(pdir / "samplesheet.with_clades.csv", "w") as fh:
        fh.write("GROUP,Species,Strain,Protein,DNA,GFF3,Short,TaxonGroup\n")
        for s in outs:
            fh.write(f"OUT,x,,p,d,g,{s},\n")
        for s, c in zip(ins, "AABBA"):
            fh.write(f"IN,y,,p,d,g,{s},{c}\n")
    with open(pdir / "strain_inventory.tsv", "w") as fh:
        fh.write("Short\tn_contigs\tn50\ttotal_length\tdedup_group\tis_representative\n")
        for k, s in enumerate(strains):
            rep = 0 if s == "i5" or (s in outs and outs.index(s) >= n_out - nonrep_out) else 1
            fh.write(f"{s}\t1\t1\t1\t{k}\t{rep}\n")
    with open(pdir / "presence_matrix.rescued.tsv", "w") as fh:
        fh.write("family\t" + "\t".join(strains) + "\n")
        with open(pdir / "frequency_table.tsv", "w") as ft:
            ft.write("family\tfrequency\tstrain_count\tbin\tfrequency_out\tstrain_count_out\tbin_out\n")
            for name, (b, n_in, n_o) in fams.items():
                cells = ["present" if k < n_o else "absent" for k in range(n_out)]
                cells += ["genome_only" if k < n_in else "absent" for k in range(4)]
                cells += ["present" if b == "nonrep_only" else "absent"]
                fh.write(name + "\t" + "\t".join(cells) + "\n")
                ft.write(f"{name}\t{n_in / 4:.4f}\t{n_in}\t{b}\t-\t-\t-\n")
    if pairs is not None:
        with open(pdir / "cooccurring_pairs.tsv", "w") as fh:
            fh.write("family_a\tfamily_b\tjaccard\tdirection_a\tasymmetry_a\tdirection_a_freq\n")
            for fa, ds, frac, df in pairs:
                fh.write(f"{fa}\tX\t0.5\t{ds}\t{frac}\t{df}\n")
    return d


def _run(run_dir, out_dir, *extra):
    subprocess.run([sys.executable, str(SCRIPT), "--run-dir", str(run_dir),
                    "--out-dir", str(out_dir), *extra], check=True, capture_output=True, text=True)
    with gzip.open(out_dir / "family_gain_loss.tsv.gz", "rt") as fh:
        return {r["family"]: r for r in csv.DictReader(fh, delimiter="\t")}


def test_rule_boundaries():
    # exactly 0.9 -> loss; exactly 0.1 -> gain; between -> ambiguous
    assert gl.polarize_by_frequency(9, 10) == "loss"
    assert gl.polarize_by_frequency(1, 10) == "gain"
    assert gl.polarize_by_frequency(2, 10) == "ambiguous"
    assert gl.polarize_by_frequency(8, 10) == "ambiguous"
    assert gl.polarize_by_frequency(0, 0) == "ambiguous"
    assert gl.polarize_direction(9, 10) == "ambiguous"
    assert gl.polarize_direction(10, 10) == "loss"
    assert gl.polarize_direction(0, 10) == "gain"
    with pytest.raises(ValueError):
        gl.validate_thresholds(0.1, 0.9)


def test_ten_strain_outgroup(tmp_path):
    fams = {
        "fGain": ("shell", 2, 0), "fGain01": ("cloud", 1, 1),     # 0.1 boundary -> gain
        "fLoss": ("shell", 3, 10), "fLoss09": ("soft_core", 3, 9),  # 0.9 boundary -> loss
        "fAmb": ("singleton", 1, 5), "fCore": ("core", 4, 0),
        "fOut": ("outgroup_only", 0, 10), "fNon": ("nonrep_only", 0, 0),
    }
    pairs = [("fGain", "gain", "0.0000", "gain"), ("fLoss09", "ambiguous", "0.9000", "loss"),
             ("fLoss09", "ambiguous", "0.9000", "loss")]
    run = _write_run(tmp_path / "run", 10, fams, pairs)
    t = _run(run, tmp_path / "o")
    d = {f: r["direction"] for f, r in t.items()}
    assert d == {"fGain": "gain", "fGain01": "gain", "fLoss": "loss", "fLoss09": "loss",
                 "fAmb": "ambiguous", "fCore": "conserved", "fOut": "absent", "fNon": "gain"}
    assert t["fLoss09"]["rule_strict"] == "ambiguous"
    assert t["fCore"]["rule_freq"] == "gain"
    assert t["fGain"]["ingroup_rep_count"] == "2"          # genome_only counts as present
    assert t["fNon"]["ingroup_nonrep_count"] == "1"
    val = (tmp_path / "o" / "validation.tsv").read_text()
    assert "direction_a_freq\tfamily_mismatches\t0\t0" in val
    assert "unique_family_a\t2\t2" in val
    scopes = (tmp_path / "o" / "summary_scopes.tsv").read_text()
    # selected = soft_core+shell+cloud+singleton: 2 gain, 2 loss, 1 ambiguous
    assert "selected:soft_core+shell+cloud+singleton\t5\t2\t2\t1\t" in scopes
    assert "selected+nonrep_only\t6\t3\t2\t1\t" in scopes
    clade = (tmp_path / "o" / "clade_counts.stratum.tsv").read_text().splitlines()
    assert clade[1].split("\t")[2:4] == ["A", "2"]


def test_two_strain_outgroup(tmp_path):
    fams = {"f0": ("shell", 2, 0), "f1": ("shell", 2, 1), "f2": ("cloud", 1, 2)}
    run = _write_run(tmp_path / "run", 2, fams)
    t = _run(run, tmp_path / "o")
    assert [t[f]["direction"] for f in ("f0", "f1", "f2")] == ["gain", "ambiguous", "loss"]
    assert all(t[f]["rule_freq"] == t[f]["rule_strict"] for f in t)
    dist = (tmp_path / "o" / "outgroup_fraction_distribution.tsv").read_text()
    assert "0/2\t1\tgain" in dist and "1/2\t1\tambiguous" in dist and "2/2\t1\tloss" in dist
    assert "not_available" in (tmp_path / "o" / "validation.tsv").read_text()


def test_outgroup_denominator_all_vs_representatives(tmp_path):
    # 10 outgroup strains, last 2 non-representative; family in 8 reps only:
    # all -> 8/10 ambiguous; representatives -> 8/8 loss
    run = _write_run(tmp_path / "run", 10, {"f": ("shell", 2, 8)}, nonrep_out=2)
    assert _run(run, tmp_path / "a")["f"]["direction"] == "ambiguous"
    assert _run(run, tmp_path / "b", "--outgroup-strains", "representatives")["f"]["direction"] == "loss"


def test_deterministic_output(tmp_path):
    run = _write_run(tmp_path / "run", 3, {"a": ("shell", 2, 0), "b": ("cloud", 1, 3)})
    digests = []
    for k in (1, 2):
        out = tmp_path / f"o{k}"
        _run(run, out)
        digests.append({p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in sorted(out.iterdir()) if p.name != "PROVENANCE.yaml"})
    assert digests[0] == digests[1]
