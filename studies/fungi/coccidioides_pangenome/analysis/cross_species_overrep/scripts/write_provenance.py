#!/usr/bin/env python3
"""Write PROVENANCE.yaml for analysis/cross_species_overrep (inputs sha256, outputs sha256,
commands, tool versions, parameters, SLURM job ids)."""
from __future__ import annotations

import hashlib
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
AF = NII / "studies/fungi/Afumigatus_pangenome"
CP = NII / "studies/fungi/coccidioides_pangenome"
OUT = CP / "analysis/cross_species_overrep"
PG = {
    "full_v072": AF / "results/full_v072/config/pangenome",
    "immitis_in_posadasii_out_v5": CP / "results/immitis_in_posadasii_out_v5/config_immitis_in_posadasii_out/pangenome",
    "posadasii_in_immitis_out_v2": CP / "results/posadasii_in_immitis_out_v2/config_posadasii_in_RS_WA211_out/pangenome",
}
STUDY = {"full_v072": AF, "immitis_in_posadasii_out_v5": CP, "posadasii_in_immitis_out_v2": CP}


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rel(p: Path) -> str:
    return str(p.relative_to(NII))


def main():
    inputs = []
    for run, pg in PG.items():
        st = STUDY[run]
        for p in [pg / "frequency_table.tsv", pg / "significant_islands.tsv", pg / "pfam.domtblout",
                  pg / "island_pfam_enrichment.tsv", pg / "cluster/tier1_rep_seq.fasta",
                  st / "analysis/gain_loss" / run / "family_gain_loss.tsv.gz",
                  st / "analysis/gain_loss_enrichment" / run / "family_pfam.tsv.gz"]:
            inputs.append(dict(path=rel(p), sha256=sha(p)))
        for sc in ("shell_soft_core", "all_accessory"):
            p = st / "analysis/gain_loss_enrichment" / run / sc / "enrichment_domains.tsv.gz"
            inputs.append(dict(path=rel(p), sha256=sha(p)))
    for p in [AF / "analysis/starship_island_overlap/full_v072/annotated_families.tsv.gz",
              Path("/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm.dat"),
              Path("/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-C")]:
        inputs.append(dict(path=str(p) if p.is_absolute() and not str(p).startswith(str(NII)) else rel(p),
                           sha256=sha(p)))
    outputs = [dict(path=rel(p), sha256=sha(p)) for p in sorted(OUT.rglob("*"))
               if p.is_file() and p.name != "PROVENANCE.yaml" and "__pycache__" not in p.parts]
    mm = subprocess.run([str(NII / ".pixi/envs/default/bin/mmseqs"), "version"],
                        capture_output=True, text=True).stdout.strip()
    import scipy
    import numpy
    import matplotlib
    import pandas
    rec = dict(
        created=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        author="Claude Code agent for Jason Stajich (not committed)",
        purpose="Pfam domains / homologous families over-represented in A. fumigatus full_v072 "
                "and in both Coccidioides runs (immitis_in_posadasii_out_v5, posadasii_in_immitis_out_v2)",
        derived_by=dict(
            scripts=["scripts/run_mmseqs_search.sbatch", "scripts/xsp_domains.py",
                     "scripts/xsp_families.py", "scripts/xsp_figures.py", "scripts/xsp_summary.py",
                     "scripts/write_provenance.py"],
            commands=[
                "sbatch scripts/run_mmseqs_search.sbatch  # SLURM job 30041523, partition stajichlab, milan",
                "python3 scripts/xsp_domains.py   # login node",
                "python3 scripts/xsp_families.py  # login node",
                "python3 scripts/xsp_figures.py   # login node",
                "python3 scripts/xsp_summary.py   # login node",
                "python3 scripts/write_provenance.py"],
            python=sys.version.split()[0],
            python_path=str(NII / ".pixi/envs/default/bin/python3"),
            packages=dict(scipy=scipy.__version__, numpy=numpy.__version__,
                          pandas=pandas.__version__, matplotlib=matplotlib.__version__),
            mmseqs_version=mm,
            mmseqs_parameters="easy-search -e 1e-5 -c 0.5 --cov-mode 0 (query and target coverage "
                              ">= 0.5), default -s 5.7, default --max-seqs 300, --threads 32; "
                              "searches af_vs_v5, v5_vs_af, af_vs_v2, v2_vs_af of tier1_rep_seq.fasta",
            best_hit_rule="highest bits per query; ties by lower evalue then target id",
            rbh_rule="best(Af family) = C and best(C) = Af family",
            slurm_job_ids=[30041523],
        ),
        definitions=dict(
            A="pipeline island test recomputed from pfam.domtblout (hmmscan --domE 1e-3, "
              "i-Evalue <= 1e-3, shell+cloud representatives), one-sided Fisher, BH over terms "
              "in >= 1 island family; counts match island_pfam_enrichment.tsv exactly; clan level new",
            B="bin/pangenome_gain_loss_enrichment.py outputs (hypergeometric, BH per test set and "
              "vocab, K >= 5), scopes shell_soft_core and all_accessory, gain and loss, Pfam and clan",
            C="new: hypergeometric one-sided, BH, K >= 5; accessory (soft_core+shell+cloud+singleton) "
              "vs universe core+soft_core+shell+cloud+singleton families with >= 1 Pfam (cut_ga); "
              "shell+cloud vs universe core+shell+cloud with >= 1 Pfam",
            tiers="1: q<0.05 in Af, v5 and v2; 2: q<0.05 in Af and exactly one of v5/v2; "
                  "3: p<0.05 in all three and not tier 1/2; only terms tested in >= 2 runs",
        ),
        license="derived from NII pangenome outputs; Pfam 38.2 annotation (CC0)",
        inputs=inputs,
        outputs=outputs,
    )
    with open(OUT / "PROVENANCE.yaml", "w") as fo:
        yaml.safe_dump(rec, fo, sort_keys=False, width=100)


if __name__ == "__main__":
    main()
