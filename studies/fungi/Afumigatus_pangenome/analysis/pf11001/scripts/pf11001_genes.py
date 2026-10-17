#!/usr/bin/env python3
"""PF11001 check, stage 2: gene-level tables from the single-profile
hmmsearch (--cut_ga) over every proteome in hmmsearch/manifest.tsv.

Writes (Afumigatus = AF tables dir, Coccidioides = CP tables dir):
  genome_counts.<set>.tsv          one row per proteome: n proteins, n PF11001 proteins
  survey_by_taxon_group.tsv        survey proteomes summarised by TaxonGroup
  survey_genomes.tsv               survey proteomes, one row each
  genes.afum_all_strains.tsv.gz    every A. fumigatus-study hit protein
  genes.cocci_all_strains.tsv.gz   every Coccidioides-study hit protein
  af293_pf11001_proteins.tsv       Af293 reference table
  cocci_reference_pf11001_proteins.tsv   RS + UCSF_Cp_Silveira table
  hits_vs_families.<run>.tsv       how hit proteins map onto pangenome families
  tandem_summary.<set>.tsv         per strain: n hits, n clusters (<=5 gene ranks apart)
  starship_captain_overlap.tsv     family-level overlap tests
"""
from __future__ import annotations

import csv
import gzip
import io
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

from scipy.stats import fisher_exact

NII = Path("/bigdata/stajichlab/jstajich/projects/NovInvenio_Investigations")
AF = NII / "studies/fungi/Afumigatus_pangenome"
CP = NII / "studies/fungi/coccidioides_pangenome"
D = AF / "analysis/pf11001"
HD = D / "hmmsearch"
AFT = D / "tables"
CPT = CP / "analysis/pf11001/tables"
RUNS = {
    "full_v072": AF / "results/full_v072/config/pangenome",
    "immitis_in_posadasii_out_v5": CP / "results/immitis_in_posadasii_out_v5/config_immitis_in_posadasii_out/pangenome",
    "posadasii_in_immitis_out_v2": CP / "results/posadasii_in_immitis_out_v2/config_posadasii_in_RS_WA211_out/pangenome",
}
GL = {
    "full_v072": AF / "analysis/gain_loss/full_v072/family_gain_loss.tsv.gz",
    "immitis_in_posadasii_out_v5": CP / "analysis/gain_loss/immitis_in_posadasii_out_v5/family_gain_loss.tsv.gz",
    "posadasii_in_immitis_out_v2": CP / "analysis/gain_loss/posadasii_in_immitis_out_v2/family_gain_loss.tsv.gz",
}
SET_RUNS = {"afum": ["full_v072"], "cocci": ["immitis_in_posadasii_out_v5", "posadasii_in_immitis_out_v2"]}
STARSHIP = AF / "analysis/starship_island_overlap"
CAPTAIN = {"afum": AF / "results/captain_gene/DUF3435_vs_study.tblout",
           "cocci": CP / "analysis/starship_island_overlap/inputs/DUF3435_vs_study.tblout.zst"}
TANDEM_RANKS = 5


def zread(path):
    path = Path(path)
    if path.suffix == ".zst":
        proc = subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE, text=True)
        return proc.stdout
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def write_tsv(path, rows, header=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if header is None:
        header = []
        for r in rows:
            for k in r:
                if k not in header:
                    header.append(k)
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "wt") as fo:
        fo.write("\t".join(header) + "\n")
        for r in rows:
            fo.write("\t".join(str(r.get(h, "")) for h in header) + "\n")


def fasta_lengths(path):
    out = {}
    name = None
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = n
                name = line[1:].split(None, 1)[0]
                n = 0
            else:
                n += len(line.strip())
    if name:
        out[name] = n
    return out


def main():
    man = {}
    with open(HD / "manifest.tsv", newline="") as fh:
        for i, r in enumerate(csv.DictReader(fh, delimiter="\t"), 1):
            man[str(i)] = r
    nseq = {}
    with open(HD / "seq_counts.tsv") as fh:
        next(fh)
        for line in fh:
            idx, s, short, n = line.split()
            nseq[idx] = int(n)
    # hits: (idx, pid) -> list of domain rows
    hits = defaultdict(list)
    for line in zread(HD / "all_proteomes.PF11001.domtbl.idx.zst"):
        p = line.rstrip("\n").split(None, 23)
        idx, pid = p[0], p[1]
        hits[(idx, pid)].append(dict(tlen=int(p[3]), full_E=p[7], full_score=float(p[8]),
                                     c_E=p[12], i_E=p[13], dom_score=float(p[14]),
                                     hmm_from=int(p[16]), hmm_to=int(p[17]),
                                     env_from=int(p[20]), env_to=int(p[21]),
                                     desc=p[23] if len(p) > 23 else ""))
    # other Pfam domains (optional)
    other = defaultdict(list)
    scan = HD / "hit_proteins_vs_PfamA.cutga.domtbl.zst"
    have_scan = scan.exists()
    if have_scan:
        for line in zread(scan):
            if line.startswith("#"):
                continue
            p = line.split()
            other[p[3]].append((int(p[19]), p[0], p[1].split(".")[0]))
    lengths = {}
    if (HD / "hit_protein_lengths.tsv").exists():
        with open(HD / "hit_protein_lengths.tsv") as fh:
            next(fh)
            for line in fh:
                k, v = line.split()
                lengths[k] = int(v)

    # ---------------- per-genome counts ----------------
    per_genome = defaultdict(lambda: [0, 0])
    for (idx, pid), doms in hits.items():
        per_genome[idx][0] += 1
        per_genome[idx][1] += len(doms)
    grows = defaultdict(list)
    for idx, r in man.items():
        grows[r["set"]].append(dict(set=r["set"], study=r["study"], group=r["group"], species=r["species"],
                                    short=r["short"], taxon_group=r["taxon_group"],
                                    n_proteins=nseq.get(idx, "NA"),
                                    n_pf11001_proteins=per_genome[idx][0],
                                    n_pf11001_domains=per_genome[idx][1],
                                    pep=Path(r["pep_path"]).name))
    write_tsv(AFT / "genome_counts.afum.tsv", grows["afum"])
    write_tsv(CPT / "genome_counts.cocci.tsv", grows["cocci"])
    write_tsv(AFT / "survey_genomes.tsv", grows["survey"])
    by_tg = defaultdict(list)
    for r in grows["survey"]:
        by_tg[(r["study"].split("/")[-1], r["taxon_group"])].append(r)
    srows = []
    for (study, tg), rs in sorted(by_tg.items()):
        c = [r["n_pf11001_proteins"] for r in rs]
        srows.append(dict(study=study, taxon_group=tg, n_proteomes=len(rs),
                          n_with_ge1_hit=sum(x > 0 for x in c), min_hits=min(c),
                          median_hits=median(c), max_hits=max(c),
                          species_with_hits=";".join(f"{r['short']}={r['n_pf11001_proteins']}" for r in rs)))
    write_tsv(AFT / "survey_by_taxon_group.tsv", srows)
    for s in ("afum", "cocci"):
        c = [r["n_pf11001_proteins"] for r in grows[s]]
        print(f"{s}: proteomes={len(c)} with_hit={sum(x>0 for x in c)} min={min(c)} "
              f"median={median(c)} max={max(c)} total={sum(c)}", file=sys.stderr)
        for grp in ("IN", "OUT"):
            cg = [r["n_pf11001_proteins"] for r in grows[s] if r["group"] == grp]
            if cg:
                print(f"  {grp}: n={len(cg)} min={min(cg)} median={median(cg)} max={max(cg)} "
                      f"total={sum(cg)} dist={sorted(Counter(cg).items())}", file=sys.stderr)

    # ---------------- per-run family context ----------------
    ctx = {}
    for run, pg in RUNS.items():
        m2r = {}
        with open(pg / "cluster/tier1_cluster.tsv") as fh:
            for line in fh:
                rep, mem = line.rstrip("\n").split("\t")
                m2r[mem] = rep
        freq = {}
        with open(pg / "frequency_table.tsv", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                freq[r["family"]] = r
        fam_isl = defaultdict(list)
        with open(pg / "significant_islands.tsv", newline="") as fh:
            for i, r in enumerate(csv.DictReader(fh, delimiter="\t")):
                for f in r["member_families"].split(","):
                    fam_isl[f].append((i, int(r["n_strains"])))
        gl = {}
        with zread(GL[run]) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                gl[r["family"]] = r["direction"]
        bys = defaultdict(list)
        for line in zread(pg / "gene_positions.tsv.zst"):
            p = line.rstrip("\n").split("\t")
            if p[0] == "Short":
                continue
            bys[(p[0], p[2])].append((int(p[3]), int(p[4]), p[1]))
        ctx[run] = dict(m2r=m2r, freq=freq, fam_isl=fam_isl, gl=gl, bys=bys)
    # extra Af293 coordinates (unplaced Swiss-Prot entries)
    extra = {}
    with open(STARSHIP / "inputs/af293_unplaced_swissprot_coords.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            extra[(r["Short"], r["protein_id"])] = (r["contig"], int(r["start"]), int(r["end"]))

    # captain genes (DUF3435, full-seq E <= 1e-5, as in starship_island_overlap.py)
    captains = {}
    for s, path in CAPTAIN.items():
        cs = set()
        for line in zread(path):
            if line.startswith("#"):
                continue
            p = line.split()
            if float(p[4]) <= 1e-5:
                short, pid = p[0].split("|", 1)
                cs.add((short, pid))
        captains[s] = cs
    print(f"captain genes: afum={len(captains['afum'])} cocci={len(captains['cocci'])}", file=sys.stderr)

    needed = {(man[idx]["short"], pid) for (idx, pid) in hits if man[idx]["set"] in ("afum", "cocci")}
    for cs in captains.values():
        needed |= cs

    # gene ranks per strain/contig; rank dict kept only for hit and captain genes
    def ranks_for(run, add_extra=False):
        bys = ctx[run].pop("bys")
        if add_extra:
            placed = {(short, pid) for (short, c), lst in bys.items() if short == "Asfu_Af293" for (_, _, pid) in lst}
            for (short, pid), (c, s, e) in extra.items():
                if (short, pid) not in placed:
                    bys[(short, c)].append((s, e, pid))
        rk = {}
        for (short, c), lst in bys.items():
            lst.sort()
            for i, (s, e, pid) in enumerate(lst):
                if (short, pid) in needed:
                    rk[(short, pid)] = (c, s, e, i)
        return rk, bys
    rank_af, bys_af = ranks_for("full_v072", add_extra=True)
    rank_cv5, bys_cv5 = ranks_for("immitis_in_posadasii_out_v5")
    rank_cv2, bys_cv2 = ranks_for("posadasii_in_immitis_out_v2")

    def nearest_captain(short, pid, rk, bys, cset):
        if (short, pid) not in rk:
            return "", ""
        c, s, e, i = rk[(short, pid)]
        best = None
        for (s2, e2, pid2) in bys.get((short, c), []):
            if (short, pid2) in cset:
                d = max(0, max(s, s2) - min(e, e2))
                j = rk[(short, pid2)][3]
                if best is None or d < best[0]:
                    best = (d, abs(j - i))
        return (best[0], best[1]) if best else ("none_on_contig", "")

    # Starship elements (Af293)
    elements = []
    with open(STARSHIP / "inputs/af293_starships_table_s7.elements.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            elements.append(r)

    def in_element(short, c, s, e):
        mid = (s + e) / 2
        hitl = [f"{r['element_id']}({r['label']})" for r in elements
                if r["Short"] == short and r["contig"] == c and int(r["start"]) <= mid <= int(r["end"])]
        return ",".join(hitl)

    # reference contig lengths
    clen = {}
    for short, fa in (("Asfu_Af293", AF / "data_dir/dna/UP000002530_330879.dna.fa"),
                      ("RS", CP / "data_dir/dna/RS.dna.fa"),
                      ("UCSF_Cp_Silveira", CP / "data_dir/dna/UCSF_Cp_Silveira.dna.fa")):
        for c, n in fasta_lengths(fa).items():
            clen[(short, c)] = n
    chrname = {}
    with open(AF / "data_dir/dna/UP000002530_330879.dna.fa") as fh:
        for line in fh:
            if line.startswith(">") and "chromosome" in line:
                c = line[1:].split()[0]
                chrname[c] = "chr" + line.split("chromosome")[1].split(",")[0].strip()

    # Af293 UniProt annotation
    af_ann = {}
    with gzip.open(AF / "annotations/UP000002530.tsv.gz", "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            af_ann[r["accession"]] = r
    # cocci reference GFF3 products
    cp_prod = {}
    for short in ("RS", "UCSF_Cp_Silveira"):
        with open(CP / f"data_dir/gff3/{short}.gff3") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) == 9 and p[2] == "mRNA":
                    a = dict(kv.split("=", 1) for kv in p[8].strip(";").split(";") if "=" in kv)
                    cp_prod[(short, a.get("ID"))] = a.get("product", "")

    # ---------------- gene rows ----------------
    gene_rows = defaultdict(list)
    for (idx, pid), doms in sorted(hits.items(), key=lambda x: (int(x[0][0]), x[0][1])):
        r = man[idx]
        if r["set"] not in ("afum", "cocci"):
            continue
        s = r["set"]
        short = r["short"]
        key = f"{short}|{pid}"
        best = max(doms, key=lambda d: d["dom_score"])
        od = sorted(other.get(key, []))
        others = [f"{acc}:{nm}" for _, nm, acc in od if acc != "PF11001"]
        row = dict(set=s, group=r["group"], species=r["species"], short=short, protein_id=pid,
                   protein_length=lengths.get(key, doms[0]["tlen"]),
                   n_pf11001_domains=len(doms),
                   pf11001_env=";".join(f"{d['env_from']}-{d['env_to']}" for d in doms),
                   pf11001_hmm=";".join(f"{d['hmm_from']}-{d['hmm_to']}" for d in doms),
                   best_dom_score=best["dom_score"], best_i_Evalue=best["i_E"],
                   full_seq_score=best["full_score"], full_seq_Evalue=best["full_E"],
                   other_pfam_cutga=(",".join(dict.fromkeys(others)) if have_scan else "not_run"),
                   pfam_architecture=("|".join(nm for _, nm, _ in od) if have_scan else "not_run"),
                   fasta_description=best["desc"])
        rk, bys = ((rank_af, bys_af) if s == "afum" else (rank_cv5, bys_cv5))
        loc = rk.get((short, pid))
        if loc is None and s == "cocci":
            loc = rank_cv2.get((short, pid))
            rk, bys = rank_cv2, bys_cv2
        if loc:
            c, st, en, i = loc
            row.update(contig=c, start=st, end=en, gene_rank_on_contig=i)
            cap = nearest_captain(short, pid, rk, bys, captains[s])
            row.update(nearest_DUF3435_captain_bp=cap[0], nearest_DUF3435_captain_gene_ranks=cap[1])
        else:
            row.update(contig="NA", start="", end="", gene_rank_on_contig="")
        for run in SET_RUNS[s]:
            cx = ctx[run]
            fam = cx["m2r"].get(key)
            tag = run.split("_")[0] if s == "cocci" else ""
            pre = f"{tag}_" if tag else ""
            if fam is None:
                row[f"{pre}family"] = "not_in_clusters"
                continue
            fr = cx["freq"].get(fam, {})
            isl = cx["fam_isl"].get(fam, [])
            row[f"{pre}family"] = fam
            row[f"{pre}family_bin"] = fr.get("bin", "NA")
            row[f"{pre}family_strain_count"] = fr.get("strain_count", "NA")
            row[f"{pre}family_direction"] = cx["gl"].get(fam, "NA")
            row[f"{pre}family_n_islands"] = len(isl)
            row[f"{pre}family_n_multi_strain_islands"] = sum(n >= 2 for _, n in isl)
            row[f"{pre}family_island_idx"] = ",".join(str(i) for i, _ in isl[:30]) + (",..." if len(isl) > 30 else "")
        gene_rows[s].append(row)

    # tandem clusters per strain
    for s, rows in gene_rows.items():
        bystrain = defaultdict(list)
        for row in rows:
            if row["contig"] != "NA":
                bystrain[(row["short"], row["contig"])].append(row)
        cl_id = 0
        for k, lst in bystrain.items():
            lst.sort(key=lambda x: x["gene_rank_on_contig"])
            prev = None
            for row in lst:
                if prev is None or row["gene_rank_on_contig"] - prev["gene_rank_on_contig"] > TANDEM_RANKS:
                    cl_id += 1
                row["pf11001_cluster"] = f"C{cl_id}"
                prev = row
        sizes = Counter(row.get("pf11001_cluster") for row in rows if row.get("pf11001_cluster"))
        for row in rows:
            row["pf11001_cluster_size"] = sizes.get(row.get("pf11001_cluster"), "")
        tsum = []
        strains = defaultdict(list)
        for row in rows:
            strains[row["short"]].append(row)
        for short, lst in sorted(strains.items()):
            cls = {r.get("pf11001_cluster") for r in lst if r.get("pf11001_cluster")}
            tsum.append(dict(short=short, group=lst[0]["group"], n_hits=len(lst),
                             n_placed=sum(r["contig"] != "NA" for r in lst), n_clusters=len(cls),
                             n_hits_in_multi_gene_clusters=sum(1 for r in lst if (r.get("pf11001_cluster_size") or 0) >= 2)))
        out = AFT if s == "afum" else CPT
        write_tsv(out / f"tandem_summary.{s}.tsv", tsum)
        allc = Counter(r["pf11001_cluster_size"] for r in rows if r.get("pf11001_cluster_size"))
        print(f"{s}: genes by cluster size (<= {TANDEM_RANKS} gene ranks): {sorted(allc.items())}", file=sys.stderr)

    # ---------------- neighbourhood: Pfam of genes within +/-5 gene ranks ----------------
    # Neighbour Pfam comes from the neighbour's family representative
    # (our --cut_ga hmmsearch family_pfam.tsv.gz), so it is a family-level label.
    FP = {"full_v072": AF / "analysis/gain_loss_enrichment/full_v072/family_pfam.tsv.gz",
          "immitis_in_posadasii_out_v5": CP / "analysis/gain_loss_enrichment/immitis_in_posadasii_out_v5/family_pfam.tsv.gz"}
    MARK = {"APH": {"PF01636"}, "DUF3723": {"PF12520", "PF29418"}, "HrmA": {"PF28515", "PF28646", "PF28647"},
            "DUF3435": {"PF11917"}}
    nb_rows = []
    for s, run, rk, bys in (("afum", "full_v072", rank_af, bys_af),
                            ("cocci", "immitis_in_posadasii_out_v5", rank_cv5, bys_cv5)):
        fam_pf = {}
        with zread(FP[run]) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["pfam_ids"]:
                    fam_pf[r["accession"]] = set(r["pfam_ids"].split("|"))
        m2r = ctx[run]["m2r"]
        # mark genes per (short, contig) order
        mark_lists = {}
        tot = Counter()
        n_genes = 0
        for (short, c), lst in bys.items():
            marks = []
            for (st, en, pid) in lst:
                pf = fam_pf.get(m2r.get(f"{short}|{pid}"), set())
                marks.append({k for k, v in MARK.items() if pf & v} | ({"PF11001"} if "PF11001" in pf else set()))
            mark_lists[(short, c)] = marks
            # background: fraction of all genes with each mark within +/-5 (excluding self)
            n = len(marks)
            for i in range(n):
                near = set()
                for j in range(max(0, i - TANDEM_RANKS), min(n, i + TANDEM_RANKS + 1)):
                    if j != i:
                        near |= marks[j]
                for k in near:
                    tot[k] += 1
            n_genes += n
        hit_near = Counter()
        n_hit = 0
        for row in gene_rows[s]:
            loc = rk.get((row["short"], row["protein_id"]))
            if loc is None:
                row["neighbour_marks_pm5"] = "NA"
                continue
            c, st, en, i = loc
            marks = mark_lists[(row["short"], c)]
            near = []
            for j in range(max(0, i - TANDEM_RANKS), min(len(marks), i + TANDEM_RANKS + 1)):
                if j != i:
                    near.extend(f"{k}@{j - i:+d}" for k in sorted(marks[j]))
            row["neighbour_marks_pm5"] = ",".join(near) or "none"
            n_hit += 1
            for k in {x.split("@")[0] for x in near}:
                hit_near[k] += 1
        for k in list(MARK) + ["PF11001"]:
            a, n1 = hit_near[k], n_hit
            b, n0 = tot[k], n_genes
            orr, p = fisher_exact([[a, n1 - a], [b, n0 - b]], alternative="greater")
            nb_rows.append(dict(set=s, run_for_positions_and_families=run, mark=k,
                                n_pf11001_genes_placed=n1, n_with_mark_within_5=a,
                                frac=round(a / n1, 4) if n1 else "",
                                n_all_genes=n0, n_all_genes_with_mark_within_5=b,
                                frac_all=round(b / n0, 4), fisher_p_greater=f"{p:.3g}",
                                note="all-genes row includes the PF11001 genes; marks are family-rep Pfam (cut_ga)"))
    write_tsv(AFT / "neighbourhood_marks.tsv", nb_rows)

    write_tsv(AFT / "genes.afum_all_strains.tsv.gz", gene_rows["afum"])
    write_tsv(CPT / "genes.cocci_all_strains.tsv.gz", gene_rows["cocci"])

    # ---------------- reference tables ----------------
    af_rows = []
    for row in gene_rows["afum"]:
        if row["short"] != "Asfu_Af293":
            continue
        acc = row["protein_id"].split("|")[1]
        a = af_ann.get(acc, {})
        r2 = dict(protein_id=row["protein_id"], uniprot=acc, locus_tag=a.get("gene_name", ""),
                  uniprot_description=a.get("description", ""), uniprot_pfam=a.get("pfam_names", ""),
                  chromosome=chrname.get(row["contig"], row["contig"]))
        r2.update({k: v for k, v in row.items() if k not in ("set", "group", "species", "short", "protein_id")})
        if row["contig"] != "NA":
            L = clen.get(("Asfu_Af293", row["contig"]))
            r2["chrom_length"] = L
            r2["dist_to_nearest_chrom_end_bp"] = min(row["start"], L - row["end"]) if L else ""
            r2["starship_element_s7"] = in_element("Asfu_Af293", row["contig"], row["start"], row["end"]) or "no"
        else:
            r2["starship_element_s7"] = "no_coordinates"
        af_rows.append(r2)
    af_rows.sort(key=lambda r: (r.get("chromosome", ""), r.get("start") or 0))
    write_tsv(AFT / "af293_pf11001_proteins.tsv", af_rows)

    cp_rows = []
    for row in gene_rows["cocci"]:
        if row["short"] not in ("RS", "UCSF_Cp_Silveira"):
            continue
        r2 = dict(strain=row["short"], species=row["species"], protein_id=row["protein_id"],
                  gff3_product=cp_prod.get((row["short"], row["protein_id"]), ""))
        r2.update({k: v for k, v in row.items() if k not in ("set", "species", "short", "protein_id")})
        if row["contig"] != "NA":
            L = clen.get((row["short"], row["contig"]))
            r2["contig_length"] = L
            r2["dist_to_nearest_contig_end_bp"] = min(row["start"], L - row["end"]) if L else ""
        cp_rows.append(r2)
    cp_rows.sort(key=lambda r: (r["strain"], r.get("contig", ""), r.get("start") or 0))
    write_tsv(CPT / "cocci_reference_pf11001_proteins.tsv", cp_rows)

    # ---------------- hits vs families ----------------
    for s, runs in SET_RUNS.items():
        for run in runs:
            cx = ctx[run]
            strains_in_run = {k.split("|", 1)[0] for k in cx["m2r"]}
            fams = Counter()
            n_hit = n_mapped = n_not = 0
            for row in gene_rows[s]:
                if row["short"] not in strains_in_run:
                    continue
                n_hit += 1
                fam = cx["m2r"].get(f"{row['short']}|{row['protein_id']}")
                if fam is None:
                    n_not += 1
                else:
                    n_mapped += 1
                    fams[fam] += 1
            # pf11001 families from stage 1
            ftab = (AFT if s == "afum" else CPT) / f"pf11001_families.{run}.tsv"
            stage1 = {}
            with open(ftab, newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    stage1[r["family"]] = r
            bins = Counter(cx["freq"].get(f, {}).get("bin", "NA") for f in fams)
            frows = []
            for f, n in fams.most_common():
                fr = cx["freq"].get(f, {})
                frows.append(dict(run=run, family=f, n_hit_proteins=n, bin=fr.get("bin", "NA"),
                                  strain_count=fr.get("strain_count", "NA"),
                                  direction=cx["gl"].get(f, "NA"),
                                  n_islands=len(cx["fam_isl"].get(f, [])),
                                  rep_has_pf11001_pipeline=stage1.get(f, {}).get("pipeline_counted", "False"),
                                  rep_has_pf11001_cutga=stage1.get(f, {}).get("cutga_hmmsearch_hit", "False")))
            write_tsv((AFT if s == "afum" else CPT) / f"hits_vs_families.{run}.tsv", frows)
            print(f"{run}: strains_in_run={len(strains_in_run)} hit_proteins={n_hit} mapped={n_mapped} "
                  f"not_in_clusters={n_not} families={len(fams)} bins={dict(bins)} "
                  f"families_rep_without_cutga_hit={sum(1 for f in fams if f not in stage1)}",
                  file=sys.stderr)

    # ---------------- Starship / captain family overlap ----------------
    ov = []
    for run, path in (("full_v072", STARSHIP / "full_v072/annotated_families.tsv.gz"),
                      ("immitis_in_posadasii_out_v5", CP / "analysis/starship_island_overlap/immitis_in_posadasii_out_v5/annotated_families.tsv.gz"),
                      ("posadasii_in_immitis_out_v2", CP / "analysis/starship_island_overlap/posadasii_in_immitis_out_v2/annotated_families.tsv.gz")):
        star, capt = set(), set()
        with gzip.open(path, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["starship_elements"]:
                    star.add(r["family"])
                if r["captain"] == "yes":
                    capt.add(r["family"])
        s = "afum" if run == "full_v072" else "cocci"
        ftab = (AFT if s == "afum" else CPT) / f"pf11001_families.{run}.tsv"
        pf = set()
        with open(ftab, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                pf.add(r["family"])
        # families that contain any PF11001 hit protein
        pf_hit = set()
        cx = ctx[run]
        for row in gene_rows[s]:
            f = cx["m2r"].get(f"{row['short']}|{row['protein_id']}")
            if f:
                pf_hit.add(f)
        freq = cx["freq"]
        universes = {
            "all_families": set(freq),
            "noncore(not core/soft_core)": {f for f, r in freq.items() if r["bin"] not in ("core", "soft_core")},
            "shell+cloud": {f for f, r in freq.items() if r["bin"] in ("shell", "cloud")},
        }
        for setname, sset in (("starship_element_family(Af293 S7)", star), ("DUF3435_captain_family", capt)):
            if not sset:
                continue
            for pfname, P in (("rep_PF11001", pf), ("any_member_PF11001_hit", pf_hit)):
                for uname, U in universes.items():
                    A = P & U
                    S = sset & U
                    a = len(A & S)
                    b = len(A - S)
                    c = len(S - A)
                    d = len(U) - a - b - c
                    orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
                    ov.append(dict(run=run, annotation=setname, pf11001_definition=pfname, universe=uname,
                                   n_universe=len(U), n_pf11001=len(A), n_annotated=len(S), n_both=a,
                                   expected_both=round(len(A) * len(S) / len(U), 3),
                                   fisher_p_greater=f"{p:.3g}", both_families=",".join(sorted(A & S))))
    write_tsv(AFT / "starship_captain_overlap.tsv", ov)


if __name__ == "__main__":
    main()
