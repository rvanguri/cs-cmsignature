#!/usr/bin/env python3
"""Where the CS-lower genes sit relative to non-failing myocardium, from published data.

No probe-based non-failing donor sections exist in this cohort, and probe-based and
polyA data are not comparable for differential expression, so the comparison against
non-failing myocardium has to come from published work: Chaffin et al., Nature 2022,
supplementary tables ST6 (DCM vs non-failing) and ST7 (HCM vs non-failing),
cardiomyocyte rows, CellBender-corrected effect sizes.

Those tables list only genes called differential, not all genes tested. Testability is
therefore one-sided: a CS-lower gene absent from the published table may be unchanged in
DCM/HCM or may simply not have been reported. The reference rate for the binomial test
is accordingly the fraction of *all* published cardiomyocyte differential genes that are
lower in disease, not 0.5 -- the question is whether the CS-lower set is enriched for
disease-lower genes beyond the base rate of the published differential set.

Inputs
  data/published/chaffin2022_cardiomyocyte_{DCM,HCM}vsNF.tsv   distilled, committed
  results/spatial_first/q1_raw_T2000_d25.tsv                   defines the CS-lower set
Rebuild the distilled tables from the original supplements with
  --from-supplements <dir containing 2021-02-03277C-ST6.DCMvsNF.xlsx and ...ST7...xlsx>

Outputs
  results/spatial_first/q1_published_nf_pergene.tsv   per-gene published effects
  results/spatial_first/q1_published_nf_summary.tsv   counts and binomial tests
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

DATA = "data/published"
RES = "results/spatial_first"
SUPP = {"DCM": "2021-02-03277C-ST6.DCMvsNF.xlsx", "HCM": "2021-02-03277C-ST7.HCMvsNF.xlsx"}


def distill(src_dir):
    """Reduce the published supplements to the cardiomyocyte rows and columns used."""
    os.makedirs(DATA, exist_ok=True)
    for tag, fn in SUPP.items():
        raw = pd.read_excel(os.path.join(src_dir, fn), sheet_name=0)
        cols = {c.replace("\n", " ").replace("_x000D_", "").strip(): c for c in raw.columns}
        cm = raw[raw["Cell Type"] == "Cardiomyocyte"]
        out = pd.DataFrame({
            "gene": cm["Gene"].values,
            "ens_id": cm["Ensembl ID"].values,
            "logFC": cm[cols["CellBender: logFC"]].values,
            "adj_p": cm[cols["CellBender: Adjusted P-Value"]].values,
            "bg_contamination_flag": cm[cols["Background Contamination Flag"]].fillna(0).values,
        })
        out.to_csv(f"{DATA}/chaffin2022_cardiomyocyte_{tag}vsNF.tsv", sep="\t", index=False)
        print(f"[distill] {tag}: {len(out)} cardiomyocyte differential genes")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-supplements", metavar="DIR",
                    help="rebuild the distilled tables from the original xlsx supplements")
    args = ap.parse_args()
    if args.from_supplements:
        distill(args.from_supplements)

    q1 = pd.read_csv(f"{RES}/q1_raw_T2000_d25.tsv", sep="\t").set_index("gene")
    lower = q1.index[(q1["adj.P.Val"] < 0.05) & (q1.logFC < 0)]

    pub, rows, per = {}, [], pd.DataFrame({"spatial_logFC": q1.loc[lower, "logFC"]})
    for tag in ("DCM", "HCM"):
        t = pd.read_csv(f"{DATA}/chaffin2022_cardiomyocyte_{tag}vsNF.tsv", sep="\t")
        # Background-contamination-flagged rows are excluded: the published flag marks
        # effects the authors attribute to ambient RNA rather than to the cell type.
        t = t[t.bg_contamination_flag == 0].set_index("gene")
        pub[tag] = t
        per[f"lfc_{tag}NF"] = t.logFC.reindex(per.index)
        per[f"fdr_{tag}NF"] = t.adj_p.reindex(per.index)

        g = per.index[per[f"lfc_{tag}NF"].notna()]
        k = int((per.loc[g, f"lfc_{tag}NF"] < 0).sum())
        base = float((t.logFC < 0).mean())
        rows.append(dict(
            comparison=f"{tag} vs non-failing", n_cs_lower=len(lower), n_testable=len(g),
            n_lower_in_disease=k, frac_lower_in_disease=round(k / len(g), 3),
            n_published_deg=len(t), base_frac_lower=round(base, 3),
            binom_p=stats.binomtest(k, len(g), base, alternative="greater").pvalue))

    per.to_csv(f"{RES}/q1_published_nf_pergene.tsv", sep="\t")
    S = pd.DataFrame(rows)
    S.to_csv(f"{RES}/q1_published_nf_summary.tsv", sep="\t", index=False)
    print(S.to_string(index=False))
    for gene in ("IDH2", "ACADVL"):
        d = {t: (None if gene not in pub[t].index else round(pub[t].loc[gene, "logFC"], 2))
             for t in pub}
        print(f"[{gene}] spatial logFC {q1.loc[gene, 'logFC']:+.2f} | published {d}")


if __name__ == "__main__":
    main()
