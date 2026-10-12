#!/usr/bin/env python3
"""ESC HF revision, section 0: universe, annotated hits and comparison with the old primary.

Run from the repo root after a0_primary_refit.R.
  results/eschf/a0_universe.tsv           new (all CM-dominant spots) vs old (>=0.5 mm) tested universe
  results/eschf/a0_hits_annotated.tsv     FDR<5% hits of the new section-level fit + existing locus-QC flags
                                          (same join and convention as 00_inputs.build_hits)
  results/eschf/a0_vs_old_primary.tsv     overlap, sign retention, Pearson r of logFC (shared genes)
  results/eschf/primary_v2/q1_hits_annotated.tsv            copy for the cascade scripts
  results/eschf/primary_v2/geneset_primary_program_qcpass_down.txt
                                          new program score genes: QC-passing lower-in-CS hits
"""
import pandas as pd

OUT, V2, RES, FDR = "results/eschf", "results/eschf/primary_v2", "results/spatial_first", 0.05
new = pd.read_csv(f"{V2}/q1_raw_allcm.tsv", sep="\t").set_index("gene")
old = pd.read_csv(f"{RES}/q1_raw_T2000_d25.tsv", sep="\t").set_index("gene")
qc = pd.read_csv(f"{RES}/q2_locus_qc.tsv", sep="\t")
dn = pd.read_csv("data/eschf/primary_v2/detection_fraction_allcm_T2000.tsv.gz", sep="\t", index_col=0)
do = pd.read_csv("data/spatial_first/detection_fraction_T2000.tsv.gz", sep="\t", index_col=0)

# ---- universe
U = pd.DataFrame(index=sorted(set(new.index) | set(old.index)))
U["in_new_universe"] = U.index.isin(new.index); U["in_old_universe"] = U.index.isin(old.index)
U["det_CS_new"] = dn.CS.reindex(U.index); U["det_COMPARATOR_new"] = dn.COMPARATOR.reindex(U.index)
U["det_CS_old"] = do.CS.reindex(U.index); U["det_COMPARATOR_old"] = do.COMPARATOR.reindex(U.index)
U.rename_axis("gene").reset_index().to_csv(f"{OUT}/a0_universe.tsv", sep="\t", index=False)

# ---- annotated hits (00_inputs.build_hits convention)
cols = ["gene", "qc_pass", "flag_length_deviation", "flag_segdup", "flag_missing_annotation", "flag_readthrough",
        "flag_source_specific_absence", "max_abs_length_dev_log2", "span_kb", "biotype"]
h = new[new["adj.P.Val"] < FDR].reset_index().merge(qc[cols], on="gene", how="left")
h["in_qc_universe"] = h.qc_pass.notna()
h["qc_pass"] = h.qc_pass.astype("boolean").fillna(False).astype(bool)
h["direction"] = (h.logFC > 0).map({True: "up", False: "down"})
h = h.rename(columns={"flag_length_deviation": "flag_locus_quant_outlier"})
h["in_old_universe"] = h.gene.isin(old.index)
h["old_primary_hit"] = h.gene.map(old["adj.P.Val"] < FDR).fillna(False).astype(bool)
h["old_logFC"] = h.gene.map(old.logFC)
h = h.sort_values("logFC")
h.to_csv(f"{OUT}/a0_hits_annotated.tsv", sep="\t", index=False)
h.drop(columns=["in_old_universe", "old_primary_hit", "old_logFC"]).to_csv(f"{V2}/q1_hits_annotated.tsv", sep="\t", index=False)
prog = sorted(h.gene[h.qc_pass & (h.direction == "down")])
open(f"{V2}/geneset_primary_program_qcpass_down.txt", "w").write("\n".join(prog) + "\n")

# ---- vs old primary
qcp = set(qc.gene[qc.qc_pass.astype(str) == "True"])
rows = []
for lab, keep in (("all", lambda g: True), ("qc_pass", lambda g: g in qcp)):
    hn = {g for g in new.index[new["adj.P.Val"] < FDR] if keep(g)}
    ho = {g for g in old.index[old["adj.P.Val"] < FDR] if keep(g)}
    shared = [g for g in new.index.intersection(old.index) if keep(g)]
    both = sorted(hn & ho)
    ho_in = [g for g in ho if g in new.index]
    rows.append(dict(gene_set=lab, n_hits_new=len(hn), n_hits_new_down=sum(new.loc[g, "logFC"] < 0 for g in hn),
                     n_hits_new_up=sum(new.loc[g, "logFC"] > 0 for g in hn), n_hits_old=len(ho),
                     n_hits_old_down=sum(old.loc[g, "logFC"] < 0 for g in ho), n_hits_old_up=sum(old.loc[g, "logFC"] > 0 for g in ho),
                     n_overlap=len(both), n_new_only=len(hn - ho), n_old_only=len(ho - hn),
                     n_old_hits_in_new_universe=len(ho_in),
                     n_old_hits_not_in_new_universe=len(ho) - len(ho_in),
                     old_hits_same_sign_in_new=int(sum((new.loc[g, "logFC"] > 0) == (old.loc[g, "logFC"] > 0) for g in ho_in)),
                     n_overlap_same_sign=int(sum((new.loc[g, "logFC"] > 0) == (old.loc[g, "logFC"] > 0) for g in both)),
                     n_shared_genes=len(shared),
                     pearson_r_logFC_shared=float(new.loc[shared, "logFC"].corr(old.loc[shared, "logFC"])),
                     pearson_r_logFC_old_hits=float(new.loc[ho_in, "logFC"].corr(old.loc[ho_in, "logFC"]))))
C = pd.DataFrame(rows)
C.to_csv(f"{OUT}/a0_vs_old_primary.tsv", sep="\t", index=False)
print(C.T.to_string()); print(len(prog), "program genes;", U.in_new_universe.sum(), U.in_old_universe.sum(),
                              (U.in_new_universe & U.in_old_universe).sum())
print(h.groupby(["direction", "qc_pass"]).size().to_string())
