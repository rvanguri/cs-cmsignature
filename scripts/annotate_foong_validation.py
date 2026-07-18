#!/usr/bin/env python3
"""Cross-cohort CM validation: combine ALL our CS-vs-mimic contrasts with independent Foong replication.

Foong et al. profiled cardiac sarcoidosis vs HCM in cardiomyocyte-dominant spatial (Visium) spots: an
INDEPENDENT CS cohort on an INDEPENDENT platform. Our earlier version annotated only CS_vs_HCM, which
(a) is our weakest contrast for the marquee genes and (b) mislabels genes like GJB7 as 'ns_ours' just
because they miss significance in that one comparison. This version scores each gene on COMBINED
evidence across every CM contrast we have, plus Foong.

Per gene we report:
  our_max_lfc            signed logFC of the contrast where the gene is most significant
  our_best_padj          that contrast's padj
  n_contrasts_sig        # of our CM contrasts where sig & |lfc|>thr (concordant direction)
  our_dir                +1/-1 consensus direction across the sig contrasts (0 if mixed)
  foong_log2FC, foong_padj
  foong_status           replicated / discordant / mixed_direction / unconfirmed_in_foong / ns_all_contrasts
  evidence_score         combined rank score (see below), higher = more robust CS-CM gene

evidence_score = n_contrasts_sig
                 + (1 if replicated in Foong else 0)
                 + expr_bonus                      # clip((mean AveExpr - expr_floor)/expr_scale, -1, 1)
                 - (2 if discordant else 0)
We reward EXPRESSION (expr_bonus), NOT raw |lfc|: a large pseudobulk logFC is usually an on/off
artifact of a near-absent gene. Genes significant+concordant across contrasts AND replicated in Foong
float to the top; single-contrast, discordant, or low-expression genes sink.

CAVEAT: the CS-vs-mimic contrasts feeding this table are cohort-confounded: CS is present
in only one study, so disease is aliased with study and these contrasts are HYPOTHESIS-GENERATING, not
batch-clean. The genuinely independent column is the Foong spatial replication, not the contrasts; rank
here is a candidate-prioritization aid, not specificity evidence. The Foong file is significant-only
(all padj<=0.05), so absence = 'not significant / undetected in Foong', not a refutation.
"""
from __future__ import annotations
import argparse, glob, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, require  # noqa: E402
LOG = log("foong_val")


def load_de(path):
    d = pd.read_csv(path, sep="\t", index_col=0)
    lfc = "logFC" if "logFC" in d.columns else "lfc"
    padj = "adj.P.Val" if "adj.P.Val" in d.columns else "padj"
    cols = {lfc: "lfc", padj: "padj"}
    if "AveExpr" in d.columns:
        cols["AveExpr"] = "aveexpr"
    out = d[list(cols)].rename(columns=cols)
    if "aveexpr" not in out.columns:
        out["aveexpr"] = np.nan
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--de_dir", required=True, help="dir with DE_CS_vs_*.tsv (e.g. de/de_cm)")
    ap.add_argument("--foong", required=True, help="Foong significant-DEG csv (gene,log2FoldChange,padj)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--lfc", type=float, default=1.0)
    ap.add_argument("--expr_floor", type=float, default=2.0,
                    help="AveExpr (log2-CPM) at which the expression bonus is neutral; below penalizes")
    ap.add_argument("--expr_scale", type=float, default=4.0)
    args = ap.parse_args()
    require(args.foong, "Foong csv")

    files = sorted(glob.glob(os.path.join(args.de_dir, "DE_CS_vs_*.tsv")))
    files = [f for f in files if "foong" not in os.path.basename(f).lower()]
    if not files:
        sys.exit("no DE_CS_vs_*.tsv in %s" % args.de_dir)
    contrasts = {os.path.basename(f).replace("DE_CS_vs_", "").replace(".tsv", ""): load_de(f) for f in files}
    LOG.info("contrasts: %s", ", ".join(contrasts))

    genes = sorted(set().union(*[d.index for d in contrasts.values()]))
    lfcs = pd.DataFrame({c: d["lfc"].reindex(genes) for c, d in contrasts.items()}, index=genes)
    padjs = pd.DataFrame({c: d["padj"].reindex(genes) for c, d in contrasts.items()}, index=genes)
    mean_aveexpr = pd.DataFrame({c: d["aveexpr"].reindex(genes)
                                 for c, d in contrasts.items()}, index=genes).mean(axis=1)

    sig = (padjs < args.fdr) & (lfcs.abs() > args.lfc)
    n_sig = sig.sum(axis=1)
    # consensus direction among the significant contrasts
    sign_in_sig = np.sign(lfcs).where(sig)
    up = (sign_in_sig > 0).sum(axis=1); dn = (sign_in_sig < 0).sum(axis=1)
    our_dir = np.where(up > dn, 1, np.where(dn > up, -1, 0))
    mean_lfc = lfcs.where(sig).mean(axis=1)
    # the single best (most significant) contrast per gene
    best_c = padjs.idxmin(axis=1)
    best_padj = padjs.min(axis=1)
    best_lfc = pd.Series([lfcs.loc[g, best_c[g]] for g in genes], index=genes)

    fo = pd.read_csv(args.foong).set_index("gene")
    fo_lfc = fo["log2FoldChange"].reindex(genes)
    fo_padj = fo["padj"].reindex(genes)
    in_foong = fo_lfc.notna()

    # status uses the CONSENSUS direction of our sig contrasts vs Foong. Genes significant but with
    # NO consensus direction across contrasts (our_dir==0) are 'mixed_direction', NOT 'discordant'.
    our_dir_s = pd.Series(our_dir, index=genes)
    mixed = (our_dir_s == 0) & (n_sig > 0)
    same = np.sign(our_dir_s) == np.sign(fo_lfc)
    status = np.where(n_sig == 0, "ns_all_contrasts",
             np.where(mixed, "mixed_direction",
             np.where(in_foong & same, "replicated",
             np.where(in_foong & ~same, "discordant", "unconfirmed_in_foong"))))

    # expression term: reward well-expressed genes, penalize low-expression fold-change spikes.
    # (We deliberately do NOT reward raw |lfc|: a pseudobulk logFC of +8 is usually an on/off
    # artifact of a near-absent gene, which is what surfaced CHRM5/GRIK1/ALK to the top before.)
    expr_bonus = np.clip((mean_aveexpr - args.expr_floor) / args.expr_scale, -1.0, 1.0).fillna(0)

    ev = (n_sig
          + (status == "replicated").astype(int)
          + expr_bonus
          - 2 * (status == "discordant").astype(int))

    out = pd.DataFrame({
        "n_contrasts_sig": n_sig,
        "our_dir": our_dir,
        "our_mean_lfc": mean_lfc.round(3),
        "mean_aveexpr": mean_aveexpr.round(2),
        "well_expressed": mean_aveexpr >= args.expr_floor,
        "best_contrast": best_c,
        "best_lfc": best_lfc.round(3),
        "best_padj": best_padj,
        "foong_log2FC": fo_lfc.round(3),
        "foong_padj": fo_padj,
        "foong_status": status,
        "evidence_score": ev.round(3),
    }, index=genes)
    out.index.name = "gene"
    for c in contrasts:                       # keep per-contrast logFC for transparency
        out[f"lfc_{c}"] = lfcs[c].round(3)
    out = out.sort_values(["evidence_score", "best_padj"], ascending=[False, True])
    out.to_csv(args.out, sep="\t")

    s = pd.Series(status, index=genes)
    LOG.info("genes: %d | contrasts: %d", len(genes), len(contrasts))
    for k in ["replicated", "discordant", "mixed_direction", "unconfirmed_in_foong", "ns_all_contrasts"]:
        LOG.info("   %-22s %d", k, int((s == k).sum()))
    hi = out[(out.n_contrasts_sig >= 2) & (out.foong_status == "replicated") & (out.our_dir > 0)]
    LOG.info("robust up-in-CS (>=2 contrasts sig AND Foong-replicated): %d genes", len(hi))
    LOG.info("top 15 by evidence_score:\n%s",
             out.head(15)[["n_contrasts_sig", "our_mean_lfc", "foong_log2FC",
                           "foong_status", "evidence_score"]].to_string())
    LOG.info("wrote %s", args.out)


if __name__ == "__main__":
    main()
