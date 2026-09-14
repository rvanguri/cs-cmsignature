# Figure 1 — legend and sourcing notes

## Legend (as drafted for the research letter)

**Figure 1. A cardiac-sarcoidosis cardiomyocyte contrast replicates in direction but not in
magnitude across platforms.**
**(a)** Spatial versus single-nucleus effect size for the 134 genes testable in both arms.
x, log2 fold change in cardiomyocyte-dominant Visium spots, cardiac sarcoidosis (CS, n = 7
sections) versus comparator cardiomyopathies (n = 6 sections: ARVC ×2, LMNA-associated
dilated cardiomyopathy ×2, Chagas, HCM); y, inverse-variance pooled log2 fold change from the
single-nucleus CS-versus-DCM and CS-versus-ARVC contrasts. Shaded quadrants are direction-
concordant. Open symbols mark genes failing locus quality control (49/134). The 18
labelled genes are direction-concordant and replicate at FDR < 5% in the single-nucleus arm
(15 reduced in CS, 3 increased). Spearman rho = +0.03 over all 134 genes (p = 0.76).
No regression or identity line is drawn: the two axes are not calibrated to each other.
**(b)** The 18 concordant, replicating genes across four contrasts. Circle area encodes
FDR, colour encodes log2 fold change on a symmetric-log diverging scale; colour scales are not
comparable between columns. Columns 1-2 are this study; columns 3-4 are the published
cardiomyocyte candidate lists of Chaffin et al. (dilated and hypertrophic cardiomyopathy versus
non-failing donor hearts, background-corrected values). Open squares mark genes absent from
those published tables — not reported, not tested null. Daggers mark locus-quantification
outliers (2 genes); asterisks mark genes without complete annotation (4 genes).

## Where the numbers come from

| element | source file |
|---|---|
| panel a points, x | `res/q1_raw_T2000_d25.tsv` (`logFC`), region T = 2000 UMI, detection 0.25 |
| panel a points, y | `res/q1_core_snrna_pergene.tsv` (`lfc_pool`, inverse-variance pooled) |
| QC flags / open symbols | `res/q2_locus_qc.tsv`; as plotted, the `flag_*` columns of `figure1_panelA_values.tsv` |
| panel b columns 1-2 | as above |
| panel b columns 3-4 | `2021-02-03277C-ST6.DCMvsNF.xlsx`, `...ST7.HCMvsNF.xlsx`, cardiomyocyte rows |
| section counts | `res/meta_T2000.tsv` (`evaluable`) |
| exact plotted values | `res/figure1_panelA_values.tsv`, `res/figure1_panelB_values.tsv` |

## Departures from the handoff, and why

(The handoff document this refers to is not in the working tree; see git tag
`pre-slim-2026-09-14`.)

1. **Spearman rho.** The handoff instructs annotating **+0.06** as the full-set value. Over all
   134 genes the value is **+0.027** (p = 0.76); **+0.062** (p = 0.53) is the down-in-CS
   subset (n = 108). The panel annotates the full-set value the handoff asked for, with the
   correct number. Within the 26 up-in-CS genes the relationship is inverse
   (rho = -0.51, p = 0.008); this is reported here rather than on the panel.
2. **QC flag description.** The handoff describes `flag_length_deviation` as an effect-size
   threshold. It is neither that nor a plain gene-length flag: per contributing dataset, mean
   expression is regressed on log gene span and a gene is flagged when its residual from that
   length-matched expectation exceeds a fixed margin in any dataset. The legend therefore calls
   it a locus-quantification outlier. The handoff's companion claim — gene length is unrelated
   to effect size — does hold.
3. **The published columns are mostly empty, and that is a property of the source.** Both
   Chaffin supplementary tables are filtered per-cell-type candidate lists, not complete
   transcriptomes: they contain non-significant rows but omit common housekeeping genes, and
   only 2/18 (DCM) and 6/18 (HCM) of the concordant genes appear at all. Absence is
   marked with an open square and means not reported. Only genes present in both published
   columns support any cross-study magnitude statement.
4. **Axis scales.** The spatial axis is linear (the data occupy enough of its range that no
   break or symlog is needed). The single-nucleus axis is symmetric-log with explicit ticks
   because one gene sits far above the rest; nothing is clipped.
5. **Leave-one-out counts.** The handoff gives a single figure for the up-in-CS count when a
   paired-section disease is dropped; the section-level and disease-level jackknife tables give
   different values depending on which section is removed. Cite
   `res/q1_leave_one_comparator_out.tsv` and `res/q1_leave_one_disease_out.tsv` directly.
6. **Alternative panel a.** `res/figure1_panelA_alt_exemplars.png` shows per-section log2 CPM
   for four exemplar genes (CS versus comparators) as a substitute if a reviewer wants
   per-sample values rather than a cross-platform scatter. It is a drop-in for panel a only.
