#!/usr/bin/env python3
"""A11 helper: which S1 CS subjects have snRNA samples in the raw GEO download on BigPurple
(<RAW>/neyazi; listing captured to data/external/neyazi_geo_snrna_sample_listing.txt).
Sample names are GSM<id>_<subject>_<region>_<rep>."""
import re, pandas as pd
L = [l.strip() for l in open("data/external/neyazi_geo_snrna_sample_listing.txt") if l.strip()]
subj = pd.Series([re.match(r"GSM\d+_(.+)_(LV|RV|SP)_\d+$", x).group(1) for x in L]).value_counts()
s1 = pd.read_excel("data/external/neyazi_table_s1.xlsx", sheet_name="Sarcoid"); s1["Subject"] = s1.Subject.astype(str)
rows = []
for s in s1.Subject:
    n = int(subj.get(s, 0))
    rows.append(dict(Subject=s, n_geo_snrna_samples=n,
                     reason="" if n else f"no snRNA sample for this subject in the GEO download ({len(L)} samples, none named {s}); "
                                         "Table S2 lists 0 snRNA samples by region but 2 by histologic class"))
pd.DataFrame(rows).to_csv("results/eschf/a11_raw_donor_check.tsv", sep="\t", index=False)
print(len(L), subj.sum(), (pd.DataFrame(rows).n_geo_snrna_samples == 0).sum())
