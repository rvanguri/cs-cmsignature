
import os, sys, numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0,"/gpfs/data/vangurilab/heartmap-cs/scripts")
from foong_spatial_figures import load_geo_visium
from cm_spatial_crossdisease import ensure_symbols, LINEAGES
from foong_regional import assign_zones, nearest_lesion_dist, LESIONAL
from geneset_robustness import score_set
import statsmodels.formula.api as smf
P="/gpfs/data/vangurilab/heartmap-cs"; RAW=P+"/raw"

GENESETS={
 "original4":["GJB7","TNNI3K","MLIP","PANK1"],
 "tnni3k":["TNNI3K"],
 "mlip_pank1":["MLIP","PANK1"],
 "cmstructural":["TNNT2","MYH7","MYH6","TTN","ACTC1","MYL2","TNNI3","TPM1","ACTN2","MYBPC3"],
 "inflammation4":["NLRC4","IL1RAP","BACH2","IL7"],
}
labels=list(GENESETS)
samples=load_geo_visium(RAW+"/foong_spatial")
print(f"loaded {len(samples)} CS sections",flush=True)

dist_rows=[]
for name,a in samples:
    a=ensure_symbols(a); a.var_names=a.var_names.astype(str).str.upper(); a.var_names_make_unique()
    umi=np.asarray(a.X.sum(1)).ravel().astype(float)          # per-spot total UMI (RAW, before normalize)
    if float(np.asarray(a.X.max()))>50:
        sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    zones=assign_zones(a,sc)
    if zones is None: print(f"  {name}: no CM markers -> skip",flush=True); continue
    for lab,genes in GENESETS.items():
        score_set(a,genes,sc,f"score_{lab}")
    if "spatial" not in a.obsm: print(f"  {name}: no spatial -> skip dist",flush=True); continue
    coords=np.asarray(a.obsm["spatial"],float); is_les=np.isin(zones,list(LESIONAL))
    d=nearest_lesion_dist(coords,is_les); pres=zones=="preserved"
    for i in np.where(pres & np.isfinite(d))[0]:
        r={"sample":name,"dist":float(d[i]),"umi":float(umi[i])}
        for lab in labels: r[lab]=float(a.obs[f"score_{lab}"].values[i])
        dist_rows.append(r)
    print(f"  {name}: preserved+finite-dist spots={int((pres&np.isfinite(d)).sum())}",flush=True)

dd=pd.DataFrame(dist_rows)
dd["patient"]=dd["sample"].astype(str).str.replace(r"-\d+$","",regex=True)
dd["dist_z"]=dd.groupby("sample")["dist"].transform(lambda x:(x-x.mean())/(x.std()+1e-9))
dd["umi_z"]=dd.groupby("sample")["umi"].transform(lambda x:(np.log10(x+1)-np.log10(x+1).mean())/(np.log10(x+1).std()+1e-9))
dd=dd.replace([np.inf,-np.inf],np.nan)
print(f"\nspots={len(dd)} patients={dd['patient'].nunique()} "
      f"corr(dist_z,umi_z)={dd[['dist_z','umi_z']].corr().iloc[0,1]:.3f}",flush=True)

rows=[]
for lab in labels:
    sub=dd.dropna(subset=["dist_z","umi_z",lab]).copy()
    rec={"geneset":lab,"n_spots":len(sub),"n_patients":sub["patient"].nunique()}
    try:
        m0=smf.mixedlm(f"{lab} ~ dist_z",sub,groups=sub["patient"].values).fit(reml=True)
        rec["beta_dist_base"]=float(m0.params["dist_z"]); rec["p_dist_base"]=float(m0.pvalues["dist_z"])
    except Exception as e: rec["beta_dist_base"]=np.nan; rec["err0"]=str(e)[:60]
    try:
        m1=smf.mixedlm(f"{lab} ~ dist_z + umi_z",sub,groups=sub["patient"].values).fit(reml=True)
        rec["beta_dist_adj"]=float(m1.params["dist_z"]); rec["p_dist_adj"]=float(m1.pvalues["dist_z"])
        rec["beta_umi_adj"]=float(m1.params["umi_z"]); rec["p_umi_adj"]=float(m1.pvalues["umi_z"])
        ci=m1.conf_int().loc["dist_z"].values; rec["dist_adj_ci_low"]=float(ci[0]); rec["dist_adj_ci_high"]=float(ci[1])
    except Exception as e: rec["beta_dist_adj"]=np.nan; rec["err1"]=str(e)[:60]
    rows.append(rec)
out=pd.DataFrame(rows)
out.to_csv("distance_umi_covariate.tsv",sep="\t",index=False)
dd.to_csv("distance_umi_perspot.tsv.gz",sep="\t",index=False)
print("\n### DISTANCE MODEL: baseline vs UMI-adjusted (beta per SD of distance) ###\n"+out.to_string(index=False),flush=True)
print("\nwrote distance_umi_covariate.tsv, distance_umi_perspot.tsv.gz",flush=True)
