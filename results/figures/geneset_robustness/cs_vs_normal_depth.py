
import os, sys, glob, numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy import sparse
sys.path.insert(0,"/gpfs/data/vangurilab/heartmap-cs/scripts")
from foong_spatial_figures import load_geo_visium
from cm_spatial_crossdisease import ensure_symbols, LINEAGES
P="/gpfs/data/vangurilab/heartmap-cs"; RAW=P+"/raw"
rng=np.random.default_rng(0)
SIG=["GJB7","TNNI3K","MLIP","PANK1"]; COMPOSITE=SIG[:]
Tgrid=[1000,1500,2000,3000]

def raw_counts_symbols(name,a,is_kuppe):
    """Return AnnData of RAW integer counts with SYMBOL var_names."""
    if is_kuppe:
        a=ensure_symbols(a)                       # maps a.var_names Ensembl->symbol, sets a.var['_ensembl']
        rc=a.raw.to_adata()                       # raw counts, Ensembl-indexed
        if "_ensembl" in a.var:
            e2s=dict(zip(a.var["_ensembl"].astype(str), a.var_names.astype(str)))
            rc.var_names=[e2s.get(str(e),str(e)) for e in rc.var_names.astype(str)]
    else:
        rc=a                                       # Foong: X is raw counts, symbols already
    rc.var_names=pd.Index(rc.var_names.astype(str)).str.upper(); rc.var_names_make_unique()
    rc.X=sparse.csr_matrix(rc.X)
    return rc

def cmdom_raw(rc):
    """normalize+log a COPY to score lineages; return raw-count subset of CM-dominant spots."""
    nn=rc.copy(); sc.pp.normalize_total(nn,target_sum=1e4); sc.pp.log1p(nn)
    lin={}
    for L,gs in LINEAGES.items():
        gg=[g for g in gs if g in nn.var_names]
        if gg:
            sc.tl.score_genes(nn,gg,score_name=f"_l_{L}",use_raw=False); lin[L]=nn.obs[f"_l_{L}"].values
    if "Cardiomyocyte" not in lin: return None
    dom=pd.DataFrame(lin).idxmax(axis=1).values=="Cardiomyocyte"
    return rc[dom].copy()

def build_cohort(objs,is_kuppe,label):
    parts=[]
    for name,a in objs:
        rc=raw_counts_symbols(name,a,is_kuppe)
        sub=cmdom_raw(rc)
        if sub is None or sub.n_obs==0: print(f"  {label}/{name}: no CM-dom -> skip",flush=True); continue
        sub.obs["section"]=name
        parts.append(sub); print(f"  {label}/{name}: {sub.n_obs} CM-dom spots",flush=True)
    R=ad.concat(parts,join="outer",fill_value=0)
    R.X=sparse.csr_matrix(R.X); R.X.data=np.rint(R.X.data)   # ensure int counts
    return R

def gene_col(R,g):
    if g not in R.var_names: return None
    j=R.var_names.get_loc(g); v=R[:,j].X
    return np.asarray(v.todense()).ravel() if sparse.issparse(v) else np.asarray(v).ravel()

def norm_log(R):
    A=R.copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A); return A

def thin(R,T):
    tot=np.asarray(R.X.sum(1)).ravel(); keep=tot>=T
    sub=R[keep].copy(); sub.X=sparse.csr_matrix(sub.X)
    totk=tot[keep]; p=np.clip(np.repeat(T/totk,np.diff(sub.X.indptr)),0,1)
    sub.X.data=rng.binomial(sub.X.data.astype(int),p).astype(np.float32); sub.X.eliminate_zeros()
    return sub, float(keep.mean()), int(keep.sum())

def analyze(R,label):
    pg=[]; comp=[]
    tot=np.asarray(R.X.sum(1)).ravel()
    # NATIVE
    A=norm_log(R)
    for g in SIG:
        c=gene_col(R,g)
        if c is None: pg.append(dict(cohort=label,T="native",gene=g,present=False,n_pos=0,n=R.n_obs,detect=0.0,mean_lognorm=0.0)); continue
        e=gene_col(A,g)
        pg.append(dict(cohort=label,T="native",gene=g,present=True,n_pos=int((c>0).sum()),n=R.n_obs,
                       detect=float((c>0).mean()),mean_lognorm=float(e.mean())))
    sc.tl.score_genes(A,[g for g in COMPOSITE if g in A.var_names],score_name="orig4",use_raw=False,random_state=0)
    comp.append(dict(cohort=label,T="native",n=R.n_obs,med_orig4=float(np.median(A.obs["orig4"])),
                     med_umi=float(np.median(tot))))
    # COMMON DEPTH
    for T in Tgrid:
        sub,frac,nk=thin(R,T)
        A2=norm_log(sub)
        for g in SIG:
            c=gene_col(sub,g)
            if c is None: pg.append(dict(cohort=label,T=T,gene=g,present=False,n_pos=0,n=nk,detect=0.0,mean_lognorm=0.0,retain=frac)); continue
            e=gene_col(A2,g)
            pg.append(dict(cohort=label,T=T,gene=g,present=True,n_pos=int((c>0).sum()),n=nk,
                           detect=float((c>0).mean()),mean_lognorm=float(e.mean()),retain=frac))
        sc.tl.score_genes(A2,[g for g in COMPOSITE if g in A2.var_names],score_name="orig4",use_raw=False,random_state=0)
        comp.append(dict(cohort=label,T=T,n=nk,med_orig4=float(np.median(A2.obs["orig4"])),
                         med_umi=float(np.median(np.asarray(sub.X.sum(1)).ravel())),retain=frac))
    return pg,comp

print("=== load CS (Foong) ===",flush=True)
cs=build_cohort(load_geo_visium(RAW+"/foong_spatial"),False,"CS")
print("=== load control (Kuppe) ===",flush=True)
kf=sorted(glob.glob(RAW+"/kuppe_control/*.h5ad"))
ct=build_cohort([(os.path.basename(f).replace(".h5ad",""),ad.read_h5ad(f)) for f in kf],True,"control")
print(f"CS CM-dom={cs.n_obs}  control CM-dom={ct.n_obs}",flush=True)

pg=[]; comp=[]
for R,lab in [(cs,"CS"),(ct,"control")]:
    a,b=analyze(R,lab); pg+=a; comp+=b
pg=pd.DataFrame(pg); comp=pd.DataFrame(comp)
pg.to_csv("cs_vs_normal_pergene.tsv",sep="\t",index=False)
comp.to_csv("cs_vs_normal_composite.tsv",sep="\t",index=False)
print("\n### PER-GENE ###\n"+pg.to_string(index=False),flush=True)
print("\n### COMPOSITE ###\n"+comp.to_string(index=False),flush=True)
print("\nwrote cs_vs_normal_pergene.tsv, cs_vs_normal_composite.tsv",flush=True)
