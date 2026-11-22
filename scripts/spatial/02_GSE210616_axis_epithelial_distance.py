#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.stats import spearmanr,rankdata,wilcoxon

ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas'); OUT=ROOT/'spatial_validation/results'
ST=Path('/mnt/stk35data/BreastCancer_Census_v1/Stage11E_SpatialTranscriptomics_v0_1/results')
A=ad.read_h5ad(ST/'GSE210616_TNBC_spatial_affinity_v0_1.h5ad')
D=np.load(ST/'GSE210616_epithelial_domain_distance_v0_1.npz')['distance_spot_units'].astype(float)
obs=A.obs.copy(); qc=obs.qc_main.to_numpy(bool); sym=np.asarray(A.var.gene_symbol.astype(str)); symmap={}
for j,s in enumerate(sym): symmap.setdefault(s.upper(),[]).append(j)
def logcpm(g):
 jj=symmap[g.upper()]; X=A.X[:,jj]
 if not sparse.issparse(X): X=sparse.csr_matrix(X)
 raw=np.asarray(X.sum(axis=1)).ravel().astype(float); lib=obs.total_counts.to_numpy(float)
 return np.log2(raw/np.maximum(lib,1)*1e6+1)
fap=logcpm('FAP'); dpp4=logcpm('DPP4')
fap_z=np.full(A.n_obs,np.nan); dpp4_z=np.full(A.n_obs,np.nan); fib_free=np.full(A.n_obs,np.nan)
# Frozen TNBC object uses 47 available fibroblast markers, including FAP. Remove FAP algebraically.
for gsm,idx0 in obs.groupby('gsm',observed=True).indices.items():
 idx=np.asarray(idx0,int); q=qc[idx]
 for e,dest in [(fap,fap_z),(dpp4,dpp4_z)]:
  mu=np.nanmean(e[idx][q]); sd=np.nanstd(e[idx][q]); sd=sd if sd>=1e-6 else 1.0
  dest[idx]=np.clip((e[idx]-mu)/sd,-4,4)
 aff=obs.iloc[idx].affinity_Fibroblast.to_numpy(float)
 fib_free[idx]=(47.0*aff-fap_z[idx])/46.0
axis=fap_z-dpp4_z

def partial_spearman(x,y,C):
 x=np.asarray(x,float); y=np.asarray(y,float); C=np.asarray(C,float); ok=np.isfinite(x)&np.isfinite(y)&np.all(np.isfinite(C),axis=1)
 if ok.sum()<20:return np.nan
 xr=rankdata(x[ok]); yr=rankdata(y[ok]); Cr=np.column_stack([rankdata(C[ok,j]) for j in range(C.shape[1])]); X=np.column_stack([np.ones(ok.sum()),Cr])
 rx=xr-X@np.linalg.lstsq(X,xr,rcond=None)[0]; ry=yr-X@np.linalg.lstsq(X,yr,rcond=None)[0]
 return float(np.corrcoef(rx,ry)[0,1]) if np.std(rx)>1e-12 and np.std(ry)>1e-12 else np.nan
rows=[]
for gsm,idx0 in obs.groupby('gsm',observed=True).indices.items():
 idx=np.asarray(idx0,int); iq=idx[qc[idx] & np.isfinite(D[idx])]
 if len(iq)<100: continue
 thr=np.nanquantile(fib_free[iq],.70); fiq=iq[fib_free[iq]>=thr]
 if len(fiq)<50: continue
 r=obs.iloc[idx[0]]
 rho=float(spearmanr(axis[fiq],D[fiq]).statistic)
 C=np.column_stack([fib_free[fiq],np.log1p(obs.iloc[fiq].total_counts.to_numpy(float)),obs.iloc[fiq].n_genes_by_counts.to_numpy(float)])
 prho=partial_spearman(axis[fiq],D[fiq],C)
 lo=np.nanquantile(axis[fiq],1/3); hi=np.nanquantile(axis[fiq],2/3)
 dlo=np.nanmedian(D[fiq][axis[fiq]<=lo]); dhi=np.nanmedian(D[fiq][axis[fiq]>=hi])
 rows.append({'gsm':gsm,'patient_id':str(r.patient_id),'section_id':str(r.section_id),'n_fibroblast_rich_spots':len(fiq),'spearman_axis_vs_distance':rho,'partial_spearman_axis_vs_distance':prho,'median_distance_low_axis':dlo,'median_distance_high_axis':dhi,'high_minus_low_distance':dhi-dlo})
S=pd.DataFrame(rows); S.to_csv(OUT/'GSE210616_axis_epithelial_distance_section.tsv',sep='\t',index=False)
P=S.groupby('patient_id').agg(n_sections=('gsm','nunique'),spearman_axis_vs_distance=('spearman_axis_vs_distance','median'),partial_spearman_axis_vs_distance=('partial_spearman_axis_vs_distance','median'),high_minus_low_distance=('high_minus_low_distance','median')).reset_index()
P.to_csv(OUT/'GSE210616_axis_epithelial_distance_patient.tsv',sep='\t',index=False)
summary=[]
for metric in ['spearman_axis_vs_distance','partial_spearman_axis_vs_distance','high_minus_low_distance']:
 v=P[metric].dropna().to_numpy(float); p=float(wilcoxon(v,zero_method='wilcox',alternative='two-sided',method='auto').pvalue)
 summary.append({'metric':metric,'n_patients':len(v),'median':float(np.median(v)),'q25':float(np.quantile(v,.25)),'q75':float(np.quantile(v,.75)),'fraction_negative':float(np.mean(v<0)),'wilcoxon_p':p})
T=pd.DataFrame(summary); T.to_csv(OUT/'GSE210616_axis_epithelial_distance_summary.tsv',sep='\t',index=False)
print(T.to_string(index=False))