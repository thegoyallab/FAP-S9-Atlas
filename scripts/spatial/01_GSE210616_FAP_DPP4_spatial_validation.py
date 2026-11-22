#!/usr/bin/env python3
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
from scipy.stats import spearmanr, rankdata, wilcoxon
from statsmodels.stats.multitest import multipletests

ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas')
SP=ROOT/'spatial_validation'; OUT=SP/'results'
ST=Path('/mnt/stk35data/BreastCancer_Census_v1/Stage11E_SpatialTranscriptomics_v0_1')
H5=ST/'results/GSE210616_TNBC_spatial_affinity_v0_1.h5ad'
SIG=ST/'results/HBCA_broad_spatial_signatures_v0_1.tsv'
MOD=ROOT/'single_cell/config/fibroblast_modules.tsv'

A=ad.read_h5ad(H5)
obs=A.obs.copy()
qc=obs.qc_main.to_numpy(bool)
symbol=np.asarray(A.var['gene_symbol'].astype(str))
symcols={}
for j,s in enumerate(symbol): symcols.setdefault(s.upper(),[]).append(j)

def expr_symbol(g):
    jj=symcols.get(g.upper())
    if not jj: return None
    X=A.X[:,jj]
    if not sparse.issparse(X): X=sparse.csr_matrix(X)
    v=np.asarray(X.sum(axis=1)).ravel().astype(float)
    lib=obs.total_counts.to_numpy(float)
    return np.log2(v/np.maximum(lib,1)*1e6+1), v

# Query genes and module genes
mods=pd.read_csv(MOD,sep='\t')
module_genes={r.module:[g.strip() for g in str(r.genes).split(',') if g.strip()] for r in mods.itertuples(index=False)}
all_needed=['FAP','DPP4']+sum(module_genes.values(),[])
expr={}; raw={}; coverage=[]
for g in dict.fromkeys(all_needed):
    z=expr_symbol(g)
    if z is not None: expr[g],raw[g]=z
for m,gs in module_genes.items():
    pr=[g for g in gs if g in expr]
    coverage.append([m,len(gs),len(pr),';'.join(pr)])
pd.DataFrame(coverage,columns=['module','n_genes','n_available','available_genes']).to_csv(OUT/'GSE210616_module_gene_coverage.tsv',sep='\t',index=False)
if 'FAP' not in expr or 'DPP4' not in expr: raise RuntimeError('FAP/DPP4 unavailable')

# FAP-free fibroblast affinity using frozen HBCA marker list and same scoring algorithm.
sig=pd.read_csv(SIG,sep='\t')
fib_markers=sig.loc[sig.lineage.eq('Fibroblast'),'gene'].astype(str).tolist()
fib_markers=[g for g in fib_markers if g.upper()!='FAP' and g.upper() in symcols]
# Cache every fibroblast-marker expression vector once; do not reslice the 36k-gene matrix per section.
for g in fib_markers:
    if g not in expr:
        z=expr_symbol(g)
        if z is not None: expr[g],raw[g]=z
fib_free=np.full(A.n_obs,np.nan,float)
module_scores={m:np.full(A.n_obs,np.nan,float) for m in module_genes}
fap_z=np.full(A.n_obs,np.nan,float); dpp4_z=np.full(A.n_obs,np.nan,float)

for gsm,idx0 in obs.groupby('gsm',observed=True).indices.items():
    idx=np.asarray(idx0,int); q=qc[idx]
    if q.sum()<20: continue
    # FAP-free fibroblast signature
    G=[]
    for g in fib_markers:
        G.append(expr[g][idx])
    G=np.vstack(G).T
    mu=np.nanmean(G[q],axis=0); sd=np.nanstd(G[q],axis=0); sd[sd<1e-6]=1
    Z=np.clip((G-mu)/sd,-4,4)
    fib_free[idx]=np.nanmean(Z,axis=1)
    # standardized FAP / DPP4 within section (for a scale-comparable axis)
    for gene,dest in [('FAP',fap_z),('DPP4',dpp4_z)]:
        e=expr[gene][idx]; mu=np.nanmean(e[q]); sd=np.nanstd(e[q]); sd=sd if sd>=1e-6 else 1
        dest[idx]=np.clip((e-mu)/sd,-4,4)
    # module scores: exact frozen FAP-project modules; standardize each gene within section
    for m,gs in module_genes.items():
        pr=[g for g in gs if g in expr]
        if not pr: continue
        MG=np.vstack([expr[g][idx] for g in pr]).T
        mu=np.nanmean(MG[q],axis=0); sd=np.nanstd(MG[q],axis=0); sd[sd<1e-6]=1
        MZ=np.clip((MG-mu)/sd,-4,4)
        module_scores[m][idx]=np.nanmean(MZ,axis=1)

axis=fap_z-dpp4_z
# FAP/DPP4 raw detection summaries per section and patient.
sec_rows=[]
for gsm,idx0 in obs.groupby('gsm',observed=True).indices.items():
    idx=np.asarray(idx0,int); q=qc[idx]; iq=idx[q]
    if len(iq)<20: continue
    r=obs.iloc[idx[0]]
    sec_rows.append({
      'gsm':gsm,'patient_id':str(r.patient_id),'section_id':str(r.section_id),'n_qc_spots':len(iq),
      'FAP_detection_fraction':float(np.mean(raw['FAP'][iq]>0)),'DPP4_detection_fraction':float(np.mean(raw['DPP4'][iq]>0)),
      'median_FAP_log2CPM1':float(np.median(expr['FAP'][iq])),'median_DPP4_log2CPM1':float(np.median(expr['DPP4'][iq])),
      'mean_FAP_log2CPM1':float(np.mean(expr['FAP'][iq])),'mean_DPP4_log2CPM1':float(np.mean(expr['DPP4'][iq]))
    })
Sdet=pd.DataFrame(sec_rows); Sdet.to_csv(OUT/'GSE210616_FAP_DPP4_section_detection.tsv',sep='\t',index=False)
Pdet=Sdet.groupby('patient_id').agg(n_sections=('gsm','nunique'),FAP_detection_fraction=('FAP_detection_fraction','median'),DPP4_detection_fraction=('DPP4_detection_fraction','median'),median_FAP_log2CPM1=('median_FAP_log2CPM1','median'),median_DPP4_log2CPM1=('median_DPP4_log2CPM1','median')).reset_index()
Pdet['FAP_minus_DPP4_detection']=Pdet.FAP_detection_fraction-Pdet.DPP4_detection_fraction
Pdet.to_csv(OUT/'GSE210616_FAP_DPP4_patient_detection.tsv',sep='\t',index=False)

# Helper partial Spearman: rank-transform y/x/covariates, residualize x/y on covariates, Pearson residual corr.
def partial_spearman(x,y,C):
    x=np.asarray(x,float); y=np.asarray(y,float); C=np.asarray(C,float)
    ok=np.isfinite(x)&np.isfinite(y)&np.all(np.isfinite(C),axis=1)
    if ok.sum()<20: return np.nan
    xr=rankdata(x[ok]); yr=rankdata(y[ok]); Cr=np.column_stack([rankdata(C[ok,j]) for j in range(C.shape[1])])
    X=np.column_stack([np.ones(ok.sum()),Cr])
    bx=np.linalg.lstsq(X,xr,rcond=None)[0]; by=np.linalg.lstsq(X,yr,rcond=None)[0]
    rx=xr-X@bx; ry=yr-X@by
    if np.std(rx)<1e-12 or np.std(ry)<1e-12: return np.nan
    return float(np.corrcoef(rx,ry)[0,1])

# Section-level tests restricted to top 30% FAP-free fibroblast-affinity spots within each section.
# Covariates: fibroblast affinity itself, epithelial affinity, log1p library size, detected genes.
rows=[]; tert=[]
for gsm,idx0 in obs.groupby('gsm',observed=True).indices.items():
    idx=np.asarray(idx0,int); q=qc[idx]; iq=idx[q]
    if len(iq)<100: continue
    threshold=np.nanquantile(fib_free[iq],0.70)
    fiq=iq[fib_free[iq]>=threshold]
    if len(fiq)<50: continue
    r=obs.iloc[idx[0]]; C=np.column_stack([fib_free[fiq],obs.iloc[fiq].affinity_Epithelial.to_numpy(float),np.log1p(obs.iloc[fiq].total_counts.to_numpy(float)),obs.iloc[fiq].n_genes_by_counts.to_numpy(float)])
    # module correlations and high-vs-low axis effect within fibroblast-rich spots
    for m,v in module_scores.items():
        rho=partial_spearman(axis[fiq],v[fiq],C)
        # tertile descriptive effect, using continuous axis thresholds within section
        lo=np.nanquantile(axis[fiq],1/3); hi=np.nanquantile(axis[fiq],2/3)
        low=v[fiq][axis[fiq]<=lo]; high=v[fiq][axis[fiq]>=hi]
        diff=float(np.nanmedian(high)-np.nanmedian(low)) if len(low) and len(high) else np.nan
        rows.append({'gsm':gsm,'patient_id':str(r.patient_id),'section_id':str(r.section_id),'module':m,'n_fibroblast_rich_spots':len(fiq),'partial_spearman_axis_module':rho,'high_minus_low_axis_module_score':diff})
    # fibroblast and epithelial associations of FAP/DPP4 separately
    for gene,e in [('FAP',expr['FAP']),('DPP4',expr['DPP4'])]:
        rho_f=partial_spearman(e[iq],fib_free[iq],np.column_stack([obs.iloc[iq].affinity_Epithelial.to_numpy(float),np.log1p(obs.iloc[iq].total_counts.to_numpy(float)),obs.iloc[iq].n_genes_by_counts.to_numpy(float)]))
        rho_e=partial_spearman(e[iq],obs.iloc[iq].affinity_Epithelial.to_numpy(float),np.column_stack([fib_free[iq],np.log1p(obs.iloc[iq].total_counts.to_numpy(float)),obs.iloc[iq].n_genes_by_counts.to_numpy(float)]))
        tert.append({'gsm':gsm,'patient_id':str(r.patient_id),'section_id':str(r.section_id),'gene':gene,'partial_rho_fibroblast_affinity':rho_f,'partial_rho_epithelial_affinity':rho_e})
S=pd.DataFrame(rows); S.to_csv(OUT/'GSE210616_FAP_DPP4_module_section_results.tsv',sep='\t',index=False)
G=pd.DataFrame(tert); G.to_csv(OUT/'GSE210616_FAP_DPP4_lineage_section_results.tsv',sep='\t',index=False)

# Collapse repeated sections within patient; signed-rank across patients. BH across modules separately for rho and high-low effect.
P=(S.groupby(['patient_id','module'],observed=True)
   .agg(n_sections=('gsm','nunique'),partial_spearman_axis_module=('partial_spearman_axis_module','median'),high_minus_low_axis_module_score=('high_minus_low_axis_module_score','median')).reset_index())
P.to_csv(OUT/'GSE210616_FAP_DPP4_module_patient_results.tsv',sep='\t',index=False)
summary=[]
for m,g in P.groupby('module',observed=True):
    for metric in ['partial_spearman_axis_module','high_minus_low_axis_module_score']:
        v=g[metric].dropna().to_numpy(float); p=np.nan
        if len(v)>=5 and np.any(v!=0): p=float(wilcoxon(v,zero_method='wilcox',alternative='two-sided',method='auto').pvalue)
        summary.append({'module':m,'metric':metric,'n_patients':len(v),'median':float(np.median(v)) if len(v) else np.nan,'q25':float(np.quantile(v,.25)) if len(v) else np.nan,'q75':float(np.quantile(v,.75)) if len(v) else np.nan,'fraction_positive':float(np.mean(v>0)) if len(v) else np.nan,'wilcoxon_p':p})
T=pd.DataFrame(summary); T['wilcoxon_q']=np.nan
for metric,idx in T.groupby('metric').groups.items():
    ids=T.loc[idx].index[T.loc[idx,'wilcoxon_p'].notna()]
    if len(ids): T.loc[ids,'wilcoxon_q']=multipletests(T.loc[ids,'wilcoxon_p'],method='fdr_bh')[1]
T.to_csv(OUT/'GSE210616_FAP_DPP4_module_patient_summary.tsv',sep='\t',index=False)

# Gene lineage association patient summary.
GP=(G.groupby(['patient_id','gene'],observed=True).agg(n_sections=('gsm','nunique'),partial_rho_fibroblast_affinity=('partial_rho_fibroblast_affinity','median'),partial_rho_epithelial_affinity=('partial_rho_epithelial_affinity','median')).reset_index())
GP.to_csv(OUT/'GSE210616_FAP_DPP4_lineage_patient_results.tsv',sep='\t',index=False)
gsum=[]
for gene,g in GP.groupby('gene',observed=True):
    for metric in ['partial_rho_fibroblast_affinity','partial_rho_epithelial_affinity']:
        v=g[metric].dropna().to_numpy(float); p=float(wilcoxon(v,zero_method='wilcox',alternative='two-sided',method='auto').pvalue) if len(v)>=5 and np.any(v!=0) else np.nan
        gsum.append({'gene':gene,'metric':metric,'n_patients':len(v),'median':float(np.median(v)),'q25':float(np.quantile(v,.25)),'q75':float(np.quantile(v,.75)),'fraction_positive':float(np.mean(v>0)),'wilcoxon_p':p})
GS=pd.DataFrame(gsum); GS['wilcoxon_q']=multipletests(GS.wilcoxon_p.fillna(1),method='fdr_bh')[1]
GS.to_csv(OUT/'GSE210616_FAP_DPP4_lineage_patient_summary.tsv',sep='\t',index=False)

# Patient-level FAP-vs-DPP4 detection paired test.
v=Pdet.FAP_detection_fraction.to_numpy()-Pdet.DPP4_detection_fraction.to_numpy()
pdet=float(wilcoxon(v,zero_method='wilcox',alternative='two-sided',method='auto').pvalue)
pd.DataFrame([{'n_patients':len(v),'median_FAP_detection':Pdet.FAP_detection_fraction.median(),'median_DPP4_detection':Pdet.DPP4_detection_fraction.median(),'median_FAP_minus_DPP4_detection':np.median(v),'fraction_FAP_gt_DPP4':np.mean(v>0),'wilcoxon_p':pdet}]).to_csv(OUT/'GSE210616_FAP_vs_DPP4_detection_summary.tsv',sep='\t',index=False)

print('patients',obs.patient_id.nunique(),'sections',obs.gsm.nunique(),'QC',int(qc.sum()))
print('\nDetection summary')
print(pd.read_csv(OUT/'GSE210616_FAP_vs_DPP4_detection_summary.tsv',sep='\t').to_string(index=False))
print('\nModule patient summary')
print(T.sort_values(['metric','wilcoxon_q']).to_string(index=False))
print('\nLineage summary')
print(GS.to_string(index=False))