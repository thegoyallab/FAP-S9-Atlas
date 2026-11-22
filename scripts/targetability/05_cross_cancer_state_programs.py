from pathlib import Path
import re, numpy as np, pandas as pd
from scipy.stats import wilcoxon, norm, binomtest
from statsmodels.stats.multitest import multipletests
ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas'); R=ROOT/'single_cell/results'; C=ROOT/'single_cell/config'; OUT=ROOT/'targetability_v1/Q2_cross_cancer'
# load genes/modules and primary fibroblast ontology universe
mods=pd.read_csv(C/'fibroblast_modules.tsv',sep='\t'); gene2mod={}
for r in mods.itertuples(index=False):
    for gene in str(r.genes).split(','):
        gene=gene.strip()
        if gene: gene2mod[gene]=r.module
module_names=list(mods.module)
fib=pd.read_csv(C/'fibroblast_cell_types.tsv',sep='\t'); primary=set(fib.loc[fib.analysis_group.eq('primary_fibroblast'),'label'].astype(str)); ex=pd.read_csv(C/'annotation_context_exclusions.tsv',sep='\t'); dm=pd.read_csv(C/'disease_context_map.tsv',sep='\t')[['label','context']]
# conservative organ-concordance rules for solid cancers; hematologic and obvious metastatic compartments excluded from primary meta-analysis
def site_ok(tg,t,d):
    q=(str(tg)+' '+str(t)).lower(); disease=str(d).lower()
    rules=[
      (r'breast', r'breast|mammary'),
      (r'lung|small cell lung|non-small cell lung|pleomorphic carcinoma', r'lung|bronch|respiratory'),
      (r'colon adenocarcinoma|colorectal cancer', r'colon|colorect|caec|large intestine'),
      (r'gastric cancer', r'stomach'),
      (r'renal.*carcinoma|clear cell renal', r'kidney|renal'),
      (r'prostatic', r'prostate'),
      (r'oropharynx squamous cell carcinoma', r'oropharynx'),
      (r'oral cavity squamous cell carcinoma', r'oral cavity'),
      (r'uveal melanoma', r'uvea|eye'),
      (r'basal cell carcinoma|melanoma|cutaneous squamous', r'skin'),
      (r'glioblastoma', r'brain|temporal lobe'),
      (r'neuroblastoma', r'adrenal gland'),
      (r'ovarian cancer|malignant ovarian serous tumor', r'ovary|fallopian tube'),
      (r'liver cancer|hepatocellular carcinoma|intrahepatic cholangiocarcinoma', r'liver'),
    ]
    return any(re.search(dp,disease,re.I) and re.search(tp,q,re.I) for dp,tp in rules)
# exact tissue donor-state-gene pseudobulk, collapsed across fibroblast annotations and assay/suspension
g=pd.read_parquet(R/'30b_FAP_DPP4_state_gene_pseudobulk.parquet')
g=g[g.cell_type.astype(str).isin(primary) & g.gene.isin(gene2mod)].copy()
for x in ex.itertuples(index=False): g=g[~((g.tissue_general.astype(str)==str(x.tissue_general))&(g.cell_type.astype(str)==str(x.cell_type)))]
g=g.merge(dm,left_on='disease',right_on='label',how='left').drop(columns='label')
# retain all noncancer contexts for exact-tissue z reference; for primary_cancer retain only organ-concordant rows
is_c=g.context.eq('primary_cancer'); g=g[(~is_c)|g.apply(lambda r:site_ok(r.tissue_general,r.tissue,r.disease),axis=1)].copy()
base=['donor_key','collection_id','tissue_general','tissue','context','state','gene','feature_id']
a=(g.groupby(base,dropna=False).agg(gene_raw_sum=('gene_raw_sum','sum'),total_library_raw=('state_total_library_raw','sum'),state_n_cells=('state_n_cells','sum')).reset_index())
a['log2_cpm']=np.log2(1+1e6*a.gene_raw_sum/a.total_library_raw.replace(0,np.nan)); a['module']=a.gene.map(gene2mod)
def zscore(v):
    sd=v.std(ddof=0); return (v-v.mean())/sd if sd and np.isfinite(sd) else v*0.0
a['gene_z_exact_tissue']=a.groupby(['tissue','gene'],dropna=False).log2_cpm.transform(zscore)
mod=(a.groupby(['donor_key','collection_id','tissue_general','tissue','context','state','module'],dropna=False).agg(module_gene_z=('gene_z_exact_tissue','mean'),module_log2_cpm=('log2_cpm','mean'),n_genes=('gene','nunique'),state_n_cells=('state_n_cells','max')).reset_index())
mod.to_parquet(OUT/'01_donor_state_module_compartment_collapsed.parquet',index=False)
# Primary cancer FAP-only vs DPP4-only per exact tissue and module
c=mod[(mod.context=='primary_cancer') & mod.state.isin(['FAP+DPP4-','FAP-DPP4+'])].copy()
rows=[]; deltas=[]
for (tg,tissue,module),z in c.groupby(['tissue_general','tissue','module']):
    p=z.pivot_table(index=['donor_key','collection_id'],columns='state',values='module_gene_z',aggfunc='first')
    if not {'FAP+DPP4-','FAP-DPP4+'}.issubset(p.columns): continue
    p=p.dropna(subset=['FAP+DPP4-','FAP-DPP4+']); n=len(p)
    if n<2: continue
    d=(p['FAP+DPP4-']-p['FAP-DPP4+']).to_numpy(float)
    try: stat,pv=wilcoxon(p['FAP+DPP4-'],p['FAP-DPP4+'],alternative='two-sided',zero_method='wilcox',method='auto')
    except ValueError: stat,pv=np.nan,1.0
    rows.append(dict(tissue_general=tg,tissue=tissue,module=module,n_paired_donors=n,n_collections=p.reset_index().collection_id.nunique(),median_delta=float(np.median(d)),mean_delta=float(np.mean(d)),sd_delta=float(np.std(d,ddof=1)) if n>1 else np.nan,se_delta=float(np.std(d,ddof=1)/np.sqrt(n)) if n>1 else np.nan,fraction_positive=float(np.mean(d>0)),p=float(pv)))
    for (donor,coll),val in zip(p.index,d): deltas.append(dict(tissue_general=tg,tissue=tissue,module=module,donor_key=donor,collection_id=coll,delta=float(val)))
r=pd.DataFrame(rows); r['q_within_tissue']=np.nan
for tissue,z in r[r.n_paired_donors>=5].groupby('tissue'):
    r.loc[z.index,'q_within_tissue']=multipletests(z.p,method='fdr_bh')[1]
r.to_csv(OUT/'02_tissue_module_FAPonly_vs_DPP4only.tsv',sep='\t',index=False); pd.DataFrame(deltas).to_parquet(OUT/'02b_donor_module_deltas.parquet',index=False)
# eligible tissue panel n>=5 paired donors for all/most modules
elig=r[r.n_paired_donors>=5].copy(); tissue_support=elig.groupby(['tissue_general','tissue']).agg(min_paired=('n_paired_donors','min'),max_paired=('n_paired_donors','max'),n_modules=('module','nunique'),n_collections=('n_collections','max')).reset_index(); tissue_support=tissue_support[tissue_support.n_modules==len(module_names)]
tissue_support.to_csv(OUT/'03_eligible_cancer_tissues.tsv',sep='\t',index=False)
elig=elig.merge(tissue_support[['tissue_general','tissue']],on=['tissue_general','tissue'],how='inner')
# Random-effects meta-analysis of mean paired deltas across tissues for each module (DerSimonian-Laird), plus tissue sign consistency
meta=[]
for module,z in elig.groupby('module'):
    z=z[np.isfinite(z.se_delta)&(z.se_delta>0)].copy(); k=len(z)
    if k<2: continue
    yi=z.mean_delta.to_numpy(); vi=(z.se_delta.to_numpy())**2; wi=1/vi; yfix=np.sum(wi*yi)/np.sum(wi); Q=np.sum(wi*(yi-yfix)**2); Cc=np.sum(wi)-np.sum(wi**2)/np.sum(wi); tau2=max(0,(Q-(k-1))/Cc) if Cc>0 else 0; wr=1/(vi+tau2); mu=np.sum(wr*yi)/np.sum(wr); se=np.sqrt(1/np.sum(wr)); zz=mu/se; pv=2*norm.sf(abs(zz)); I2=max(0,(Q-(k-1))/Q)*100 if Q>0 else 0
    med=float(np.median(z.median_delta)); npos=int((z.median_delta>0).sum()); nneg=int((z.median_delta<0).sum()); p_sign=float(binomtest(npos,npos+nneg,p=.5,alternative='two-sided').pvalue) if npos+nneg else 1
    meta.append(dict(module=module,n_tissues=k,total_paired_donors=int(z.n_paired_donors.sum()),tissues_positive=npos,tissues_negative=nneg,median_tissue_delta=med,random_effect_mean=mu,random_effect_se=se,ci95_low=mu-1.96*se,ci95_high=mu+1.96*se,p_meta=pv,tau2=tau2,I2=I2,p_sign=p_sign))
m=pd.DataFrame(meta); m['q_meta']=multipletests(m.p_meta,method='fdr_bh')[1]; m['q_sign']=multipletests(m.p_sign,method='fdr_bh')[1]; m=m.sort_values('random_effect_mean',ascending=False); m.to_csv(OUT/'04_cross_cancer_module_meta_analysis.tsv',sep='\t',index=False)
print('ELIGIBLE ORGAN-CONCORDANT SOLID-CANCER TISSUES (>=5 paired FAP-only/DPP4-only donors, all 8 modules)')
print(tissue_support.sort_values('max_paired',ascending=False).to_string(index=False))
print('\nCROSS-CANCER META')
print(m.to_string(index=False))
print('\nTISSUE EFFECT MATRIX')
mat=elig.pivot(index='tissue',columns='module',values='median_delta'); print(mat.round(3).to_string())
