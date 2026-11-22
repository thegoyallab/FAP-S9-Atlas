from pathlib import Path
import re,numpy as np,pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas'); R=ROOT/'single_cell/results'; C=ROOT/'single_cell/config'; OUT=ROOT/'targetability_v1/Q5_cellular_selectivity'
fib=pd.read_csv(C/'fibroblast_cell_types.tsv',sep='\t'); primary=set(fib.loc[fib.analysis_group.eq('primary_fibroblast'),'label'].astype(str)); dm=pd.read_csv(C/'disease_context_map.tsv',sep='\t')[['label','context']]
def site_ok(tg,d):
    q=str(tg).lower(); disease=str(d).lower()
    rules=[(r'breast',r'breast|exocrine gland'),(r'lung|small cell lung|non-small cell lung|pleomorphic carcinoma',r'lung|respiratory system'),(r'colon adenocarcinoma|colorectal cancer',r'colon|large intestine'),(r'gastric cancer',r'stomach'),(r'renal.*carcinoma|clear cell renal',r'kidney'),(r'prostatic',r'prostate'),(r'oropharynx squamous cell carcinoma|oral cavity squamous cell carcinoma',r'digestive system|mucosa'),(r'uveal melanoma',r'eye'),(r'basal cell carcinoma|melanoma|cutaneous squamous',r'skin'),(r'glioblastoma',r'brain'),(r'neuroblastoma',r'adrenal gland'),(r'ovarian cancer|malignant ovarian serous tumor',r'ovary|fallopian tube'),(r'liver cancer|hepatocellular carcinoma|intrahepatic cholangiocarcinoma',r'liver')]
    return any(re.search(dp,disease,re.I) and re.search(tp,q,re.I) for dp,tp in rules)
d=pd.read_parquet(R/'25_donor_tissue_disease_cell_type.parquet',filters=[('gene','==','FAP')]); d['collection_id']=d.donor_key.str.split('::').str[0]; d=d.merge(dm,left_on='disease',right_on='label',how='left').drop(columns='label'); d['compartment']=np.where(d.cell_type.astype(str).isin(primary),'primary_fibroblast','other_annotated_cells'); d['organ_concordant']=d.apply(lambda r: site_ok(r.tissue_general,r.disease) if r.context=='primary_cancer' else True,axis=1)
# primary analysis: healthy + organ-concordant cancer
q=d[(d.context=='healthy')|((d.context=='primary_cancer')&d.organ_concordant)].copy()
base=['donor_key','collection_id','tissue_general','disease','context','compartment']; a=q.groupby(base).agg(raw_sum=('donor_raw_sum','sum'),total_library=('donor_total_library_raw','sum'),n_cells=('donor_n_cells','sum'),n_cells_detected=('donor_n_cells_detected','sum')).reset_index(); a['log2_cpm']=np.log2(1+1e6*a.raw_sum/a.total_library.replace(0,np.nan)); a['pct_detected']=100*a.n_cells_detected/a.n_cells.replace(0,np.nan); a.to_parquet(OUT/'04_donor_primaryfibroblast_vs_other_organconcordant.parquet',index=False)
p=a.pivot_table(index=['donor_key','collection_id','tissue_general','disease','context'],columns='compartment',values=['log2_cpm','pct_detected'],aggfunc='first').reset_index(); p.columns=['_'.join([x for x in c if x]) if isinstance(c,tuple) else c for c in p.columns]
rows=[]
for (tg,ctx),g in p.groupby(['tissue_general','context']):
    for metric in ['log2_cpm','pct_detected']:
        f=f'{metric}_primary_fibroblast'; o=f'{metric}_other_annotated_cells'; z=g.dropna(subset=[f,o]); n=len(z)
        if n<5: continue
        delta=(z[f]-z[o]).to_numpy();
        try:stat,pv=wilcoxon(z[f],z[o],alternative='two-sided',zero_method='wilcox',method='auto')
        except ValueError:stat,pv=np.nan,1
        rows.append(dict(tissue_general=tg,context=ctx,metric=metric,n_paired_donors=n,n_collections=z.collection_id.nunique(),median_primary_fibroblast=float(z[f].median()),median_other=float(z[o].median()),median_paired_delta=float(np.median(delta)),fraction_primary_gt_other=float(np.mean(delta>0)),p=float(pv)))
r=pd.DataFrame(rows); r['q']=np.nan
for (ctx,m),g in r.groupby(['context','metric']): r.loc[g.index,'q']=multipletests(g.p,method='fdr_bh')[1]
r.to_csv(OUT/'05_organconcordant_tissue_cellular_selectivity.tsv',sep='\t',index=False)
# exact outside-primary labels in organ-concordant cancer and healthy, robust >=5 donors >=2 collections
other=q[q.compartment=='other_annotated_cells']; s=other.groupby(['context','cell_type']).agg(n_donors=('donor_key','nunique'),n_collections=('collection_id','nunique'),median_log2_cpm=('donor_log2_cpm','median'),median_pct_detected=('donor_pct_cells_detected','median'),q25_log2_cpm=('donor_log2_cpm',lambda v:v.quantile(.25)),q75_log2_cpm=('donor_log2_cpm',lambda v:v.quantile(.75))).reset_index(); s['robust']=(s.n_donors>=5)&(s.n_collections>=2); s.to_csv(OUT/'06_other_celltype_FAP_organconcordant_summary.tsv',sep='\t',index=False)
# selected clinically interpretable exact labels
selected=['endothelial cell','blood vessel endothelial cell','endothelial cell of artery','vein endothelial cell','capillary endothelial cell','endothelial cell of vascular tree','pericyte','vascular associated smooth muscle cell','smooth muscle cell','malignant cell','neoplastic cell','epithelial cell','B cell','T cell','macrophage','monocyte','natural killer cell']
sel=s[(s.context=='primary_cancer')&s.cell_type.isin(selected)].sort_values('median_log2_cpm',ascending=False); sel.to_csv(OUT/'07_selected_offcompartment_celltypes_cancer.tsv',sep='\t',index=False)
print('ORGAN-CONCORDANT CANCER SELECTIVITY log2CPM')
z=r[(r.context=='primary_cancer')&(r.metric=='log2_cpm')].sort_values('median_paired_delta',ascending=False); print(z[['tissue_general','n_paired_donors','n_collections','median_primary_fibroblast','median_other','median_paired_delta','fraction_primary_gt_other','q']].to_string(index=False)); print('positive',int((z.median_paired_delta>0).sum()),'of',len(z),'signif positive',int(((z.median_paired_delta>0)&(z.q<.05)).sum()))
print('\nTOP OTHER-ANNOTATED CELL TYPES, ORGAN-CONCORDANT CANCER')
print(s[(s.context=='primary_cancer')&s.robust].sort_values('median_log2_cpm',ascending=False).head(30)[['cell_type','n_donors','n_collections','median_log2_cpm','median_pct_detected']].to_string(index=False))
print('\nSELECTED OFF-COMPARTMENT LABELS')
print(sel[['cell_type','n_donors','n_collections','median_log2_cpm','median_pct_detected','robust']].to_string(index=False))
