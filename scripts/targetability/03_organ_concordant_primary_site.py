from pathlib import Path
import re, numpy as np, pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas'); R=ROOT/'single_cell/results'; C=ROOT/'single_cell/config'; OUT=ROOT/'targetability_v1'
fib=pd.read_csv(C/'fibroblast_cell_types.tsv',sep='\t'); primary=set(fib.loc[fib.analysis_group.eq('primary_fibroblast'),'label'].astype(str)); ex=pd.read_csv(C/'annotation_context_exclusions.tsv',sep='\t'); dm=pd.read_csv(C/'disease_context_map.tsv',sep='\t')[['label','context']]
# explicit site-concordance rules for robust exact-tissue comparisons; lymphoma analyzed separately
rules={
 'breast': r'breast', 'lung':r'lung|small cell lung|non-small cell lung|pleomorphic carcinoma',
 'colorectum':r'colon adenocarcinoma|colorectal cancer', 'sigmoid colon':r'colorectal cancer|colon adenocarcinoma',
 'stomach':r'gastric cancer', 'oropharynx':r'oropharynx squamous cell carcinoma',
 'kidney':r'renal.*carcinoma|clear cell renal', 'prostate gland':r'prostatic', 'skin epidermis':r'melanoma|basal cell carcinoma|cutaneous squamous',
 'liver':r'liver cancer|hepatocellular carcinoma|cholangiocarcinoma', 'uvea':r'uveal melanoma',
 'adrenal gland':r'neuroblastoma|adrenocortical',
}
def concordant(tissue,disease):
    pat=rules.get(str(tissue)); return bool(pat and re.search(pat,str(disease),flags=re.I))
# source state counts
s=pd.read_csv(R/'30_FAP_DPP4_state_by_donor.csv.gz'); s=s[s.cell_type.astype(str).isin(primary)].copy()
for x in ex.itertuples(index=False): s=s[~((s.tissue_general.astype(str)==str(x.tissue_general))&(s.cell_type.astype(str)==str(x.cell_type)))]
s=s.merge(dm,left_on='disease',right_on='label',how='left').drop(columns='label')
s['site_concordant']=s.apply(lambda r: concordant(r.tissue,r.disease) if r.context=='primary_cancer' else True,axis=1)
# exclude nonconcordant cancer rows but retain healthy
z=s[(s.context=='healthy')|((s.context=='primary_cancer')&s.site_concordant)].copy()
base=['donor_key','collection_id','tissue_general','tissue','context']
g=z.groupby(base+['state'],dropna=False).state_n_cells.sum().reset_index(); w=g.pivot_table(index=base,columns='state',values='state_n_cells',aggfunc='sum',fill_value=0).reset_index()
for c in ['FAP-DPP4-','FAP-DPP4+','FAP+DPP4-','FAP+DPP4+']:
    if c not in w:w[c]=0
w['n_cells']=w[['FAP-DPP4-','FAP-DPP4+','FAP+DPP4-','FAP+DPP4+']].sum(axis=1); w['FAP_pos_fraction']=(w['FAP+DPP4-']+w['FAP+DPP4+'])/w.n_cells; w['DPP4_pos_fraction']=(w['FAP-DPP4+']+w['FAP+DPP4+'])/w.n_cells; w['FAP_minus_DPP4_fraction']=w.FAP_pos_fraction-w.DPP4_pos_fraction
w.to_parquet(OUT/'10_donor_exacttissue_organconcordant_cancer_vs_healthy.parquet',index=False)
# disease support after concordance
cs=s[(s.context=='primary_cancer')&s.site_concordant].drop_duplicates(['donor_key','tissue','disease']).groupby(['tissue_general','tissue','disease']).donor_key.nunique().reset_index(name='n')
cs.to_csv(OUT/'10b_organconcordant_cancer_disease_support.tsv',sep='\t',index=False)
# descriptive exact-tissue medians
summ=w.groupby(['tissue_general','tissue','context']).agg(n_donors=('donor_key','nunique'),n_collections=('collection_id','nunique'),median_FAP=('FAP_pos_fraction','median'),median_DPP4=('DPP4_pos_fraction','median'),median_axis=('FAP_minus_DPP4_fraction','median')).reset_index()
piv=summ.pivot(index=['tissue_general','tissue'],columns='context'); piv.columns=[f'{a}_{b}' for a,b in piv.columns]; piv=piv.reset_index()
if {'n_donors_healthy','n_donors_primary_cancer'}.issubset(piv.columns):
    piv=piv[piv.n_donors_healthy.notna()&piv.n_donors_primary_cancer.notna()].copy(); piv['robust5each']=(piv.n_donors_healthy>=5)&(piv.n_donors_primary_cancer>=5); piv['delta_FAP']=piv.median_FAP_primary_cancer-piv.median_FAP_healthy; piv['delta_axis']=piv.median_axis_primary_cancer-piv.median_axis_healthy
else: raise SystemExit('missing contexts')
piv.to_csv(OUT/'11_organconcordant_cancer_vs_healthy_exact_tissue.tsv',sep='\t',index=False)
# paired tests organ-concordant
rows=[]
for (tg,tissue),q in w.groupby(['tissue_general','tissue']):
    for m in ['FAP_pos_fraction','FAP_minus_DPP4_fraction']:
        p=q.pivot_table(index=['donor_key','collection_id'],columns='context',values=m,aggfunc='first')
        if not {'healthy','primary_cancer'}.issubset(p.columns):continue
        p=p.dropna(subset=['healthy','primary_cancer']); n=len(p)
        if n<2:continue
        delta=(p.primary_cancer-p.healthy).to_numpy()
        try:stat,pv=wilcoxon(p.primary_cancer,p.healthy,alternative='two-sided',zero_method='wilcox',method='auto')
        except ValueError:stat,pv=np.nan,1
        rows.append(dict(tissue_general=tg,tissue=tissue,metric=m,n_pairs=n,n_collections=p.reset_index().collection_id.nunique(),median_healthy=float(np.median(p.healthy)),median_cancer=float(np.median(p.primary_cancer)),median_paired_delta=float(np.median(delta)),pairs_positive=int((delta>0).sum()),pairs_zero=int((delta==0).sum()),pairs_negative=int((delta<0).sum()),p=float(pv)))
r=pd.DataFrame(rows); r['q_primary_n5']=np.nan
for m,q in r[r.n_pairs>=5].groupby('metric'):r.loc[q.index,'q_primary_n5']=multipletests(q.p,method='fdr_bh')[1]
r.to_csv(OUT/'12_organconcordant_matched_paired_tests.tsv',sep='\t',index=False)
print('ORGAN-CONCORDANT ROBUST EXACT-TISSUE COMPARISONS >=5 EACH')
print(piv[piv.robust5each][['tissue_general','tissue','n_donors_healthy','n_donors_primary_cancer','n_collections_healthy','n_collections_primary_cancer','median_FAP_healthy','median_FAP_primary_cancer','delta_FAP','median_axis_healthy','median_axis_primary_cancer','delta_axis']].sort_values('delta_FAP',ascending=False).to_string(index=False))
print('\nORGAN-CONCORDANT MATCHED n>=5')
print(r[r.n_pairs>=5].to_string(index=False))
print('\nEXCLUDED CANCER ROW SUPPORT FROM PREVIOUS ROBUST TISSUES')
allc=s[s.context.eq('primary_cancer')].drop_duplicates(['donor_key','tissue','disease']); exc=allc[~allc.site_concordant]
print(exc.groupby(['tissue_general','tissue','disease']).donor_key.nunique().reset_index(name='n').sort_values('n',ascending=False).head(35).to_string(index=False))
