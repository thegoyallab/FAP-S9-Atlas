from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas'); R=ROOT/'single_cell/results'; C=ROOT/'single_cell/config'; OUT=ROOT/'targetability_v1'
g=pd.read_parquet(R/'30b_FAP_DPP4_state_gene_pseudobulk.parquet',filters=[('gene','in',['FAP','DPP4'])])
fib=pd.read_csv(C/'fibroblast_cell_types.tsv',sep='\t'); primary=set(fib.loc[fib.analysis_group.eq('primary_fibroblast'),'label'].astype(str))
ex=pd.read_csv(C/'annotation_context_exclusions.tsv',sep='\t'); dmap=pd.read_csv(C/'disease_context_map.tsv',sep='\t')[['label','context']]
g=g[g.cell_type.astype(str).isin(primary)].copy()
for x in ex.itertuples(index=False):
    g=g[~((g.tissue_general.astype(str)==str(x.tissue_general))&(g.cell_type.astype(str)==str(x.cell_type)))]
g=g.merge(dmap,left_on='disease',right_on='label',how='left').drop(columns='label')
base=['donor_key','collection_id','tissue_general','tissue','context','gene']
a=(g.groupby(base,dropna=False).agg(gene_raw_sum=('gene_raw_sum','sum'),total_library_raw=('state_total_library_raw','sum'),n_cells=('state_n_cells','sum'),n_measured=('n_datasets_gene_measured','sum')).reset_index())
a['log2_cpm']=np.log2(1+1e6*a.gene_raw_sum/a.total_library_raw.replace(0,np.nan))
a.to_parquet(OUT/'06_donor_exacttissue_context_FAP_DPP4_log2cpm.parquet',index=False)
# paired cancer-normal tests n>=2, exact tissue
rows=[]
for (tg,tissue,gene),z in a[a.context.isin(['healthy','primary_cancer'])].groupby(['tissue_general','tissue','gene']):
    p=z.pivot_table(index=['donor_key','collection_id'],columns='context',values='log2_cpm',aggfunc='first')
    if not {'healthy','primary_cancer'}.issubset(p.columns): continue
    p=p.dropna(subset=['healthy','primary_cancer']); n=len(p)
    if n<2: continue
    delta=(p.primary_cancer-p.healthy).to_numpy(float)
    try: stat,pv=wilcoxon(p.primary_cancer,p.healthy,alternative='two-sided',zero_method='wilcox',method='auto')
    except ValueError: stat,pv=np.nan,1.0
    rows.append(dict(tissue_general=tg,tissue=tissue,gene=gene,n_pairs=n,n_collections=p.reset_index().collection_id.nunique(),median_healthy=float(np.median(p.healthy)),median_cancer=float(np.median(p.primary_cancer)),median_paired_delta=float(np.median(delta)),pairs_positive=int((delta>0).sum()),pairs_zero=int((delta==0).sum()),pairs_negative=int((delta<0).sum()),p=float(pv)))
r=pd.DataFrame(rows); r['q_within_gene']=np.nan
for gene,z in r.groupby('gene'): r.loc[z.index,'q_within_gene']=multipletests(z.p,method='fdr_bh')[1]
r['primary_gate_n5']=r.n_pairs>=5
r.to_csv(OUT/'07_matched_tumor_normal_FAP_DPP4_log2cpm.tsv',sep='\t',index=False)
# exact tissue context medians
s=a.groupby(['tissue_general','tissue','context','gene']).agg(n_donors=('donor_key','nunique'),n_collections=('collection_id','nunique'),median_log2_cpm=('log2_cpm','median'),q25_log2_cpm=('log2_cpm',lambda v:v.quantile(.25)),q75_log2_cpm=('log2_cpm',lambda v:v.quantile(.75))).reset_index()
s.to_csv(OUT/'08_exact_tissue_context_FAP_DPP4_log2cpm_summary.tsv',sep='\t',index=False)
print('MATCHED FAP LOG2 CPM, n>=5')
print(r[(r.gene=='FAP')&(r.primary_gate_n5)].sort_values('median_paired_delta',ascending=False).to_string(index=False))
print('\nMATCHED DPP4 LOG2 CPM, n>=5')
print(r[(r.gene=='DPP4')&(r.primary_gate_n5)].sort_values('median_paired_delta',ascending=False).to_string(index=False))
print('\nLUNG CONTEXT FAP LOG2 CPM')
print(s[(s.tissue.eq('lung'))&(s.gene.eq('FAP'))&s.context.isin(['healthy','primary_cancer','fibrosis_remodeling','inflammatory_infectious'])].to_string(index=False))
