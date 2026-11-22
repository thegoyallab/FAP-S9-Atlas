#!/usr/bin/env python3
from pathlib import Path
import numpy as np, pandas as pd
ROOT=Path(__file__).resolve().parents[1]; R=ROOT/'results'; C=ROOT/'config'
# Primary fibroblast definitions/exclusions
ex=pd.read_csv(C/'annotation_context_exclusions.tsv',sep='\t')
fib=pd.read_csv(C/'fibroblast_cell_types.tsv',sep='\t')
primary=set(fib.loc[fib.analysis_group=='primary_fibroblast','label'].astype(str))
dmap=pd.read_csv(C/'disease_context_map.tsv',sep='\t')[["label","context","needs_manual_review"]]
def excluded(df):
    m=~df.cell_type.astype(str).isin(primary)
    for r in ex.itertuples(index=False):
        m |= (df.tissue_general.astype(str)==str(r.tissue_general)) & (df.cell_type.astype(str)==str(r.cell_type))
    return m
# ---- state fractions: collapse assay/suspension before context summaries ----
s=pd.read_csv(R/'30_FAP_DPP4_state_by_donor.csv.gz')
base=['donor_key','collection_id','tissue_general','tissue','disease','cell_type']
state=(s.groupby(base+['state'],dropna=False,observed=True)
       .agg(state_n_cells=('state_n_cells','sum'),n_datasets=('n_datasets','sum'))
       .reset_index())
den=state.groupby(base,dropna=False,observed=True).state_n_cells.sum().rename('donor_fibroblast_cells').reset_index()
state=state.merge(den,on=base,how='left'); state['state_fraction']=state.state_n_cells/state.donor_fibroblast_cells
wide=state.pivot_table(index=base,columns='state',values='state_fraction',aggfunc='sum',fill_value=0).reset_index()
for c in ['FAP-DPP4-','FAP-DPP4+','FAP+DPP4-','FAP+DPP4+']:
    if c not in wide: wide[c]=0.0
wide['FAP_pos_fraction']=wide['FAP+DPP4-']+wide['FAP+DPP4+']; wide['DPP4_pos_fraction']=wide['FAP-DPP4+']+wide['FAP+DPP4+']
wide['double_positive_fraction']=wide['FAP+DPP4+']; wide['positive_union_fraction']=wide['FAP+DPP4-']+wide['FAP-DPP4+']+wide['FAP+DPP4+']
wide['double_positive_share_of_positive_union']=np.where(wide.positive_union_fraction>0,wide.double_positive_fraction/wide.positive_union_fraction,np.nan)
wide['FAP_minus_DPP4_fraction']=wide.FAP_pos_fraction-wide.DPP4_pos_fraction
wide=wide.merge(dmap,left_on='disease',right_on='label',how='left').drop(columns='label')
wide['annotation_context_excluded']=excluded(wide); wide['primary_fibroblast_analysis']=~wide.annotation_context_excluded
wide.to_csv(R/'31d_FAP_DPP4_axis_by_donor_collapsed_primary.csv.gz',index=False,compression='gzip')
a=wide[wide.primary_fibroblast_analysis].copy()
metrics=['FAP_pos_fraction','DPP4_pos_fraction','double_positive_fraction','double_positive_share_of_positive_union','FAP_minus_DPP4_fraction']
rows=[]
for keys,g in a.groupby(['tissue_general','tissue','disease','context','cell_type'],dropna=False,observed=True):
    rec=dict(zip(['tissue_general','tissue','disease','context','cell_type'],keys)); rec['n_donors']=g.donor_key.nunique(); rec['n_collections']=g.collection_id.nunique()
    for m in metrics:
        rec['median_'+m]=g[m].median(); rec['q25_'+m]=g[m].quantile(.25); rec['q75_'+m]=g[m].quantile(.75)
    rows.append(rec)
ctx=pd.DataFrame(rows); ctx['reliable_3donors']=ctx.n_donors>=3; ctx['replicated_2collections']=ctx.n_collections>=2
ctx.to_csv(R/'35c_FAP_DPP4_axis_context_summary_donorcollapsed_primary.csv',index=False)
# ---- module profiles: collapse state-gene pseudobulk across assay/suspension ----
g=pd.read_parquet(R/'30b_FAP_DPP4_state_gene_pseudobulk.parquet')
g=g[~excluded(g)].copy()
gbase=base+['state','gene','feature_id']
gc=(g.groupby(gbase,dropna=False,observed=True)
    .agg(gene_raw_sum=('gene_raw_sum','sum'),state_total_library_raw=('state_total_library_raw','sum'),state_n_cells=('state_n_cells','sum'),n_datasets_gene_measured=('n_datasets_gene_measured','sum'))
    .reset_index())
gc['log2_cpm']=np.log2(1+1e6*gc.gene_raw_sum/gc.state_total_library_raw.replace(0,np.nan))
mods=pd.read_csv(C/'fibroblast_modules.tsv',sep='\t'); gene2module={}
for r in mods.itertuples(index=False):
    for gene in str(r.genes).split(','):
        gene=gene.strip()
        if gene: gene2module[gene]=r.module
gm=gc[gc.gene.isin(gene2module)].copy(); gm['module']=gm.gene.map(gene2module)
def zwithin(v):
    sd=v.std(ddof=0); return (v-v.mean())/sd if sd and np.isfinite(sd) else v*0.0
# Exact tissue is mandatory here (not tissue_general).
gm['gene_z_within_exact_tissue']=gm.groupby(['tissue','gene'],dropna=False,observed=True).log2_cpm.transform(zwithin)
module=(gm.groupby(base+['state','module'],dropna=False,observed=True)
        .agg(module_mean_log2_cpm=('log2_cpm','mean'),module_mean_gene_z=('gene_z_within_exact_tissue','mean'),n_module_genes_measured=('gene','nunique'))
        .reset_index())
module.to_parquet(R/'33c_state_module_profiles_by_donor_collapsed_exacttissue.parquet',index=False,compression='zstd')
ms=(module.groupby(['tissue_general','tissue','disease','cell_type','state','module'],dropna=False,observed=True)
    .agg(n_donors=('donor_key','nunique'),n_collections=('collection_id','nunique'),median_module_log2_cpm=('module_mean_log2_cpm','median'),median_module_gene_z=('module_mean_gene_z','median'),q25_module_gene_z=('module_mean_gene_z',lambda v:v.quantile(.25)),q75_module_gene_z=('module_mean_gene_z',lambda v:v.quantile(.75)),median_n_module_genes_measured=('n_module_genes_measured','median')).reset_index())
ms.to_csv(R/'33d_state_module_profiles_by_context_donorcollapsed_exacttissue.csv',index=False)
print('axis unique donor-context rows',len(wide),'primary',len(a),'contexts',len(ctx))
print('module donor rows',len(module),'context rows',len(ms))