from pathlib import Path
import numpy as np, pandas as pd
ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas'); OUT=ROOT/'targetability_v1'
d=pd.read_parquet(OUT/'00_donor_exacttissue_context_primary_fibroblast.parquet')
metrics=['FAP_pos_fraction','DPP4_pos_fraction','double_positive_fraction','FAP_minus_DPP4_fraction']
# exact tissue/context summaries after full fibroblast-compartment collapse
rows=[]
for keys,g in d.groupby(['tissue_general','tissue','context'],dropna=False):
    rec=dict(zip(['tissue_general','tissue','context'],keys)); rec['n_donors']=g.donor_key.nunique(); rec['n_collections']=g.collection_id.nunique(); rec['median_fibroblast_cells']=float(g.n_fibroblast_cells.median())
    for m in metrics:
        rec['median_'+m]=float(g[m].median()); rec['q25_'+m]=float(g[m].quantile(.25)); rec['q75_'+m]=float(g[m].quantile(.75))
    rows.append(rec)
s=pd.DataFrame(rows)
s['robust_n5']=s.n_donors>=5; s['replicated_2collections']=s.n_collections>=2
s.to_csv(OUT/'04_exact_tissue_context_landscape.tsv',sep='\t',index=False)
# cancer vs healthy descriptive table, exact tissue
x=s[s.context.isin(['healthy','primary_cancer'])].copy()
idx=['tissue_general','tissue']
w=x.pivot(index=idx,columns='context')
w.columns=[f'{a}_{b}' for a,b in w.columns]; w=w.reset_index()
needed=['n_donors_healthy','n_donors_primary_cancer','n_collections_healthy','n_collections_primary_cancer','median_FAP_pos_fraction_healthy','median_FAP_pos_fraction_primary_cancer','median_DPP4_pos_fraction_healthy','median_DPP4_pos_fraction_primary_cancer','median_FAP_minus_DPP4_fraction_healthy','median_FAP_minus_DPP4_fraction_primary_cancer']
for c in needed:
    if c not in w: w[c]=np.nan
w=w[(w.n_donors_healthy.notna())&(w.n_donors_primary_cancer.notna())].copy()
w['delta_FAP_median']=w.median_FAP_pos_fraction_primary_cancer-w.median_FAP_pos_fraction_healthy
w['delta_DPP4_median']=w.median_DPP4_pos_fraction_primary_cancer-w.median_DPP4_pos_fraction_healthy
w['delta_axis_median']=w.median_FAP_minus_DPP4_fraction_primary_cancer-w.median_FAP_minus_DPP4_fraction_healthy
w['robust_5_each']=(w.n_donors_healthy>=5)&(w.n_donors_primary_cancer>=5)
w['replicated_both']=(w.n_collections_healthy>=2)&(w.n_collections_primary_cancer>=2)
w.to_csv(OUT/'05_cancer_vs_healthy_exact_tissue_descriptive.tsv',sep='\t',index=False)
print('ROBUST >=5 DONORS EACH')
cols=['tissue_general','tissue','n_donors_healthy','n_donors_primary_cancer','n_collections_healthy','n_collections_primary_cancer','median_FAP_pos_fraction_healthy','median_FAP_pos_fraction_primary_cancer','delta_FAP_median','median_FAP_minus_DPP4_fraction_healthy','median_FAP_minus_DPP4_fraction_primary_cancer','delta_axis_median','replicated_both']
print(w[w.robust_5_each].sort_values('delta_FAP_median',ascending=False)[cols].to_string(index=False))
# Context overlaps with primary cancer for nonmalignant remodeling
contexts=['healthy','primary_cancer','fibrosis_remodeling','inflammatory_infectious','injury_remodeling']
y=s[s.context.isin(contexts)]
sets=y.groupby(['tissue_general','tissue']).context.apply(set)
print('\nCANCER + NONMALIGNANT CONTEXT OVERLAPS')
for key,cs in sets.items():
    if 'primary_cancer' in cs and len(cs & {'fibrosis_remodeling','inflammatory_infectious','injury_remodeling'})>0:
        z=y[(y.tissue_general==key[0])&(y.tissue==key[1])][['context','n_donors','n_collections','median_FAP_pos_fraction','median_DPP4_pos_fraction','median_FAP_minus_DPP4_fraction']]
        print('\n###',key); print(z.sort_values('context').to_string(index=False))
# top benign/remodeling contexts by FAP with n>=5
b=s[s.context.isin(['healthy','fibrosis_remodeling','inflammatory_infectious','injury_remodeling']) & (s.n_donors>=5)].copy()
print('\nTOP NONMALIGNANT FAP CONTEXTS n>=5')
print(b.sort_values('median_FAP_pos_fraction',ascending=False)[['tissue_general','tissue','context','n_donors','n_collections','median_FAP_pos_fraction','median_DPP4_pos_fraction','median_FAP_minus_DPP4_fraction']].head(35).to_string(index=False))
