#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parents[1]; R=ROOT/'results'
d=pd.read_parquet(R/'33c_state_module_profiles_by_donor_collapsed_exacttissue.parquet')
comparisons=[('FAP+DPP4-','FAP-DPP4+','FAP_only_minus_DPP4_only'),('FAP+DPP4-','FAP-DPP4-','FAP_only_minus_double_negative'),('FAP-DPP4+','FAP-DPP4-','DPP4_only_minus_double_negative'),('FAP+DPP4+','FAP-DPP4-','double_positive_minus_double_negative')]
rows=[]; groupcols=['tissue_general','tissue','disease','cell_type','module']
for keys,g in d.groupby(groupcols,dropna=False,observed=True):
    pv=g.pivot(index=['donor_key','collection_id'],columns='state',values='module_mean_gene_z')
    for A,B,name in comparisons:
        if A not in pv or B not in pv: continue
        q=pv[[A,B]].dropna();
        if len(q)<2: continue
        diff=q[A]-q[B]; p=np.nan
        if len(q)>=5 and np.any(diff!=0):
            try:p=float(wilcoxon(diff,zero_method='wilcox',alternative='two-sided',method='auto').pvalue)
            except Exception: pass
        rows.append({**dict(zip(groupcols,keys)),'comparison':name,'state_A':A,'state_B':B,'n_paired_donors':len(q),'n_collections':q.reset_index().collection_id.nunique(),'median_paired_difference_gene_z':diff.median(),'q25_paired_difference_gene_z':diff.quantile(.25),'q75_paired_difference_gene_z':diff.quantile(.75),'fraction_pairs_A_gt_B':float((diff>0).mean()),'wilcoxon_p':p})
out=pd.DataFrame(rows); out['wilcoxon_q']=np.nan
for comp,idx in out.groupby('comparison').groups.items():
    ids=out.loc[idx].index[out.loc[idx,'wilcoxon_p'].notna()]
    if len(ids): out.loc[ids,'wilcoxon_q']=multipletests(out.loc[ids,'wilcoxon_p'],method='fdr_bh')[1]
out.to_csv(R/'38b_paired_state_module_contrasts_donorcollapsed_exacttissue.csv',index=False)
print('rows',len(out),'finite p',out.wilcoxon_p.notna().sum(),'q<.05',(out.wilcoxon_q<.05).sum())