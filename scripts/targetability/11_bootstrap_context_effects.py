from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import mannwhitneyu,kruskal
from statsmodels.stats.multitest import multipletests
OUT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas/targetability_v1/Q4_malignant_vs_benign')
a=pd.read_parquet(OUT/'01_same_collection_lung_FAPonly_module_profiles.parquet')
mods=['ECM_matrix','TGFb_LRRC15','myofibroblast','COL15A1_like','inflammatory','detoxification','antigen_presentation','PI16_homeostatic']
core=['ECM_matrix','TGFb_LRRC15','myofibroblast','COL15A1_like','inflammatory','detoxification']
contexts=['healthy','fibrosis_remodeling','inflammatory_infectious','primary_cancer']
# donor-wide core score
w=a.pivot_table(index=['donor_key','context'],columns='module',values='module_gene_z',aggfunc='first').reset_index(); w['core_score']=w[core].mean(axis=1); w.to_parquet(OUT/'05_same_collection_lung_donor_core_scores.parquet',index=False)
# core omnibus/pairwise
H,pkw=kruskal(*[w[w.context==c].core_score for c in contexts]); rows=[]
c=w[w.context=='primary_cancer'].core_score.to_numpy()
rng=np.random.default_rng(20261122)
for refctx in contexts[:-1]:
    ref=w[w.context==refctx].core_score.to_numpy(); U,p=mannwhitneyu(c,ref,alternative='two-sided'); diffs=(c[:,None]-ref[None,:]).ravel(); hl=float(np.median(diffs)); boots=[]
    for i in range(5000):
        cb=rng.choice(c,len(c),replace=True); rb=rng.choice(ref,len(ref),replace=True); boots.append(np.median((cb[:,None]-rb[None,:]).ravel()))
    lo,hi=np.quantile(boots,[.025,.975]); rows.append(dict(contrast=f'cancer_minus_{refctx}',n_cancer=len(c),n_reference=len(ref),median_cancer=np.median(c),median_reference=np.median(ref),HL_shift=hl,ci95_low=lo,ci95_high=hi,p=p,prob_superiority=U/(len(c)*len(ref))))
corecmp=pd.DataFrame(rows); corecmp['q']=multipletests(corecmp.p,method='fdr_bh')[1]; corecmp.to_csv(OUT/'06_core_score_context_contrasts.tsv',sep='\t',index=False)
# bootstrap HL CIs per module for cancer vs fibrosis/inflammation/healthy
rows=[]
for m in mods:
    can=a[(a.context=='primary_cancer')&(a.module==m)].module_gene_z.to_numpy()
    for refctx in contexts[:-1]:
        ref=a[(a.context==refctx)&(a.module==m)].module_gene_z.to_numpy(); U,p=mannwhitneyu(can,ref,alternative='two-sided'); hl=np.median((can[:,None]-ref[None,:]).ravel()); boots=[]
        for i in range(3000):
            cb=rng.choice(can,len(can),replace=True); rb=rng.choice(ref,len(ref),replace=True); boots.append(np.median((cb[:,None]-rb[None,:]).ravel()))
        lo,hi=np.quantile(boots,[.025,.975]); rows.append(dict(module=m,contrast=f'cancer_minus_{refctx}',HL_shift=hl,ci95_low=lo,ci95_high=hi,p=p,prob_superiority=U/(len(can)*len(ref))))
r=pd.DataFrame(rows); r['q_within_contrast']=np.nan
for cc,g in r.groupby('contrast'): r.loc[g.index,'q_within_contrast']=multipletests(g.p,method='fdr_bh')[1]
r.to_csv(OUT/'07_module_context_effects_bootstrap.tsv',sep='\t',index=False)
print('CORE KW',H,pkw)
print(corecmp.to_string(index=False))
print('\nMODULE CANCER-FIBROSIS')
print(r[r.contrast=='cancer_minus_fibrosis_remodeling'].sort_values('HL_shift',ascending=False).to_string(index=False))
print('\nMODULE CANCER-INFLAMMATORY')
print(r[r.contrast=='cancer_minus_inflammatory_infectious'].sort_values('HL_shift',ascending=False).to_string(index=False))
