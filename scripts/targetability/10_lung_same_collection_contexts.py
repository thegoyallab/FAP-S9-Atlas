from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import kruskal,mannwhitneyu,spearmanr
from statsmodels.stats.multitest import multipletests
ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas'); R=ROOT/'single_cell/results'; Q2=ROOT/'targetability_v1/Q2_cross_cancer'; OUT=ROOT/'targetability_v1/Q4_malignant_vs_benign'
COLL='6f6d381a-7701-4781-935c-db10d30de293'; contexts=['healthy','fibrosis_remodeling','inflammatory_infectious','primary_cancer']; mods=['ECM_matrix','TGFb_LRRC15','myofibroblast','COL15A1_like','inflammatory','detoxification','antigen_presentation','PI16_homeostatic']
d=pd.read_parquet(Q2/'01_donor_state_module_compartment_collapsed.parquet')
x=d[(d.collection_id==COLL)&(d.tissue=='lung')&(d.state=='FAP+DPP4-')&d.context.isin(contexts)].copy()
# donor x module one value; export
a=x.groupby(['donor_key','context','module']).module_gene_z.median().reset_index(); a.to_parquet(OUT/'01_same_collection_lung_FAPonly_module_profiles.parquet',index=False)
# disease labels for these donors
mapd=pd.read_csv(R/'30_FAP_DPP4_state_by_donor.csv.gz',usecols=['donor_key','collection_id','tissue','disease']).drop_duplicates(); labels=a[['donor_key','context']].drop_duplicates().merge(mapd[(mapd.collection_id==COLL)&(mapd.tissue=='lung')],on='donor_key',how='left').drop_duplicates(); labels.to_csv(OUT/'01b_same_collection_lung_disease_labels.tsv',sep='\t',index=False)
print('DISEASE LABELS')
print(labels.groupby(['context','disease']).donor_key.nunique().to_string())
# omnibus + pairwise cancer vs 3 contexts
rows=[]; kw=[]
for m in mods:
    groups=[a[(a.module==m)&(a.context==c)].module_gene_z.dropna().to_numpy() for c in contexts]
    H,pkw=kruskal(*groups); kw.append(dict(module=m,H=H,p=pkw))
    cancer=groups[-1]
    for c,ref in zip(contexts[:-1],groups[:-1]):
        U,p=mannwhitneyu(cancer,ref,alternative='two-sided',method='auto'); diffs=(cancer[:,None]-ref[None,:]).ravel(); rows.append(dict(module=m,contrast=f'cancer_minus_{c}',n_cancer=len(cancer),n_reference=len(ref),median_cancer=float(np.median(cancer)),median_reference=float(np.median(ref)),HL_shift=float(np.median(diffs)),probability_superiority=float(U/(len(cancer)*len(ref))),p=float(p)))
kw=pd.DataFrame(kw); kw['q']=multipletests(kw.p,method='fdr_bh')[1]; kw.to_csv(OUT/'02_module_omnibus_KW.tsv',sep='\t',index=False)
r=pd.DataFrame(rows); r['q_within_contrast']=np.nan
for c,g in r.groupby('contrast'): r.loc[g.index,'q_within_contrast']=multipletests(g.p,method='fdr_bh')[1]
r.to_csv(OUT/'03_cancer_vs_benign_module_contrasts.tsv',sep='\t',index=False)
# median profiles, centroid correlations
med=a.groupby(['context','module']).module_gene_z.median().unstack().reindex(contexts)[mods]; med.to_csv(OUT/'04_context_median_module_profiles.tsv',sep='\t')
cor=[]
for c in contexts[:-1]:
    rho,p=spearmanr(med.loc['primary_cancer'],med.loc[c]); cor.append(dict(reference_context=c,rho_module_profile=rho,p=p,euclidean_distance=float(np.linalg.norm(med.loc['primary_cancer'].to_numpy()-med.loc[c].to_numpy()))))
pd.DataFrame(cor).to_csv(OUT/'04b_cancer_context_profile_similarity.tsv',sep='\t',index=False)
print('\nOMNIBUS')
print(kw.to_string(index=False))
print('\nPAIRWISE')
print(r.sort_values(['contrast','HL_shift'],ascending=[True,False]).to_string(index=False))
print('\nMEDIAN PROFILES')
print(med.round(3).to_string())
print('\nCANCER PROFILE SIMILARITY')
print(pd.DataFrame(cor).to_string(index=False))
