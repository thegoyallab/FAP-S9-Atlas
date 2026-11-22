from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas'); OUT=ROOT/'targetability_v1/Q6_integrated_targetability'
q1=pd.read_csv(ROOT/'targetability_v1/11_organconcordant_cancer_vs_healthy_exact_tissue.tsv',sep='\t')
q2=pd.read_csv(ROOT/'targetability_v1/Q2_cross_cancer/05_canonical_cancer_module_effects.tsv',sep='\t')
q5=pd.read_csv(ROOT/'targetability_v1/Q5_cellular_selectivity/05_organconcordant_tissue_cellular_selectivity.tsv',sep='\t')
core=['ECM_matrix','TGFb_LRRC15','myofibroblast','COL15A1_like','inflammatory','detoxification']
corefx=q2[q2.module.isin(core)].groupby('cancer_group').agg(core_remodeling_effect=('median_delta','mean'),min_paired_state_donors=('n_paired_donors','min'),max_state_collections=('n_collections','max'),core_modules_positive=('median_delta',lambda v:int((v>0).sum()))).reset_index()
# mappings deliberately explicit
maps=[
 ('breast cancer','breast','breast'),('lung cancer','lung','lung'),('colorectal cancer','colorectum','large intestine'),('prostate cancer','prostate gland','prostate gland'),('adrenal neuroblastoma','adrenal gland','adrenal gland'),('gastric cancer','stomach','stomach'),('oropharyngeal cancer','oropharynx','digestive system')]
rows=[]
for cg,tissue,tg in maps:
    a=q1[q1.tissue.eq(tissue)]
    b=corefx[corefx.cancer_group.eq(cg)]
    c=q5[(q5.context=='primary_cancer')&(q5.metric=='log2_cpm')&(q5.tissue_general.eq(tg))]
    if len(a)!=1 or len(b)!=1 or len(c)!=1: print('MISSING',cg,len(a),len(b),len(c)); continue
    a=a.iloc[0]; b=b.iloc[0]; c=c.iloc[0]
    rows.append(dict(cancer_group=cg,exact_tissue=tissue,tissue_general=tg,
        healthy_FAP_fraction=a.median_FAP_healthy,cancer_FAP_fraction=a.median_FAP_primary_cancer,cancer_minus_healthy_FAP_fraction=a.delta_FAP,
        healthy_donors=int(a.n_donors_healthy),cancer_donors=int(a.n_donors_primary_cancer),healthy_collections=int(a.n_collections_healthy),cancer_collections=int(a.n_collections_primary_cancer),
        core_remodeling_effect=b.core_remodeling_effect,state_paired_donors=int(b.min_paired_state_donors),state_collections=int(b.max_state_collections),core_modules_positive=int(b.core_modules_positive),
        fibroblast_selectivity_log2cpm=c.median_paired_delta,cellular_selectivity_donors=int(c.n_paired_donors),cellular_selectivity_collections=int(c.n_collections),cellular_selectivity_q=c.q))
d=pd.DataFrame(rows); d.to_csv(OUT/'01_integrated_targetability_components.tsv',sep='\t',index=False)
print(d.sort_values('cancer_minus_healthy_FAP_fraction',ascending=False).to_string(index=False))
print('\nCORRELATIONS AMONG COMPONENTS')
for x,y in [('cancer_minus_healthy_FAP_fraction','core_remodeling_effect'),('cancer_minus_healthy_FAP_fraction','fibroblast_selectivity_log2cpm'),('core_remodeling_effect','fibroblast_selectivity_log2cpm')]:
 from scipy.stats import spearmanr
 rho,p=spearmanr(d[x],d[y]); print(x,'vs',y,'rho',rho,'p',p)
