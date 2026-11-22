from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import wilcoxon,norm,binomtest
from statsmodels.stats.multitest import multipletests
OUT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas/targetability_v1/Q2_cross_cancer')
d=pd.read_parquet(OUT/'02b_donor_module_deltas.parquet')
# only exact tissues that passed >=5 paired donors for all 8 modules in first pass
support=pd.read_csv(OUT/'03_eligible_cancer_tissues.tsv',sep='\t'); d=d.merge(support[['tissue_general','tissue']],on=['tissue_general','tissue'],how='inner')
def canon(t):
    if t in {'ovary','left ovary','right ovary'}: return 'ovarian cancer'
    return {'breast':'breast cancer','lung':'lung cancer','colorectum':'colorectal cancer','prostate gland':'prostate cancer','adrenal gland':'adrenal neuroblastoma','stomach':'gastric cancer','oropharynx':'oropharyngeal cancer','oral cavity':'oral cavity cancer'}.get(t,t)
d['cancer_group']=d.tissue.map(canon)
# one donor per cancer group/module; median if donor contributes multiple exact tissue labels
u=d.groupby(['cancer_group','module','donor_key','collection_id']).delta.median().reset_index()
rows=[]
for (cg,m),g in u.groupby(['cancer_group','module']):
    x=g.delta.to_numpy(float); n=len(x)
    if n<5: continue
    try:stat,p=wilcoxon(x,alternative='two-sided',zero_method='wilcox',method='auto')
    except ValueError:stat,p=np.nan,1
    rows.append(dict(cancer_group=cg,module=m,n_paired_donors=n,n_collections=g.collection_id.nunique(),median_delta=float(np.median(x)),mean_delta=float(np.mean(x)),sd_delta=float(np.std(x,ddof=1)),se_delta=float(np.std(x,ddof=1)/np.sqrt(n)),fraction_positive=float(np.mean(x>0)),p=float(p)))
r=pd.DataFrame(rows); r['q_within_cancer']=np.nan
for cg,g in r.groupby('cancer_group'): r.loc[g.index,'q_within_cancer']=multipletests(g.p,method='fdr_bh')[1]
r.to_csv(OUT/'05_canonical_cancer_module_effects.tsv',sep='\t',index=False)
# keep cancer groups with all 8 modules
nmods=r.groupby('cancer_group').module.nunique(); keep=nmods[nmods==8].index; r=r[r.cancer_group.isin(keep)].copy()
meta=[]
for m,g in r.groupby('module'):
    g=g[np.isfinite(g.se_delta)&(g.se_delta>0)]; yi=g.mean_delta.to_numpy(); vi=g.se_delta.to_numpy()**2; k=len(g); wi=1/vi; yfix=(wi*yi).sum()/wi.sum(); Q=(wi*(yi-yfix)**2).sum(); C=wi.sum()-(wi**2).sum()/wi.sum(); tau2=max(0,(Q-(k-1))/C) if C>0 else 0; wr=1/(vi+tau2); mu=(wr*yi).sum()/wr.sum(); se=(1/wr.sum())**.5; p=2*norm.sf(abs(mu/se)); I2=max(0,(Q-(k-1))/Q)*100 if Q>0 else 0; med=float(g.median_delta.median()); pos=int((g.median_delta>0).sum()); neg=int((g.median_delta<0).sum()); ps=float(binomtest(pos,pos+neg,.5).pvalue) if pos+neg else 1
    meta.append(dict(module=m,n_cancer_groups=k,total_unique_group_donors=int(g.n_paired_donors.sum()),groups_positive=pos,groups_negative=neg,median_group_delta=med,random_effect_mean=mu,se=se,ci95_low=mu-1.96*se,ci95_high=mu+1.96*se,p_meta=p,tau2=tau2,I2=I2,p_sign=ps))
meta=pd.DataFrame(meta); meta['q_meta']=multipletests(meta.p_meta,method='fdr_bh')[1]; meta['q_sign']=multipletests(meta.p_sign,method='fdr_bh')[1]; meta=meta.sort_values('random_effect_mean',ascending=False); meta.to_csv(OUT/'06_canonical_cross_cancer_meta.tsv',sep='\t',index=False)
print('CANONICAL CANCER GROUPS')
print(r.groupby('cancer_group').agg(n_paired=('n_paired_donors','max'),n_collections=('n_collections','max')).sort_values('n_paired',ascending=False).to_string())
print('\nMETA')
print(meta.to_string(index=False))
print('\nEFFECT MATRIX')
print(r.pivot(index='cancer_group',columns='module',values='median_delta').round(3).to_string())
