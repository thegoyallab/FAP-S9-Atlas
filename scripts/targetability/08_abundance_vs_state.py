from pathlib import Path
import re,numpy as np,pandas as pd
from scipy.stats import spearmanr,norm
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas'); R=ROOT/'single_cell/results'; C=ROOT/'single_cell/config'; Q2=ROOT/'targetability_v1/Q2_cross_cancer'; OUT=ROOT/'targetability_v1/Q3_abundance_vs_state'
# canonical group map matching Q2 eligible groups
def cancer_group(tg,t,d):
    q=(str(tg)+' '+str(t)).lower(); disease=str(d).lower()
    if re.search(r'breast',disease) and re.search(r'breast|mammary',q): return 'breast cancer'
    if re.search(r'lung|small cell lung|non-small cell lung|pleomorphic carcinoma',disease) and re.search(r'lung|bronch|respiratory',q): return 'lung cancer'
    if re.search(r'colon adenocarcinoma|colorectal cancer',disease) and re.search(r'colon|colorect|caec|large intestine',q): return 'colorectal cancer'
    if re.search(r'gastric cancer',disease) and 'stomach' in q: return 'gastric cancer'
    if re.search(r'prostatic',disease) and 'prostate' in q: return 'prostate cancer'
    if re.search(r'neuroblastoma',disease) and 'adrenal gland' in q: return 'adrenal neuroblastoma'
    if re.search(r'ovarian cancer|malignant ovarian serous tumor',disease) and re.search(r'ovary|fallopian tube',q): return 'ovarian cancer'
    if re.search(r'oropharynx squamous cell carcinoma',disease) and 'oropharynx' in q: return 'oropharyngeal cancer'
    if re.search(r'oral cavity squamous cell carcinoma',disease) and 'oral cavity' in q: return 'oral cavity cancer'
    return None
eligible=['breast cancer','lung cancer','colorectal cancer','ovarian cancer','prostate cancer','adrenal neuroblastoma','gastric cancer','oropharyngeal cancer','oral cavity cancer']
# donor FAP-only module scores from Q2; one row per donor/group after averaging anatomical variants
mod=pd.read_parquet(Q2/'01_donor_state_module_compartment_collapsed.parquet'); mod=mod[(mod.context=='primary_cancer')&(mod.state=='FAP+DPP4-')].copy()
# Need disease for group assignment: recover donor/tissue disease map from state counts
mapd=pd.read_csv(R/'30_FAP_DPP4_state_by_donor.csv.gz',usecols=['donor_key','tissue_general','tissue','disease']).drop_duplicates()
mod=mod.merge(mapd,on=['donor_key','tissue_general','tissue'],how='left'); mod['cancer_group']=mod.apply(lambda r:cancer_group(r.tissue_general,r.tissue,r.disease),axis=1); mod=mod[mod.cancer_group.isin(eligible)]
modu=mod.groupby(['cancer_group','donor_key','collection_id','module']).module_gene_z.median().reset_index(); wide=modu.pivot_table(index=['cancer_group','donor_key','collection_id'],columns='module',values='module_gene_z',aggfunc='first').reset_index()
core=['ECM_matrix','TGFb_LRRC15','myofibroblast','COL15A1_like','inflammatory','detoxification']; accessory=['antigen_presentation','PI16_homeostatic']; wide['core_remodeling_score']=wide[core].mean(axis=1); wide['accessory_score']=wide[accessory].mean(axis=1)
# donor FAP state fractions from raw state counts, aggregated to cancer group
a=pd.read_csv(R/'30_FAP_DPP4_state_by_donor.csv.gz'); fib=pd.read_csv(C/'fibroblast_cell_types.tsv',sep='\t'); primary=set(fib.loc[fib.analysis_group.eq('primary_fibroblast'),'label'].astype(str)); ex=pd.read_csv(C/'annotation_context_exclusions.tsv',sep='\t'); dm=pd.read_csv(C/'disease_context_map.tsv',sep='\t')[['label','context']]
a=a[a.cell_type.astype(str).isin(primary)].copy()
for x in ex.itertuples(index=False): a=a[~((a.tissue_general.astype(str)==str(x.tissue_general))&(a.cell_type.astype(str)==str(x.cell_type)))]
a=a.merge(dm,left_on='disease',right_on='label',how='left').drop(columns='label'); a=a[a.context.eq('primary_cancer')]; a['cancer_group']=a.apply(lambda r:cancer_group(r.tissue_general,r.tissue,r.disease),axis=1); a=a[a.cancer_group.isin(eligible)]
g=a.groupby(['cancer_group','donor_key','collection_id','state']).state_n_cells.sum().reset_index(); st=g.pivot_table(index=['cancer_group','donor_key','collection_id'],columns='state',values='state_n_cells',aggfunc='sum',fill_value=0).reset_index()
for c in ['FAP-DPP4-','FAP-DPP4+','FAP+DPP4-','FAP+DPP4+']:
    if c not in st: st[c]=0
st['n_cells']=st[['FAP-DPP4-','FAP-DPP4+','FAP+DPP4-','FAP+DPP4+']].sum(axis=1); st['FAP_fraction']=(st['FAP+DPP4-']+st['FAP+DPP4+'])/st.n_cells; st['DPP4_fraction']=(st['FAP-DPP4+']+st['FAP+DPP4+'])/st.n_cells; st['FAP_minus_DPP4']=st.FAP_fraction-st.DPP4_fraction
# donor FAP expression from raw state-gene table across same cancer group
e=pd.read_parquet(R/'30b_FAP_DPP4_state_gene_pseudobulk.parquet',filters=[('gene','==','FAP')]); e=e[e.cell_type.astype(str).isin(primary)].copy()
for x in ex.itertuples(index=False): e=e[~((e.tissue_general.astype(str)==str(x.tissue_general))&(e.cell_type.astype(str)==str(x.cell_type)))]
e=e.merge(dm,left_on='disease',right_on='label',how='left').drop(columns='label'); e=e[e.context.eq('primary_cancer')]; e['cancer_group']=e.apply(lambda r:cancer_group(r.tissue_general,r.tissue,r.disease),axis=1); e=e[e.cancer_group.isin(eligible)]
e=e.groupby(['cancer_group','donor_key','collection_id']).agg(FAP_raw=('gene_raw_sum','sum'),lib_raw=('state_total_library_raw','sum')).reset_index(); e['FAP_log2_cpm']=np.log2(1+1e6*e.FAP_raw/e.lib_raw.replace(0,np.nan))
# merge
d=wide.merge(st,on=['cancer_group','donor_key','collection_id'],how='inner').merge(e[['cancer_group','donor_key','collection_id','FAP_log2_cpm']],on=['cancer_group','donor_key','collection_id'],how='inner')
# within-cancer standardized abundance/state for pooled cancer-type-adjusted correlation
for col in ['FAP_fraction','FAP_log2_cpm','core_remodeling_score','accessory_score']:
    d[col+'_within_z']=d.groupby('cancer_group')[col].transform(lambda v:(v-v.mean())/v.std(ddof=0) if v.std(ddof=0)>0 else 0)
# per-cancer correlations
rows=[]
for cg,gp in d.groupby('cancer_group'):
    for x in ['FAP_fraction','FAP_log2_cpm']:
      for y in ['core_remodeling_score','accessory_score']:
        if len(gp)>=5:
            rho,p=spearmanr(gp[x],gp[y]); rows.append(dict(cancer_group=cg,x=x,y=y,n=len(gp),rho=rho,p=p))
corr=pd.DataFrame(rows); corr.to_csv(OUT/'01_per_cancer_abundance_state_correlations.tsv',sep='\t',index=False)
# pooled within-cancer standardized correlation
summary=[]
for x in ['FAP_fraction','FAP_log2_cpm']:
  for y in ['core_remodeling_score','accessory_score']:
    rho,p=spearmanr(d[x+'_within_z'],d[y+'_within_z']); summary.append(dict(x=x,y=y,n=len(d),rho_within_cancer_adjusted=rho,p=p))
sumdf=pd.DataFrame(summary); sumdf.to_csv(OUT/'02_pooled_within_cancer_correlations.tsv',sep='\t',index=False)
# PCA of eight module scores to quantify state dimensionality
mods8=core+accessory; complete=d.dropna(subset=mods8).copy(); X=StandardScaler().fit_transform(complete[mods8]); pca=PCA(n_components=min(8,len(mods8))).fit(X); scores=pca.transform(X); complete['PC1']=scores[:,0]; complete['PC2']=scores[:,1]; complete.to_parquet(OUT/'03_donor_abundance_state_matrix.parquet',index=False); pd.DataFrame({'PC':np.arange(1,len(pca.explained_variance_ratio_)+1),'variance_ratio':pca.explained_variance_ratio_}).to_csv(OUT/'03b_PCA_variance.tsv',sep='\t',index=False); pd.DataFrame(pca.components_.T,index=mods8,columns=[f'PC{i}' for i in range(1,len(mods8)+1)]).to_csv(OUT/'03c_PCA_loadings.tsv',sep='\t')
print('DONORS',len(d)); print(d.groupby('cancer_group').size().sort_values(ascending=False).to_string())
print('\nPOOLED WITHIN-CANCER CORRELATIONS'); print(sumdf.to_string(index=False))
print('\nPER-CANCER FAP LOG2CPM vs CORE'); print(corr[(corr.x=='FAP_log2_cpm')&(corr.y=='core_remodeling_score')].sort_values('rho').to_string(index=False))
print('\nPCA variance',pca.explained_variance_ratio_[:4])
print('PC1 loadings'); print(pd.Series(pca.components_[0],index=mods8).sort_values(ascending=False).to_string())
print('PC2 loadings'); print(pd.Series(pca.components_[1],index=mods8).sort_values(ascending=False).to_string())
