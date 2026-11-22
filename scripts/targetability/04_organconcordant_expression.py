from pathlib import Path
import re,numpy as np,pandas as pd
ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas'); R=ROOT/'single_cell/results'; C=ROOT/'single_cell/config'; OUT=ROOT/'targetability_v1'
rules={'breast':r'breast','lung':r'lung|small cell lung|non-small cell lung|pleomorphic carcinoma','colorectum':r'colon adenocarcinoma|colorectal cancer','sigmoid colon':r'colorectal cancer|colon adenocarcinoma','stomach':r'gastric cancer','oropharynx':r'oropharynx squamous cell carcinoma','kidney':r'renal.*carcinoma|clear cell renal','prostate gland':r'prostatic','skin epidermis':r'melanoma|basal cell carcinoma|cutaneous squamous','liver':r'liver cancer|hepatocellular carcinoma|cholangiocarcinoma','uvea':r'uveal melanoma','adrenal gland':r'neuroblastoma|adrenocortical'}
def ok(t,d):
 p=rules.get(str(t)); return bool(p and re.search(p,str(d),re.I))
g=pd.read_parquet(R/'30b_FAP_DPP4_state_gene_pseudobulk.parquet',filters=[('gene','in',['FAP','DPP4'])]); fib=pd.read_csv(C/'fibroblast_cell_types.tsv',sep='\t'); primary=set(fib.loc[fib.analysis_group.eq('primary_fibroblast'),'label'].astype(str)); ex=pd.read_csv(C/'annotation_context_exclusions.tsv',sep='\t'); dm=pd.read_csv(C/'disease_context_map.tsv',sep='\t')[['label','context']]
g=g[g.cell_type.astype(str).isin(primary)].copy()
for x in ex.itertuples(index=False): g=g[~((g.tissue_general.astype(str)==str(x.tissue_general))&(g.cell_type.astype(str)==str(x.cell_type)))]
g=g.merge(dm,left_on='disease',right_on='label',how='left').drop(columns='label'); g=g[(g.context!='primary_cancer')|g.apply(lambda r:ok(r.tissue,r.disease),axis=1)]
base=['donor_key','collection_id','tissue_general','tissue','context','gene']; a=g.groupby(base).agg(gene_raw_sum=('gene_raw_sum','sum'),total_library_raw=('state_total_library_raw','sum')).reset_index(); a['log2_cpm']=np.log2(1+1e6*a.gene_raw_sum/a.total_library_raw.replace(0,np.nan)); a.to_parquet(OUT/'10c_donor_exacttissue_organconcordant_FAP_DPP4_log2cpm.parquet',index=False)
