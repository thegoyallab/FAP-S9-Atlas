from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
OUT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas/targetability_v1/Q2_cross_cancer')
r=pd.read_csv(OUT/'05_canonical_cancer_module_effects.tsv',sep='\t'); m=pd.read_csv(OUT/'06_canonical_cross_cancer_meta.tsv',sep='\t')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'svg.fonttype':'none','pdf.fonttype':42})
mods=['ECM_matrix','TGFb_LRRC15','myofibroblast','COL15A1_like','inflammatory','detoxification','antigen_presentation','PI16_homeostatic']
mlab={'ECM_matrix':'ECM/matrix','TGFb_LRRC15':'TGF-β/LRRC15','myofibroblast':'Myofibroblast','COL15A1_like':'COL15A1-like','inflammatory':'Inflammatory','detoxification':'Detoxification','antigen_presentation':'Antigen presentation','PI16_homeostatic':'PI16/homeostatic'}
order=['breast cancer','lung cancer','colorectal cancer','ovarian cancer','prostate cancer','adrenal neuroblastoma','gastric cancer','oropharyngeal cancer','oral cavity cancer']
clab={'breast cancer':'Breast','lung cancer':'Lung','colorectal cancer':'Colorectal','ovarian cancer':'Ovarian','prostate cancer':'Prostate','adrenal neuroblastoma':'Adrenal neuroblastoma','gastric cancer':'Gastric','oropharyngeal cancer':'Oropharyngeal','oral cavity cancer':'Oral cavity'}
mat=r.pivot(index='cancer_group',columns='module',values='median_delta').reindex(order)[mods]
fig=plt.figure(figsize=(7.08,7.25)); gs=GridSpec(2,2,figure=fig,height_ratios=[1.35,1],width_ratios=[1.45,.75],hspace=.42,wspace=.38)
# A heatmap
ax=fig.add_subplot(gs[0,:]); vmax=max(1.6,float(np.nanmax(np.abs(mat.values)))); im=ax.imshow(mat.values,aspect='auto',cmap='RdBu_r',vmin=-vmax,vmax=vmax)
ax.set_xticks(range(len(mods)),[mlab[x] for x in mods],rotation=35,ha='right'); ax.set_yticks(range(len(order)),[clab[x] for x in order])
for i in range(mat.shape[0]):
 for j in range(mat.shape[1]):
  v=mat.iloc[i,j]; ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=6.3,color='white' if abs(v)>.75 else 'black')
cb=fig.colorbar(im,ax=ax,pad=.012,fraction=.025); cb.set_label('Median paired FAP-only − DPP4-only\nmodule score (gene-z units)')
ax.set_title('A  A conserved FAP-dominant fibroblast program recurs across solid cancers',loc='left',fontweight='bold')
# B forest meta
ax=fig.add_subplot(gs[1,0]); mm=m.set_index('module').loc[mods].iloc[::-1].reset_index(); y=np.arange(len(mm));
ax.axvline(0,color='#888',lw=.8); ax.errorbar(mm.random_effect_mean,y,xerr=[mm.random_effect_mean-mm.ci95_low,mm.ci95_high-mm.random_effect_mean],fmt='o',ms=4,color='#C23B3B',ecolor='#555',elinewidth=1,capsize=2)
ax.set_yticks(y,[mlab[x] for x in mm.module]); ax.set_xlabel('Random-effects mean paired difference (gene-z units)'); ax.set_title('B  Cross-cancer random-effects meta-analysis',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
# C support
ax=fig.add_subplot(gs[1,1]); support=r.groupby('cancer_group').agg(n=('n_paired_donors','max'),collections=('n_collections','max')).reindex(order).iloc[::-1]; yy=np.arange(len(support)); ax.barh(yy,support.n,color='#666',alpha=.75); ax.set_yticks(yy,[clab[x] for x in support.index]); ax.set_xlabel('Paired donors'); ax.set_title('C  Cancer-group support',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
for i,(n,c) in enumerate(zip(support.n,support.collections)): ax.text(n+2,i,f'{int(c)} coll.',va='center',fontsize=6.5)
ax.set_xlim(0,max(support.n)*1.22)
fig.suptitle('Within-donor FAP-only versus DPP4-only contrasts identify a conserved cancer fibroblast program',fontsize=10.5,fontweight='bold',y=.995)
fig.subplots_adjust(top=.94,bottom=.09,left=.18,right=.95)
for ext in ['svg','pdf','png']: fig.savefig(OUT/f'Figure_CrossCancer_FAP_Program_Q2_FINAL.{ext}',dpi=400 if ext=='png' else None,bbox_inches='tight')
print(OUT/'Figure_CrossCancer_FAP_Program_Q2_FINAL.svg')
