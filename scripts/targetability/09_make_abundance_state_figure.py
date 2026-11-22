from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import spearmanr
OUT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas/targetability_v1/Q3_abundance_vs_state')
d=pd.read_parquet(OUT/'03_donor_abundance_state_matrix.parquet')
corr=pd.read_csv(OUT/'01_per_cancer_abundance_state_correlations.tsv',sep='\t')
var=pd.read_csv(OUT/'03b_PCA_variance.tsv',sep='\t'); load=pd.read_csv(OUT/'03c_PCA_loadings.tsv',sep='\t',index_col=0)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':8.4,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'svg.fonttype':'none','pdf.fonttype':42})
colors={'breast cancer':'#C84C4C','lung cancer':'#4C78A8','colorectal cancer':'#59A14F','ovarian cancer':'#B279A2','prostate cancer':'#F28E2B','gastric cancer':'#E15759','adrenal neuroblastoma':'#76B7B2','oral cavity cancer':'#EDC949','oropharyngeal cancer':'#9C755F'}
short={'breast cancer':'Breast','lung cancer':'Lung','colorectal cancer':'Colorectal','ovarian cancer':'Ovarian','prostate cancer':'Prostate','gastric cancer':'Gastric','adrenal neuroblastoma':'Adrenal NB','oral cavity cancer':'Oral cavity','oropharyngeal cancer':'Oropharyngeal'}
fig=plt.figure(figsize=(7.08,7.2)); gs=GridSpec(2,2,figure=fig,hspace=.42,wspace=.36)
# A abundance vs core within-cancer z
ax=fig.add_subplot(gs[0,0])
for cg,g in d.groupby('cancer_group'):
    ax.scatter(g['FAP_log2_cpm_within_z'],g['core_remodeling_score_within_z'],s=13,alpha=.5,color=colors.get(cg,'gray'),label=short.get(cg,cg),linewidth=0)
rho,p=spearmanr(d['FAP_log2_cpm_within_z'],d['core_remodeling_score_within_z'])
coef=np.polyfit(d['FAP_log2_cpm_within_z'],d['core_remodeling_score_within_z'],1); xx=np.linspace(d['FAP_log2_cpm_within_z'].min(),d['FAP_log2_cpm_within_z'].max(),100); ax.plot(xx,np.polyval(coef,xx),color='black',lw=1.1)
ax.text(.03,.96,f'ρ = {rho:.3f}\nP = {p:.3g}',transform=ax.transAxes,va='top',fontsize=7.5)
ax.set_xlabel('FAP abundance within cancer type (z)'); ax.set_ylabel('Core remodeling score within cancer type (z)'); ax.set_title('A  FAP abundance is a weak state proxy',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
# B FAP fraction vs core
ax=fig.add_subplot(gs[0,1])
for cg,g in d.groupby('cancer_group'):
    ax.scatter(g['FAP_fraction_within_z'],g['core_remodeling_score_within_z'],s=13,alpha=.5,color=colors.get(cg,'gray'),linewidth=0)
rho2,p2=spearmanr(d['FAP_fraction_within_z'],d['core_remodeling_score_within_z']); coef=np.polyfit(d['FAP_fraction_within_z'],d['core_remodeling_score_within_z'],1); xx=np.linspace(d['FAP_fraction_within_z'].min(),d['FAP_fraction_within_z'].max(),100); ax.plot(xx,np.polyval(coef,xx),color='black',lw=1.1); ax.text(.03,.96,f'ρ = {rho2:.3f}\nP = {p2:.3g}',transform=ax.transAxes,va='top',fontsize=7.5)
ax.set_xlabel('FAP+ fibroblast fraction within cancer type (z)'); ax.set_ylabel('Core remodeling score within cancer type (z)'); ax.set_title('B  FAP detection is a weak state proxy',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
# C per-cancer correlations FAP abundance vs core
ax=fig.add_subplot(gs[1,0]); c=corr[(corr.x=='FAP_log2_cpm')&(corr.y=='core_remodeling_score')].copy().sort_values('rho'); yy=np.arange(len(c)); ax.axvline(0,color='#888',lw=.8); ax.scatter(c.rho,yy,s=28,color=[colors.get(x,'gray') for x in c.cancer_group]); ax.set_yticks(yy,[short.get(x,x) for x in c.cancer_group]); ax.set_xlabel('Spearman ρ: FAP abundance vs core state'); ax.set_xlim(-.65,.65); ax.set_title('C  Abundance–state coupling differs by cancer',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
for i,(rho,pv) in enumerate(zip(c.rho,c.p)): ax.text(rho+(.025 if rho>=0 else -.025),i,f'{rho:.2f}',ha='left' if rho>=0 else 'right',va='center',fontsize=6.5)
# D PCA loadings showing multidimensional state
ax=fig.add_subplot(gs[1,1]); mods=load.index.tolist(); x=load['PC1']; y=load['PC2']; ax.axhline(0,color='#bbb',lw=.6); ax.axvline(0,color='#bbb',lw=.6); ax.scatter(x,y,s=32,color='#555')
offsets={'ECM_matrix':(.012,.035),'myofibroblast':(.012,-.035),'TGFb_LRRC15':(.014,.015),'antigen_presentation':(.014,.025),'inflammatory':(.014,.015),'COL15A1_like':(.014,-.015),'PI16_homeostatic':(.014,.018),'detoxification':(.014,-.018)}
for m,xx0,yy0 in zip(mods,x,y):
    dx,dy=offsets.get(m,(.015,.015)); ax.text(xx0+dx,yy0+dy,m.replace('_',' '),fontsize=6.2)
ax.set_xlabel(f'PC1 loading ({var.loc[var.PC==1,"variance_ratio"].iloc[0]*100:.1f}% variance)'); ax.set_ylabel(f'PC2 loading ({var.loc[var.PC==2,"variance_ratio"].iloc[0]*100:.1f}% variance)'); ax.set_title('D  FAP-positive state is multidimensional',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
# legend
handles=[]; labels=[]
for cg in colors:
    handles.append(plt.Line2D([0],[0],marker='o',linestyle='',markersize=4,color=colors[cg])); labels.append(short[cg])
fig.legend(handles,labels,loc='lower center',ncol=5,frameon=False,bbox_to_anchor=(.5,.005),fontsize=6.8)
fig.suptitle('FAP abundance does not uniquely specify the biological state of FAP-positive cancer stroma',fontsize=10.5,fontweight='bold',y=.995)
fig.subplots_adjust(top=.94,bottom=.13,left=.17,right=.98)
for ext in ['svg','pdf','png']: fig.savefig(OUT/f'Figure_Abundance_vs_State_Q3_FINAL.{ext}',dpi=400 if ext=='png' else None,bbox_inches='tight')
print(OUT/'Figure_Abundance_vs_State_Q3_FINAL.svg')
