from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import spearmanr
ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas'); Q1=ROOT/'targetability_v1'; Q6=ROOT/'targetability_v1/Q6_integrated_targetability'; OUT=Q6
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.4,'axes.titlesize':8.3,'axes.labelsize':7.5,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'svg.fonttype':'none','pdf.fonttype':42})
healthy='#3973AC'; cancer='#C23B3B'; fibrosis='#7A5195'; infl='#E07A1F'; gray='#777777'
primary=['colorectum','lung','oropharynx','kidney']; labels={'colorectum':'Colorectum','lung':'Lung','oropharynx':'Oropharynx','kidney':'Kidney'}
w=pd.read_parquet(Q1/'10_donor_exacttissue_organconcordant_cancer_vs_healthy.parquet')
ex=pd.read_parquet(Q1/'10c_donor_exacttissue_organconcordant_FAP_DPP4_log2cpm.parquet'); ex=ex[(ex.gene=='FAP')&ex.tissue.isin(primary)&ex.context.isin(['healthy','primary_cancer'])]
site=pd.read_csv(Q1/'11_organconcordant_cancer_vs_healthy_exact_tissue.tsv',sep='\t'); site=site[site.robust5each].sort_values('delta_FAP')
allfrac=pd.read_parquet(Q1/'00_donor_exacttissue_context_primary_fibroblast.parquet'); exall=pd.read_parquet(Q1/'10c_donor_exacttissue_organconcordant_FAP_DPP4_log2cpm.parquet'); exall_full=pd.read_parquet(Q1/'06_donor_exacttissue_context_FAP_DPP4_log2cpm.parquet')
integ=pd.read_csv(Q6/'01_integrated_targetability_components.tsv',sep='\t'); short={'oropharyngeal cancer':'Oropharyngeal','lung cancer':'Lung','colorectal cancer':'Colorectal','adrenal neuroblastoma':'Adrenal NB','breast cancer':'Breast','gastric cancer':'Gastric','prostate cancer':'Prostate'}
fig=plt.figure(figsize=(7.08,9.35)); gs=GridSpec(4,2,figure=fig,height_ratios=[1.02,1.15,1.02,1.15],hspace=.55,wspace=.42)
# A matched fraction
ax=fig.add_subplot(gs[0,0])
for j,t in enumerate(primary):
 q=w[(w.tissue==t)&w.context.isin(['healthy','primary_cancer'])].pivot_table(index=['donor_key','collection_id'],columns='context',values='FAP_pos_fraction',aggfunc='first')
 if not {'healthy','primary_cancer'}.issubset(q.columns): continue
 q=q.dropna(); x0=j*2.0; x1=x0+.62
 for _,r in q.iterrows(): ax.plot([x0,x1],[r.healthy,r.primary_cancer],color='#C9C9C9',lw=.55,zorder=1)
 ax.scatter(np.full(len(q),x0),q.healthy,s=10,color=healthy,edgecolor='white',linewidth=.2,zorder=2); ax.scatter(np.full(len(q),x1),q.primary_cancer,s=10,color=cancer,edgecolor='white',linewidth=.2,zorder=2)
 ax.plot([x0-.1,x0+.1],[q.healthy.median()]*2,color='black',lw=1); ax.plot([x1-.1,x1+.1],[q.primary_cancer.median()]*2,color='black',lw=1)
ax.set_xticks([j*2+.31 for j in range(4)],[labels[t] for t in primary],rotation=20,ha='right'); ax.set_ylabel('FAP+ fibroblast fraction'); ax.set_ylim(-.04,1.02); ax.set_title('A  Matched tumor–normal FAP detection',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
# B matched abundance
ax=fig.add_subplot(gs[0,1])
for j,t in enumerate(primary):
 q=ex[ex.tissue==t].pivot_table(index=['donor_key','collection_id'],columns='context',values='log2_cpm',aggfunc='first')
 if not {'healthy','primary_cancer'}.issubset(q.columns): continue
 q=q.dropna(); x0=j*2.0; x1=x0+.62
 for _,r in q.iterrows(): ax.plot([x0,x1],[r.healthy,r.primary_cancer],color='#C9C9C9',lw=.55,zorder=1)
 ax.scatter(np.full(len(q),x0),q.healthy,s=10,color=healthy,edgecolor='white',linewidth=.2,zorder=2); ax.scatter(np.full(len(q),x1),q.primary_cancer,s=10,color=cancer,edgecolor='white',linewidth=.2,zorder=2)
 ax.plot([x0-.1,x0+.1],[q.healthy.median()]*2,color='black',lw=1); ax.plot([x1-.1,x1+.1],[q.primary_cancer.median()]*2,color='black',lw=1)
ax.set_xticks([j*2+.31 for j in range(4)],[labels[t] for t in primary],rotation=20,ha='right'); ax.set_ylabel('FAP log$_2$(CPM + 1)'); ax.set_title('B  Matched tumor–normal FAP abundance',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
# C broader landscape full width
ax=fig.add_subplot(gs[1,:]); yy=np.arange(len(site)); ax.hlines(yy,site.median_FAP_healthy,site.median_FAP_primary_cancer,color='#B7B7B7',lw=.9); ax.scatter(site.median_FAP_healthy,yy,s=22,color=healthy,label='Healthy'); ax.scatter(site.median_FAP_primary_cancer,yy,s=22,color=cancer,label='Organ-concordant cancer'); ax.set_yticks(yy,[x.replace('_',' ').title() for x in site.tissue]); ax.set_xlabel('Median FAP+ fraction across primary fibroblasts'); ax.set_xlim(-.03,1); ax.set_title('C  Exact-tissue cancer–healthy FAP landscape (≥5 donors per group)',loc='left',fontweight='bold'); ax.legend(frameon=False,ncol=2,loc='lower right',fontsize=6.5); ax.spines[['top','right']].set_visible(False)
# D/E lung benign context
ctx=['healthy','fibrosis_remodeling','inflammatory_infectious','primary_cancer']; clab=['Healthy','Fibrosis','Inflammatory/\ninfectious','Lung cancer']; colors=[healthy,fibrosis,infl,cancer]
ax=fig.add_subplot(gs[2,0])
for k,cx in enumerate(ctx):
 v=(w[(w.tissue=='lung')&(w.context=='primary_cancer')].FAP_pos_fraction if cx=='primary_cancer' else allfrac[(allfrac.tissue=='lung')&(allfrac.context==cx)].FAP_pos_fraction).dropna().to_numpy(); ax.boxplot([v],positions=[k],widths=.48,patch_artist=True,showfliers=False,medianprops={'color':'black','linewidth':1},boxprops={'facecolor':colors[k],'alpha':.25,'edgecolor':colors[k]},whiskerprops={'color':colors[k]},capprops={'color':colors[k]}); rng=np.random.default_rng(100+k); ax.scatter(k+rng.normal(0,.06,len(v)),v,s=4.5,color=colors[k],alpha=.3,linewidth=0)
ax.set_xticks(range(4),clab,rotation=15,ha='right'); ax.set_ylabel('FAP+ fibroblast fraction'); ax.set_ylim(-.04,1.02); ax.set_title('D  Lung FAP detection across contexts',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
ax=fig.add_subplot(gs[2,1])
for k,cx in enumerate(ctx):
 src=exall if cx=='primary_cancer' else exall_full; v=src[(src.tissue=='lung')&(src.gene=='FAP')&(src.context==cx)].log2_cpm.dropna().to_numpy(); ax.boxplot([v],positions=[k],widths=.48,patch_artist=True,showfliers=False,medianprops={'color':'black','linewidth':1},boxprops={'facecolor':colors[k],'alpha':.25,'edgecolor':colors[k]},whiskerprops={'color':colors[k]},capprops={'color':colors[k]}); rng=np.random.default_rng(200+k); ax.scatter(k+rng.normal(0,.06,len(v)),v,s=4.5,color=colors[k],alpha=.3,linewidth=0)
ax.set_xticks(range(4),clab,rotation=15,ha='right'); ax.set_ylabel('FAP log$_2$(CPM + 1)'); ax.set_title('E  Lung FAP abundance across contexts',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
# F integrated map full width
ax=fig.add_subplot(gs[3,:]); sc=ax.scatter(integ.cancer_minus_healthy_FAP_fraction,integ.core_remodeling_effect,s=28+1.8*integ.state_paired_donors,c=integ.fibroblast_selectivity_log2cpm,cmap='viridis',edgecolor='black',linewidth=.45,alpha=.9)
offsets={'oropharyngeal cancer':(-.10,.012),'prostate cancer':(.010,-.030),'lung cancer':(.010,.028),'colorectal cancer':(.020,-.022),'adrenal neuroblastoma':(.010,-.035),'breast cancer':(.010,.010),'gastric cancer':(.010,.010)}
for _,r in integ.iterrows():
 dx,dy=offsets.get(r.cancer_group,(.01,.01)); ax.text(r.cancer_minus_healthy_FAP_fraction+dx,r.core_remodeling_effect+dy,short[r.cancer_group],fontsize=6.5)
rho,p=spearmanr(integ.cancer_minus_healthy_FAP_fraction,integ.core_remodeling_effect); ax.text(.62,.96,f'Across 7 groups: ρ={rho:.2f}, P={p:.3f}',transform=ax.transAxes,va='top',fontsize=6.5); ax.set_xlabel('Cancer − healthy median FAP+ fibroblast fraction'); ax.set_ylabel('Mean conserved six-module\nFAP-only remodeling effect'); ax.set_title('F  Target enrichment, remodeling strength and cellular selectivity are distinct dimensions',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False); cb=fig.colorbar(sc,ax=ax,pad=.012,fraction=.025); cb.set_label('Fibroblast FAP selectivity\n(log2 CPM difference vs other cells)',fontsize=6.5); cb.ax.tick_params(labelsize=6)
for n in [10,50,100]: ax.scatter([],[],s=28+1.8*n,color='lightgray',edgecolor='black',label=f'{n} donors')
ax.legend(title='State-program support',frameon=False,loc='lower right',fontsize=6,title_fontsize=6.2,ncol=1)
fig.suptitle('FAP targetability is organ-dependent and separates target abundance from stromal biology',fontsize=10.2,fontweight='bold',y=.997)
fig.subplots_adjust(top=.965,bottom=.055,left=.18,right=.96)
for ext in ['svg','pdf','png']: fig.savefig(OUT/f'Figure7_FAP_Targetability_FINAL.{ext}',dpi=400 if ext=='png' else None,bbox_inches='tight')
print(OUT/'Figure7_FAP_Targetability_FINAL.svg')
