from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
OUT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas/targetability_v1/Q4_malignant_vs_benign')
a=pd.read_parquet(OUT/'01_same_collection_lung_FAPonly_module_profiles.parquet'); core=pd.read_parquet(OUT/'05_same_collection_lung_donor_core_scores.parquet'); eff=pd.read_csv(OUT/'07_module_context_effects_bootstrap.tsv',sep='\t')
med=pd.read_csv(OUT/'04_context_median_module_profiles.tsv',sep='\t',index_col=0)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':8.7,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'svg.fonttype':'none','pdf.fonttype':42})
mods=['ECM_matrix','TGFb_LRRC15','myofibroblast','COL15A1_like','inflammatory','detoxification','antigen_presentation','PI16_homeostatic']; mlab={'ECM_matrix':'ECM/matrix','TGFb_LRRC15':'TGF-β/LRRC15','myofibroblast':'Myofibroblast','COL15A1_like':'COL15A1-like','inflammatory':'Inflammatory','detoxification':'Detoxification','antigen_presentation':'Antigen presentation','PI16_homeostatic':'PI16/homeostatic'}
ctx=['healthy','fibrosis_remodeling','inflammatory_infectious','primary_cancer']; clab={'healthy':'Healthy','fibrosis_remodeling':'Fibrosis/remodeling','inflammatory_infectious':'Inflammatory/infectious','primary_cancer':'Lung cancer'}; colors={'healthy':'#4C78A8','fibrosis_remodeling':'#7A5195','inflammatory_infectious':'#F28E2B','primary_cancer':'#C23B3B'}
fig=plt.figure(figsize=(7.08,7.3)); gs=GridSpec(2,2,figure=fig,height_ratios=[1.08,1],hspace=.43,wspace=.45)
# A heatmap
ax=fig.add_subplot(gs[0,:]); arr=med.reindex(ctx)[mods].values; vmax=max(1.0,np.nanmax(np.abs(arr))); im=ax.imshow(arr,aspect='auto',cmap='RdBu_r',vmin=-vmax,vmax=vmax)
ax.set_xticks(range(len(mods)),[mlab[m] for m in mods],rotation=32,ha='right'); ax.set_yticks(range(len(ctx)),[clab[c] for c in ctx]);
for i in range(len(ctx)):
 for j in range(len(mods)):
  v=arr[i,j]; ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=6.5,color='white' if abs(v)>.65 else 'black')
cb=fig.colorbar(im,ax=ax,pad=.012,fraction=.025); cb.set_label('Median FAP-only module score\n(exact-tissue gene-z units)'); ax.set_title('A  Same-study lung FAP-only states differ across malignant and benign remodeling',loc='left',fontweight='bold')
# B core score distributions
ax=fig.add_subplot(gs[1,0]); pos=np.arange(4)
for i,c in enumerate(ctx):
    v=core[core.context==c].core_score.to_numpy(); ax.boxplot([v],positions=[i],widths=.5,patch_artist=True,showfliers=False,medianprops={'color':'black','lw':1.1},boxprops={'facecolor':colors[c],'alpha':.25,'edgecolor':colors[c]},whiskerprops={'color':colors[c]},capprops={'color':colors[c]}); rng=np.random.default_rng(500+i); ax.scatter(i+rng.normal(0,.06,len(v)),v,s=7,color=colors[c],alpha=.45,linewidth=0); ax.text(i,ax.get_ylim()[1] if False else v.max()+.08,f'n={len(v)}',ha='center',fontsize=6.5)
ax.axhline(0,color='#bbb',lw=.6); ax.set_xticks(pos,['Healthy','Fibrosis','Inflammatory','Cancer'],rotation=20,ha='right'); ax.set_ylabel('Six-module remodeling-core score'); ax.set_title('B  Cancer carries the strongest remodeling core',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
# helper forest
def forest(ax,contrast,title):
    z=eff[eff.contrast==contrast].set_index('module').loc[mods].iloc[::-1].reset_index(); y=np.arange(len(z)); sig=z.q_within_contrast<.05; cols=np.where(sig,'#C23B3B','#777777'); ax.axvline(0,color='#999',lw=.7)
    for i,r in z.iterrows(): ax.errorbar(r.HL_shift,i,xerr=[[r.HL_shift-r.ci95_low],[r.ci95_high-r.HL_shift]],fmt='o',ms=3.7,color=cols[i],ecolor=cols[i],elinewidth=1,capsize=2)
    ax.set_yticks(y,[mlab[m] for m in z.module]); ax.set_xlabel('Hodges–Lehmann shift\n(cancer − reference)'); ax.set_title(title,loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
# create inset style two forests? replace B width? We have only one slot right, split vertically using inset axes
container=fig.add_subplot(gs[1,1]); container.axis('off'); b=container.get_position();
ax1=fig.add_axes([b.x0,b.y0+b.height*.58,b.width,b.height*.42]); forest(ax1,'cancer_minus_fibrosis_remodeling','C  Cancer versus fibrosis'); ax1.set_xlabel(''); ax1.tick_params(axis='x',labelbottom=False)
ax2=fig.add_axes([b.x0,b.y0,b.width,b.height*.40]); forest(ax2,'cancer_minus_inflammatory_infectious','D  Cancer versus inflammation')
fig.suptitle('High FAP in benign lung remodeling does not reproduce the full cancer-associated fibroblast state',fontsize=10.4,fontweight='bold',y=.995)
fig.subplots_adjust(top=.94,bottom=.09,left=.18,right=.96)
for ext in ['svg','pdf','png']: fig.savefig(OUT/f'Figure_Malignant_vs_Benign_Q4_FINAL.{ext}',dpi=400 if ext=='png' else None,bbox_inches='tight')
print(OUT/'Figure_Malignant_vs_Benign_Q4_FINAL.svg')
