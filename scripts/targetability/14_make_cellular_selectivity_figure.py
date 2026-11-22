from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
OUT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas/targetability_v1/Q5_cellular_selectivity')
r=pd.read_csv(OUT/'05_organconcordant_tissue_cellular_selectivity.tsv',sep='\t'); s=pd.read_csv(OUT/'07_selected_offcompartment_celltypes_cancer.tsv',sep='\t')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':8.8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'svg.fonttype':'none','pdf.fonttype':42})
# cancer tissue selectivity
z=r[(r.context=='primary_cancer')&(r.metric=='log2_cpm')].sort_values('median_paired_delta');
# selected cell labels order and display
labels=['endothelial cell of artery','endothelial cell of vascular tree','endothelial cell','vein endothelial cell','blood vessel endothelial cell','pericyte','vascular associated smooth muscle cell','malignant cell','epithelial cell','macrophage','monocyte','B cell','T cell','natural killer cell']
disp={'endothelial cell of artery':'Arterial endothelial','endothelial cell of vascular tree':'Vascular-tree endothelial','endothelial cell':'Endothelial','vein endothelial cell':'Venous endothelial','blood vessel endothelial cell':'Blood-vessel endothelial','pericyte':'Pericyte','vascular associated smooth muscle cell':'Vascular smooth muscle','malignant cell':'Malignant cell','epithelial cell':'Epithelial cell','macrophage':'Macrophage','monocyte':'Monocyte','B cell':'B cell','T cell':'T cell','natural killer cell':'NK cell'}
q=s[s.cell_type.isin(labels)].set_index('cell_type').reindex(labels).reset_index()
fig=plt.figure(figsize=(7.08,7.0)); gs=GridSpec(2,2,figure=fig,height_ratios=[1.25,1],width_ratios=[1.2,.8],hspace=.42,wspace=.38)
# A dumbbell tissue contexts
ax=fig.add_subplot(gs[0,:]); y=np.arange(len(z)); ax.hlines(y,z.median_other,z.median_primary_fibroblast,color='#B5B5B5',lw=1); ax.scatter(z.median_other,y,s=25,color='#777',label='All other annotated cells'); ax.scatter(z.median_primary_fibroblast,y,s=25,color='#C23B3B',label='Primary fibroblast compartment'); ax.set_yticks(y,[x.replace('_',' ').title() for x in z.tissue_general]); ax.set_xlabel('Median donor FAP log$_2$(CPM + 1)'); ax.set_title('A  FAP is concentrated in primary fibroblasts across most organ-concordant cancers',loc='left',fontweight='bold'); ax.legend(frameon=False,ncol=2,loc='lower right'); ax.spines[['top','right']].set_visible(False)
# B paired difference vs support
ax=fig.add_subplot(gs[1,0]); zz=z.sort_values('median_paired_delta'); yy=np.arange(len(zz)); sig=(zz.q<.05)&(zz.median_paired_delta>0); ax.axvline(0,color='#888',lw=.7); ax.scatter(zz.median_paired_delta,yy,s=np.clip(zz.n_paired_donors,10,80),c=np.where(sig,'#C23B3B','#777'),alpha=.85); ax.set_yticks(yy,[x.replace('_',' ').title() for x in zz.tissue_general]); ax.set_xlabel('Median paired FAP abundance difference\n(primary fibroblast − other cells)'); ax.set_title('B  Donor-paired cellular selectivity',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
# C selected exact off-compartment labels, bars + detection markers
ax=fig.add_subplot(gs[1,1]); q=q.iloc[::-1]; yy=np.arange(len(q)); ax.barh(yy,q.median_log2_cpm,color='#4C78A8',alpha=.75); ax.set_yticks(yy,[disp[x] for x in q.cell_type]); ax.set_xlabel('Median FAP log$_2$(CPM + 1)'); ax.set_title('C  Reproducible off-compartment signals',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
# annotate detection fraction for top vascular/malignant only
for i,row in q.iterrows():
    pass
for pos,(_,row) in enumerate(q.iterrows()):
    if row.median_log2_cpm>0:
        ax.text(row.median_log2_cpm+.1,pos,f'{row.median_pct_detected:.1f}%+',va='center',fontsize=6.2)
ax.set_xlim(0,max(6.4,float(q.median_log2_cpm.max()+1.2)))
fig.suptitle('FAP is predominantly fibroblast-associated in cancer but is not fibroblast-exclusive',fontsize=10.5,fontweight='bold',y=.995)
fig.subplots_adjust(top=.94,bottom=.09,left=.24,right=.98)
for ext in ['svg','pdf','png']: fig.savefig(OUT/f'Figure_Cellular_Selectivity_Q5_FINAL.{ext}',dpi=400 if ext=='png' else None,bbox_inches='tight')
print(OUT/'Figure_Cellular_Selectivity_Q5_FINAL.svg')
