from pathlib import Path
import pandas as pd,numpy as np,matplotlib.pyplot as plt
ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas'); R=ROOT/'single_cell/results'; F=ROOT/'single_cell/figures'
GENES=['FAP','DPP4','DPP8','DPP9','DPP6','DPP10','PREP','PREPL','APEH']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':.8,'svg.fonttype':'none','pdf.fonttype':42})
g=pd.read_csv(R/'01_global_cell_type.csv.gz'); t=pd.read_csv(R/'04_tissue_cell_type.csv.gz')
gr=g[g.robust_support==True].copy(); tr=t[t.robust_support==True].copy()
counts=pd.DataFrame({'gene':GENES,'robust_global_cell_types':[gr[gr.gene==x].cell_type.nunique() for x in GENES],'robust_tissue_cell_strata':[tr[tr.gene==x].shape[0] for x in GENES]})
# top 3 robust global cell types per gene; union in gene order, de-duplicated
rows=[]
for gene in GENES:
 q=gr[gr.gene==gene].sort_values('median_donor_log2_cpm',ascending=False).head(3)
 for ct in q.cell_type:
  if ct not in rows: rows.append(ct)
mat=gr[gr.cell_type.isin(rows)].pivot_table(index='cell_type',columns='gene',values='median_donor_log2_cpm',aggfunc='first').reindex(index=rows,columns=GENES)
# use 0 only where a robust combination is absent, visually gray masked instead of biological zero
arr=mat.to_numpy(float); mask=np.isnan(arr)
fig=plt.figure(figsize=(9.6,8.3)); gs=fig.add_gridspec(2,1,height_ratios=[.85,2.8],hspace=.35)
ax=fig.add_subplot(gs[0]); x=np.arange(len(GENES)); w=.38
ax.bar(x-w/2,counts.robust_global_cell_types,w,label='Robust global cell types',color='#4c78a8')
# scale tissue-cell strata to secondary axis because counts are much larger
ax2=ax.twinx(); ax2.bar(x+w/2,counts.robust_tissue_cell_strata,w,label='Robust tissue×cell-type strata',color='#f58518')
ax.set_xticks(x,GENES); ax.set_ylabel('No. robust global cell types'); ax2.set_ylabel('No. robust tissue×cell-type strata')
ax.set_title('A  Broad robust support across the nine-gene S9 atlas',loc='left',fontweight='bold',fontsize=11)
ax.spines['top'].set_visible(False); ax2.spines['top'].set_visible(False)
handles=[plt.Rectangle((0,0),1,1,color='#4c78a8'),plt.Rectangle((0,0),1,1,color='#f58518')]
ax.legend(handles,['Global cell types','Tissue×cell-type strata'],frameon=False,ncol=2,loc='upper left')

ax=fig.add_subplot(gs[1]); cmap=plt.cm.viridis.copy(); cmap.set_bad('#eeeeee')
im=ax.imshow(np.ma.masked_invalid(arr),aspect='auto',cmap=cmap,vmin=0,vmax=np.nanpercentile(arr,98))
ax.set_xticks(np.arange(len(GENES)),GENES,fontweight='bold'); ax.set_yticks(np.arange(len(rows)),rows,fontsize=7.4)
ax.set_title('B  Top robust cellular niches for each S9-family gene',loc='left',fontweight='bold',fontsize=11,pad=8)
cb=fig.colorbar(im,ax=ax,fraction=.025,pad=.02); cb.set_label('Median donor log2 CPM')
# mark each gene's own top three cells with a small dot at the corresponding heatmap cell
for j,gene in enumerate(GENES):
 tops=gr[gr.gene==gene].sort_values('median_donor_log2_cpm',ascending=False).head(3).cell_type.tolist()
 for ct in tops:
  if ct in rows:
   ax.scatter(j,rows.index(ct),s=13,facecolors='none',edgecolors='white',linewidth=.8)
ax.text(0,-.075,'Rows are the union of the top three robust global cell types for each gene; white circles mark the selecting gene–cell-type combinations. Gray indicates no robust combination.',transform=ax.transAxes,fontsize=7.8,color='#444444')
fig.suptitle('Donor-aware CELLxGENE atlas reveals distinct cellular specialization across the human S9 family',x=.03,ha='left',fontsize=13,fontweight='bold')
fig.subplots_adjust(left=.30,right=.92,top=.93,bottom=.08)
for ext in ['svg','pdf']: fig.savefig(F/f'Figure1_Human_S9_global_atlas_FINAL.{ext}',bbox_inches='tight')
fig.savefig(F/'Figure1_Human_S9_global_atlas_FINAL.png',dpi=600,bbox_inches='tight'); plt.close(fig)
counts.to_csv(ROOT/'manuscript/TABLE_S9_ROBUST_ATLAS_SUPPORT_FINAL.csv',index=False)
mat.to_csv(ROOT/'manuscript/TABLE_S9_TOP_ROBUST_CELLTYPE_EXPRESSION_FINAL.csv')
print(counts.to_string(index=False)); print('rows',len(rows))