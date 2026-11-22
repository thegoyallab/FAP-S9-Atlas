from pathlib import Path
import pandas as pd, numpy as np, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas'); M=ROOT/'manuscript'; R=ROOT/'single_cell/results'; F=ROOT/'single_cell/figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.linewidth':.8,'svg.fonttype':'none','pdf.fonttype':42})

# ---------- Figure 2 final249: all replicated contexts ----------
x=pd.read_csv(M/'TABLE_FAP_DPP4_REPLICATED_TISSUE_CONTEXTS_FINAL249.csv')
fig=plt.figure(figsize=(9.5,7.2)); gs=fig.add_gridspec(2,2,height_ratios=[1.15,1],hspace=.38,wspace=.42)
ax=fig.add_subplot(gs[0,:])
vals=x.median_FAP_minus_DPP4_fraction.to_numpy(); lim=max(abs(vals.min()),abs(vals.max()))
sc=ax.scatter(x.median_DPP4_pos_fraction,x.median_FAP_pos_fraction,c=vals,cmap='coolwarm',vmin=-lim,vmax=lim,s=24+2.4*np.sqrt(x.n_donors)*8,alpha=.84,edgecolor='white',linewidth=.45)
mx=max(x.median_DPP4_pos_fraction.max(),x.median_FAP_pos_fraction.max())*1.03
ax.plot([0,mx],[0,mx],ls='--',lw=.8,color='#999999'); ax.set_xlim(-.01,mx); ax.set_ylim(-.01,mx)
ax.set_xlabel('Median DPP4-positive donor fraction'); ax.set_ylabel('Median FAP-positive donor fraction'); ax.set_title('A  All replicated exact-tissue fibroblast contexts',loc='left',fontweight='bold',fontsize=11)
ax.spines[['top','right']].set_visible(False)
cb=fig.colorbar(sc,ax=ax,fraction=.027,pad=.02); cb.set_label('Median FAP − DPP4 fraction')
# predefined biologically important labels
sel=[('breast','normal','fibroblast','healthy breast'),('breast','breast cancer','fibroblast','breast cancer'),('lung','normal','fibroblast','healthy lung'),('lung','lung adenocarcinoma','bronchus fibroblast of lung','lung adenocarcinoma'),('lung','squamous cell lung carcinoma','bronchus fibroblast of lung','lung SCC')]
offsets={'healthy breast':(6,4),'breast cancer':(6,6),'healthy lung':(6,10),'lung adenocarcinoma':(6,-10),'lung SCC':(6,6)}
for tissue,dis,ct,label in sel:
    q=x[(x.tissue==tissue)&(x.disease==dis)&(x.cell_type==ct)]
    if len(q):
        r=q.iloc[0]; ax.annotate(label,(r.median_DPP4_pos_fraction,r.median_FAP_pos_fraction),xytext=offsets.get(label,(4,4)),textcoords='offset points',fontsize=7.2)
# extremes among >=5 donors: 8 FAP and 8 DPP4 dominant
q=x[x.n_donors>=5].copy(); top=q.nlargest(8,'median_FAP_minus_DPP4_fraction'); bot=q.nsmallest(8,'median_FAP_minus_DPP4_fraction'); e=pd.concat([bot,top]).drop_duplicates().sort_values('median_FAP_minus_DPP4_fraction');
def short(r):
    t=str(r.tissue).replace('mammary gland connective tissue','mammary CT').replace('upper outer quadrant of breast','breast UOQ').replace('skin of external ear','external-ear skin').replace('corneo-scleral junction','corneoscleral junction')
    d='' if r.disease=='normal' else ' | '+str(r.disease).replace('squamous cell lung carcinoma','lung SCC').replace('non-specific interstitial pneumonia','NSIP')
    c='' if r.cell_type=='fibroblast' else ' | '+str(r.cell_type).replace('alveolar type 1 fibroblast cell','AT1 fibroblast').replace('bronchus fibroblast of lung','bronchus fibroblast')
    return t+d+c
ax=fig.add_subplot(gs[1,0]); yy=np.arange(len(e)); colors=np.where(e.median_FAP_minus_DPP4_fraction>=0,'#D55E00','#0072B2')
ax.barh(yy,e.median_FAP_minus_DPP4_fraction,color=colors,height=.68); ax.axvline(0,color='#555555',lw=.8); ax.set_yticks(yy,[short(r) for r in e.itertuples()],fontsize=6.7); ax.set_xlabel('Median FAP − DPP4 fraction'); ax.set_title('B  Strongest replicated context biases',loc='left',fontweight='bold',fontsize=11); ax.spines[['top','right']].set_visible(False)
# DP distribution
ax=fig.add_subplot(gs[1,1]); dp=x.median_double_positive_fraction.dropna().to_numpy(); ax.hist(dp,bins=18,color='#777777',edgecolor='white',linewidth=.5)
med=np.median(dp); p90=np.quantile(dp,.90); ax.axvline(med,color='#333333',lw=1.2,label=f'median={med:.3f}'); ax.axvline(p90,color='#777777',lw=1.0,ls='--',label=f'90th percentile={p90:.3f}')
ax.set_xlabel('Median FAP+DPP4+ donor fraction'); ax.set_ylabel('No. replicated contexts'); ax.set_title('C  Double-positive fraction across contexts',loc='left',fontweight='bold',fontsize=11); ax.legend(frameon=False,fontsize=7.5); ax.spines[['top','right']].set_visible(False)
fig.suptitle('FAP and DPP4 partition across a broad, tissue-dependent fibroblast landscape',x=.03,ha='left',fontsize=13,fontweight='bold'); fig.subplots_adjust(left=.20,right=.95,top=.92,bottom=.08)
for ext in ['svg','pdf']: fig.savefig(F/f'Figure2_FAP_DPP4_replicated_contexts_FINAL249.{ext}',bbox_inches='tight')
fig.savefig(F/'Figure2_FAP_DPP4_replicated_contexts_FINAL249.png',dpi=600,bbox_inches='tight'); plt.close(fig)

# ---------- Figure 3 final249: breast health/cancer ----------
b=pd.read_csv(R/'31d_FAP_DPP4_axis_by_donor_collapsed_primary.csv.gz'); b=b[(b.tissue=='breast')&(b.cell_type=='fibroblast')]
cs=pd.read_csv(M/'TABLE_BREAST_HEALTH_CANCER_COLLECTION_MEDIANS_FINAL249.csv')
bp=pd.read_csv(M/'TABLE_BREAST_FAP_ONLY_vs_DPP4_ONLY_FINAL249.csv')
fig=plt.figure(figsize=(9.6,7.5)); gs=fig.add_gridspec(2,2,hspace=.42,wspace=.36)
# A disease-label medians (n>=4)
ax=fig.add_subplot(gs[0,0]); ds=(b.groupby('disease').agg(n=('donor_key','nunique'),collections=('collection_id','nunique'),FAP=('FAP_pos_fraction','median'),DPP4=('DPP4_pos_fraction','median'),DP=('double_positive_fraction','median')).reset_index()); ds=ds[ds.n>=4].copy(); order=['normal','breast cancer','estrogen-receptor positive breast cancer','triple-negative breast carcinoma','invasive ductal breast carcinoma','HER2 positive breast carcinoma']; ds=ds.set_index('disease').reindex([d for d in order if d in set(ds.disease)]).reset_index(); yy=np.arange(len(ds))
ax.scatter(ds.DPP4,yy,s=42,label='DPP4+',color='#0072B2'); ax.scatter(ds.FAP,yy,s=42,label='FAP+',color='#D55E00')
for i,r in ds.iterrows(): ax.plot([r.DPP4,r.FAP],[i,i],color='#bbbbbb',lw=.8,zorder=0); ax.text(max(r.FAP,r.DPP4)+.012,i,f'n={int(r.n)}',va='center',fontsize=7)
lab=[d.replace('estrogen-receptor positive breast cancer','ER+ breast cancer').replace('triple-negative breast carcinoma','TNBC').replace('invasive ductal breast carcinoma','IDC').replace('HER2 positive breast carcinoma','HER2+ breast cancer') for d in ds.disease]
ax.set_yticks(yy,lab,fontsize=7.5); ax.set_xlabel('Median donor-positive fraction'); ax.set_title('A  Breast fibroblast FAP/DPP4 detection',loc='left',fontweight='bold',fontsize=10.5); ax.legend(frameon=False,ncol=2,fontsize=7.5); ax.spines[['top','right']].set_visible(False)
# B collection medians normal vs generic cancer
ax=fig.add_subplot(gs[0,1]); mapx={'normal':0,'breast cancer':1}; jitter={'normal':-.04,'breast cancer':.04}
for dis,c in [('normal','#4c78a8'),('breast cancer','#e45756')]:
 q=cs[cs.disease==dis]; ax.scatter(np.repeat(mapx[dis],len(q)),q.median_axis,s=35+2*q.n_donors,c=c,alpha=.85,edgecolor='white',linewidth=.5)
 ax.hlines(q.median_axis.median(),mapx[dis]-.18,mapx[dis]+.18,color='#222222',lw=2)
ax.set_xticks([0,1],['Healthy\n(4 collections)','Breast cancer\n(3 collections)']); ax.set_ylabel('Collection median FAP − DPP4 fraction'); ax.set_title('B  Collection-level context sensitivity',loc='left',fontweight='bold',fontsize=10.5); ax.spines[['top','right']].set_visible(False); ax.text(.02,.02,'Descriptive only; no pooled disease-effect claim',transform=ax.transAxes,fontsize=7,color='#555555')
# C paired module heatmap normal vs cancer
ax=fig.add_subplot(gs[1,0]); mods=['PI16_homeostatic','COL15A1_like','ECM_matrix','myofibroblast','TGFb_LRRC15','inflammatory','antigen_presentation','detoxification']; diseases=['normal','breast cancer']; mat=pd.DataFrame(np.nan,index=diseases,columns=mods); qm=mat.copy(); nm=mat.copy()
for r in bp.itertuples(): mat.loc[r.disease,r.module]=r.median_paired_difference_gene_z; qm.loc[r.disease,r.module]=r.wilcoxon_q; nm.loc[r.disease,r.module]=r.n_paired_donors
im=ax.imshow(mat.to_numpy(float),aspect='auto',vmin=-.9,vmax=.9,cmap='coolwarm'); ax.set_xticks(np.arange(len(mods)),[m.replace('_','\n') for m in mods],rotation=35,ha='right',fontsize=6.8); ax.set_yticks([0,1],['Healthy breast\n(n=142 paired)','Breast cancer\n(n=44 paired)'],fontsize=8)
for i in range(2):
 for j in range(len(mods)):
  v=mat.iloc[i,j]
  if np.isfinite(v): ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=7.0,color='white' if abs(v)>.42 else 'black',fontweight='bold')
ax.set_title('C  FAP-only minus DPP4-only module scores (all q<0.05)',loc='left',fontweight='bold',fontsize=10.5); cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.02); cb.set_label('Median paired gene-z difference',fontsize=7.5)
# D collection-specific cancer effects
ax=fig.add_subplot(gs[1,1]); d=pd.read_parquet(R/'33c_state_module_profiles_by_donor_collapsed_exacttissue.parquet'); d=d[(d.tissue=='breast')&(d.cell_type=='fibroblast')&(d.disease=='breast cancer')]
rows=[]
for (col,mod),g in d.groupby(['collection_id','module']):
 p=g.pivot(index='donor_key',columns='state',values='module_mean_gene_z')
 if 'FAP+DPP4-' not in p or 'FAP-DPP4+' not in p: continue
 q=p[['FAP+DPP4-','FAP-DPP4+']].dropna()
 if len(q)>=3: rows.append((col,mod,len(q),(q['FAP+DPP4-']-q['FAP-DPP4+']).median()))
co=pd.DataFrame(rows,columns=['collection_id','module','n','effect']); cmap={c:plt.cm.tab10(i) for i,c in enumerate(sorted(co.collection_id.unique()))}
xx=np.arange(len(mods))
for col,g in co.groupby('collection_id'):
 vals=[g.loc[g.module==m,'effect'].iloc[0] if (g.module==m).any() else np.nan for m in mods]; ax.plot(xx,vals,marker='o',lw=1.0,ms=4,label=f'collection {list(sorted(co.collection_id.unique())).index(col)+1}',color=cmap[col])
ax.axhline(0,color='#666666',lw=.8); ax.set_xticks(xx,[m.replace('_','\n') for m in mods],rotation=35,ha='right',fontsize=6.8); ax.set_ylabel('Collection median paired difference'); ax.set_title('D  Cancer-state contrast across 3 collections',loc='left',fontweight='bold',fontsize=10.5); ax.legend(frameon=False,fontsize=7,ncol=1); ax.spines[['top','right']].set_visible(False)
fig.suptitle('FAP-only fibroblasts show a broad remodeling-associated program in breast cancer',x=.03,ha='left',fontsize=13,fontweight='bold'); fig.subplots_adjust(left=.20,right=.95,top=.92,bottom=.14)
for ext in ['svg','pdf']: fig.savefig(F/f'Figure3_Breast_FAP_DPP4_health_cancer_FINAL249.{ext}',bbox_inches='tight')
fig.savefig(F/'Figure3_Breast_FAP_DPP4_health_cancer_FINAL249.png',dpi=600,bbox_inches='tight'); plt.close(fig)
print('generated Figure2/3 final249')