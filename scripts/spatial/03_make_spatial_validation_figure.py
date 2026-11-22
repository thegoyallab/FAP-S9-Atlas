from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
ROOT=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas'); SP=ROOT/'spatial_validation'; R=SP/'results'; F=SP/'figures'
ST=Path('/mnt/stk35data/BreastCancer_Census_v1/Stage11E_SpatialTranscriptomics_v0_1/results')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.linewidth':.8,'svg.fonttype':'none','pdf.fonttype':42})

Pdet=pd.read_csv(R/'GSE210616_FAP_DPP4_patient_detection.tsv',sep='\t')
G=pd.read_csv(R/'GSE210616_FAP_DPP4_lineage_patient_results.tsv',sep='\t')
T=pd.read_csv(R/'GSE210616_FAP_DPP4_module_patient_results.tsv',sep='\t')
Ts=pd.read_csv(R/'GSE210616_FAP_DPP4_module_patient_summary.tsv',sep='\t')
Dist=pd.read_csv(R/'GSE210616_axis_epithelial_distance_patient.tsv',sep='\t')
Ds=pd.read_csv(R/'GSE210616_axis_epithelial_distance_summary.tsv',sep='\t')

fig=plt.figure(figsize=(10.2,8.8)); gs=fig.add_gridspec(3,2,height_ratios=[1,1.15,.95],hspace=.48,wspace=.34)
# A detection fractions
ax=fig.add_subplot(gs[0,0]); y1=Pdet.FAP_detection_fraction.to_numpy(); y2=Pdet.DPP4_detection_fraction.to_numpy(); x=np.array([0,1])
for a,b in zip(y1,y2): ax.plot(x,[a,b],color='#cfcfcf',lw=.7,zorder=1)
ax.scatter(np.zeros(len(y1)),y1,s=24,color='#D55E00',zorder=2,label='FAP'); ax.scatter(np.ones(len(y2)),y2,s=24,color='#0072B2',zorder=2,label='DPP4')
ax.set_xticks([0,1],['FAP','DPP4']); ax.set_ylabel('Patient median detection fraction'); ax.set_title('A  FAP is detected more broadly than DPP4 in TNBC',loc='left',fontweight='bold',fontsize=10.5); ax.spines[['top','right']].set_visible(False)
ax.text(.03,.96,'22 patients / 43 sections\npaired Wilcoxon P=0.00094',transform=ax.transAxes,va='top',fontsize=7.4)
# B fibroblast affinity
ax=fig.add_subplot(gs[0,1]); p=G.pivot(index='patient_id',columns='gene',values='partial_rho_fibroblast_affinity').dropna(); y1=p.FAP.to_numpy(); y2=p.DPP4.to_numpy()
for a,b in zip(y1,y2): ax.plot(x,[a,b],color='#cfcfcf',lw=.7,zorder=1)
ax.scatter(np.zeros(len(y1)),y1,s=24,color='#D55E00',zorder=2); ax.scatter(np.ones(len(y2)),y2,s=24,color='#0072B2',zorder=2)
ax.axhline(0,color='#999999',lw=.7); ax.set_xticks([0,1],['FAP','DPP4']); ax.set_ylabel('Partial ρ with FAP-free fibroblast affinity'); ax.set_title('B  FAP is more strongly localized to fibroblast-rich spots',loc='left',fontweight='bold',fontsize=10.5); ax.spines[['top','right']].set_visible(False)
ax.text(.98,.96,'median paired Δρ=+0.151\n22/22 patients FAP>DPP4\nP=4.77×10⁻⁷',transform=ax.transAxes,ha='right',va='top',fontsize=7.4)
# C module correlations
ax=fig.add_subplot(gs[1,:]); modules=['ECM_matrix','myofibroblast','TGFb_LRRC15','COL15A1_like','detoxification','antigen_presentation','PI16_homeostatic','inflammatory']; labels=['ECM/matrix','Myofibroblast','TGF-β/LRRC15','COL15A1-like','Detoxification','Antigen presentation','PI16/homeostatic','Inflammatory']
vals=[]
for m in modules: vals.append(T[(T.module==m)].partial_spearman_axis_module.dropna().to_numpy())
parts=ax.boxplot(vals,positions=np.arange(len(modules)),widths=.55,patch_artist=True,showfliers=False,medianprops={'color':'black','linewidth':1})
for b in parts['boxes']: b.set_facecolor('#e6e6e6'); b.set_edgecolor('#777777')
for i,v in enumerate(vals):
    # deterministic horizontal offsets
    off=np.linspace(-.14,.14,len(v)); ax.scatter(np.full(len(v),i)+off,v,s=13,alpha=.75,color='#555555',zorder=3)
ax.axhline(0,color='#999999',lw=.7); ax.set_xticks(np.arange(len(modules)),labels,rotation=28,ha='right',fontsize=7.5); ax.set_ylabel('Patient partial Spearman ρ\nFAP–DPP4 axis vs module'); ax.set_title('C  FAP-dominant fibroblast regions carry remodeling programs',loc='left',fontweight='bold',fontsize=10.5); ax.spines[['top','right']].set_visible(False)
for i,m in enumerate(modules):
    r=Ts[(Ts.module==m)&(Ts.metric=='partial_spearman_axis_module')].iloc[0]; txt=f"q={r.wilcoxon_q:.1g}" if r.wilcoxon_q<.05 else 'ns'; ax.text(i,ax.get_ylim()[1]*.93,txt,ha='center',va='top',fontsize=6.8,fontweight='bold' if r.wilcoxon_q<.05 else 'normal')
# D independent histology annotations
ax=fig.add_subplot(gs[2,0]); cats=['stromal','mixed/invasive','invasive','DCIS #1','DCIS #2']; fap=json.load(open(ST/'FAP_SpatialGene_v0_4.json'))['validation']['GSE243275_published_annotation']['published_annotation_expression']; dpp=json.load(open(ST/'DPP4_SpatialGene_v0_4.json'))['validation']['GSE243275_published_annotation']['published_annotation_expression']; fd={r['published_annotation']:r for r in fap}; dd={r['published_annotation']:r for r in dpp}; xx=np.arange(len(cats)); w=.36
ax.bar(xx-w/2,[fd[c]['detection_fraction'] for c in cats],w,label='FAP',color='#D55E00'); ax.bar(xx+w/2,[dd[c]['detection_fraction'] for c in cats],w,label='DPP4',color='#0072B2')
ax.set_xticks(xx,['Stromal','Mixed/\ninvasive','Invasive','DCIS 1','DCIS 2'],fontsize=7.3); ax.set_ylabel('Detection fraction'); ax.set_title('D  Published histology: stromal FAP enrichment',loc='left',fontweight='bold',fontsize=10.0); ax.legend(frameon=False,fontsize=7.5); ax.spines[['top','right']].set_visible(False); ax.text(.03,.95,'GSE243275; 1 annotated section\nFAP stromal median CPM=216',transform=ax.transAxes,va='top',fontsize=7.2)
# E interface/distance result
ax=fig.add_subplot(gs[2,1]); v=Dist.partial_spearman_axis_vs_distance.to_numpy(); off=np.linspace(-.12,.12,len(v)); ax.scatter(off,v,s=24,color='#666666'); ax.boxplot([v],positions=[0],widths=.32,showfliers=False,medianprops={'color':'black','linewidth':1.2},boxprops={'color':'#555555'},whiskerprops={'color':'#777777'},capprops={'color':'#777777'}); ax.axhline(0,color='#999999',lw=.7); ax.set_xlim(-.45,.45); ax.set_xticks([0],['22 TNBC patients']); ax.set_ylabel('Partial ρ: FAP–DPP4 axis vs\nepithelial-domain distance'); ax.set_title('E  FAP-dominant axis is not interface-enriched',loc='left',fontweight='bold',fontsize=10.0); ax.spines[['top','right']].set_visible(False); r=Ds[Ds.metric=='partial_spearman_axis_vs_distance'].iloc[0]; ax.text(.03,.95,f"median ρ={r['median']:.3f}\nP={r.wilcoxon_p:.4g}\npositive = farther from epithelial-rich domain",transform=ax.transAxes,va='top',fontsize=7.2)
fig.suptitle('Spatial transcriptomics independently validates FAP-dominant remodeling in breast cancer',x=.03,ha='left',fontsize=13,fontweight='bold')
fig.subplots_adjust(left=.11,right=.96,top=.93,bottom=.09)
for ext in ['svg','pdf']: fig.savefig(F/f'Figure4_Spatial_FAP_DPP4_TNBC_FINAL.{ext}',bbox_inches='tight')
fig.savefig(F/'Figure4_Spatial_FAP_DPP4_TNBC_FINAL.png',dpi=600,bbox_inches='tight'); plt.close(fig)
print('generated')