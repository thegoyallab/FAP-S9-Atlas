from pathlib import Path
import pandas as pd,numpy as np,matplotlib.pyplot as plt
EV=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas/evolution'); R=EV/'results'; F=EV/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':.8,'svg.fonttype':'none','pdf.fonttype':42})
def save(fig,stem):
 for ext in ['svg','pdf']: fig.savefig(F/f'{stem}.{ext}',bbox_inches='tight')
 fig.savefig(F/f'{stem}.png',dpi=600,bbox_inches='tight'); plt.close(fig)
# D
a=pd.read_csv(R/'85b_FAP_A657_ancestral_probabilities.tsv',sep='\t'); order=['all_vertebrate_FAPDPP4','gar_teleost','teleost_four','FAP_DPP4_predup','DPP4_lungfish_tetrapod','FAP_lungfish_tetrapod']; a=a.set_index('biological_node').loc[order].reset_index(); labels=['All vertebrate FAP/DPP4','Gar + teleost','Teleost duplicated ancestor','Pre-duplication FAP/DPP4','Lungfish/tetrapod DPP4','Lungfish/tetrapod FAP']
fig,ax=plt.subplots(figsize=(7.0,4.0)); y=np.arange(len(a)); left=np.zeros(len(a)); cols=[('p_D','D','#4c78a8'),('p_A','A','#e45756'),('p_S','S','#72b7b2'),('p_G','G','#f2cf5b')]
for col,lab,c in cols:
 vals=a[col].to_numpy(float); ax.barh(y,vals,left=left,height=.58,color=c,edgecolor='white',linewidth=.4)
 for i,v in enumerate(vals):
  if v>=.06: ax.text(left[i]+v/2,i,lab,ha='center',va='center',fontsize=8,color='white' if lab in ['D','A'] else '#222222',fontweight='bold')
 left+=vals
ax.set_yticks(y,labels); ax.set_xlim(0,1.08); ax.set_xlabel('Posterior probability'); ax.set_title('FAP A657-equivalent ancestral-state reconstruction',loc='left',fontweight='bold',fontsize=11,pad=7); ax.spines[['top','right']].set_visible(False)
for i,r in a.iterrows(): ax.text(1.01,i,f'{r.ML_state}  {r.p_max:.5f}',ha='left',va='center',fontsize=7.7)
fig.subplots_adjust(left=.29,right=.91,top=.88,bottom=.16); save(fig,'Panel_D_FAP_A657_ASR_FINAL')
# F
d=pd.read_csv(R/'87b_DPP6_DPP10_S630eq_ancestral_probabilities.tsv',sep='\t'); order=['full_DPP6_DPP10_family','vertebrate_DPP6_DPP10','DPP6_crown','DPP10_crown']; d=d.set_index('biological_node').loc[order].reset_index(); labels=['Deep DPP6/DPP10 precursor','Vertebrate DPP6/DPP10 ancestor','DPP6 crown','DPP10 crown']
fig,ax=plt.subplots(figsize=(7.0,3.8)); y=np.arange(len(d)); left=np.zeros(len(d)); cols=[('p_S','S','#72b7b2'),('p_A','A','#e45756'),('p_D','D','#4c78a8'),('p_G','G','#f2cf5b')]
for col,lab,c in cols:
 vals=d[col].to_numpy(float); ax.barh(y,vals,left=left,height=.58,color=c,edgecolor='white',linewidth=.4)
 for i,v in enumerate(vals):
  if v>=.06: ax.text(left[i]+v/2,i,lab,ha='center',va='center',fontsize=8,color='white' if lab in ['D','A'] else '#222222',fontweight='bold')
 left+=vals
ax.set_yticks(y,labels); ax.set_xlim(0,1.08); ax.set_xlabel('Posterior probability'); ax.set_title('DPP6/DPP10 catalytic-serine-equivalent ancestral states',loc='left',fontweight='bold',fontsize=11,pad=7); ax.spines[['top','right']].set_visible(False)
for i,r in d.iterrows(): ax.text(1.01,i,f'{r.ML_state}  {r.p_max:.5f}',ha='left',va='center',fontsize=7.7)
fig.subplots_adjust(left=.34,right=.91,top=.86,bottom=.17); save(fig,'Panel_F_DPP6_DPP10_ASR_FINAL')