from pathlib import Path
import re,csv,toytree
import toyplot.svg,toyplot.pdf,toyplot.png
EV=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas/evolution'); OUT=EV/'figures'
COL={'S9A':'#0072B2','S9B':'#D55E00','S9C':'#009E73','S9D':'#CC79A7','other':'#666666'}
meta={}
with (EV/'results/51_deep_S9_curated_registry.tsv').open(newline='') as fh:
 for r in csv.DictReader(fh,delimiter='\t'):
  fam=(r.get('Protein families') or ''); name=r['Entry Name'].strip(); entry=r['Entry'].strip(); group=r['query_group'].strip()
  if 'S9A' in fam: cl='S9A'
  elif 'S9B' in fam: cl='S9B'
  elif 'S9C' in fam: cl='S9C'
  elif 'S9D' in fam: cl='S9D'
  elif name.startswith(('DPP4_','SEPR_','DPP6_','DPP10_','DPP8_','DPP9_')): cl='S9B'
  elif name.startswith(('PPCE_','PPCEL_')): cl='S9A'
  elif name.startswith(('ACPH_','APEH_')): cl='S9C'
  else: cl='other'
  shortgrp={'Bacteria':'Bac','Archaea':'Arc','Arabidopsis':'Plant','Yeast':'Yeast','Human':'Human'}.get(group,group[:6])
  meta[entry]=(name,shortgrp,cl)
tree=toytree.tree(EV/'results/51_deep_S9_curated_ML.treefile').unroot(); tips=[]; colors=[]
for n in tree[:tree.ntips]:
 m=re.search(r'\|([^|]+)\|([^:]+)$',n.name); entry=m.group(1) if m else n.name; fallback=m.group(2) if m else n.name
 name,grp,cl=meta.get(entry,(fallback,'','other')); tips.append(f'{name} [{grp}; {cl}]'); colors.append(COL[cl])
labs=[]
for n in tree:
 if n.is_leaf() or not n.name: labs.append(''); continue
 m=re.fullmatch(r'([0-9.]+)/([0-9.]+)',str(n.name)); labs.append(f'{float(m.group(1)):g}/{float(m.group(2)):g}' if m and float(m.group(1))>=90 and float(m.group(2))>=95 else '')
c,a,m=tree.draw(width=1850,height=1050,layout='r',use_edge_lengths=True,tip_labels=tips,tip_labels_colors=colors,tip_labels_style={'font-size':'17px','font-weight':'500'},node_labels=labs,node_labels_style={'font-size':'11px','fill':'#444444'},node_sizes=0,edge_widths=2.0,edge_style={'stroke':'#555555','stroke-linecap':'round'},padding=90,scale_bar=True,edge_type='p'); c.style={'background-color':'white'}
stem='Panel_A_Deep_S9_ML_FINAL'
toyplot.svg.render(c,str(OUT/f'{stem}.svg')); toyplot.pdf.render(c,str(OUT/f'{stem}.pdf')); toyplot.png.render(c,str(OUT/f'{stem}.png'),scale=3)