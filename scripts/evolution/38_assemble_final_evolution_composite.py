from pathlib import Path
from PIL import Image
import io,base64,html
import cairosvg

EV=Path('/mnt/stk35data/BreastCancer_Census_v1/FAP_Prolyl_Peptidase_Atlas/evolution')
F=EV/'figures'
OUT=EV/'figure_package_FINAL_EVOLUTION_1086'/'composite'
OUT.mkdir(parents=True,exist_ok=True)

panels={
 'A':F/'Panel_A_Deep_S9_ML_FINAL.png',
 'B':F/'Panel_B_S9_duplication_history_FINAL.png',
 'C':F/'Panel_C_FAP_DPP4_synteny_FINAL.png',
 'D':F/'Panel_D_FAP_A657_ASR_FINAL.png',
 'E':F/'Panel_E_FAP_branchsite_structure_FINAL_compact.png',
 'F':F/'Panel_F_DPP6_DPP10_ASR_FINAL.png',
 'G':F/'S9_sister_pair_expression_divergence_final1086.png',
}
# Canvas 180 x 240 mm with 20 SVG units/mm.
W,H=3600,4800
cells={
 'A':(120,80,3360,1770),
 'B':(120,1910,1600,820),
 'C':(1880,1910,1600,820),
 'D':(120,2810,1600,780),
 'F':(1880,2810,1600,780),
 'E':(120,3675,1600,900),
 'G':(1880,3675,1600,900),
}

def crop_and_b64(path):
    im=Image.open(path).convert('RGBA')
    # Crop only nearly-white outer margins; retain all plotted gray/text elements.
    rgb=im.convert('RGB')
    import numpy as np
    a=np.asarray(rgb)
    mask=np.any(a<247,axis=2)
    ys,xs=np.where(mask)
    if len(xs):
        pad=18
        x0=max(0,int(xs.min())-pad); x1=min(im.width,int(xs.max())+pad+1)
        y0=max(0,int(ys.min())-pad); y1=min(im.height,int(ys.max())+pad+1)
        im=im.crop((x0,y0,x1,y1))
    bio=io.BytesIO(); im.save(bio,format='PNG',optimize=True)
    return base64.b64encode(bio.getvalue()).decode(), im.size

items=[]
for letter,path in panels.items():
    b64,(iw,ih)=crop_and_b64(path)
    x,y,cw,ch=cells[letter]
    scale=min(cw/iw,ch/ih)
    rw,rh=iw*scale,ih*scale
    ix=x+(cw-rw)/2; iy=y+(ch-rh)/2
    items.append(f'<image x="{ix:.2f}" y="{iy:.2f}" width="{rw:.2f}" height="{rh:.2f}" href="data:image/png;base64,{b64}"/>')
    # Large panel letter outside the image content but inside cell.
    items.append(f'<text x="{x-88}" y="{y+82}" font-family="Arial, Helvetica, sans-serif" font-size="82" font-weight="700" fill="#111111">{letter}</text>')

svg=f'''<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="240mm" viewBox="0 0 {W} {H}">
<rect x="0" y="0" width="{W}" height="{H}" fill="white"/>
{''.join(items)}
</svg>'''
svg_path=OUT/'Evolution_main_figure_FINAL_1086.svg'
svg_path.write_text(svg)
# Vector container with embedded high-resolution panel rasters.
cairosvg.svg2pdf(bytestring=svg.encode(),write_to=str(OUT/'Evolution_main_figure_FINAL_1086.pdf'),output_width=180/25.4*72,output_height=240/25.4*72)
# 600 dpi raster submission image.
png_w=round(180/25.4*600); png_h=round(240/25.4*600)
cairosvg.svg2png(bytestring=svg.encode(),write_to=str(OUT/'Evolution_main_figure_FINAL_1086_600dpi.png'),output_width=png_w,output_height=png_h)
# TIFF LZW with explicit 600 dpi metadata.
im=Image.open(OUT/'Evolution_main_figure_FINAL_1086_600dpi.png').convert('RGB')
im.save(OUT/'Evolution_main_figure_FINAL_1086_600dpi.tif',compression='tiff_lzw',dpi=(600,600))
# Lightweight preview.
im.thumbnail((1800,2400),Image.Resampling.LANCZOS)
im.save(OUT/'Evolution_main_figure_FINAL_1086_preview.jpg',quality=90,dpi=(150,150))
print(svg_path)
print('png',png_w,png_h)