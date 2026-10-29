# compose_fig7_hpa.py — assemble HPA IHC figure (landscape 2 rows x 5 gene columns)
from PIL import Image, ImageDraw, ImageFont
import pathlib

D = pathlib.Path(r'D:/yuproject/前列腺癌/08_QC_Revision/data/external/hpa')
FIG = pathlib.Path(r'D:/yuproject/前列腺癌/08_QC_Revision/figures')

PICK = {  # gene -> (normal_file, cancer_file, ab, normal_id, tumour_id)
    'AMACR': ('AMACR_normal_0.jpg', 'AMACR_cancer_2.jpg', 'ab1809',  '4880_A_1_5',  '4790_A_7_3'),
    'HPN':   ('HPN_normal_0.jpg',   'HPN_cancer_1.jpg',   'ab6804',  '19524_A_1_5', '19521_A_7_3'),
    'MYC':   ('MYC_normal_0.jpg',   'MYC_cancer_2.jpg',   'ab10307', '30810_A_1_5', '30656_A_7_5'),
    'EZH2':  ('EZH2_normal_0.jpg',  'EZH2_cancer_1.jpg',  'ab9589',  '23664_A_1_5', '23661_A_7_2'),
    'FASN':  ('FASN_normal_0.jpg',  'FASN_cancer_1.jpg',  'ab15417', '155365_A_1_5','155362_A_7_3'),
}
GENES = list(PICK)
CELL, GAP, LBL, HEAD, PROV = 760, 26, 800, 150, 48
W = LBL + CELL * len(GENES) + GAP * (len(GENES) + 1)
H = HEAD + (CELL + PROV + GAP) * 2 + 100

img = Image.new('RGB', (W, H), 'white')
d = ImageDraw.Draw(img)
try:
    fb = ImageFont.truetype(r'C:/Windows/Fonts/arialbd.ttf', 60)
    fr = ImageFont.truetype(r'C:/Windows/Fonts/arialbd.ttf', 50)
    fs = ImageFont.truetype(r'C:/Windows/Fonts/arial.ttf', 34)
    ft = ImageFont.truetype(r'C:/Windows/Fonts/arial.ttf', 32)
except Exception:
    fb = fr = fs = ft = ImageFont.load_default()

# column headers: gene names
for j, g in enumerate(GENES):
    x = LBL + GAP + j * (CELL + GAP)
    d.text((x + CELL // 2, HEAD // 2), g, font=fb, fill='black', anchor='mm')

rows = [('Benign prostate', 'normal', 'A'), ('Prostate adenocarcinoma', 'cancer', 'B')]
for r, (title, key, tag) in enumerate(rows):
    y = HEAD + r * (CELL + PROV + GAP)
    ymid = y + CELL // 2
    d.text((44, ymid), tag, font=fb, fill='black', anchor='lm')
    d.text((LBL - 30, ymid), title, font=fr, fill='black', anchor='rm')
    for j, g in enumerate(GENES):
        x = LBL + GAP + j * (CELL + GAP)
        im = Image.open(D / PICK[g][0 if key == 'normal' else 1]).convert('RGB')
        im = im.resize((CELL, CELL), Image.LANCZOS)
        img.paste(im, (x, y))
        prov = f"{PICK[g][2]} \u00b7 {PICK[g][3 if key == 'normal' else 4]}"
        d.text((x + CELL // 2, y + CELL + 8), prov, font=fs, fill='#595959', anchor='ma')

foot = ('Source: Human Protein Atlas (www.proteinatlas.org), v24, CC BY-SA 3.0. '
        'PCA3 (lncRNA) and HOXC6 are not profiled by IHC in prostate tissue on the atlas; '
        'full image identifiers are listed in Table S17.')
d.text((LBL + GAP, H - 60), foot, font=ft, fill='#404040', anchor='la')

out = FIG / 'Fig7_HPA_protein.png'
img.save(out, dpi=(600, 600))
print('saved', out, img.size, '->', round(W / 600, 2), 'x', round(H / 600, 2), 'in @600dpi')
