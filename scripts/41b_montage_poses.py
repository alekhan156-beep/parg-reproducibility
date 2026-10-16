import sys, os
from PIL import Image, ImageDraw, ImageFont
SRC, OUT = sys.argv[1], sys.argv[2]
FD = '/usr/share/fonts/truetype/dejavu'
f_t  = ImageFont.truetype(f'{FD}/DejaVuSans-Bold.ttf', 40)
f_l  = ImageFont.truetype(f'{FD}/DejaVuSans-Bold.ttf', 30)
f_s  = ImageFont.truetype(f'{FD}/DejaVuSans.ttf', 26)
f_sb = ImageFont.truetype(f'{FD}/DejaVuSans-Bold.ttf', 26)
f_c  = ImageFont.truetype(f'{FD}/DejaVuSans.ttf', 24)

GREEN, ORANGE, BLUE, INK, GREY = (34,110,34), (230,126,0), (0,114,178), (26,26,26), (110,110,110)
# tag, dry, wet, note
P = [('6OAL_M0P',3.66,0.57,''), ('7KG7_WDJ',4.59,0.74,''),
     ('6O9X_M0S',3.74,0.62,''), ('6OA1_M0P',3.68,0.69,''),
     ('7KG8_WDM',3.36,0.62,''), ('7KFP_WD7',2.73,1.94,''),
     ('7KG6_WDG',3.63,6.42,'water makes it worse'),
     ('6HMK_7JC',4.38,5.32,'water-displacing — fails either way')]

PW, GAP, HDR, LAB, PAD = 700, 18, 150, 112, 26
src0 = Image.open(f'{SRC}/{P[0][0]}.png'); PH = int(PW * src0.size[1] / src0.size[0])
W = PAD*2 + PW*4 + GAP*3
H = HDR + (PH + LAB)*2 + GAP + PAD
im = Image.new('RGB', (W, H), 'white'); dr = ImageDraw.Draw(im)

dr.text((PAD, 24), 'Bridging waters change which pose is selected, not which are available',
        font=f_t, fill=INK)
x = PAD; y = 86
for col, txt in ((GREEN,'crystal ligand'), (ORANGE,'top pose, dry receptor'),
                 (BLUE,'top pose, wet receptor')):
    dr.rectangle([x, y+6, x+30, y+26], fill=col); x += 40
    dr.text((x, y), txt, font=f_s, fill=INK); x += int(dr.textlength(txt, font=f_s)) + 46
dr.ellipse([x, y+6, x+20, y+26], fill=(139,26,26)); x += 30
dr.text((x, y), 'retained bridging waters', font=f_s, fill=INK)

for i, (tag, d, w, note) in enumerate(P):
    r, c = divmod(i, 4)
    px = PAD + c*(PW+GAP); py = HDR + r*(PH+LAB+GAP)
    im.paste(Image.open(f'{SRC}/{tag}.png').resize((PW, PH), Image.LANCZOS), (px, py))
    dr.rectangle([px, py, px+PW, py+PH], outline=(224,224,224), width=1)
    ly = py + PH + 8
    dr.text((px+2, ly), tag, font=f_l, fill=INK)
    tx = px + 2; ty = ly + 38
    dr.text((tx, ty), 'dry ', font=f_s, fill=GREY)
    tx += dr.textlength('dry ', font=f_s)
    s = f'{d:.2f}'; dr.text((tx, ty), s, font=f_sb, fill=ORANGE)
    tx += dr.textlength(s, font=f_sb) + 16
    dr.text((tx, ty), 'wet ', font=f_s, fill=GREY)
    tx += dr.textlength('wet ', font=f_s)
    s = f'{w:.2f}'; dr.text((tx, ty), s, font=f_sb, fill=BLUE if w < 2 else ORANGE)
    tx += dr.textlength(s, font=f_sb) + 10
    dr.text((tx, ty), 'Å', font=f_s, fill=GREY)
    if note:
        dr.text((px+2, ty+36), note, font=f_c, fill=GREY)

im.save(f'{OUT}/fig5_poses.png', dpi=(300,300))
im.convert('RGB').save(f'{OUT}/fig5_poses.pdf', 'PDF', resolution=300)
print('wrote', f'{OUT}/fig5_poses.png', im.size)
