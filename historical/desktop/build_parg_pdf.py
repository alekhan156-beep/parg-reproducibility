from pathlib import Path
import re,html,shutil
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.enums import TA_LEFT,TA_CENTER
root=Path(__file__).parent; out=root/'outputs'; text=(out/'PARG_Methodology.md').read_text(encoding='utf-8')
for name,file in [('Body','times.ttf'),('Bold','timesbd.ttf')]: pdfmetrics.registerFont(TTFont(name,str(Path('C:/Windows/Fonts')/file)))
body=ParagraphStyle('Body',fontName='Body',fontSize=10.5,leading=13,spaceAfter=6,allowWidows=0,allowOrphans=0)
title=ParagraphStyle('Title',parent=body,fontName='Bold',fontSize=21,leading=24,spaceAfter=14,keepWithNext=True)
h1=ParagraphStyle('H1',parent=body,fontName='Bold',fontSize=13,leading=16,spaceBefore=10,spaceAfter=6,keepWithNext=True)
h2=ParagraphStyle('H2',parent=h1,fontSize=11.5,leading=14,spaceBefore=8)
cell=ParagraphStyle('Cell',parent=body,fontSize=8.5,leading=10,spaceAfter=0)
head=ParagraphStyle('CellHead',parent=cell,fontName='Bold')
def esc(s):
    return html.escape(s).replace('−','-').replace('–','-')
story=[]; lines=text.splitlines(); i=0; width=A4[0]-108
while i<len(lines):
    s=lines[i].strip()
    if not s: i+=1; continue
    if s.startswith('|'):
        rows=[]
        while i<len(lines) and lines[i].strip().startswith('|'):
            vals=[x.strip() for x in lines[i].strip().strip('|').split('|')]
            if not all(re.match(r'^[-: ]+$',v) for v in vals): rows.append(vals)
            i+=1
        n=len(rows[0]); weights=[1]*n
        if n==6: weights=[.5,1.1,1.1,1.1,1.1,1.35]
        elif rows[0][0]=='Software': weights=[1.4,2.1,2.8]
        elif rows[0][-1]=='Cases': weights=[1.4,4.4,.5]
        elif n==3: weights=[1.3,3.4,1.6]
        data=[[Paragraph(esc(v),head if j==0 else cell) for v in row] for j,row in enumerate(rows)]
        t=Table(data,colWidths=[width*w/sum(weights) for w in weights],repeatRows=1,hAlign='LEFT')
        t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#E7E6E6')),('GRID',(0,0),(-1,-1),.4,colors.HexColor('#D9D9D9')),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6)]))
        story.extend([t,Spacer(1,8)]); continue
    if s.startswith('# '): story.append(Paragraph(esc(s[2:]),title))
    elif s.startswith('## '): story.append(Paragraph(esc(s[3:]),h1))
    elif s.startswith('### '): story.append(Paragraph(esc(s[4:]),h2))
    else:
        v=esc(s)
        v=re.sub(r'https://[^\s]+',lambda m:'<link href="'+m.group(0)+'" color="#000000">'+m.group(0)+'</link>',v)
        story.append(Paragraph(v,body))
    i+=1
def footer(c,d):
    c.setFont('Body',9); c.drawCentredString(A4[0]/2,29,f'Page {d.page}')
path=out/'PARG_Methodology.pdf'
doc=SimpleDocTemplate(str(path),pagesize=A4,leftMargin=54,rightMargin=54,topMargin=48,bottomMargin=48,title='Computational methodology for PARG ligand modelling and docking',author='Muhammad Ali Khan')
doc.build(story,onFirstPage=footer,onLaterPages=footer)
from pypdf import PdfReader
reader=PdfReader(path)
assert len(reader.pages)>5
alltext='\n'.join(p.extract_text() for p in reader.pages)
assert '{{STATUS}}' not in alltext
assert all(x in alltext for x in ['0.445215','518','H121','650 IC50','Appendix C','References'])
shutil.copy2(path,Path(r'\\wsl.localhost\Ubuntu-24.04\home\aliwhoami\parg\runs\20261016_132012')/path.name)
print('PAGES',len(reader.pages)); print('CHARACTERS',len(alltext))
