import sys, os
from pymol import cmd
ROOT, REC, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)
WAT = [1211, 1232, 1307, 1317, 1402]
PANELS = ['6OAL_M0P','7KG7_WDJ','6O9X_M0S','6OA1_M0P',
          '7KG8_WDM','7KFP_WD7','7KG6_WDG','6HMK_7JC']
cmd.feedback('disable','all','everything')
cmd.set('ray_opaque_background',1); cmd.set('ambient',0.35)
cmd.set('stick_radius',0.13); cmd.set('sphere_scale',0.22)
cmd.set('ray_shadows',0); cmd.set('antialias',2); cmd.set('depth_cue',0)
cmd.set('orthoscopic',1)

# establish one view that contains every ligand in every panel, so the panels
# are directly comparable and nothing clips at the edges
cmd.reinitialize()
for tag in PANELS:
    cmd.load(f'{ROOT}/structures/xdock/{tag}.pdb', f'v_{tag}_x')
    for rec in ('dry','wet'):
        cmd.load(f'{ROOT}/results/tmp_poses/{tag}_{rec}1.pdb', f'v_{tag}_{rec}')
cmd.orient('v_*'); cmd.zoom('v_*', 2.0)
view = cmd.get_view()
cmd.reinitialize()

for tag in PANELS:
    cmd.reinitialize(); cmd.bg_color('white')
    cmd.set('ray_opaque_background',1); cmd.set('ambient',0.35)
    cmd.set('stick_radius',0.13); cmd.set('sphere_scale',0.22)
    cmd.set('ray_shadows',0); cmd.set('antialias',2); cmd.set('depth_cue',0)
    cmd.set('orthoscopic',1)
    cmd.load(REC,'rec'); cmd.remove('rec and not chain A')
    cmd.load(f'{ROOT}/structures/xdock/{tag}.pdb','xtal')
    cmd.load(f'{ROOT}/results/tmp_poses/{tag}_dry1.pdb','dryp')
    cmd.load(f'{ROOT}/results/tmp_poses/{tag}_wet1.pdb','wetp')
    cmd.remove('hydro')
    cmd.hide('everything')
    # pocket
    cmd.show('sticks','rec and polymer and byres (polymer within 5.5 of xtal)')
    cmd.color('grey90','rec and polymer')
    cmd.set('stick_radius',0.055,'rec')
    cmd.set('stick_transparency',0.45,'rec')
    # retained waters
    wsel = 'rec and resn HOH and (' + ' or '.join(f'resi {r}' for r in WAT) + ')'
    cmd.show('spheres', wsel); cmd.color('firebrick', wsel)
    # ligands
    for obj,col in (('xtal','forest'),('dryp','orange'),('wetp','marine')):
        cmd.show('sticks',obj); cmd.color(col,obj)
        cmd.set('stick_radius',0.145,obj)
    # one flat colour per pose: the point of the panel is telling them apart
    for obj,col in (('xtal','forest'),('dryp','orange'),('wetp','marine')):
        cmd.color(col, obj)
    cmd.set_view(view)
    cmd.png(f'{OUT}/{tag}.png', width=1100, height=900, dpi=200, ray=1)
    print('rendered', tag)
