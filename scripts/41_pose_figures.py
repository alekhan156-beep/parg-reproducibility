#!/usr/bin/env python3
"""
Pose figures: crystal ligand vs top-ranked dry pose vs top-ranked wet pose.

The pose coordinates only exist on your machine, so this has to run there.

    cd ~/parg
    pymol -cq scripts/41_pose_figures.py

Needs PyMOL (conda install -c conda-forge pymol-open-source).
Reads the per-mode PDBs that 40_allposes.py already wrote into
results/tmp_poses/ , so run that first (you have).

Writes figures/poses/<tag>.png  — one panel per ligand:
    green  = crystal ligand (truth)
    orange = top-ranked pose, dry receptor
    blue   = top-ranked pose, wet receptor
    red spheres = the five retained bridging waters
"""
import os, glob, sys

try:
    from pymol import cmd
except ImportError:
    sys.exit('PyMOL not importable. Run with:  pymol -cq scripts/41_pose_figures.py\n'
             'Install:  conda install -c conda-forge pymol-open-source')

OUT = 'figures/poses'
os.makedirs(OUT, exist_ok=True)

# the six where waters rescued the pose, plus the one where they hurt,
# plus a water-displacing control that fails in both
LIGANDS = ['6O9X_M0S', '6OA1_M0P', '6OAL_M0P', '7KFP_WD7', '7KG7_WDJ',
           '7KG8_WDM', '7KG6_WDG', '6HMK_7JC']

# find the receptor however it was saved
REC = None
for pat in ('structures/rec_wet.pdb', 'structures/rec_wet.pdbqt',
            'rec_wet.pdb', 'rec_wet.pdbqt', 'structures/6OAK*.pdb'):
    hit = glob.glob(pat)
    if hit:
        REC = hit[0]; break
print(f'receptor: {REC or "NOT FOUND - panels will show ligands only"}')

made, skipped = [], []
for tag in LIGANDS:
    ref  = f'structures/xdock/{tag}.pdb'
    dry  = f'results/tmp_poses/{tag}_dry1.pdb'
    wet  = f'results/tmp_poses/{tag}_wet1.pdb'
    have = [p for p in (ref, dry, wet) if os.path.exists(p)]
    if ref not in have:
        skipped.append(f'{tag}: no crystal reference at {ref}'); continue

    cmd.reinitialize()
    cmd.bg_color('white')
    cmd.load(ref, 'xtal')
    cmd.show_as('sticks', 'xtal'); cmd.color('forest', 'xtal')

    if os.path.exists(dry):
        cmd.load(dry, 'dry'); cmd.show_as('sticks', 'dry'); cmd.color('orange', 'dry')
    if os.path.exists(wet):
        cmd.load(wet, 'wet'); cmd.show_as('sticks', 'wet'); cmd.color('marine', 'wet')

    if REC:
        cmd.load(REC, 'rec')
        cmd.hide('everything', 'rec')
        cmd.show('sticks', 'rec and byres (polymer within 5 of xtal)')
        cmd.color('grey70', 'rec and polymer')
        cmd.set('stick_transparency', 0.55, 'rec')
        # retained waters, if the receptor file carries them
        cmd.show('spheres', 'rec and resn HOH')
        cmd.color('red', 'rec and resn HOH')
        cmd.set('sphere_scale', 0.25, 'rec and resn HOH')

    cmd.hide('everything', 'hydro')
    cmd.orient('xtal')
    cmd.zoom('xtal', 4)
    cmd.set('ray_opaque_background', 1)
    cmd.set('ambient', 0.3)
    png = f'{OUT}/{tag}.png'
    cmd.png(png, width=1400, height=1100, dpi=300, ray=1)
    made.append(png)
    print(f'  wrote {png}')

print(f'\n{len(made)} panel(s) written to {OUT}/')
for s in skipped:
    print('  skipped', s)
if not made:
    print('\nNothing rendered. Check that results/tmp_poses/ has *_dry1.pdb / *_wet1.pdb\n'
          '(run scripts/40_allposes.py first) and that structures/xdock/<tag>.pdb exists.')
