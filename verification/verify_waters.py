"""Independently verify section 12: distance from each superposed ligand to the
nearest retained bridging water in 6OAK."""
from pymol import cmd
import numpy as np

WAT = [1211, 1232, 1307, 1317, 1402]
LIG = {'6OAK':'M0M','6O9X':'M0S','6OA1':'M0P','6OAL':'M0P','7KFP':'WD7',
       '7KG7':'WDJ','7KG8':'WDM','7KG6':'WDG','6HMK':'7JC'}
DISPLACING = {'6HMK'}          # per sec 12, the 5LHB/6HMx series

cmd.feedback('disable','all','everything')
cmd.load('pdb/6OAK.pdb','ref')
cmd.remove('ref and not chain A')

wsel = ' or '.join(f'(resi {r})' for r in WAT)
wat = cmd.get_coords(f'ref and resn HOH and ({wsel})')
print(f'retained waters found in 6OAK chain A: '
      f'{0 if wat is None else len(wat)} of {len(WAT)}')
if wat is None or len(wat) == 0:
    wat = cmd.get_coords('ref and resn HOH')
    print(f'  (falling back to all {len(wat)} waters in chain A)')

print(f"\n{'PDB':6} {'lig':5} {'min dist to retained water':>27}   class")
print('-'*62)
rows=[]
for pdb, lig in LIG.items():
    if pdb == '6OAK':
        cmd.create('mob','ref'); 
    else:
        cmd.load(f'pdb/{pdb}.pdb','mob')
        cmd.remove('mob and not chain A')
        cmd.align('mob and name CA','ref and name CA', cycles=5)
    lc = cmd.get_coords(f'mob and resn {lig} and not elem H')
    if lc is None or len(lc)==0:
        print(f'{pdb:6} {lig:5}   ligand not found in chain A'); cmd.delete('mob'); continue
    d = np.sqrt(((lc[:,None,:]-wat[None,:,:])**2).sum(-1)).min()
    cls = 'water-DISPLACING' if d < 2.0 else 'water-bridging'
    rows.append((pdb,lig,d,cls))
    print(f'{pdb:6} {lig:5} {d:>22.2f} A   {cls}')
    cmd.delete('mob')

br=[r[2] for r in rows if r[3]=='water-bridging']
dp=[r[2] for r in rows if r[3]=='water-DISPLACING']
print(f'\nbridging   n={len(br)}  range {min(br):.2f}-{max(br):.2f} A'
      f'   (record says 2.59-2.78)')
if dp: print(f'displacing n={len(dp)}  range {min(dp):.2f}-{max(dp):.2f} A'
             f'   (record says 0.82-1.27)')
