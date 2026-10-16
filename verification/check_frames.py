import sys, os
import numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
ROOT = sys.argv[1]
def load(p):
    m = Chem.MolFromPDBFile(p, removeHs=True, sanitize=False)
    if m is None: return None
    return m.GetConformer().GetPositions()
print(f"{'ligand':12} {'crystal centroid':>26} {'dry1 off':>9} {'wet1 off':>9}   expect dry/wet")
EXP={'6O9X_M0S':'3.74/0.62','6OA1_M0P':'3.68/0.69','6OAL_M0P':'3.66/0.57',
     '7KFP_WD7':'2.73/1.94','7KG7_WDJ':'4.59/0.74','7KG8_WDM':'3.36/0.62',
     '7KG6_WDG':'3.63/6.42','6HMK_7JC':'4.38/5.32'}
for tag in EXP:
    r = load(f'{ROOT}/structures/xdock/{tag}.pdb')
    if r is None: print(f'{tag:12}  reference unreadable'); continue
    c = r.mean(0)
    out=[]
    for rec in ('dry','wet'):
        p = load(f'{ROOT}/results/tmp_poses/{tag}_{rec}1.pdb')
        out.append(np.linalg.norm(p.mean(0)-c) if p is not None else float('nan'))
    print(f'{tag:12} ({c[0]:7.2f},{c[1]:6.2f},{c[2]:6.2f}) {out[0]:9.2f} {out[1]:9.2f}   {EXP[tag]}')
