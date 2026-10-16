from pymol import cmd
import numpy as np
cmd.feedback('disable','all','everything')
WAT=[1211,1232,1307,1317,1402]
cmd.load('pdb/6OAK.pdb','ref'); cmd.remove('ref and not chain A')
wsel=' or '.join(f'(resi {r})' for r in WAT)
wat=cmd.get_coords(f'ref and resn HOH and ({wsel})')
wres=[a.resi for a in cmd.get_model(f'ref and resn HOH and ({wsel})').atom]

for pdb,lig in (('7KFP','WD7'),('7KG8','WDM'),('6HMK','7JC')):
    cmd.load(f'pdb/{pdb}.pdb','m'); cmd.remove('m and not chain A')
    r=cmd.align('m and name CA','ref and name CA',cycles=5)
    print(f'\n=== {pdb} / {lig} ===')
    print(f'  align RMSD {r[0]:.3f} A over {r[1]} atoms')
    # copies and altlocs
    mdl=cmd.get_model(f'm and resn {lig}')
    print(f'  {lig} atoms={len(mdl.atom)}  chains={set(a.chain for a in mdl.atom)} '
          f'resi={sorted(set(a.resi for a in mdl.atom))} '
          f'altlocs={sorted(set(a.alt for a in mdl.atom))}')
    lc=cmd.get_coords(f'm and resn {lig} and not elem H')
    names=[a.name for a in cmd.get_model(f'm and resn {lig} and not elem H').atom]
    D=np.sqrt(((lc[:,None,:]-wat[None,:,:])**2).sum(-1))
    i,j=np.unravel_index(D.argmin(),D.shape)
    print(f'  closest contact: ligand atom {names[i]} -- HOH {wres[j]}  {D[i,j]:.2f} A')
    print(f'  per-water minima: ' +
          '  '.join(f'{wres[k]}:{D[:,k].min():.2f}' for k in range(len(wres))))
    cmd.delete('m')
