exec(open('scripts/24_ensemble.py').read().split('\ntags =')[0])
import os, subprocess, csv
import numpy as np

os.makedirs('results/box16', exist_ok=True)
tags = [l.strip() for l in open('results/xdock_targets.txt') if l.strip()]
print(f"{'ligand':14} {'A16':>7} {'B16':>7} {'pick':>5} {'ens':>7}")
rows = []
for tag in tags:
    ref, rel = parse(f'structures/xdock/{tag}.pdb')
    c = ref.mean(0)
    res = {}
    for rec, rp in (('A','structures/rec_wet.pdbqt'), ('B','structures/rec_B.pdbqt')):
        dst = f'results/box16/{tag}_{rec}.pdbqt'
        if not os.path.exists(dst):
            subprocess.run(['./vina_1.2.5_linux_x86_64','--receptor',rp,
                '--ligand',f'structures/xdock/{tag}.pdbqt',
                '--center_x',f'{c[0]:.2f}','--center_y',f'{c[1]:.2f}',
                '--center_z',f'{c[2]:.2f}',
                '--size_x','16','--size_y','16','--size_z','16',
                '--exhaustiveness','16','--num_modes','9','--seed','42',
                '--out',dst], capture_output=True)
        p, pe = parse(dst, model=1, pdbqt=True)
        res[rec] = (floor_rmsd(p,pe,ref,rel) if len(p)==len(ref) else np.nan, aff(dst)[0])
    pk = 'A' if res['A'][1] <= res['B'][1] else 'B'
    rows.append((tag, res['A'][0], res['B'][0], pk, res[pk][0], min(res['A'][0],res['B'][0])))
    print(f"{tag:14} {res['A'][0]:7.2f} {res['B'][0]:7.2f} {pk:>5} {res[pk][0]:7.2f}", flush=True)

a=np.array([r[1] for r in rows]); b=np.array([r[2] for r in rows])
e=np.array([r[4] for r in rows]); o=np.array([r[5] for r in rows]); n=len(rows)
print(f"\n24A box, A        : 7/{n}  (previous)")
print(f"16A box, A        : {(a<2).sum():2d}/{n}")
print(f"16A box, B        : {(b<2).sum():2d}/{n}")
print(f"16A box, ensemble : {(e<2).sum():2d}/{n}")
print(f"16A box, oracle   : {(o<2).sum():2d}/{n}")
csv.writer(open('results/box16.csv','w',newline='')).writerows(
    [['ligand','A','B','pick','ens','oracle']] + [list(r) for r in rows])
