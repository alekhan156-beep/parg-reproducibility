import os, subprocess, csv
import numpy as np
from scipy.optimize import linear_sum_assignment

def parse(path, model=1, pdbqt=False):
    xyz, els, cur = [], [], 0
    for l in open(path):
        if l.startswith('MODEL'): cur = int(l.split()[1])
        if l.startswith(('ATOM','HETATM')):
            if cur and cur != model: continue
            raw = l[76:78].strip()
            if pdbqt:
                t = raw.upper()
                el = ('H' if t.startswith('H') else 'C' if t in ('A','C')
                      else 'N' if t in ('N','NA','NS') else 'O' if t in ('O','OA','OS')
                      else 'S' if t in ('S','SA') else t.capitalize())
            else:
                el = raw.capitalize()
            if el == 'H': continue
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            els.append(el)
    return np.array(xyz), els

def floor_rmsd(p, pe, r, re_):
    c = ((p[:,None,:] - r[None,:,:])**2).sum(-1)
    for i,x in enumerate(pe):
        for j,y in enumerate(re_):
            if x != y: c[i,j] = 1e6
    i,j = linear_sum_assignment(c)
    return float(np.sqrt(c[i,j].mean()))

def aff(path):
    return [float(l.split()[3]) for l in open(path)
            if l.startswith('REMARK VINA RESULT')]

tags = [l.strip() for l in open('results/xdock_targets.txt') if l.strip()]
A = {}
for r in csv.DictReader(open('results/crossdock.csv')):
    if r['receptor'] == 'wet':
        A[r['ligand']] = (float(r['top1_rmsd']), float(r['top1_kcal']))

print(f"{'ligand':14} {'A_rmsd':>7} {'A_kcal':>7} {'B_rmsd':>7} {'B_kcal':>7} {'pick':>5} {'ens':>6}")
print("-" * 62)
rows = []
for tag in tags:
    ref, rel = parse(f'structures/xdock/{tag}.pdb')
    dst = f'results/xdock/{tag}_B.pdbqt'
    if not os.path.exists(dst):
        subprocess.run(['./vina_1.2.5_linux_x86_64',
            '--receptor', 'structures/rec_B.pdbqt',
            '--ligand', f'structures/xdock/{tag}.pdbqt',
            '--center_x','27.20','--center_y','2.84','--center_z','0.17',
            '--size_x','24','--size_y','24','--size_z','24',
            '--exhaustiveness','16','--num_modes','9','--seed','42',
            '--out', dst], capture_output=True)
    pose, pel = parse(dst, model=1, pdbqt=True)
    br = floor_rmsd(pose, pel, ref, rel) if len(pose)==len(ref) else np.nan
    bk = aff(dst)[0]
    ar, ak = A[tag]
    pick = 'A' if ak <= bk else 'B'
    ens = ar if pick == 'A' else br
    rows.append((tag, ar, ak, br, bk, pick, ens, min(ar, br)))
    print(f"{tag:14} {ar:7.2f} {ak:7.2f} {br:7.2f} {bk:7.2f} {pick:>5} {ens:6.2f}", flush=True)

a = np.array([r[1] for r in rows]); b = np.array([r[3] for r in rows])
e = np.array([r[6] for r in rows]); o = np.array([r[7] for r in rows])
n = len(rows)
print("\n" + "=" * 62)
print(f"state A (6OAK) alone      top-1 under 2A: {(a<2).sum():2d}/{n}")
print(f"state B (7KG1) alone      top-1 under 2A: {(b<2).sum():2d}/{n}")
print(f"ensemble, score-selected  top-1 under 2A: {(e<2).sum():2d}/{n}")
print(f"oracle (upper bound)      top-1 under 2A: {(o<2).sum():2d}/{n}")
print("=" * 62)

with open('results/ensemble.csv','w',newline='') as f:
    w = csv.writer(f)
    w.writerow(['ligand','A_rmsd','A_kcal','B_rmsd','B_kcal','pick','ensemble','oracle'])
    w.writerows(rows)
print("wrote results/ensemble.csv")
