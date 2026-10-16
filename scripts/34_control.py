import os, random, subprocess, csv, time
import numpy as np

VINA = './vina_1.2.5_linux_x86_64'
REC = 'structures/rec_wet.pdbqt'
N = 100
os.makedirs('results/control', exist_ok=True)
os.makedirs('ligands', exist_ok=True)

top = []
for line in open('results/screen_ranked.smi'):
    p = line.split()
    if len(p) >= 2: top.append((p[0], p[1]))
top = top[:N]
topids = {z for _, z in top}

lib = []
for line in open('data/zinc_instock.smi'):
    p = line.split()
    if len(p) >= 2 and p[1] not in topids: lib.append((p[0], p[1]))
random.seed(42)
rnd = random.sample(lib, N)
print(f"{len(top)} top-ranked vs {len(rnd)} random, from {len(lib):,} candidates")

jobs = []
for i in range(N):
    jobs.append(('top', top[i][0], top[i][1]))
    jobs.append(('random', rnd[i][0], rnd[i][1]))
print(f"{len(jobs)} docking jobs, interleaved\n")

def aff(p):
    for l in open(p):
        if l.startswith('REMARK VINA RESULT'): return float(l.split()[3])
    return None

rows, t0 = [], time.time()
for k, (grp, smi, zid) in enumerate(jobs, 1):
    lig = f"ligands/{zid}.pdbqt"
    if not os.path.exists(lig):
        subprocess.run(['obabel', f'-:{smi}', '-opdbqt', '-O', lig, '--gen3d', '-h'],
                       capture_output=True)
    if not os.path.exists(lig) or os.path.getsize(lig) == 0:
        print(f"  [{k}] {zid} {grp:7} 3D failed"); continue
    out = f"results/control/{zid}.pdbqt"
    if not os.path.exists(out):
        subprocess.run([VINA, '--receptor', REC, '--ligand', lig,
            '--center_x','27.20','--center_y','2.84','--center_z','0.17',
            '--size_x','16','--size_y','16','--size_z','16',
            '--exhaustiveness','16','--num_modes','9','--seed','42',
            '--out', out], capture_output=True)
    a = aff(out) if os.path.exists(out) else None
    if a is None:
        print(f"  [{k}] {zid} {grp:7} dock failed"); continue
    rows.append({'group': grp, 'zinc_id': zid, 'smiles': smi, 'vina_kcal': a})
    el = (time.time()-t0)/60
    print(f"  [{k}/{len(jobs)}] {grp:7} {zid:14} {a:7.2f}   "
          f"{el:5.1f} min, eta {el/k*(len(jobs)-k):5.1f} min", flush=True)
    with open('results/control.csv','w',newline='') as f:
        w = csv.DictWriter(f, fieldnames=['group','zinc_id','smiles','vina_kcal'])
        w.writeheader(); w.writerows(rows)

from scipy.stats import mannwhitneyu
t = np.array([r['vina_kcal'] for r in rows if r['group']=='top'])
r_ = np.array([r['vina_kcal'] for r in rows if r['group']=='random'])
print("\n" + "="*58)
print(f"{'':10} {'n':>4} {'mean':>8} {'median':>8} {'best':>8} {'<-9':>6}")
for lab, v in (('top-ranked', t), ('random', r_)):
    print(f"{lab:10} {len(v):4d} {v.mean():8.2f} {np.median(v):8.2f} "
          f"{v.min():8.2f} {(v < -9).sum():6d}")
d = r_.mean() - t.mean()
u, p = mannwhitneyu(t, r_, alternative='less')
print(f"\ntop-ranked are {d:+.2f} kcal/mol better on average")
print(f"Mann-Whitney U p = {p:.2g}")
print("VERDICT:", "model beat random" if p < 0.05 else
      "NO significant difference - model added nothing on this library")
print("="*58)
