import os, subprocess, csv, time
import numpy as np
from scipy.optimize import linear_sum_assignment

VINA = './vina_1.2.5_linux_x86_64'
XD, OUT = 'structures/xdock', 'results/xdock'
os.makedirs(OUT, exist_ok=True)

def parse(path, model=1, pdbqt=False):
    xyz, els, cur = [], [], 0
    for line in open(path):
        if line.startswith('MODEL'):
            cur = int(line.split()[1])
        if line.startswith(('ATOM', 'HETATM')):
            if cur and cur != model:
                continue
            raw = line[76:78].strip()
            if pdbqt:
                t = raw.upper()
                el = ('H' if t.startswith('H') else 'C' if t in ('A', 'C')
                      else 'N' if t in ('N', 'NA', 'NS') else 'O' if t in ('O', 'OA', 'OS')
                      else 'S' if t in ('S', 'SA') else t.capitalize())
            else:
                el = raw.capitalize()
            if el == 'H':
                continue
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            els.append(el)
    return np.array(xyz), els

def floor_rmsd(pose, pel, ref, rel):
    cost = ((pose[:, None, :] - ref[None, :, :]) ** 2).sum(-1)
    for i, a in enumerate(pel):
        for j, b in enumerate(rel):
            if a != b:
                cost[i, j] = 1e6
    r, c = linear_sum_assignment(cost)
    return float(np.sqrt(cost[r, c].mean()))

def affinities(path):
    return [float(l.split()[3]) for l in open(path) if l.startswith('REMARK VINA RESULT')]

tags = [l.strip() for l in open('results/xdock_targets.txt') if l.strip()]
print(f"{len(tags)} ligands x 2 receptors = {len(tags) * 2} runs\n")
print(f"{'ligand':14} {'dry#1':>7} {'wet#1':>7} {'dryBest':>8} {'wetBest':>8} {'elapsed':>8}")
print("-" * 58)

records, t0 = [], time.time()
for tag in tags:
    ref, rel = parse(f'{XD}/{tag}.pdb')
    row = {}
    for rec in ('dry', 'wet'):
        dst = f'{OUT}/{tag}_{rec}.pdbqt'
        if not os.path.exists(dst):
            subprocess.run([VINA,
                '--receptor', f'structures/rec_{rec}.pdbqt',
                '--ligand', f'{XD}/{tag}.pdbqt',
                '--center_x', '27.20', '--center_y', '2.84', '--center_z', '0.17',
                '--size_x', '24', '--size_y', '24', '--size_z', '24',
                '--exhaustiveness', '16', '--num_modes', '9', '--seed', '42',
                '--out', dst], capture_output=True)
        aff = affinities(dst)
        vals = []
        for m in range(1, 10):
            pose, pel = parse(dst, model=m, pdbqt=True)
            if len(pose) == 0:
                break
            vals.append(floor_rmsd(pose, pel, ref, rel) if len(pose) == len(ref) else np.nan)
        best = int(np.nanargmin(vals)) + 1
        row[rec] = vals
        records.append({'ligand': tag, 'receptor': rec, 'n_heavy': len(ref),
                        'top1_rmsd': round(vals[0], 3),
                        'top1_kcal': round(aff[0], 2) if aff else '',
                        'best_rmsd': round(vals[best - 1], 3),
                        'best_rank': best, 'n_modes': len(vals)})
    print(f"{tag:14} {row['dry'][0]:7.2f} {row['wet'][0]:7.2f} "
          f"{np.nanmin(row['dry']):8.2f} {np.nanmin(row['wet']):8.2f} "
          f"{(time.time()-t0)/60:7.1f}m", flush=True)

with open('results/crossdock.csv', 'w', newline='') as f:
    wr = csv.DictWriter(f, fieldnames=list(records[0].keys()))
    wr.writeheader()
    wr.writerows(records)

print("\n" + "=" * 58)
for rec in ('dry', 'wet'):
    s = [r for r in records if r['receptor'] == rec]
    t1 = np.array([r['top1_rmsd'] for r in s])
    bb = np.array([r['best_rmsd'] for r in s])
    rk = np.array([r['best_rank'] for r in s])
    print(f"{rec.upper():4}  top-1 under 2A: {(t1 < 2).sum():2d}/{len(t1)}"
          f"   median top-1: {np.median(t1):.2f} A"
          f"   right pose sampled: {(bb < 2).sum():2d}/{len(bb)}"
          f"   ranked #1: {(rk == 1).sum():2d}/{len(rk)}")
print("=" * 58)
