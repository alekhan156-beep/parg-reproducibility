import sys
import numpy as np
from scipy.optimize import linear_sum_assignment

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

def affinities(path):
    return [float(l.split()[3]) for l in open(path) if l.startswith('REMARK VINA RESULT')]

ref, ref_el = parse('structures/6OAK_ligand.pdb')
print(f"crystal reference: {len(ref)} heavy atoms\n")

for path in sys.argv[1:]:
    print(f"=== {path} ===")
    aff = affinities(path)
    rows = []
    for mode in range(1, 10):
        pose, pose_el = parse(path, model=mode, pdbqt=True)
        if len(pose) == 0:
            break
        cen = np.linalg.norm(pose.mean(0) - ref.mean(0))
        cost = ((pose[:, None, :] - ref[None, :, :]) ** 2).sum(-1)
        for i, a in enumerate(pose_el):
            for j, b in enumerate(ref_el):
                if a != b:
                    cost[i, j] = 1e6
        r, c = linear_sum_assignment(cost)
        rmsd = np.sqrt(cost[r, c].mean())
        rows.append((mode, cen, rmsd))
        kcal = aff[mode - 1] if mode <= len(aff) else float('nan')
        print(f"  mode {mode}: {kcal:7.2f} kcal/mol   centroid {cen:5.2f} A   RMSD floor {rmsd:5.2f} A")
    best = min(rows, key=lambda x: x[2])
    print(f"  --> closest to crystal: mode {best[0]} at {best[2]:.2f} A (Vina ranked it #{best[0]} of {len(rows)})\n")
