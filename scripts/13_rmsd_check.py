import numpy as np
from collections import Counter
from scipy.optimize import linear_sum_assignment

def parse(path, model=1, pdbqt=False):
    xyz, els, names, cur = [], [], [], 0
    for line in open(path):
        if line.startswith('MODEL'):
            cur = int(line.split()[1])
        if line.startswith(('ATOM', 'HETATM')):
            if cur and cur != model:
                continue
            raw = line[76:78].strip()
            name = line[12:16].strip()
            if pdbqt:
                t = raw.upper()
                if t.startswith('H'): el = 'H'
                elif t in ('A', 'C'): el = 'C'
                elif t in ('N', 'NA', 'NS'): el = 'N'
                elif t in ('O', 'OA', 'OS'): el = 'O'
                elif t in ('S', 'SA'): el = 'S'
                elif t == 'CL': el = 'Cl'
                elif t == 'BR': el = 'Br'
                else: el = t.capitalize()
            else:
                el = raw.capitalize() if raw else name[0]
            if el == 'H':
                continue
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            els.append(el); names.append(name)
    return np.array(xyz), els, names

ref, ref_el, ref_nm = parse('structures/6OAK_ligand.pdb')
print(f"crystal : {len(ref)} heavy atoms, centroid {ref.mean(0).round(2)}")
print(f"  comp  : {dict(Counter(ref_el))}")
print(f"  order : {' '.join(ref_nm)}\n")

for mode in range(1, 4):
    pose, pose_el, pose_nm = parse('results/redock_dry.pdbqt', model=mode, pdbqt=True)
    if len(pose) == 0:
        break
    cen = np.linalg.norm(pose.mean(0) - ref.mean(0))
    naive = np.sqrt(((pose - ref) ** 2).sum(1).mean()) if len(pose) == len(ref) else float('nan')
    cost = ((pose[:, None, :] - ref[None, :, :]) ** 2).sum(-1)
    for i, a in enumerate(pose_el):
        for j, b in enumerate(ref_el):
            if a != b:
                cost[i, j] = 1e6
    r, c = linear_sum_assignment(cost)
    matched = np.sqrt(cost[r, c].mean())
    print(f"mode {mode}: centroid offset {cen:5.2f} A | in-order RMSD {naive:5.2f} A | matched RMSD {matched:5.2f} A")
    if mode == 1:
        print(f"  comp  : {dict(Counter(pose_el))}")
        print(f"  order : {' '.join(pose_nm)}")
