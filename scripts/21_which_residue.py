import csv
import numpy as np
from collections import defaultdict

BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}
CATALYTIC = {'737', '740', '754', '755', '756', '795', '902'}

def prot_atoms(path):
    d = {}
    for line in open(path):
        if line.startswith('ATOM') and line[21] == 'A' and line[16] in (' ', 'A'):
            if line[76:78].strip() == 'H':
                continue
            b = line[60:66].strip()
            d[(line[22:27], line[12:16].strip())] = (
                line[17:20].strip(),
                np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                float(b) if b else 0.0)
    return d

def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - Pc).T @ (Q - Qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, Qc - R @ Pc

ref = prot_atoms('structures/pdb/6OAK.pdb')
m0m = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                for l in open('structures/6OAK_ligand.pdb')
                if l.startswith(('ATOM', 'HETATM')) and l[76:78].strip() != 'H'])
site = {k[0] for k, v in ref.items() if np.linalg.norm(m0m - v[1], axis=1).min() < 6.0}

found = {}
for r in csv.DictReader(open('results/crossdock.csv')):
    if r['receptor'] == 'wet':
        found[r['ligand']] = float(r['best_rmsd']) < 2.0

shifts = defaultdict(dict)
for tag in sorted(found):
    tgt = prot_atoms(f"structures/pdb/{tag.split('_')[0]}.pdb")
    keys = [k for k in ref if k[1] == 'CA' and k in tgt]
    R, t = kabsch(np.array([tgt[k][1] for k in keys]),
                  np.array([ref[k][1] for k in keys]))
    for (rk, an), (rn, xyz, b) in tgt.items():
        if rk not in site or an in BACKBONE or (rk, an) not in ref:
            continue
        if ref[(rk, an)][0] != rn:
            continue
        d = np.linalg.norm((xyz @ R.T) + t - ref[(rk, an)][1])
        key = f"{rn}{rk.strip()}"
        shifts[tag][key] = max(shifts[tag].get(key, 0), d)

print("top 3 moving side chains per ligand\n")
for tag in sorted(found, key=lambda x: found[x]):
    top = sorted(shifts[tag].items(), key=lambda kv: -kv[1])[:3]
    mark = "Y" if found[tag] else "n"
    print(f"{tag:14} {mark}  " + "   ".join(f"{k}:{v:.1f}" for k, v in top))

print("\n" + "=" * 62)
print(f"{'residue':12} {'found':>8} {'notfound':>9} {'gap':>7} {'refB':>7}")
print("-" * 62)
allres = sorted({k for s in shifts.values() for k in s})
table = []
for res in allres:
    f = [shifts[t][res] for t in shifts if found[t] and res in shifts[t]]
    n = [shifts[t][res] for t in shifts if not found[t] and res in shifts[t]]
    if not f or not n:
        continue
    num = ''.join(c for c in res if c.isdigit())
    bs = [v[2] for k, v in ref.items() if k[0].strip() == num and k[1] not in BACKBONE]
    table.append((np.mean(n) - np.mean(f), res, np.mean(f), np.mean(n),
                  np.mean(bs) if bs else 0.0, num))
for gap, res, f, n, b, num in sorted(table, reverse=True)[:8]:
    star = " *catalytic" if num in CATALYTIC else ""
    print(f"{res:12} {f:8.2f} {n:9.2f} {gap:7.2f} {b:7.1f}{star}")
