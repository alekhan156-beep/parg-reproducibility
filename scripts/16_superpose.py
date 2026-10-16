import os
import numpy as np
from collections import defaultdict

BOX = np.array([27.20, 2.84, 0.17])
SKIPRES = {'HOH','GOL','EDO','SO4','PO4','CL','NA','MG','ZN','CA','ACT','DMS',
           'PEG','PG4','TRS','MES','EPE','CME','FMT','IOD','K','BR','NO3','CIT',
           'MSE','CSO','CSS','BME','DTV','PPI','CAC','DTT','6WK'}

def ca_atoms(path):
    out = {}
    for line in open(path):
        if line.startswith('ATOM') and line[21] == 'A' and line[12:16].strip() == 'CA':
            if line[16] not in (' ', 'A'):
                continue
            out[line[22:27]] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    return out

def het_groups(path):
    g = defaultdict(lambda: {'xyz': [], 'el': []})
    for line in open(path):
        if line.startswith('HETATM') and line[21] == 'A':
            res = line[17:20].strip()
            if res in SKIPRES or line[16] not in (' ', 'A'):
                continue
            el = line[76:78].strip().capitalize()
            if el == 'H':
                continue
            k = (res, line[22:27].strip())
            g[k]['xyz'].append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            g[k]['el'].append(el)
    return g

def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, Qc - R @ Pc

def write_pdb(path, xyz, els, resname):
    with open(path, 'w') as f:
        for i, (c, e) in enumerate(zip(xyz, els), 1):
            nm = f"{e}{i}"
            f.write(f"HETATM{i:5d} {nm:<4s} {resname:>3s} A1001    "
                    f"{c[0]:8.3f}{c[1]:8.3f}{c[2]:8.3f}  1.00  0.00          {e:>2s}\n")
        f.write("END\n")

os.makedirs('structures/xdock', exist_ok=True)
ref_ca = ca_atoms('structures/pdb/6OAK.pdb') if os.path.exists('structures/pdb/6OAK.pdb') \
         else ca_atoms('structures/6OAK_receptor_dry.pdb')
print(f"reference 6OAK: {len(ref_ca)} CA atoms\n")

print(f"{'PDB':6} {'lig':5} {'n':>4} {'type':12} {'CA-rmsd':>8} {'nCA':>5} {'to-box':>8}")
print("-" * 56)
rows = []
for fn in sorted(os.listdir('structures/pdb')):
    if not fn.endswith('.pdb'):
        continue
    pid = fn[:-4]
    if pid == '6OAK':
        continue
    path = f'structures/pdb/{fn}'
    tgt_ca = ca_atoms(path)
    shared = sorted(set(ref_ca) & set(tgt_ca))
    if len(shared) < 100:
        print(f"{pid:6} too few shared CA ({len(shared)})")
        continue
    P = np.array([tgt_ca[k] for k in shared])
    Q = np.array([ref_ca[k] for k in shared])
    R, t = kabsch(P, Q)
    ca_rmsd = np.sqrt((((P @ R.T) + t - Q) ** 2).sum(1).mean())
    for (res, num), d in sorted(het_groups(path).items(), key=lambda k: -len(k[1]['xyz'])):
        n = len(d['xyz'])
        if n < 15 or n > 50:
            continue
        xyz = (np.array(d['xyz']) @ R.T) + t
        kind = 'substrate' if 'P' in d['el'] else 'inhibitor'
        dist = np.linalg.norm(xyz.mean(0) - BOX)
        print(f"{pid:6} {res:5} {n:4d} {kind:12} {ca_rmsd:8.2f} {len(shared):5d} {dist:8.2f}")
        write_pdb(f'structures/xdock/{pid}_{res}.pdb', xyz, d['el'], res)
        rows.append((pid, res, kind, dist))

print(f"\nwrote {len(rows)} superposed ligands to structures/xdock/")
inh = [r for r in rows if r[2] == 'inhibitor' and r[3] < 8]
print(f"inhibitors in the active site (<8 A from box centre): {len(inh)}")
