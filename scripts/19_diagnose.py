import os, csv
import numpy as np

BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}

def prot_atoms(path):
    d = {}
    for line in open(path):
        if line.startswith('ATOM') and line[21] == 'A' and line[16] in (' ', 'A'):
            if line[76:78].strip() == 'H':
                continue
            d[(line[22:27], line[12:16].strip())] = (
                line[17:20].strip(),
                np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    return d

def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, Qc - R @ Pc

def lig_coords(path):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in open(path)
                     if l.startswith(('ATOM', 'HETATM')) and l[76:78].strip() != 'H'])

ref = prot_atoms('structures/pdb/6OAK.pdb')
m0m = lig_coords('structures/6OAK_ligand.pdb')

waters = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                   for l in open('structures/6OAK_receptor_wet.pdb')
                   if l[17:20].strip() == 'HOH' and l[12:16].strip() == 'O'])
print(f"kept waters: {len(waters)}")

site = {k[0] for k, (rn, xyz) in ref.items()
        if np.linalg.norm(m0m - xyz, axis=1).min() < 6.0}
print(f"binding-site residues within 6 A of M0M: {len(site)}\n")

wet = {}
for r in csv.DictReader(open('results/crossdock.csv')):
    if r['receptor'] == 'wet':
        wet[r['ligand']] = (float(r['top1_rmsd']), float(r['best_rmsd']))

print(f"{'ligand':14} {'ok':>3} {'sideRMSD':>9} {'nAt':>4} {'maxShift':>9} {'H2Oclash':>9}")
print("-" * 54)
rows = []
for tag in sorted(wet):
    pid = tag.split('_')[0]
    tgt = prot_atoms(f'structures/pdb/{pid}.pdb')
    shared = sorted({k[0] for k in ref} & {k[0] for k in tgt})
    P = np.array([tgt[(k, 'CA')][1] for k in shared if (k, 'CA') in tgt and (k, 'CA') in ref])
    Q = np.array([ref[(k, 'CA')][1] for k in shared if (k, 'CA') in tgt and (k, 'CA') in ref])
    R, t = kabsch(P, Q)

    devs = []
    for (rk, an), (rn, xyz) in tgt.items():
        if rk not in site or an in BACKBONE:
            continue
        if (rk, an) not in ref or ref[(rk, an)][0] != rn:
            continue
        devs.append(np.linalg.norm((xyz @ R.T) + t - ref[(rk, an)][1]))
    devs = np.array(devs) if devs else np.array([np.nan])

    lig = lig_coords(f'structures/xdock/{tag}.pdb')
    clash = np.linalg.norm(lig[:, None, :] - waters[None, :, :], axis=-1).min() if len(waters) else np.nan

    ok = "Y" if wet[tag][1] < 2.0 else "n"
    print(f"{tag:14} {ok:>3} {np.sqrt((devs**2).mean()):9.2f} {len(devs):4d} "
          f"{devs.max():9.2f} {clash:9.2f}")
    rows.append((ok, np.sqrt((devs**2).mean()), devs.max(), clash))

for lab, want in (("found   ", "Y"), ("not found", "n")):
    s = [r for r in rows if r[0] == want]
    if s:
        print(f"\n{lab} (n={len(s)}): mean side-chain RMSD {np.mean([r[1] for r in s]):.2f} A"
              f"   mean max shift {np.mean([r[2] for r in s]):.2f} A"
              f"   mean closest water {np.mean([r[3] for r in s]):.2f} A")
