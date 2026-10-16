import numpy as np

def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P-Pc).T @ (Q-Qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1., 1., d]) @ U.T
    return R, Qc - R @ Pc

def ca(p):
    d = {}
    for l in open(p):
        if l[:4]=='ATOM' and l[21]=='A' and l[12:16].strip()=='CA' and l[16] in ' A':
            d[l[22:27]] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
    return d

a, b = ca('structures/pdb/6OAK.pdb'), ca('structures/pdb/7KG1.pdb')
k = sorted(set(a) & set(b))
R, t = kabsch(np.array([b[x] for x in k]), np.array([a[x] for x in k]))

wat = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                for l in open('structures/pdb/7KG1.pdb')
                if l[:6]=='HETATM' and l[17:20].strip()=='HOH'
                and l[12:16].strip()=='O' and l[21]=='A'])
wat = (wat @ R.T) + t
print(f"7KG1 crystallographic waters: {len(wat)}")

site  = np.array([28.70,  2.42,  3.30])
decoy = np.array([29.44, -0.38, -5.78])
for nm, p in (('true site ', site), ('decoy site', decoy)):
    d = np.linalg.norm(wat - p, axis=1)
    print(f"{nm}  within 5A: {(d<5).sum():3d}    within 8A: {(d<8).sum():3d}")

near = set()
for l in open('structures/pdb/7KG1.pdb'):
    if l[:4]=='ATOM' and l[21]=='A' and l[16] in ' A' and l[76:78].strip()!='H':
        c = (np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) @ R.T) + t
        if np.linalg.norm(c - decoy) < 6:
            near.add(f"{l[17:20].strip()}{l[22:26].strip()}")
print("residues lining decoy site:", " ".join(sorted(near)))
