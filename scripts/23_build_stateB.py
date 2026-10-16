import numpy as np

def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P-Pc).T @ (Q-Qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1., 1., d]) @ U.T
    return R, Qc - R @ Pc

def ca(path):
    d = {}
    for l in open(path):
        if l[:4]=='ATOM' and l[21]=='A' and l[12:16].strip()=='CA' and l[16] in ' A':
            d[l[22:27]] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
    return d

def xyz(l): return np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
def put(l, c): return l[:30] + f"{c[0]:8.3f}{c[1]:8.3f}{c[2]:8.3f}" + l[54:]

SRC, LIG = 'structures/pdb/7KG1.pdb', 'WDD'
a, b = ca('structures/pdb/6OAK.pdb'), ca(SRC)
k = sorted(set(a) & set(b))
R, t = kabsch(np.array([b[x] for x in k]), np.array([a[x] for x in k]))
print(f"{len(k)} CA matched")

lines = [l for l in open(SRC) if l[:6] in ('ATOM  ','HETATM')
         and l[21]=='A' and l[16] in ' A']
lig = np.array([xyz(l) for l in lines
                if l[17:20].strip()==LIG and l[76:78].strip()!='H'])
pol = np.array([xyz(l) for l in lines
                if l[:4]=='ATOM' and l[76:78].strip() in ('N','O')])
print(f"ligand {len(lig)} atoms, polar protein atoms {len(pol)}")

keep = set()
for i, l in enumerate(lines):
    if l[17:20].strip()=='HOH' and l[12:16].strip()=='O':
        c = xyz(l)
        if (np.linalg.norm(lig-c, axis=1).min() < 3.5
                and np.linalg.norm(pol-c, axis=1).min() < 3.5):
            keep.add(i)
print(f"bridging waters: {len(keep)}")

with open('structures/rec_B_wet.pdb','w') as f:
    for i, l in enumerate(lines):
        if l[:4]=='ATOM' or i in keep:
            f.write(put(l, (xyz(l) @ R.T) + t))
    f.write('END\n')
print("wrote structures/rec_B_wet.pdb")
