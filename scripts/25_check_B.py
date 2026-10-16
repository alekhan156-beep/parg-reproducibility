import numpy as np

def ca(path):
    d = {}
    for l in open(path):
        if l[:4]=='ATOM' and l[21]=='A' and l[12:16].strip()=='CA' and l[16] in ' A':
            d[l[22:27]] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return d

def atoms(path):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in open(path) if l[:6] in ('ATOM  ','HETATM')])

a = ca('structures/pdb/6OAK.pdb')
b = ca('structures/6OAK_receptor_dry.pdb')
k = sorted(set(a) & set(b))
d = np.array([np.linalg.norm(a[x]-b[x]) for x in k])
print(f"6OAK download vs receptor_dry : {len(k)} CA  mean {d.mean():.3f} A  max {d.max():.3f} A")

c = ca('structures/rec_B_wet.pdb')
k2 = sorted(set(a) & set(c))
d2 = np.array([np.linalg.norm(a[x]-c[x]) for x in k2])
print(f"rec_B vs 6OAK frame          : {len(k2)} CA  mean {d2.mean():.2f} A  max {d2.max():.2f} A")

lig = atoms('structures/xdock/7KG1_WDD.pdb')
rec = atoms('structures/rec_B_wet.pdb')
dist = np.linalg.norm(lig[:,None,:] - rec[None,:,:], axis=-1)
print(f"WDD reference to rec_B       : min {dist.min():.2f} A   contacts<4A {(dist<4).sum()}")
print(f"WDD centroid {lig.mean(0).round(2)}   box centre [27.2  2.84  0.17]")
print(f"rec_B.pdb atoms   {len(rec)}")
print(f"rec_B.pdbqt atoms {sum(1 for l in open('structures/rec_B.pdbqt') if l[:6] in ('ATOM  ','HETATM'))}")
