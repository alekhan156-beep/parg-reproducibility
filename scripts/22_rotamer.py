import csv, numpy as np

def dih(p):
    b0, b1, b2 = p[0]-p[1], p[2]-p[1], p[3]-p[2]
    b1 = b1/np.linalg.norm(b1)
    v = b0 - (b0@b1)*b1; w = b2 - (b2@b1)*b1
    return np.degrees(np.arctan2(np.cross(b1,v)@w, v@w))

def chi1b(path, n):
    a = {}
    for l in open(path):
        if l[:4]=='ATOM' and l[21]=='A' and l[22:26].strip()==str(n) and l[16] in ' A':
            a[l[12:16].strip()] = (np.array([float(l[30:38]), float(l[38:46]),
                                   float(l[46:54])]), float(l[60:66]))
    k = ('N','CA','CB','CG')
    if not all(x in a for x in k): return None, None
    b = np.mean([v[1] for x,v in a.items() if x not in ('N','CA','C','O')])
    return dih([a[x][0] for x in k]), b

st = {}
for r in csv.DictReader(open('results/crossdock.csv')):
    if r['receptor']=='wet':
        st[r['ligand'].split('_')[0]] = 'Y' if float(r['best_rmsd'])<2 else 'n'
st['6OAK'] = 'ref'

print("PDB     grp   Y795chi1  Y795B   F902chi1   F902B")
g = {'Y':[], 'n':[]}
for p in sorted(st):
    c1, b1 = chi1b(f'structures/pdb/{p}.pdb', 795)
    c2, b2 = chi1b(f'structures/pdb/{p}.pdb', 902)
    if c1 is None or c2 is None:
        print(f"{p:7} {st[p]:4}  incomplete"); continue
    print(f"{p:7} {st[p]:4} {c1:9.1f} {b1:6.1f} {c2:10.1f} {b2:7.1f}")
    if st[p] in g: g[st[p]].append((c1, c2))

for k in ('Y','n'):
    a = np.array(g[k])
    print(f"{k}: Y795 {a[:,0].mean():7.1f} +/- {a[:,0].std():5.1f}    "
          f"F902 {a[:,1].mean():7.1f} +/- {a[:,1].std():5.1f}")
