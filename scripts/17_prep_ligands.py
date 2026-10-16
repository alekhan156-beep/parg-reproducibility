import os, subprocess
import numpy as np

XD = 'structures/xdock'
BOX = np.array([27.20, 2.84, 0.17])

def parse(path, pdbqt=False):
    xyz, els = [], []
    for line in open(path):
        if line.startswith(('ATOM', 'HETATM')):
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

tags = sorted(f[:-4] for f in os.listdir(XD) if f.endswith('.pdb'))
print(f"{'ligand':14} {'pdb':>4} {'pdbqt':>6} {'torsions':>9}  status")
print("-" * 52)
ok = []
for tag in tags:
    src = f'{XD}/{tag}.pdb'
    ref, ref_el = parse(src)
    if 'P' in ref_el:
        continue
    if np.linalg.norm(ref.mean(0) - BOX) > 8:
        continue
    dst = f'{XD}/{tag}.pdbqt'
    r = subprocess.run(['obabel', '-ipdb', src, '-opdbqt', '-O', dst, '-h'],
                       capture_output=True, text=True)
    if not os.path.exists(dst):
        print(f"{tag:14} {len(ref):4d} {'-':>6} {'-':>9}  CONVERSION FAILED")
        continue
    pq, _ = parse(dst, pdbqt=True)
    tors = sum(1 for l in open(dst) if l.startswith('BRANCH'))
    status = "ok" if len(pq) == len(ref) else f"ATOM MISMATCH"
    if len(pq) == len(ref):
        ok.append(tag)
    print(f"{tag:14} {len(ref):4d} {len(pq):6d} {tors:9d}  {status}")

print(f"\n{len(ok)} of {len(tags)} ligands ready to dock")
open('results/xdock_targets.txt', 'w').write('\n'.join(ok) + '\n')
