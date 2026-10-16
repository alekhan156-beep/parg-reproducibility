#!/usr/bin/env python3
"""
True RMSD for ALL 9 Vina poses per ligand-receptor pair.

Script 38 corrected only the top-ranked pose. The sampling analysis in
section 9 of LAB_RECORD.md ("near-native pose sampled 7/17", "ranked first
7/7 wet vs 2/7 dry") was computed across all 9 modes using the old
optimal-assignment lower bounds, so it is not yet on the same footing as the
corrected top-1 numbers. This recomputes it.

Same method as script 38: PDB Chemical Component Dictionary ideal structure
matched onto both molecules as a generic graph, minimum over all matches,
no superposition. Hydrogens deleted explicitly (removeHs is a no-op when
sanitize=False).

Run from ~/parg :  python scripts/40_allposes.py
Writes results/crossdock_allposes.csv
"""
import os, csv, glob, re, subprocess, urllib.request
import numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

os.makedirs('data/ccd', exist_ok=True)
os.makedirs('results/tmp_poses', exist_ok=True)
MAXM = 300
THR = 2.0


def template(code):
    f = f'data/ccd/{code}_ideal.sdf'
    if not os.path.exists(f):
        urllib.request.urlretrieve(
            f'https://files.rcsb.org/ligands/download/{code}_ideal.sdf', f)
    return Chem.MolFromMolFile(f, removeHs=True)


def generic_query(t):
    q = Chem.Mol(t)
    try:
        Chem.Kekulize(q, clearAromaticFlags=True)
    except Exception:
        pass
    qp = Chem.AdjustQueryParameters.NoAdjustments()
    qp.makeBondsGeneric = True
    return Chem.AdjustQueryProperties(q, qp)


def load_pdb(path):
    if not path or not os.path.exists(path):
        return None
    m = Chem.MolFromPDBFile(path, removeHs=True, sanitize=False)
    if m is None or m.GetNumAtoms() == 0:
        return None
    hs = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() == 1]
    if hs:
        e = Chem.RWMol(m)
        for i in sorted(hs, reverse=True):
            e.RemoveAtom(i)
        m = e.GetMol()
    if m.GetNumAtoms() == 0:
        return None
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    return m


def rmsd(ref, pose, q):
    mr = ref.GetSubstructMatches(q, uniquify=False, maxMatches=MAXM)
    mp = pose.GetSubstructMatches(q, uniquify=False, maxMatches=MAXM)
    if not mr or not mp:
        return None
    R = ref.GetConformer().GetPositions()
    P = pose.GetConformer().GetPositions()
    Ra = np.array([R[list(x)] for x in mr])
    Pa = np.array([P[list(y)] for y in mp])
    d = ((Pa[:, None, :, :] - Ra[None, :, :, :]) ** 2).sum(-1).mean(-1)
    return float(np.sqrt(d.min()))


def affinities(pdbqt):
    """Vina writes modes in rank order; pull the score of each for the record."""
    out = []
    for line in open(pdbqt):
        m = re.match(r'REMARK VINA RESULT:\s+(-?\d+\.\d+)', line)
        if m:
            out.append(float(m.group(1)))
    return out


tags = sorted({r['ligand'] for r in csv.DictReader(open('results/crossdock.csv'))})
missing  = [f'results/xdock/{t}_{r}.pdbqt' for t in tags for r in ('dry', 'wet')
            if not os.path.exists(f'results/xdock/{t}_{r}.pdbqt')]
missing += [f'structures/xdock/{t}.pdb' for t in tags
            if not os.path.exists(f'structures/xdock/{t}.pdb')]
if missing:
    print(f'PREFLIGHT FAILED - {len(missing)} input file(s) not found:')
    for m in missing[:10]:
        print('   ', m)
    if len(missing) > 10:
        print(f'    ... and {len(missing) - 10} more')
    raise SystemExit(1)
_ob = subprocess.run(['obabel', '-V'], capture_output=True, text=True)
if _ob.returncode != 0:
    print('PREFLIGHT FAILED - obabel not on PATH'); raise SystemExit(1)
print(f'preflight OK: {len(tags)} ligands, {len(tags) * 2} pose files, '
      f'{_ob.stdout.strip()}\n')

rows = []
hdr = f"{'ligand':14} {'rec':4} {'top1':>6} {'best':>6} {'rank':>5} {'n':>3}"
print(hdr); print('-' * len(hdr))

for tag in tags:
    code = tag.split('_')[1]
    tmpl = template(code)
    if tmpl is None:
        print(f'{tag:14}  template unreadable'); continue
    q = generic_query(tmpl)
    ref = load_pdb(f'structures/xdock/{tag}.pdb')
    if ref is None:
        print(f'{tag:14}  reference unreadable'); continue

    for rec in ('dry', 'wet'):
        src = f'results/xdock/{tag}_{rec}.pdbqt'
        if not os.path.exists(src):
            continue
        base = f'results/tmp_poses/{tag}_{rec}'
        for old in glob.glob(f'{base}*.pdb'):
            os.remove(old)
        # -m splits the multi-model PDBQT into one PDB per mode, in rank order
        subprocess.run(['obabel', '-ipdbqt', src, '-opdb', '-O', f'{base}.pdb',
                        '-m'], capture_output=True)
        files = sorted(glob.glob(f'{base}*.pdb'),
                       key=lambda p: int(re.search(r'(\d+)\.pdb$', p).group(1))
                       if re.search(r'(\d+)\.pdb$', p) else 0)
        aff = affinities(src)
        vals = []
        for i, f in enumerate(files, 1):
            pose = load_pdb(f)
            v = rmsd(ref, pose, q) if pose is not None else None
            vals.append(v)
            rows.append({'ligand': tag, 'receptor': rec, 'rank': i,
                         'affinity': aff[i-1] if i-1 < len(aff) else '',
                         'true_rmsd': round(v, 3) if v is not None else ''})
        good = [v for v in vals if v is not None]
        if not good:
            print(f'{tag:14} {rec:4}   no pose computed'); continue
        best = min(good); brank = vals.index(best) + 1
        t1 = vals[0] if vals[0] is not None else float('nan')
        print(f'{tag:14} {rec:4} {t1:6.2f} {best:6.2f} {brank:5d} {len(good):3d}')

with open('results/crossdock_allposes.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['ligand', 'receptor', 'rank',
                                      'affinity', 'true_rmsd'])
    w.writeheader(); w.writerows(rows)

print('\n' + '=' * 64)
print(f"{'receptor':9} {'top-1 <2A':>10} {'any pose <2A':>13} {'near-native ranked 1st':>24}")
for rec in ('dry', 'wet'):
    per = {}
    for r in rows:
        if r['receptor'] == rec and r['true_rmsd'] != '':
            per.setdefault(r['ligand'], {})[r['rank']] = r['true_rmsd']
    t1 = sum(1 for v in per.values() if v.get(1, 99) < THR)
    anyp = sum(1 for v in per.values() if min(v.values()) < THR)
    first = sum(1 for v in per.values()
                if min(v.values()) < THR and v.get(1, 99) < THR)
    print(f'{rec:9} {t1:7d}/{len(per)} {anyp:10d}/{len(per)} {first:19d}/{anyp}')
print('=' * 64)
print('wrote results/crossdock_allposes.csv')
