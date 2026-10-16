#!/usr/bin/env python3
"""
Cross-docking RMSD, symmetry-corrected, no alignment.

Why this replaces scripts 13/14 and 37:
  * script 13/14 used element-matched optimal assignment, which ignores
    connectivity. It is a LOWER BOUND on the true RMSD, not the RMSD.
  * script 37 used rdMolAlign.CalcRMS, which is correct but requires
    perceiving bond orders in both molecules. It failed for 12 of 17 ligands.
    Diagnosed cause: RDKit's MolFromPDBFile(removeHs=True) does not remove
    hydrogens when sanitize=False. The 12 failures are exactly the ligands
    deposited with explicit hydrogens (e.g. 5LHB/6WJ: 55 PDB records, 31
    heavy atoms), so the reference entered as a 55-atom molecule and could
    not be matched to a 31-atom docked pose. The 5 that succeeded are
    exactly those deposited without hydrogens. Fixed here by deleting
    atomic-number-1 atoms explicitly after loading.
  * RMSD needs only coordinates, elements and connectivity. This script
    matches the PDB Chemical Component Dictionary ideal structure onto both
    molecules as a GENERIC graph (bond orders ignored, aromaticity cleared),
    then takes the minimum over all matches. Enumerating all matches is what
    makes it symmetry-corrected; it also resolves multiple crystallographic
    copies of the ligand by selecting the nearest one.
  * Molecules are NOT superposed. GetBestRMS would superpose them, which
    measures shape similarity rather than placement accuracy and would
    spuriously improve every number.

Validated on three synthetic cases: identity -> 0.0000, 1.0 A rigid
translation -> 1.0000, duplicated reference copy -> 1.0000.

Run from ~/parg :  python scripts/38_true_rmsd.py
Writes results/crossdock_truermsd.csv
"""
import os, csv, subprocess, urllib.request
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS
RDLogger.DisableLog('rdApp.*')

os.makedirs('data/ccd', exist_ok=True)
os.makedirs('results/tmp_pdb', exist_ok=True)
MAXM = 300          # cap on substructure matches enumerated per molecule


def template(code):
    f = f'data/ccd/{code}_ideal.sdf'
    if not os.path.exists(f):
        urllib.request.urlretrieve(
            f'https://files.rcsb.org/ligands/download/{code}_ideal.sdf', f)
    return Chem.MolFromMolFile(f, removeHs=True)


def generic_query(t):
    """Template as a connectivity-only query: elements must match, bond
    orders and aromatic flags must not."""
    q = Chem.Mol(t)
    try:
        Chem.Kekulize(q, clearAromaticFlags=True)
    except Exception:
        pass
    qp = Chem.AdjustQueryParameters.NoAdjustments()
    qp.makeBondsGeneric = True
    return Chem.AdjustQueryProperties(q, qp)


def load_pdb(path):
    """Heavy atoms + coordinates, no sanitization (bond orders are irrelevant
    here, and sanitization is what was failing). Hydrogens are deleted
    explicitly because removeHs=True is a no-op when sanitize=False."""
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


def pose_to_pdb(pdbqt, out):
    subprocess.run(['obabel', '-ipdbqt', pdbqt, '-opdb', '-O', out,
                    '-f', '1', '-l', '1'], capture_output=True)
    return out if os.path.exists(out) and os.path.getsize(out) > 0 else None


def best_rmsd(ref, pose, q):
    """min over all template matches in ref x all template matches in pose.
    Returns (rmsd, n_atoms_used, how)."""
    mr = ref.GetSubstructMatches(q, uniquify=False, maxMatches=MAXM)
    mp = pose.GetSubstructMatches(q, uniquify=False, maxMatches=MAXM)
    how = 'template'
    if not mr or not mp:
        # Fall back to the largest substructure common to ref and pose.
        # Used when the deposited ligand is incompletely modelled.
        res = rdFMCS.FindMCS([ref, pose], timeout=20,
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareAny,
                             ringMatchesRingOnly=False, completeRingsOnly=False)
        if res.canceled or res.numAtoms < 4:
            return (None, 0, f'no match (ref {ref.GetNumAtoms()}at '
                    f'pose {pose.GetNumAtoms()}at query {q.GetNumAtoms()}at)')
        qq = Chem.MolFromSmarts(res.smartsString)
        mr = ref.GetSubstructMatches(qq, uniquify=False, maxMatches=MAXM)
        mp = pose.GetSubstructMatches(qq, uniquify=False, maxMatches=MAXM)
        how = f'MCS {res.numAtoms}at'
        if not mr or not mp:
            return None, 0, 'no match even by MCS'
    R = ref.GetConformer().GetPositions()
    P = pose.GetConformer().GetPositions()
    Ra = np.array([R[list(x)] for x in mr])            # (nr, k, 3)
    Pa = np.array([P[list(y)] for y in mp])            # (np, k, 3)
    d = ((Pa[:, None, :, :] - Ra[None, :, :, :]) ** 2).sum(-1).mean(-1)
    return float(np.sqrt(d.min())), Ra.shape[1], how


rows = list(csv.DictReader(open('results/crossdock.csv')))
tags = sorted({r['ligand'] for r in rows})
old = {(r['ligand'], r['receptor']): float(r['top1_rmsd']) for r in rows}

hdr = f"{'ligand':14} {'rec':4} {'lower bd':>9} {'true':>8} {'shift':>7}  {'atoms':>5} {'method':<12}"
print(hdr); print('-' * len(hdr))
out_rows, problems = [], []

for tag in tags:
    code = tag.split('_')[1]
    try:
        tmpl = template(code)
    except Exception as e:
        problems.append(f'{tag}: template download failed ({e})'); continue
    if tmpl is None:
        problems.append(f'{tag}: template unreadable'); continue
    q = generic_query(tmpl)
    ref = load_pdb(f'structures/xdock/{tag}.pdb')
    if ref is None:
        problems.append(f'{tag}: reference PDB unreadable'); continue

    for rec in ('dry', 'wet'):
        src = f'results/xdock/{tag}_{rec}.pdbqt'
        if not os.path.exists(src):
            continue
        pose = load_pdb(pose_to_pdb(src, f'results/tmp_pdb/{tag}_{rec}.pdb'))
        if pose is None:
            problems.append(f'{tag} {rec}: pose unreadable'); continue
        rms, nat, how = best_rmsd(ref, pose, q)
        if rms is None:
            problems.append(f'{tag} {rec}: {how}'); continue
        lb = old.get((tag, rec), float('nan'))
        print(f'{tag:14} {rec:4} {lb:9.2f} {rms:8.2f} {rms-lb:+7.2f}  {nat:5d} {how:<12}')
        out_rows.append({'ligand': tag, 'receptor': rec,
                         'lower_bound_rmsd': round(lb, 3),
                         'true_rmsd': round(rms, 3),
                         'n_atoms': nat, 'method': how,
                         'tmpl_atoms': tmpl.GetNumAtoms(),
                         'ref_atoms': ref.GetNumAtoms(),
                         'pose_atoms': pose.GetNumAtoms()})

with open('results/crossdock_truermsd.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['ligand', 'receptor', 'lower_bound_rmsd',
                                      'true_rmsd', 'n_atoms', 'method',
                                      'tmpl_atoms', 'ref_atoms', 'pose_atoms'])
    w.writeheader(); w.writerows(out_rows)

print('\n' + '=' * len(hdr))
for rec in ('dry', 'wet'):
    s = [r for r in out_rows if r['receptor'] == rec]
    if not s:
        continue
    lb = sum(1 for r in s if r['lower_bound_rmsd'] < 2)
    tr = sum(1 for r in s if r['true_rmsd'] < 2)
    print(f'{rec:4}  under 2 A:   lower bound {lb}/{len(s)}      TRUE {tr}/{len(s)}')
print(f'\ncomputed {len(out_rows)} of {2*len(tags)} ligand-receptor pairs')
if problems:
    print('\nnot computed:')
    for p in problems:
        print('  ' + p)
print('wrote results/crossdock_truermsd.csv')
print('=' * len(hdr))
