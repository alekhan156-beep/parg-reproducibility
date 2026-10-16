#!/usr/bin/env python3
"""
Recompute cross-docking RMSD with proper symmetry-corrected values.

The values in results/crossdock.csv were computed by element-matched optimal
assignment, which ignores connectivity and is therefore a LOWER BOUND on the
true RMSD. This replaces them with RDKit rdMolAlign.CalcRMS, which is
symmetry-aware and does NOT align the molecules (GetBestRMS would align them,
which inflates apparent accuracy and is wrong for pose evaluation).

Bond orders come from each ligand's ideal SDF in the PDB Chemical Component
Dictionary, applied to both poses via AssignBondOrdersFromTemplate.

Run from ~/parg :  python scripts/37_true_rmsd.py
Writes results/crossdock_truermsd.csv
"""
import os, csv, subprocess, urllib.request
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign
RDLogger.DisableLog('rdApp.*')

os.makedirs('data/ccd', exist_ok=True)
os.makedirs('results/tmp_pdb', exist_ok=True)

def template(code):
    """Ideal-geometry reference molecule with correct bond orders."""
    f = f'data/ccd/{code}_ideal.sdf'
    if not os.path.exists(f):
        urllib.request.urlretrieve(
            f'https://files.rcsb.org/ligands/download/{code}_ideal.sdf', f)
    m = Chem.MolFromMolFile(f, removeHs=True)
    return m

def with_bonds(path, tmpl):
    m = Chem.MolFromPDBFile(path, removeHs=True, sanitize=False)
    if m is None:
        return None
    try:
        return AllChem.AssignBondOrdersFromTemplate(tmpl, m)
    except Exception:
        return None

def pose_to_pdb(pdbqt, out):
    subprocess.run(['obabel', '-ipdbqt', pdbqt, '-opdb', '-O', out,
                    '-f', '1', '-l', '1'], capture_output=True)
    return out if os.path.exists(out) and os.path.getsize(out) > 0 else None

rows = list(csv.DictReader(open('results/crossdock.csv')))
tags = sorted({r['ligand'] for r in rows})
old = {(r['ligand'], r['receptor']): float(r['top1_rmsd']) for r in rows}

print(f"{'ligand':14} {'rec':4} {'lower bound':>12} {'true RMSD':>10} {'shift':>7}")
print("-" * 54)
out_rows, fails = [], []
for tag in tags:
    code = tag.split('_')[1]
    try:
        tmpl = template(code)
    except Exception as e:
        print(f"{tag:14}  template download failed ({e})"); fails.append(tag); continue
    if tmpl is None:
        print(f"{tag:14}  template unreadable"); fails.append(tag); continue

    ref = with_bonds(f'structures/xdock/{tag}.pdb', tmpl)
    if ref is None:
        print(f"{tag:14}  reference bond assignment failed"); fails.append(tag); continue

    for rec in ('dry', 'wet'):
        src = f'results/xdock/{tag}_{rec}.pdbqt'
        if not os.path.exists(src):
            continue
        p = pose_to_pdb(src, f'results/tmp_pdb/{tag}_{rec}.pdb')
        pose = with_bonds(p, tmpl) if p else None
        if pose is None:
            print(f"{tag:14} {rec:4}  pose bond assignment failed"); continue
        try:
            rms = rdMolAlign.CalcRMS(pose, ref)
        except Exception as e:
            print(f"{tag:14} {rec:4}  CalcRMS failed ({e})"); continue
        lb = old.get((tag, rec), float('nan'))
        print(f"{tag:14} {rec:4} {lb:12.2f} {rms:10.2f} {rms-lb:+7.2f}")
        out_rows.append({'ligand': tag, 'receptor': rec,
                         'lower_bound_rmsd': round(lb, 3),
                         'true_rmsd': round(rms, 3)})

with open('results/crossdock_truermsd.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['ligand', 'receptor',
                                      'lower_bound_rmsd', 'true_rmsd'])
    w.writeheader(); w.writerows(out_rows)

print("\n" + "=" * 54)
for rec in ('dry', 'wet'):
    s = [r for r in out_rows if r['receptor'] == rec]
    lb = sum(1 for r in s if r['lower_bound_rmsd'] < 2)
    tr = sum(1 for r in s if r['true_rmsd'] < 2)
    print(f"{rec:4}  under 2 A:  lower bound {lb}/{len(s)}   true {tr}/{len(s)}")
if fails:
    print(f"\nfailed entirely: {', '.join(fails)}")
print(f"wrote results/crossdock_truermsd.csv")
print("=" * 54)
