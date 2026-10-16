#!/usr/bin/env python3
"""
PARG virtual screen - stage 2 (structure-based).

Docks the top-ranked compounds from stage 1 into the validated wet receptor.

Usage
-----
  python dock_hits.py --hits screen_ranked.smi --receptor rec_wet.pdbqt \
                      --vina ./vina_1.2.5_linux_x86_64 --n 500

Box is 16 A centred on the 6OAK M0M site. Do not enlarge it: a 24 A box spans
the whole poly(ADP-ribose) substrate groove and lets Vina park small ligands in
the wrong subpocket.
"""
import argparse, csv, os, subprocess, sys, time

CENTER = ('27.20', '2.84', '0.17')
SIZE = ('16', '16', '16')


def affinity(path):
    for line in open(path):
        if line.startswith('REMARK VINA RESULT'):
            return float(line.split()[3])
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hits', required=True, help='.smi from stage 1')
    ap.add_argument('--receptor', required=True, help='rec_wet.pdbqt')
    ap.add_argument('--vina', default='./vina_1.2.5_linux_x86_64')
    ap.add_argument('--n', type=int, default=500)
    ap.add_argument('--exhaustiveness', type=int, default=16)
    ap.add_argument('--outdir', default='docked')
    ap.add_argument('--out', default='docking_results.csv')
    a = ap.parse_args()

    if not os.path.exists(a.receptor):
        sys.exit(f"receptor not found: {a.receptor}")
    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs('ligands', exist_ok=True)

    hits = []
    for line in open(a.hits):
        p = line.split()
        if len(p) >= 2:
            hits.append((p[0], p[1]))
    hits = hits[:a.n]
    print(f"docking {len(hits)} compounds, exhaustiveness {a.exhaustiveness}\n")

    rows, t0 = [], time.time()
    for k, (smi, zid) in enumerate(hits, 1):
        lig = f"ligands/{zid}.pdbqt"
        if not os.path.exists(lig):
            r = subprocess.run(['obabel', f'-:{smi}', '-opdbqt', '-O', lig,
                                '--gen3d', '-h'], capture_output=True, text=True)
            if not os.path.exists(lig) or os.path.getsize(lig) == 0:
                print(f"  [{k}] {zid}  3D generation failed, skipped")
                continue

        out = f"{a.outdir}/{zid}.pdbqt"
        if not os.path.exists(out):
            subprocess.run([a.vina, '--receptor', a.receptor, '--ligand', lig,
                            '--center_x', CENTER[0], '--center_y', CENTER[1],
                            '--center_z', CENTER[2],
                            '--size_x', SIZE[0], '--size_y', SIZE[1],
                            '--size_z', SIZE[2],
                            '--exhaustiveness', str(a.exhaustiveness),
                            '--num_modes', '9', '--seed', '42',
                            '--out', out], capture_output=True)
        aff = affinity(out) if os.path.exists(out) else None
        if aff is None:
            print(f"  [{k}] {zid}  docking failed, skipped")
            continue
        rows.append({'zinc_id': zid, 'smiles': smi, 'vina_kcal': aff})

        el = (time.time() - t0) / 60
        print(f"  [{k}/{len(hits)}] {zid}  {aff:7.2f} kcal/mol   "
              f"{el:5.1f} min elapsed, eta {el/k*(len(hits)-k):5.1f} min", flush=True)

        with open(a.out, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=['zinc_id', 'smiles', 'vina_kcal'])
            w.writeheader()
            w.writerows(sorted(rows, key=lambda r: r['vina_kcal']))

    print(f"\ndocked {len(rows)} compounds -> {a.out}")
    if rows:
        best = sorted(rows, key=lambda r: r['vina_kcal'])[:10]
        print("\ntop 10 by Vina score:")
        for r in best:
            print(f"  {r['zinc_id']:20} {r['vina_kcal']:7.2f}")
        print("\nNOTE: Vina score is NOT a reliable proxy for pose correctness.")
        print("      Inspect poses before trusting this ranking.")


if __name__ == '__main__':
    main()
