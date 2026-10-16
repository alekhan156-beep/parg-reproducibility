#!/usr/bin/env python3
"""
PARG virtual screen - stage 1 (ligand-based).

Scores a large SMILES library with a pre-trained LightGBM model and writes the
top-ranked compounds for downstream docking.

Memory-safe: the library is processed in chunks, so a 3-million-molecule run
needs ~1 GB of RAM regardless of library size.

Usage
-----
  python screen_parg.py --library zinc_instock.smi --model lgbm_final.pkl

Options
-------
  --top N        how many top compounds to write out      (default 5000)
  --jobs N       parallel worker processes                (default: all cores)
  --chunk N      molecules featurised per batch           (default 50000)
  --fp-only      fingerprints only, skip RDKit descriptors (~50x faster,
                 only valid if the model was trained fingerprint-only)
"""
import argparse, pickle, sys, time
from multiprocessing import Pool

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors

RDLogger.DisableLog('rdApp.*')
DESC = [f for _, f in Descriptors._descList]

FP_ONLY = False


def featurise(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    fp = np.frombuffer(
        AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048).ToBitString().encode(),
        'u1') - ord('0')
    if FP_ONLY:
        return fp.astype(np.float32)
    try:
        d = np.array([f(m) for f in DESC], dtype=np.float64)
    except Exception:
        return None
    d = np.nan_to_num(d, nan=0.0, posinf=0.0, neginf=0.0)
    return np.concatenate([fp, d]).astype(np.float32)


def init(fp_only):
    global FP_ONLY
    FP_ONLY = fp_only


def chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--library', required=True)
    ap.add_argument('--model', required=True)
    ap.add_argument('--out', default='screen_ranked.csv')
    ap.add_argument('--top', type=int, default=5000)
    ap.add_argument('--jobs', type=int, default=0)
    ap.add_argument('--chunk', type=int, default=50000)
    ap.add_argument('--fp-only', action='store_true')
    a = ap.parse_args()

    blob = pickle.load(open(a.model, 'rb'))
    model = blob['model'] if isinstance(blob, dict) else blob
    want = getattr(model, 'n_features_in_', None)
    have = 2048 if a.fp_only else 2048 + len(DESC)
    if want and want != have:
        sys.exit(f"ERROR: model expects {want} features, this run produces {have}.\n"
                 f"       {'Drop --fp-only.' if a.fp_only else 'Add --fp-only.'}")
    print(f"model loaded: {want} features")

    lib = []
    for line in open(a.library):
        p = line.split()
        if len(p) >= 2 and p[0].lower() != 'smiles':
            lib.append((p[0], p[1]))
    print(f"library: {len(lib):,} molecules")

    jobs = a.jobs or None
    t0 = time.time()
    smis, ids, scores = [], [], []
    done = bad = 0

    with Pool(jobs, initializer=init, initargs=(a.fp_only,)) as pool:
        for batch in chunks(lib, a.chunk):
            feats = pool.map(featurise, [s for s, _ in batch], chunksize=200)
            keep = [i for i, f in enumerate(feats) if f is not None]
            bad += len(batch) - len(keep)
            if keep:
                X = np.vstack([feats[i] for i in keep])
                p = model.predict(X)
                for j, i in enumerate(keep):
                    smis.append(batch[i][0]); ids.append(batch[i][1]); scores.append(p[j])
            done += len(batch)
            el = time.time() - t0
            rate = done / el
            print(f"  {done:>9,} / {len(lib):,}   {rate:7.0f} mol/s   "
                  f"eta {(len(lib)-done)/rate/60:6.1f} min", flush=True)

    scores = np.array(scores)
    order = np.argsort(-scores)
    print(f"\nscored {len(scores):,} molecules in {(time.time()-t0)/60:.1f} min "
          f"({bad:,} unparseable)")
    print(f"predicted pChEMBL  max {scores.max():.2f}  "
          f"mean {scores.mean():.2f}  min {scores.min():.2f}")

    with open(a.out, 'w') as f:
        f.write('rank,zinc_id,smiles,predicted_pchembl\n')
        for r, i in enumerate(order[:a.top], 1):
            f.write(f"{r},{ids[i]},{smis[i]},{scores[i]:.4f}\n")
    print(f"wrote top {min(a.top, len(order)):,} -> {a.out}")

    with open(a.out.replace('.csv', '.smi'), 'w') as f:
        for i in order[:a.top]:
            f.write(f"{smis[i]} {ids[i]}\n")
    print(f"wrote {a.out.replace('.csv', '.smi')} for docking")


if __name__ == '__main__':
    main()
