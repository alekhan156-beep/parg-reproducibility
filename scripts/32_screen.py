import pickle, time
import numpy as np
from multiprocessing import Pool
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors
RDLogger.DisableLog('rdApp.*')
DESC = [f for _, f in Descriptors._descList]

def feat(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    fp = np.array(AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048))
    try:
        d = np.array([f(m) for f in DESC], dtype=np.float64)
    except Exception:
        return None
    d = np.nan_to_num(d, nan=0, posinf=0, neginf=0)
    return np.concatenate([fp, d]).astype(np.float32)

if __name__ == '__main__':
    blob = pickle.load(open('models/lgbm_final.pkl', 'rb'))
    model = blob['model'] if isinstance(blob, dict) else blob
    lib = []
    for line in open('data/zinc_instock.smi'):
        p = line.split()
        if len(p) >= 2: lib.append((p[0], p[1]))
    print(f"{len(lib):,} molecules, model wants {model.n_features_in_} features")

    smis, ids, sc = [], [], []
    t0 = time.time()
    with Pool(8) as pool:
        for i in range(0, len(lib), 50000):
            b = lib[i:i+50000]
            fs = pool.map(feat, [s for s, _ in b], chunksize=200)
            k = [j for j, f in enumerate(fs) if f is not None]
            if k:
                pr = model.predict(np.vstack([fs[j] for j in k]))
                for n, j in enumerate(k):
                    smis.append(b[j][0]); ids.append(b[j][1]); sc.append(pr[n])
            d = i + len(b); el = time.time() - t0
            print(f"  {d:>9,}/{len(lib):,}  {d/el:6.0f} mol/s  "
                  f"eta {(len(lib)-d)/(d/el)/60:5.1f} min", flush=True)

    sc = np.array(sc); o = np.argsort(-sc)
    print(f"\nscored {len(sc):,} in {(time.time()-t0)/60:.1f} min    "
          f"max {sc.max():.2f}  mean {sc.mean():.2f}  min {sc.min():.2f}")
    with open('results/screen_ranked.csv', 'w') as f:
        f.write('rank,zinc_id,smiles,predicted_pchembl\n')
        for r, j in enumerate(o[:5000], 1):
            f.write(f"{r},{ids[j]},{smis[j]},{sc[j]:.4f}\n")
    with open('results/screen_ranked.smi', 'w') as f:
        for j in o[:5000]: f.write(f"{smis[j]} {ids[j]}\n")
    print("wrote results/screen_ranked.csv and results/screen_ranked.smi")
