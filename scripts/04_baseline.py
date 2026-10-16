import pandas as pd, numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr

df = pd.read_csv('data/parg_split.csv')
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
df['fp'] = df.smiles.apply(lambda s: gen.GetFingerprint(Chem.MolFromSmiles(s)))

for split in ['split_random', 'split_scaffold']:
    tr = df[df[split] == 'train'].reset_index(drop=True)
    te = df[df[split] == 'test'].reset_index(drop=True)
    tr_fps = list(tr.fp)
    tr_y   = tr.label.values

    preds = []
    for f in te.fp:
        sims = np.array(DataStructs.BulkTanimotoSimilarity(f, tr_fps))
        top = np.argsort(sims)[-5:]
        w = sims[top]
        preds.append(np.average(tr_y[top], weights=w) if w.sum() > 0 else tr_y.mean())
    preds = np.array(preds)

    rmse = np.sqrt(mean_squared_error(te.label, preds))
    r, _ = pearsonr(te.label, preds)
    dummy = np.sqrt(mean_squared_error(te.label, np.full(len(te), tr_y.mean())))
    print(f"{split:16s} RMSE {rmse:.3f}  R2 {r2_score(te.label, preds):6.3f}  r {r:.3f}   (mean-only RMSE {dummy:.3f})")
