import pandas as pd, numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, Descriptors
from sklearn.metrics import mean_squared_error, r2_score
import lightgbm as lgb
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

df = pd.read_csv('data/parg_split.csv')
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def featurise(smi):
    m = Chem.MolFromSmiles(smi)
    fp = np.array(gen.GetFingerprint(m))
    desc = np.array(list(Descriptors.CalcMolDescriptors(m).values()))
    return np.concatenate([fp, desc])

X = np.vstack([featurise(s) for s in df.smiles])
X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
y = df.label.values
print(f"features: {X.shape}")

for split in ['split_random', 'split_scaffold']:
    tr = (df[split] == 'train').values
    te = ~tr
    rmses, r2s = [], []
    for seed in range(5):
        m = lgb.LGBMRegressor(n_estimators=500, learning_rate=0.05,
                              num_leaves=31, random_state=seed, verbose=-1)
        m.fit(X[tr], y[tr])
        p = m.predict(X[te])
        rmses.append(np.sqrt(mean_squared_error(y[te], p)))
        r2s.append(r2_score(y[te], p))
    print(f"{split:16s} RMSE {np.mean(rmses):.3f} +/- {np.std(rmses):.3f}"
          f"   R2 {np.mean(r2s):.3f}")
