import pickle, os
import numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from lightgbm import LGBMRegressor

df = pd.read_csv('data/parg_curated.csv')
print("columns:", list(df.columns))
sc = 'smiles' if 'smiles' in df.columns else [c for c in df.columns if 'smile' in c.lower()][0]
yc = 'label'  if 'label'  in df.columns else 'pchembl'
print(f"using smiles='{sc}'  target='{yc}'  n={len(df)}")

def feats(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    fp = np.array(AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048))
    d = np.array([f(m) for _, f in Descriptors._descList])
    return np.nan_to_num(np.concatenate([fp, d]), nan=0, posinf=0, neginf=0)

X, y = [], []
for s, p in zip(df[sc], df[yc]):
    f = feats(s)
    if f is not None:
        X.append(f); y.append(p)
X, y = np.array(X), np.array(y)
print(f"features {X.shape}")

m = LGBMRegressor(n_estimators=500, learning_rate=0.05, num_leaves=31,
                  random_state=0, verbose=-1).fit(X, y)
os.makedirs('models', exist_ok=True)
pickle.dump({'model': m, 'n_desc': len(Descriptors._descList)},
            open('models/lgbm_final.pkl', 'wb'))
print(f"saved models/lgbm_final.pkl  (trained on {len(y)} compounds)")
