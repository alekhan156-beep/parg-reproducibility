import pandas as pd, numpy as np, glob
from sklearn.metrics import mean_squared_error, r2_score

truth = pd.read_csv('data/parg_split.csv')[['smiles','label']].rename(columns={'label':'y'})
files = sorted(glob.glob('models/chemprop_seed*/**/test_predictions.csv', recursive=True))
print(f"found {len(files)} file(s)")

rmses, r2s = [], []
for f in files:
    d = pd.read_csv(f).rename(columns={'label':'pred'})
    m = d.merge(truth, on='smiles')
    r = np.sqrt(mean_squared_error(m.y, m.pred))
    r2 = r2_score(m.y, m.pred)
    rmses.append(r); r2s.append(r2)
    print(f"  {f.split('/')[1]:18s} n={len(m):3d}  RMSE {r:.3f}  R2 {r2:.3f}")

if rmses:
    print()
    print(f"Chemprop  RMSE {np.mean(rmses):.3f} +/- {np.std(rmses):.3f}   R2 {np.mean(r2s):.3f}")
    print("Baseline  RMSE 0.509")
    print("LightGBM  RMSE 0.445")
