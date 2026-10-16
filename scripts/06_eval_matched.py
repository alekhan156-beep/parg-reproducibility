import pandas as pd, numpy as np, glob
from sklearn.metrics import mean_squared_error, r2_score

truth = pd.read_csv('data/parg_split.csv')[['smiles','label']].rename(columns={'label':'y'})
files = glob.glob('models/chemprop_matched/**/test_predictions.csv', recursive=True)
print(f"found {len(files)} prediction file(s)")

for f in sorted(files):
    d = pd.read_csv(f).rename(columns={'label':'pred'})
    m = d.merge(truth, on='smiles')
    r = np.sqrt(mean_squared_error(m.y, m.pred))
    print(f"  n={len(m):3d}  RMSE {r:.3f}  R2 {r2_score(m.y, m.pred):.3f}")

print()
print("--- scaffold split, same 260 test compounds ---")
print("Mean predictor      RMSE 0.801   R2 0.000")
print("Tanimoto baseline   RMSE 0.509   R2 0.561")
print("LightGBM            RMSE 0.445   R2 0.664")
