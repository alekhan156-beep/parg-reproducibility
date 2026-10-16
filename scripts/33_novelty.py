import numpy as np, pandas as pd
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import AllChem
RDLogger.DisableLog('rdApp.*')

def fp(s):
    m = Chem.MolFromSmiles(s)
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) if m else None

tr = pd.read_csv('data/parg_curated.csv')
sc = 'smiles' if 'smiles' in tr.columns else [c for c in tr.columns if 'smile' in c.lower()][0]
known = [f for f in (fp(s) for s in tr[sc]) if f is not None]
print(f"{len(known)} known PARG compounds")

hits = pd.read_csv('results/screen_ranked.csv')
print(f"{len(hits)} screened hits\n")

sims = []
for s in hits['smiles']:
    f = fp(s)
    sims.append(max(DataStructs.BulkTanimotoSimilarity(f, known)) if f else np.nan)
hits['max_sim_to_known'] = sims

for lo, hi, lab in ((0.0,0.3,'very novel  '), (0.3,0.4,'novel       '),
                    (0.4,0.6,'related     '), (0.6,1.01,'close analogue')):
    n = ((hits['max_sim_to_known'] >= lo) & (hits['max_sim_to_known'] < hi)).sum()
    print(f"  {lab} (Tanimoto {lo:.1f}-{hi:.1f}): {n:5d}  ({n/len(hits)*100:4.1f}%)")

hits.to_csv('results/screen_ranked_novelty.csv', index=False)
nov = hits[hits['max_sim_to_known'] < 0.4].head(20)
print(f"\ntop 20 novel hits (predicted pChEMBL, similarity to nearest known):")
for _, r in nov.iterrows():
    print(f"  {r['zinc_id']:18} {r['predicted_pchembl']:.2f}   sim {r['max_sim_to_known']:.2f}")
print(f"\nwrote results/screen_ranked_novelty.csv")
