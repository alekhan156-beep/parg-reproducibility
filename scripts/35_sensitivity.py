
import random

import numpy as np, pandas as pd

from rdkit import Chem, RDLogger

from rdkit.Chem import Descriptors, rdMolDescriptors

from scipy.stats import mannwhitneyu, fisher_exact

RDLogger.DisableLog('rdApp.*')

N = 100



top = []

for line in open('results/screen_ranked.smi'):

    p = line.split()

    if len(p) >= 2: top.append((p[0], p[1]))

top = top[:N]

topids = {z for _, z in top}



lib = []

for line in open('data/zinc_instock.smi'):

    p = line.split()

    if len(p) >= 2 and p[1] not in topids: lib.append((p[0], p[1]))

random.seed(42)

rnd = random.sample(lib, N)

del lib



df = pd.read_csv('results/control.csv')

done = set(df['zinc_id'].astype(str))

miss_t = [(s, z) for s, z in top if z not in done]

miss_r = [(s, z) for s, z in rnd if z not in done]

print(f"attempted: {N} top, {N} random")

print(f"docked   : {(df['group']=='top').sum()} top, {(df['group']=='random').sum()} random")

print(f"FAILED   : {len(miss_t)} top, {len(miss_r)} random\n")




def props(s):

    m = Chem.MolFromSmiles(s)

    if m is None: return None

    return (Descriptors.MolWt(m), m.GetNumHeavyAtoms(),

            rdMolDescriptors.CalcNumRotatableBonds(m),

            rdMolDescriptors.CalcNumRings(m), Descriptors.MolLogP(m))



print("failed compounds:")

print(f"  {'zinc_id':16} {'MW':>7} {'heavy':>6} {'rotB':>5} {'rings':>6} {'logP':>6}")

fp_ = []

for s, z in miss_t + miss_r:

    p = props(s)

    if p: fp_.append(p); print(f"  {z:16} {p[0]:7.1f} {p[1]:6d} {p[2]:5d} {p[3]:6d} {p[4]:6.2f}")



ok = [props(s) for s in df['smiles']]

ok = np.array([o for o in ok if o])

if fp_:

    f = np.array(fp_)

    print(f"\n  {'':16} {'MW':>7} {'heavy':>6} {'rotB':>5} {'rings':>6} {'logP':>6}")

    print(f"  {'failed  mean':16} {f[:,0].mean():7.1f} {f[:,1].mean():6.1f} "

          f"{f[:,2].mean():5.1f} {f[:,3].mean():6.1f} {f[:,4].mean():6.2f}")

    print(f"  {'docked  mean':16} {ok[:,0].mean():7.1f} {ok[:,1].mean():6.1f} "

          f"{ok[:,2].mean():5.1f} {ok[:,3].mean():6.1f} {ok[:,4].mean():6.2f}")




t = df[df['group']=='top']['vina_kcal'].values

r = df[df['group']=='random']['vina_kcal'].values

worst = max(t.max(), r.max())



print("\n" + "="*60)

print(f"{'scenario':34} {'diff':>7} {'p':>9}  verdict")

print("-"*60)



def rep(lab, a, b):

    d = b.mean() - a.mean()

    p = mannwhitneyu(a, b, alternative='less').pvalue

    print(f"{lab:34} {d:+7.3f} {p:9.3g}  "

          f"{'holds' if p < 0.05 else 'FAILS'}")



rep("as reported (92 v 100)", t, r)

t_worst = np.concatenate([t, np.full(len(miss_t), worst)])

rep(f"failures imputed at {worst:.2f} (worst)", t_worst, r)

t_best = np.concatenate([t, np.full(len(miss_t), t.min())])

rep(f"failures imputed at {t.min():.2f} (best)", t_best, r)



print("\nenrichment below -9 kcal/mol:")

a, b = (t < -9).sum(), (r < -9).sum()

odds, pf = fisher_exact([[a, len(t)-a], [b, len(r)-b]], alternative='greater')

print(f"  top {a}/{len(t)} = {a/len(t)*100:.1f}%   "

      f"random {b}/{len(r)} = {b/len(r)*100:.1f}%")

print(f"  {a/len(t)/(b/len(r)):.2f}-fold enrichment, Fisher exact p = {pf:.3g}")

print("="*60)

