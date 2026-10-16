#!/usr/bin/env python3
"""
Figure 6 — bridging waters act on ranking, not on sampling.

Panel A: how often a near-native pose (<2 A) appears anywhere among the 9 modes.
Panel B: where that pose is ranked, for the ligands sampled in BOTH receptors.

Run from ~/parg :  python scripts/42_fig6.py
Reads  results/crossdock_allposes.csv  (written by 40_allposes.py).
Falls back to the values printed by that script if the CSV isn't present.
Writes figures/fig6_sampling_vs_ranking.png / .pdf
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BLUE, VERM, GREY, INK = '#0072B2', '#D55E00', '#9a9a9a', '#1a1a1a'
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9,
    'axes.edgecolor': '#cccccc', 'axes.linewidth': 0.8,
    'axes.labelcolor': INK, 'text.color': INK,
    'xtick.color': '#555555', 'ytick.color': '#555555',
    'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.facecolor': 'white', 'axes.facecolor': 'white'})
os.makedirs('figures', exist_ok=True)
THR, NLIG = 2.0, 17

SRC = 'results/crossdock_allposes.csv'
if os.path.exists(SRC):
    per = {}
    for r in csv.DictReader(open(SRC)):
        if r['true_rmsd'] == '':
            continue
        per.setdefault(r['receptor'], {}).setdefault(r['ligand'], {})[int(r['rank'])] = float(r['true_rmsd'])
    best = {rec: {lg: min(v.values()) for lg, v in d.items()} for rec, d in per.items()}
    rank = {rec: {lg: min(v, key=v.get) for lg, v in d.items()} for rec, d in per.items()}
    samp = {rec: {lg for lg, b in best[rec].items() if b < THR} for rec in ('dry', 'wet')}
    print(f'read {SRC}')
else:                                   # values as printed by 40_allposes.py
    samp = {'dry': {'6O9X_M0S','6OA1_M0P','6OA3_M0M','6OAL_M0P','7KFP_WD7','7KG8_WDM'},
            'wet': {'6O9X_M0S','6OA1_M0P','6OA3_M0M','6OAL_M0P','7KFP_WD7','7KG7_WDJ','7KG8_WDM'}}
    rank = {'dry': {'6O9X_M0S':6,'6OA1_M0P':5,'6OA3_M0M':5,'6OAL_M0P':6,'7KFP_WD7':4,'7KG8_WDM':3},
            'wet': {'6O9X_M0S':1,'6OA1_M0P':1,'6OA3_M0M':2,'6OAL_M0P':1,'7KFP_WD7':1,'7KG7_WDJ':1,'7KG8_WDM':1}}
    print('CSV not found - using the values printed by 40_allposes.py')

common = sorted(samp['dry'] & samp['wet'])
dr = np.array([rank['dry'][k] for k in common], float)
wr = np.array([rank['wet'][k] for k in common], float)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.6, 4.3),
                               gridspec_kw={'width_ratios': [1.0, 1.45]})

# ---- A: sampling ---------------------------------------------------------
nd, nw = len(samp['dry']), len(samp['wet'])
axA.barh([1, 0], [nd, nw], height=0.46, color=[VERM, BLUE],
         edgecolor='white', linewidth=0.8, zorder=3)
for y, n in ((1, nd), (0, nw)):
    axA.text(n + 0.3, y, f'{n}/{NLIG}', va='center', ha='left',
             fontsize=10.5, fontweight='bold', color=INK, zorder=4)
axA.set_yticks([1, 0]); axA.set_yticklabels(['dry', 'wet'], fontsize=9.5)
axA.set_xlim(0, NLIG + 1.5); axA.set_xticks([0, 4, 8, 12, 16])
axA.set_xlabel('ligands with a near-native pose\namong the 9 modes', fontsize=8.8)
axA.set_title('A   Sampling — unchanged', loc='left',
              fontsize=10, fontweight='bold', pad=12)
axA.set_ylim(-0.95, 1.6)
axA.text(0, -0.78, 'McNemar  p = 1.0   (no detectable difference)',
         fontsize=8.4, color='#555555', ha='left')

# ---- B: ranking ----------------------------------------------------------
def spread(vals):
    """vertical offsets so points sharing a rank are individually visible"""
    out, seen = [], {}
    counts = {}
    for v in vals:
        counts[v] = counts.get(v, 0) + 1
    for v in vals:
        i = seen.get(v, 0); seen[v] = i + 1
        out.append(v + (i - (counts[v] - 1) / 2) * 0.30)
    return np.array(out)

x0, x1 = 0.0, 1.0
dy, wy = spread(list(dr)), spread(list(wr))
axB.axhspan(0.45, 1.55, color='#f2f7fa', zorder=0)
for i in range(len(common)):
    axB.plot([x0, x1], [dy[i], wy[i]], color='#d9d9d9', lw=1.5,
             zorder=2, solid_capstyle='round')
axB.scatter([x0]*len(dy), dy, s=46, marker='o', color=VERM, zorder=4,
            edgecolor='white', linewidth=0.8)
axB.scatter([x1]*len(wy), wy, s=50, marker='s', color=BLUE, zorder=4,
            edgecolor='white', linewidth=0.8)
for i, k in enumerate(common):
    axB.text(x1 + 0.08, wy[i], k.split('_')[0], va='center', ha='left',
             fontsize=8.0, color='#444444')
axB.text(x0 - 0.10, 1.0, 'ranked\nfirst', va='center', ha='right',
         fontsize=8.0, color=BLUE, fontweight='bold', linespacing=1.3)
axB.set_xlim(-0.52, 1.72); axB.set_ylim(7.2, 0.15)
axB.set_xticks([x0, x1]); axB.set_xticklabels(['dry', 'wet'], fontsize=9.5)
axB.set_yticks(range(1, 8))
axB.set_ylabel('rank of the near-native pose among the 9 modes', fontsize=8.8)
axB.set_title('B   Ranking — transformed', loc='left',
              fontsize=10, fontweight='bold', pad=12)
axB.set_xlabel(f'mean rank {dr.mean():.2f} → {wr.mean():.2f}     '
               f'better in {int((wr<dr).sum())}/{len(common)}     '
               f'Wilcoxon p = 0.031', fontsize=8.4, color='#555555', labelpad=10)
axB.spines['bottom'].set_visible(False)
axB.tick_params(axis='x', length=0, pad=6)

fig.suptitle('Retaining bridging waters does not change which poses are found — '
             'only which one is chosen',
             fontsize=10.5, fontweight='bold', y=0.985)
fig.tight_layout(rect=[0, 0.02, 1, 0.95])
for ext in ('png', 'pdf'):
    fig.savefig(f'figures/fig6_sampling_vs_ranking.{ext}', dpi=300,
                bbox_inches='tight', facecolor='white')
print(f'sampled  dry {nd}/{NLIG}  wet {nw}/{NLIG}')
print(f'ranking  n={len(common)}  mean {dr.mean():.2f} -> {wr.mean():.2f}')
print('wrote figures/fig6_sampling_vs_ranking.png / .pdf')
