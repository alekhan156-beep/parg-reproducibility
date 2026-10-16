#!/usr/bin/env python3
"""
Regenerate all paper figures from the project's own result files.

Run from ~/parg :   python scripts/36_figures.py
Writes PNG (300 dpi) and PDF into figures/.

Reads:  results/crossdock.csv, results/control.csv, structures/pdb/*.pdb
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

def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(f'figures/{name}.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  wrote figures/{name}.png / .pdf")

# ---------------------------------------------------------------- figure 1
models = ['Mean predictor', 'Tanimoto 5-NN\n(ECFP4)',
          'Chemprop D-MPNN', 'LightGBM\n(ECFP4 + RDKit)']
rmse, err = [0.801, 0.509, 0.504, 0.445], [0, 0, 0.004, 0]
y = np.arange(len(models))[::-1]
fig, ax = plt.subplots(figsize=(6.0, 2.9))
ax.barh(y, rmse, height=0.52, color=BLUE, zorder=3)
ax.errorbar(rmse, y, xerr=err, fmt='none', ecolor=INK, elinewidth=1.1, capsize=3, zorder=4)
ax.axvline(0.509, color=VERM, ls='--', lw=1.2, zorder=2)
ax.text(0.515, 3.42, 'similarity baseline', color=VERM, fontsize=8, va='center')
for yi, v in zip(y, rmse):
    ax.text(v + 0.012, yi, f'{v:.3f}', va='center', fontsize=8.5)
ax.set_yticks(y); ax.set_yticklabels(models)
ax.set_xlabel('RMSE on held-out scaffold split (pChEMBL units)')
ax.set_xlim(0, 0.92)
ax.xaxis.grid(True, color='#eeeeee', lw=0.8, zorder=0); ax.set_axisbelow(True)
ax.set_title('Gradient boosting outperforms the graph neural network',
             loc='left', fontsize=10, fontweight='bold', pad=8)
save(fig, 'fig1_model_comparison')

# ---------------------------------------------------------------- figure 2
rows = list(csv.DictReader(open('results/crossdock.csv')))
dry = {r['ligand']: float(r['top1_rmsd']) for r in rows if r['receptor'] == 'dry'}
wet = {r['ligand']: float(r['top1_rmsd']) for r in rows if r['receptor'] == 'wet'}
lg = sorted(wet, key=lambda k: -wet[k])
d = np.array([dry[k] for k in lg]); w = np.array([wet[k] for k in lg])
yy = np.arange(len(lg))
fig, ax = plt.subplots(figsize=(6.2, 4.6))
ax.axvspan(0, 2.0, color='#f2f7fa', zorder=0)
ax.axvline(2.0, color=GREY, ls='--', lw=1.0, zorder=2)
for i in yy:
    ax.plot([d[i], w[i]], [i, i], color='#d8d8d8', lw=1.6, zorder=2, solid_capstyle='round')
ax.scatter(d, yy, s=42, marker='o', color=VERM, zorder=4, label='dry receptor',
           edgecolor='white', linewidth=0.8)
ax.scatter(w, yy, s=46, marker='s', color=BLUE, zorder=4,
           label='wet receptor (5 waters)', edgecolor='white', linewidth=0.8)
ax.set_yticks(yy); ax.set_yticklabels(lg, fontsize=7.5, family='monospace')
ax.set_xlabel('Top-ranked pose RMSD vs crystal structure (Å)')
ax.set_xlim(0, max(d.max(), w.max()) + 0.8); ax.set_ylim(-0.8, len(lg) - 0.2)
ax.xaxis.grid(True, color='#f0f0f0', lw=0.8, zorder=1); ax.set_axisbelow(True)
ax.legend(loc='upper right', frameon=False, fontsize=8.5)
ax.set_title(f"Bridging waters: {(w < 2).sum()}/{len(w)} correct vs {(d < 2).sum()}/{len(d)}",
             loc='left', fontsize=10, fontweight='bold', pad=8)
save(fig, 'fig2_dry_vs_wet')

# ---------------------------------------------------------------- figure 3
def dih(p):
    b0, b1, b2 = p[0]-p[1], p[2]-p[1], p[3]-p[2]
    b1 = b1/np.linalg.norm(b1)
    v = b0 - (b0@b1)*b1; w_ = b2 - (b2@b1)*b1
    return np.degrees(np.arctan2(np.cross(b1, v)@w_, v@w_))

def chi1(path, n):
    a = {}
    for l in open(path):
        if l[:4]=='ATOM' and l[21]=='A' and l[22:26].strip()==str(n) and l[16] in ' A':
            a[l[12:16].strip()] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    k = ('N','CA','CB','CG')
    return dih([a[x] for x in k]) if all(x in a for x in k) else None

status = {}
for r in rows:
    if r['receptor'] == 'wet':
        status[r['ligand'].split('_')[0]] = 'Y' if float(r['best_rmsd']) < 2.0 else 'n'
status['6OAK'] = 'ref'

pts = []
for pid, st in sorted(status.items()):
    p = f'structures/pdb/{pid}.pdb'
    if not os.path.exists(p): continue
    c1, c2 = chi1(p, 795), chi1(p, 902)
    if c1 is not None and c2 is not None:
        pts.append((pid, st, c1, c2))

if pts:
    g = np.array([p[1] for p in pts])
    x = np.array([p[2] for p in pts]); yv = np.array([p[3] for p in pts])
    fig, ax = plt.subplots(figsize=(5.8, 4.6))
    m = g == 'n'
    ax.scatter(x[m], yv[m], s=58, marker='X', color=VERM, zorder=4,
               edgecolor='white', linewidth=0.8, label='docking failed')
    m = g == 'Y'
    ax.scatter(x[m], yv[m], s=52, marker='o', color=BLUE, zorder=4,
               edgecolor='white', linewidth=0.8, label='docked < 2 Å')
    m = g == 'ref'
    ax.scatter(x[m], yv[m], s=95, marker='D', facecolor='none', edgecolor=INK,
               linewidth=1.4, zorder=5, label='6OAK (reference)')
    ax.set_xlabel('Tyr795  $\\chi_1$  (degrees)')
    ax.set_ylabel('Phe902  $\\chi_1$  (degrees)')
    ax.grid(True, color='#f2f2f2', lw=0.8); ax.set_axisbelow(True)
    ax.legend(loc='lower left', frameon=False, fontsize=8.5)
    ax.set_title('Two coupled rotamer states partition the inhibitors',
                 loc='left', fontsize=10, fontweight='bold', pad=8)
    for lab, xs, ys in (('State A', x[g!='n'].mean(), yv[g!='n'].mean()+28),
                        ('State B', x[g=='n'].mean(), yv[g=='n'].mean()-26)):
        ax.annotate(lab, xy=(xs, ys), fontsize=9, ha='center', color='#555555')
    save(fig, 'fig3_rotamer_states')
else:
    print("  fig3 skipped - structures/pdb/*.pdb not found")

# ---------------------------------------------------------------- figure 4
cr = list(csv.DictReader(open('results/control.csv')))
t = np.array([float(r['vina_kcal']) for r in cr if r['group'] == 'top'])
rr = np.array([float(r['vina_kcal']) for r in cr if r['group'] == 'random'])
bins = np.linspace(min(t.min(), rr.min()) - 0.3, max(t.max(), rr.max()) + 0.3, 26)

fig, ax = plt.subplots(figsize=(6.2, 3.6))
ax.hist(rr, bins=bins, color=VERM, alpha=0.55, label=f'random (n={len(rr)})', zorder=3)
ax.hist(t, bins=bins, color=BLUE, alpha=0.55, label=f'model top-ranked (n={len(t)})', zorder=3)
ax.axvline(rr.mean(), color=VERM, ls='--', lw=1.4, zorder=4)
ax.axvline(t.mean(), color=BLUE, ls='--', lw=1.4, zorder=4)
ax.set_xlabel('Vina binding score (kcal/mol)')
ax.set_ylabel('compounds')
ax.yaxis.grid(True, color='#f0f0f0', lw=0.8, zorder=0); ax.set_axisbelow(True)
ax.legend(loc='upper left', frameon=False, fontsize=8.5)
ax.text(0.98, 0.94, f"means  {t.mean():.2f}  vs  {rr.mean():.2f}\n"
        f"8 top-ranked failed to dock",
        transform=ax.transAxes, ha='right', va='top', fontsize=8.5, color='#555555')
ax.set_title('Model-ranked compounds vs random selection',
             loc='left', fontsize=10, fontweight='bold', pad=8)
save(fig, 'fig4_random_control')

print("\nall figures written to figures/")
