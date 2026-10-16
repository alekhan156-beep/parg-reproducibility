#!/usr/bin/env python3
"""
Figure 2, rebuilt on the corrected (symmetry-corrected, non-superposed) RMSD
values from script 38, and grouped by the inhibitor classes assigned in
section 8 from crystallographic water contacts.

Supersedes the fig2 block of scripts/36_figures.py, which read the
optimal-assignment lower bounds in results/crossdock.csv.

Run from ~/parg :   python scripts/39_fig2.py
Reads   results/crossdock_truermsd.csv
Writes  figures/fig2_dry_vs_wet.png / .pdf
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

# assigned in section 8 from crystal water contacts, before these RMSD values
# were computed -- the split is a priori, not fitted to this outcome
DISPLACING = {'5LHB_6WJ', '6HMK_7JC', '6HML_73L', '6HMM_7JB', '6HMN_70J'}

rows = list(csv.DictReader(open('results/crossdock_truermsd.csv')))
dry = {r['ligand']: float(r['true_rmsd']) for r in rows if r['receptor'] == 'dry'}
wet = {r['ligand']: float(r['true_rmsd']) for r in rows if r['receptor'] == 'wet'}

brid = sorted([k for k in wet if k not in DISPLACING], key=lambda k: -wet[k])
disp = sorted([k for k in wet if k in DISPLACING], key=lambda k: -wet[k])
lg = disp + brid                      # displacing at the bottom of the axis
gap = len(disp) - 0.5                 # y position of the class divider

d = np.array([dry[k] for k in lg]); w = np.array([wet[k] for k in lg])
yy = np.arange(len(lg))

fig, ax = plt.subplots(figsize=(6.6, 5.0))
ax.axvspan(0, 2.0, color='#f2f7fa', zorder=0)
ax.axvline(2.0, color=GREY, ls='--', lw=1.0, zorder=2)
for i in yy:
    ax.plot([d[i], w[i]], [i, i], color='#d8d8d8', lw=1.6,
            zorder=2, solid_capstyle='round')
ax.scatter(d, yy, s=42, marker='o', color=VERM, zorder=4,
           label='dry receptor', edgecolor='white', linewidth=0.8)
ax.scatter(w, yy, s=46, marker='s', color=BLUE, zorder=4,
           label='wet receptor (5 waters)', edgecolor='white', linewidth=0.8)

xmax = max(d.max(), w.max()) + 1.9
ax.axhline(gap, color='#cccccc', lw=0.9, ls='-', zorder=3)
ax.text(xmax - 0.15, len(disp) + len(brid) / 2 - 0.5,
        f'water-bridging\n{(w[len(disp):] < 2).sum()}/{len(brid)} correct',
        ha='right', va='center', fontsize=8.5, color=BLUE, fontweight='bold',
        linespacing=1.5)
ax.text(xmax - 0.15, (len(disp) - 1) / 2,
        f'water-displacing\n{(w[:len(disp)] < 2).sum()}/{len(disp)} correct',
        ha='right', va='center', fontsize=8.5, color='#777777',
        fontweight='bold', linespacing=1.5)

ax.set_yticks(yy); ax.set_yticklabels(lg, fontsize=7.5, family='monospace')
for t, k in zip(ax.get_yticklabels(), lg):
    t.set_color('#888888' if k in DISPLACING else INK)
ax.set_xlabel('Top-ranked pose RMSD vs crystal structure (Å),'
              ' symmetry-corrected, not superposed')
ax.set_xlim(0, xmax); ax.set_ylim(-0.8, len(lg) - 0.2)
ax.xaxis.grid(True, color='#f0f0f0', lw=0.8, zorder=1); ax.set_axisbelow(True)
ax.legend(loc='lower left', frameon=False, fontsize=8.5, ncol=2,
          bbox_to_anchor=(0.0, 1.005), handletextpad=0.4, columnspacing=1.6)
ax.set_title(f'Bridging waters are required for correct placement: '
             f'{(w < 2).sum()}/{len(w)} vs {(d < 2).sum()}/{len(d)}',
             loc='left', fontsize=10, fontweight='bold', pad=26)
for ext in ('png', 'pdf'):
    fig.savefig(f'figures/fig2_dry_vs_wet.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)
print(f'  wet {(w<2).sum()}/{len(w)}   dry {(d<2).sum()}/{len(d)}')
print('  wrote figures/fig2_dry_vs_wet.png / .pdf')
