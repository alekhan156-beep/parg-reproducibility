#!/usr/bin/env python3
"""
Figure 2 - random control: affinity distributions, both groups.
Reads results/control.csv (written by 34_control.py).
Writes figures/fig2_random_control.png / .pdf
"""
import os, sys, csv
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

SRC = sys.argv[1] if len(sys.argv) > 1 else 'results/control.csv'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'figures'
BLUE, VERM, GREY, INK = '#0072B2', '#D55E00', '#9a9a9a', '#1a1a1a'
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9,
    'axes.edgecolor': '#cccccc', 'axes.linewidth': 0.8,
    'axes.labelcolor': INK, 'text.color': INK,
    'xtick.color': '#555555', 'ytick.color': '#555555',
    'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.facecolor': 'white', 'axes.facecolor': 'white'})
os.makedirs(OUT, exist_ok=True)

g = {}
for r in csv.DictReader(open(SRC)):
    g.setdefault(r['group'], []).append(float(r['vina_kcal']))
top, rnd = np.array(g['top']), np.array(g['random'])
p2 = stats.mannwhitneyu(top, rnd, alternative='two-sided').pvalue
p1 = stats.mannwhitneyu(top, rnd, alternative='less').pvalue
imp = stats.mannwhitneyu(np.concatenate([top, np.full(8, -6.42)]), rnd,
                         alternative='two-sided').pvalue

fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.0, 5.1), sharex=True,
                             gridspec_kw={'height_ratios': [3.1, 1], 'hspace': 0.10})
lo, hi = min(top.min(), rnd.min()) - .35, max(top.max(), rnd.max()) + .35
bins = np.linspace(lo, hi, 28)
ax.hist(rnd, bins=bins, color=GREY, alpha=.60, label=f'random (n = {len(rnd)})', zorder=3)
ax.hist(top, bins=bins, color=BLUE, alpha=.62, label=f'model top-ranked (n = {len(top)})', zorder=4)
ax.axvline(-9, color=VERM, ls='--', lw=1.1, zorder=5)
ax.text(-9.12, ax.get_ylim()[1]*.60, 'enrichment\nthreshold', fontsize=8,
        color=VERM, va='top', ha='right', linespacing=1.35)
ax.set_ylabel('compounds', fontsize=9)
ax.legend(frameon=False, fontsize=8.6, loc='upper right')
ax.set_title('Model-ranked compounds do not outperform random selection',
             loc='left', fontsize=11, fontweight='bold', pad=12)
ax.text(0.0, 1.005, '', transform=ax.transAxes)

rng = np.random.default_rng(0)
for i, (v, c, lab) in enumerate(((rnd, GREY, 'random'), (top, BLUE, 'top-ranked'))):
    y = i * .95
    bx.scatter(v, np.full(len(v), y) + rng.normal(0, .055, len(v)), s=11,
               color=c, alpha=.62, zorder=3, edgecolor='none')
    bx.scatter([v.mean()], [y], marker='|', s=430, color=INK, zorder=5, lw=1.7)
    bx.text(hi + .04, y, f'mean {v.mean():.2f}', fontsize=8, va='center', color='#555555')
bx.axvline(-9, color=VERM, ls='--', lw=1.1, zorder=2)
bx.set_yticks([0, .95]); bx.set_yticklabels(['random', 'top-ranked'], fontsize=8.8)
bx.set_ylim(-.55, 1.45); bx.set_xlim(lo, hi + 1.05)
ax.set_xlim(lo, hi + 1.05)
bx.spines['left'].set_visible(False); bx.tick_params(axis='y', length=0)
bx.set_xlabel('AutoDock Vina best affinity (kcal/mol)', fontsize=9, labelpad=8)

best = min(top.min(), rnd.min())
bx.annotate('best of all 192\ncame from random',
            xy=(best, 0), xytext=(best + .22, .52), fontsize=7.9, color='#555555',
            ha='left', va='center', linespacing=1.35,
            arrowprops=dict(arrowstyle='-', color='#aaaaaa', lw=.8,
                            connectionstyle='arc3,rad=0.2'))

fig.text(0.013, -0.045,
    f'Mann–Whitney  p = {p2:.3f} two-sided  ({p1:.3f} one-sided, as reported in §16).   '
    f'With the eight failed top-ranked compounds imputed at the worst observed\n'
    f'affinity the difference falls to +0.067 kcal/mol, p = {imp:.3f}.   '
    f'Enrichment below −9 kcal/mol is 21.7% vs 15.0%, Fisher p = 0.264 two-sided — not significant.',
    fontsize=8.1, color='#555555', ha='left', va='top')

for ext in ('png', 'pdf'):
    fig.savefig(f'{OUT}/fig2_random_control.{ext}', dpi=300,
                bbox_inches='tight', facecolor='white')
print(f'top n={len(top)} mean {top.mean():.2f} | random n={len(rnd)} mean {rnd.mean():.2f}')
print(f'two-sided p = {p2:.4f}   one-sided p = {p1:.4f}   imputed p = {imp:.4f}')
print(f'wrote {OUT}/fig2_random_control.png / .pdf')
