import os, csv
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
os.makedirs('figures', exist_ok=True)
r = list(csv.DictReader(open('results/control.csv')))
t = np.array([float(x['vina_kcal']) for x in r if x['group'] == 'top'])
q = np.array([float(x['vina_kcal']) for x in r if x['group'] == 'random'])
b = np.linspace(min(t.min(), q.min()) - .3, max(t.max(), q.max()) + .3, 26)
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.facecolor': 'white'})
fig, ax = plt.subplots(figsize=(6.2, 3.6))
ax.hist(q, bins=b, color='#D55E00', alpha=.55, label='random (n=%d)' % len(q))
ax.hist(t, bins=b, color='#0072B2', alpha=.55, label='model top-ranked (n=%d)' % len(t))
ax.axvline(q.mean(), color='#D55E00', ls='--', lw=1.4)
ax.axvline(t.mean(), color='#0072B2', ls='--', lw=1.4)
ax.set_xlabel('Vina binding score (kcal/mol)')
ax.set_ylabel('compounds')
ax.legend(loc='upper left', frameon=False, fontsize=8.5)
ax.text(.98, .94, 'means %.2f vs %.2f\n8 top-ranked failed to dock' % (t.mean(), q.mean()),
        transform=ax.transAxes, ha='right', va='top', fontsize=8.5, color='#555555')
ax.set_title('Model-ranked compounds vs random selection',
             loc='left', fontsize=10, fontweight='bold', pad=8)
for e in ('png', 'pdf'):
    fig.savefig('figures/fig4_random_control.' + e, dpi=300, bbox_inches='tight')
print('wrote figures/fig4_random_control.png and .pdf')
