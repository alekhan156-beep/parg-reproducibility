import sys, csv
import numpy as np
from scipy import stats
rows = list(csv.DictReader(open(sys.argv[1])))
g = {}
for r in rows:
    g.setdefault(r['group'], []).append(float(r['vina_kcal']))
top, rnd = np.array(g['top']), np.array(g['random'])
print(f"{'':12}{'n':>4}{'mean':>8}{'median':>9}{'best':>8}{'<-9':>6}   record says")
print(f"{'top-ranked':12}{len(top):4d}{top.mean():8.2f}{np.median(top):9.2f}"
      f"{top.min():8.2f}{(top<-9).sum():6d}   92 / -8.41 / -8.42 / -10.21 / 20")
print(f"{'random':12}{len(rnd):4d}{rnd.mean():8.2f}{np.median(rnd):9.2f}"
      f"{rnd.min():8.2f}{(rnd<-9).sum():6d}   100 / -8.18 / -8.17 / -10.60 / 15")
U = stats.mannwhitneyu(top, rnd)
print(f"\ndifference        {rnd.mean()-top.mean():+.3f} kcal/mol   (record +0.227)")
print(f"Mann-Whitney U    p = {U.pvalue:.4f}            (record 0.0154)")
ct = [[int((top<-9).sum()), int((top>=-9).sum())],
      [int((rnd<-9).sum()), int((rnd>=-9).sum())]]
print(f"Fisher enrichment p = {stats.fisher_exact(ct)[1]:.4f}            (record 0.154)")
print(f"enrichment        {(top<-9).mean()/(rnd<-9).mean():.2f}-fold        (record 1.45)")
print(f"\nbest of all 192   {min(top.min(), rnd.min()):.2f}  from the "
      f"{'RANDOM' if rnd.min() < top.min() else 'TOP'} group   (record: random)")
for lab, v in (('as run', None), ('failures at -6.42', -6.42), ('failures at -10.21', -10.21)):
    t = top if v is None else np.concatenate([top, np.full(8, v)])
    p = stats.mannwhitneyu(t, rnd).pvalue
    print(f"  {lab:20} diff {rnd.mean()-t.mean():+.3f}   p = {p:.4f}")
