from scipy import stats
import numpy as np, csv
D = {  # ligand: (dry_lb, dry_true, wet_lb, wet_true)   all 17, from results/crossdock_truermsd.csv
'5LHB_6WJ':(4.977,5.454,5.254,5.977),'6HMK_7JC':(3.88,4.38,4.81,5.32),
'6HML_73L':(5.50,6.54,5.67,6.20),'6HMM_7JB':(6.10,7.67,5.35,7.60),
'6HMN_70J':(5.08,7.54,7.82,8.27),'6O9X_M0S':(2.56,3.74,0.62,0.62),
'6O9Y_M0V':(7.04,7.55,6.47,6.62),'6OA0_M1D':(5.62,7.01,4.71,7.10),
'6OA1_M0P':(2.58,3.68,0.69,0.69),'6OA3_M0M':(1.91,3.68,1.99,2.03),
'6OAL_M0P':(3.46,3.66,0.57,0.57),'7KFP_WD7':(1.96,2.73,1.88,1.94),
'7KG0_WDA':(5.66,9.08,6.35,6.43),'7KG1_WDD':(6.26,6.87,4.18,4.67),
'7KG6_WDG':(3.55,3.63,5.98,6.42),'7KG7_WDJ':(3.93,4.59,0.74,0.74),
'7KG8_WDM':(3.11,3.36,0.62,0.62)}
DISP = {'5LHB_6WJ','6HMK_7JC','6HML_73L','6HMM_7JB','6HMN_70J'}
np.save('data.npy', D, allow_pickle=True)

dry = np.array([v[1] for v in D.values()]); wet = np.array([v[3] for v in D.values()])
lbd = np.array([v[0] for v in D.values()]); lbw = np.array([v[2] for v in D.values()])
names = list(D.keys()); isd = np.array([n in DISP for n in names])

print(f'n = {len(D)} ligands, 34 pairs\n')
print('OVERALL  (threshold 2.0 A)')
print(f'  dry   lower bound {(lbd<2).sum()}/17   TRUE {(dry<2).sum()}/17')
print(f'  wet   lower bound {(lbw<2).sum()}/17   TRUE {(wet<2).sum()}/17')
b = int(((wet<2)&~(dry<2)).sum()); c = int(((dry<2)&~(wet<2)).sum())
print(f'  McNemar exact (paired): {b} wet-only, {c} dry-only -> p = {stats.binomtest(b,b+c,0.5).pvalue:.4f}')
print(f'  Wilcoxon signed-rank on paired RMSD          -> p = {stats.wilcoxon(wet,dry).pvalue:.4f}')
print(f'  median  dry {np.median(dry):.2f}   wet {np.median(wet):.2f}  (NOT the right summary - see split)\n')

print('BY INHIBITOR CLASS  (classes assigned in sec 8 from crystal water contacts)')
for lab, mask in (('water-bridging', ~isd), ('water-displacing', isd)):
    d_, w_ = dry[mask], wet[mask]; g = w_-d_
    print(f'  {lab:17} n={mask.sum():2d}   dry {(d_<2).sum()}/{mask.sum()}   wet {(w_<2).sum()}/{mask.sum()}'
          f'   mean wet-dry {g.mean():+.2f} A   improved {(g<0).sum()}/{mask.sum()}')
gb, gd = (wet-dry)[~isd], (wet-dry)[isd]
print(f'  Mann-Whitney on water effect, bridging vs displacing -> p = {stats.mannwhitneyu(gb,gd).pvalue:.4f}')
bb, cc = int(((wet<2)&~(dry<2))[~isd].sum()), int(((dry<2)&~(wet<2))[~isd].sum())
print(f'  McNemar within bridging class: {bb} vs {cc} -> p = {stats.binomtest(bb,bb+cc,0.5).pvalue:.4f}\n')

print('BOUND TIGHTNESS  (validates that the old numbers flattered failures)')
sh = np.concatenate([wet-lbw, dry-lbd]); ok = np.concatenate([wet<2, dry<2])
print(f'  correct poses   n={ok.sum():2d}  mean shift {sh[ok].mean():+.3f} A  max {sh[ok].max():.2f}')
print(f'  incorrect poses n={(~ok).sum():2d}  mean shift {sh[~ok].mean():+.3f} A  max {sh[~ok].max():.2f}')
print(f'  Mann-Whitney -> p = {stats.mannwhitneyu(sh[ok],sh[~ok]).pvalue:.5f}')
print(f'  every shift is >= 0 : {bool((sh>=-0.005).all())}  (optimal assignment is a true lower bound)')
