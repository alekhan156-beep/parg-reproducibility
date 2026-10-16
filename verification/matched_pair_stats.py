from scipy import stats
import numpy as np
# ligand: (dry_best, dry_rank, wet_best, wet_rank)   transcribed from the run
D = {
'5LHB_6WJ':(None,None,4.81,8),          # dry row scrolled off; top-1 5.45, not <2
'6HMK_7JC':(2.42,4,3.99,7), '6HML_73L':(6.18,2,6.20,1),
'6HMM_7JB':(7.31,2,7.60,1), '6HMN_70J':(5.70,4,6.60,4),
'6O9X_M0S':(1.05,6,0.62,1), '6O9Y_M0V':(4.03,8,4.41,6),
'6OA0_M1D':(4.37,9,4.54,3), '6OA1_M0P':(0.94,5,0.69,1),
'6OA3_M0M':(1.10,5,1.13,2), '6OAL_M0P':(0.53,6,0.57,1),
'7KFP_WD7':(1.79,4,1.94,1), '7KG0_WDA':(6.43,3,6.43,1),
'7KG1_WDD':(4.37,5,4.10,3), '7KG6_WDG':(3.63,1,3.56,6),
'7KG7_WDJ':(2.53,8,0.74,1), '7KG8_WDM':(0.67,3,0.62,1)}
T=2.0
sd={k for k,v in D.items() if v[0] is not None and v[0]<T}
sw={k for k,v in D.items() if v[2] is not None and v[2]<T}
print(f'SAMPLING (near-native pose present among the 9 modes)')
print(f'  dry {len(sd)}/17   wet {len(sw)}/17')
print(f'  dry set subset of wet set: {sd <= sw}')
print(f'  wet-only: {sorted(sw-sd)}   dry-only: {sorted(sd-sw)}')
b,c=len(sw-sd),len(sd-sw)
print(f'  McNemar exact -> p = {stats.binomtest(b,b+c,0.5).pvalue:.4f}  (no sampling difference)\n')

common=sorted(sd&sw)
print(f'RANKING of the near-native pose, {len(common)} ligands sampled in BOTH')
print(f"  {'ligand':12} {'dry rank':>9} {'wet rank':>9}")
for k in common: print(f'  {k:12} {D[k][1]:9d} {D[k][3]:9d}')
dr=np.array([D[k][1] for k in common]); wr=np.array([D[k][3] for k in common])
print(f'  mean rank  dry {dr.mean():.2f}   wet {wr.mean():.2f}')
print(f'  wet ranks better in {(wr<dr).sum()}/{len(common)}, worse in {(wr>dr).sum()}')
print(f'  Wilcoxon signed-rank on rank -> p = {stats.wilcoxon(wr,dr).pvalue:.4f}')
f1d=sum(D[k][1]==1 for k in common); f1w=sum(D[k][3]==1 for k in common)
bb=sum(1 for k in common if D[k][3]==1 and D[k][1]!=1)
cc=sum(1 for k in common if D[k][1]==1 and D[k][3]!=1)
print(f'  ranked 1st: dry {f1d}/{len(common)}  wet {f1w}/{len(common)}')
print(f'  McNemar on ranked-1st -> p = {stats.binomtest(bb,bb+cc,0.5).pvalue:.4f}\n')
print('AS THE SCRIPT REPORTED (unmatched denominators)')
print(f'  dry 0/{len(sd)} ranked 1st     wet {sum(D[k][3]==1 for k in sw)}/{len(sw)} ranked 1st')
