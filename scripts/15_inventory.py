import os, urllib.request
import numpy as np
from collections import defaultdict

IDS = ['4A0D','4B1G','4B1H','4B1I','4B1J','5A7R','5LHB','6HH6','6HMK','6HML',
       '6HMM','6HMN','6O9X','6O9Y','6OA0','6OA1','6OA3','6OAL','7KFP','7KG0',
       '7KG1','7KG6','7KG7','7KG8']

SKIP = {'HOH','GOL','EDO','SO4','PO4','CL','NA','MG','ZN','CA','ACT','DMS',
        'PEG','PG4','TRS','MES','EPE','CME','FMT','IOD','K','BR','NO3','CIT'}

os.makedirs('structures/pdb', exist_ok=True)
ok = []
for pid in IDS:
    path = f'structures/pdb/{pid}.pdb'
    if not os.path.exists(path):
        try:
            urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pid}.pdb', path)
        except Exception as e:
            print(f'{pid}: download failed - {e}')
            continue
    ok.append(pid)
print(f"downloaded {len(ok)} structures\n")

print(f"{'PDB':6} {'lig':5} {'resnum':>7} {'heavy':>6}")
print("-" * 30)
for pid in ok:
    ligs = defaultdict(list)
    for line in open(f'structures/pdb/{pid}.pdb'):
        if line.startswith('HETATM') and line[21] == 'A':
            res = line[17:20].strip()
            if res in SKIP or line[76:78].strip() == 'H':
                continue
            ligs[(res, line[22:27].strip())].append(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not ligs:
        print(f"{pid:6} (none)")
    for (res, num), xyz in sorted(ligs.items(), key=lambda k: -len(k[1])):
        print(f"{pid:6} {res:5} {num:>7} {len(xyz):6}")
