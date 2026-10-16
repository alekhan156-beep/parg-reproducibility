import pandas as pd, numpy as np
from collections import defaultdict

df = pd.read_csv('data/parg_curated.csv')
rng = np.random.default_rng(42)

# scaffold split: whole families to one side
groups = defaultdict(list)
for i, s in enumerate(df.scaffold):
    groups[s].append(i)

train_idx, test_idx = [], []
for g in sorted(groups.values(), key=len, reverse=True):
    if len(train_idx) < 0.7 * len(df):
        train_idx += g
    else:
        test_idx += g

df['split_scaffold'] = 'test'
df.loc[train_idx, 'split_scaffold'] = 'train'

# random split, for comparison
perm = rng.permutation(len(df))
cut = int(0.7 * len(df))
df['split_random'] = 'test'
df.loc[perm[:cut], 'split_random'] = 'train'

# temporal split
df['split_time'] = np.where(df.year <= 2019, 'train', 'test')

df.to_csv('data/parg_split.csv', index=False)

for col in ['split_scaffold', 'split_random', 'split_time']:
    tr = (df[col] == 'train').sum()
    te = (df[col] == 'test').sum()
    print(f"{col:16s} train {tr:4d}  test {te:4d}")
