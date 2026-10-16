import numpy as np

def heavy_coords(path, model=1):
    out, cur = [], 0
    for line in open(path):
        if line.startswith('MODEL'):
            cur = int(line.split()[1])
        if line.startswith(('ATOM', 'HETATM')):
            if cur and cur != model:
                continue
            name = line[12:16].strip()
            if name.startswith('H') or line[76:78].strip() == 'H':
                continue
            out.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(out)

ref = heavy_coords('structures/6OAK_ligand.pdb')
print(f"crystal pose: {len(ref)} heavy atoms")

for mode in range(1, 4):
    pose = heavy_coords('results/redock_dry.pdbqt', model=mode)
    if len(pose) == 0:
        break
    if len(pose) != len(ref):
        print(f"mode {mode}: {len(pose)} atoms vs {len(ref)} in crystal - cannot compare")
        continue
    rmsd = np.sqrt(((pose - ref) ** 2).sum(axis=1).mean())
    verdict = "PASS" if rmsd < 2.0 else ("marginal" if rmsd < 3.0 else "FAIL")
    print(f"mode {mode}: RMSD {rmsd:.2f} A   {verdict}")
