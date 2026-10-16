import subprocess
import numpy as np

def counts(p):
    n = h = 0
    for l in open(p):
        if l[:6] in ('ATOM  ','HETATM'):
            n += 1
            if l[76:78].strip().upper().startswith('H'): h += 1
    return n, h

for p in ['structures/6OAK_receptor_wet.pdb', 'structures/rec_wet.pdbqt',
          'structures/rec_B_wet.pdb', 'structures/rec_B.pdbqt']:
    n, h = counts(p)
    print(f"{p:36} atoms {n:6d}   H {h:5d}")

def cen(path, model=1):
    xs, cur = [], 0
    for l in open(path):
        if l.startswith('MODEL'): cur = int(l.split()[1])
        if l[:6] in ('ATOM  ','HETATM') and (not cur or cur == model):
            if l[76:78].strip().upper().startswith('H'): continue
            xs.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(xs).mean(0)

ref = cen('structures/xdock/7KG1_WDD.pdb')
dA = cen('results/xdock/7KG1_WDD_wet.pdbqt')
dB = cen('results/xdock/7KG1_WDD_B.pdbqt')
print(f"\nWDD crystal centroid   {ref.round(2)}")
print(f"docked into A, mode 1  {dA.round(2)}   offset {np.linalg.norm(dA-ref):.2f} A")
print(f"docked into B, mode 1  {dB.round(2)}   offset {np.linalg.norm(dB-ref):.2f} A")

print("\ntight-box self-redock into rec_B ...")
subprocess.run(['./vina_1.2.5_linux_x86_64',
    '--receptor','structures/rec_B.pdbqt',
    '--ligand','structures/xdock/7KG1_WDD.pdbqt',
    '--center_x',f'{ref[0]:.2f}','--center_y',f'{ref[1]:.2f}','--center_z',f'{ref[2]:.2f}',
    '--size_x','16','--size_y','16','--size_z','16',
    '--exhaustiveness','32','--num_modes','9','--seed','42',
    '--out','results/selfredock_B.pdbqt'], capture_output=True)
t = cen('results/selfredock_B.pdbqt')
print(f"tight-box mode 1       {t.round(2)}   offset {np.linalg.norm(t-ref):.2f} A")
