from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np

KEEP = {1211, 1232, 1307, 1317, 1402}
p = PDBParser(QUIET=True)
st = p.get_structure('parg', 'structures/6OAK.pdb')
ch = st[0]['A']

cme = [r for r in ch if r.get_resname() == 'CME']
lig = [r for r in ch if r.get_resname() == 'M0M'][0]
lc = np.array([a.coord for a in lig]).mean(axis=0)
for r in cme:
    d = min(np.linalg.norm(a.coord - lc) for a in r)
    print(f"CME at residue {r.id[1]}: {d:.1f} A from ligand centroid")

class Rec(Select):
    def __init__(self, waters): self.waters = waters
    def accept_chain(self, c):   return c.id == 'A'
    def accept_residue(self, r):
        if r.get_resname() == 'M0M': return 0
        if r.id[0] == 'W':           return 1 if r.id[1] in self.waters else 0
        return 1

io = PDBIO(); io.set_structure(st)
io.save('structures/6OAK_receptor_dry.pdb', Rec(set()))
io.save('structures/6OAK_receptor_wet.pdb', Rec(KEEP))
print("\nwrote 6OAK_receptor_dry.pdb  (no waters)")
print("wrote 6OAK_receptor_wet.pdb  (5 bridging waters kept)")
