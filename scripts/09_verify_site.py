from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np

p = PDBParser(QUIET=True)
st = p.get_structure('parg', 'structures/6OAK.pdb')
ch = st[0]['A']

print("Catalytic residue check (expect D, N, E, E, Y, F):")
for num, expect in [(737,'ASP'), (740,'ASN'), (755,'GLU'),
                    (756,'GLU'), (795,'TYR'), (902,'PHE')]:
    try:
        r = ch[num]
        ok = "OK " if r.get_resname() == expect else "MISMATCH"
        print(f"  {num:4d}  found {r.get_resname():4s} expected {expect:4s}  {ok}")
    except KeyError:
        print(f"  {num:4d}  MISSING from structure")

nums = sorted(r.id[1] for r in ch if r.id[0] == ' ')
gaps = [(a, b) for a, b in zip(nums, nums[1:]) if b - a > 1]
print(f"\nchain breaks: {len(gaps)}")
for a, b in gaps:
    print(f"  gap after {a}, resumes at {b}  ({b-a-1} missing)")

lig = [r for r in ch if r.get_resname() == 'M0M']
if lig:
    coords = np.array([a.coord for a in lig[0]])
    c = coords.mean(axis=0)
    print(f"\nM0M ligand: {len(coords)} atoms")
    print(f"  centroid  x={c[0]:.2f}  y={c[1]:.2f}  z={c[2]:.2f}")
    print(f"  box 24 A cube centred there covers the site")

    class LigOnly(Select):
        def accept_residue(self, r):
            return r.get_resname() == 'M0M'
    io = PDBIO(); io.set_structure(st)
    io.save('structures/6OAK_ligand.pdb', LigOnly())
    print("  saved structures/6OAK_ligand.pdb")
