from Bio.PDB import PDBParser, NeighborSearch
import numpy as np

p = PDBParser(QUIET=True)
st = p.get_structure('parg', 'structures/6OAK.pdb')
ch = st[0]['A']

lig = [r for r in ch if r.get_resname() == 'M0M'][0]
lig_atoms = list(lig)
waters = [r for r in ch if r.id[0] == 'W']
prot = [a for r in ch if r.id[0] == ' ' for a in r if a.element in ('N','O')]
ns = NeighborSearch(prot)

print(f"{len(waters)} waters total; checking which bridge M0M to protein\n")
bridging = []
for w in waters:
    wa = list(w)[0]
    dl = min(np.linalg.norm(wa.coord - a.coord) for a in lig_atoms)
    if dl > 3.5:
        continue
    near = ns.search(wa.coord, 3.5)
    if not near:
        continue
    partners = sorted({f"{a.get_parent().get_resname()}{a.get_parent().id[1]}" for a in near})
    bridging.append(w.id[1])
    print(f"  HOH {w.id[1]:4d}  {dl:.2f} A to ligand  ->  {', '.join(partners)}")

print(f"\n{len(bridging)} bridging waters: {bridging}")
np.save('structures/bridging_waters.npy', np.array(bridging))
