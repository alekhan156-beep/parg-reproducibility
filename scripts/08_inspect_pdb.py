from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from collections import Counter

p = PDBParser(QUIET=True)
st = p.get_structure('parg', 'structures/6OAK.pdb')

for model in st:
    for chain in model:
        res = [r for r in chain if r.id[0] == ' ']
        het = [r for r in chain if r.id[0] not in (' ', 'W')]
        wat = [r for r in chain if r.id[0] == 'W']
        if not res and not het:
            continue
        print(f"chain {chain.id}: {len(res)} residues, {len(het)} hetero, {len(wat)} waters")
        if res:
            print(f"   numbering runs {res[0].id[1]} to {res[-1].id[1]}")
        if het:
            print(f"   hetero groups: {Counter(r.get_resname() for r in het)}")
    break
