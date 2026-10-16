import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

df = pd.read_csv('data/parg_raw.csv')
print(f"Start: {len(df)}")

df = df[df['data_validity_comment'].isna()]
df = df[df['standard_relation'] == '=']
df = df[df['canonical_smiles'].notna()]
df = df[df['pchembl_value'].notna()]
df['pchembl_value'] = pd.to_numeric(df['pchembl_value'], errors='coerce')
df = df[df['pchembl_value'].notna()]
print(f"After validity filters: {len(df)}")

chooser   = rdMolStandardize.LargestFragmentChooser()
uncharger = rdMolStandardize.Uncharger()

def standardise(smi):
    try:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        m = chooser.choose(m)
        m = uncharger.uncharge(m)
        if m.GetNumHeavyAtoms() > 60:
            return None
        return Chem.MolToSmiles(m)
    except Exception:
        return None

df['clean_smiles'] = df['canonical_smiles'].apply(standardise)
df = df[df['clean_smiles'].notna()]
df['inchikey'] = df['clean_smiles'].apply(
    lambda s: Chem.MolToInchiKey(Chem.MolFromSmiles(s)))
print(f"After standardisation: {len(df)}")

agg = df.groupby('inchikey').agg(
    smiles      = ('clean_smiles', 'first'),
    pchembl     = ('pchembl_value', 'median'),
    pchembl_std = ('pchembl_value', 'std'),
    n_records   = ('pchembl_value', 'size'),
    assay_type  = ('standard_type', 'first'),
    assay_id    = ('assay_chembl_id', 'first'),
    year        = ('document_year', 'min'),
).reset_index()
print(f"Unique compounds: {len(agg)}")

agg['label'] = agg['pchembl']

agg['scaffold'] = agg['smiles'].apply(
    lambda s: MurckoScaffold.MurckoScaffoldSmiles(smiles=s, includeChirality=False))

agg.to_csv('data/parg_curated.csv', index=False)

print()
print(f"FINAL: {len(agg)} compounds")
print(f"  actives   (pChEMBL >= 6): {(agg.label == 1).sum()}")
print(f"  inactives (pChEMBL <= 5): {(agg.label == 0).sum()}")
print(f"  unique scaffolds: {agg.scaffold.nunique()}")
