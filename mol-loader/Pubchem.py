"""
Encapsulate reading and formatting of PubChem
"""

import glob
import gzip
import json
import os, os.path
import rdkit.Chem

class Pubchem:
    def __init__(self, path):
        '''
        Encapsulate PubChem files found at a particular path.

        Args:
        path - a file, directory or glob pattern for the SDF files of PubChem
        '''
        self.glob(path)

    def __iter__(self):
        for name in self.files:
            print("Reading {0}".format(name))
            f = gzip.open(name)
            # Leave molecules untouched on reading (sanitize, removeHs, strictParsing)
            # We'll touch them up later.
            mols = rdkit.Chem.ForwardSDMolSupplier(f, sanitize=False, removeHs=False, strictParsing=False)
            molnum = 0
            for mol in mols:
                molnum += 1
                if mol is None:
                    print("Failed to read molecule {0} from {1}".format(molnum, name))
                    continue
                yield self.mol_to_dict(mol)
            f.close()

    SYNONYM_PROPS = [
        'PUBCHEM_IUPAC_OPENEYE_NAME',
        'PUBCHEM_IUPAC_CAS_NAME',
        'PUBCHEM_IUPAC_NAME',
        'PUBCHEM_IUPAC_SYSTEMATIC_NAME',
        'PUBCHEM_IUPAC_TRADITIONAL_NAME'
    ]
    SMILES_PROPS = [
        'PUBCHEM_OPENEYE_CAN_SMILES',
        'PUBCHEM_OPENEYE_ISO_SMILES'
    ]

    def mol_to_dict(self, mol):
        """
        Capture all the information from an rdkit.Chem.Mol into a dictionary
        """
        props = set(mol.GetPropNames())
        # All the properties defined in the molecule, as-is
        d = {p.lower() : self.cast(mol.GetProp(p)) for p in props}
        d['id'] = 'https://pubchem.ncbi.nlm.nih.gov/compound/' + mol.GetProp('PUBCHEM_COMPOUND_CID')

        smiles = [mol.GetProp(p) for p in self.SMILES_PROPS if p in props]
        smiles.append(rdkit.Chem.MolToSmiles(mol)) # Original SMILES
        # Attempt normalization; this may fail if there are ... oddities ...
        # in the molecule.
        rdkit_mol = rdkit.Chem.MolFromSmiles(rdkit.Chem.MolToSmiles(mol)) # Canonicalizes
        if rdkit_mol is not None:
            d['rdkit_smiles'] = rdkit.Chem.MolToSmiles(rdkit_mol)
            smiles.append(d['rdkit_smiles'])
        if smiles:
            d['smiles'] = list(set(smiles))

        synonyms = [mol.GetProp(p) for p in self.SYNONYM_PROPS if p in props]
        if synonyms:
            d['synonyms'] = list(set(synonyms))

        if 'PUBCHEM_COMPOUND_CANONICALIZED' in props:
            canon = mol.GetProp('PUBCHEM_COMPOUND_CANONICALIZED')
            if canon == 1:
                d['pubchem_compound_canonicalized'] = True
            elif canon == 0:
                d['pubchem_compound_canonicalized'] = False
            else:
                del d['pubchem_compound_canonicalized']

        return d

    def glob(self, path):
        path = os.path.abspath(path)
        if (os.path.isdir(path)):
            self.files = [d for d in [
                    os.path.join(path, f) for f in os.listdir(path)
                ] if os.path.isfile(d)]
        else:
            self.files = glob.glob(path)

    def cast(self, val):
        try:
            return float(val)
        except ValueError:
            pass
        return val

def test_glob():
    p = Pubchem(".")
    print p.files
    p = Pubchem("Pubchem.py")
    print p.files
    p = Pubchem("*.py")
    print p.files

def test_iter():
    p = Pubchem(path="/Users/nik/Data/PubChem/Compound_000000001_000025000.sdf.gz")
    mols = [m for m in p]
    for n in range(5):
        print(json.dumps(mols[n]))

if __name__ == "__main__":
    test_glob()
    test_iter()
