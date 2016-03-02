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
                    print("Failed to read molecule {0} from {1}", (molnum, name))
                    continue
                yield self.mol_to_dict(mol)
            f.close()

    def mol_to_dict(self, mol):
        """
        Capture all the information from an rdkit.Chem.Mol into a dictionary
        """
        # All the properties defined in the molecule, as-is
        d = {p.lower() : mol.GetProp(p) for p in mol.GetPropNames()}
        # Original SMILES
        d['pubchem_smiles'] = rdkit.Chem.MolToSmiles(mol)
        d['id'] = 'https://pubchem.ncbi.nlm.nih.gov/compound/' + mol.GetProp('PUBCHEM_COMPOUND_CID')
        # Attempt normalization; this may fail if there are ... oddities ...
        # in the molecule.
        mol = rdkit.Chem.MolFromSmiles(d['pubchem_smiles']) # Canonicalizes
        if mol is not None:
            d['rdkit_smiles'] = rdkit.Chem.MolToSmiles(mol)
        return d

    def glob(self, path):
        path = os.path.abspath(path)
        if (os.path.isdir(path)):
            self.files = [d for d in [
                    os.path.join(path, f) for f in os.listdir(path)
                ] if os.path.isfile(d)]
        else:
            self.files = glob.glob(path)

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
