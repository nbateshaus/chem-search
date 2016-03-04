"""
Encapsulate reading and formatting of an SDF collection
"""

import glob
import gzip
import os, os.path
import rdkit.Chem

class Sdf:
    '''
    Encapsulate SDF files found at a particular path.

    Subclasses must provide method mol_to_dict
    '''

    def __init__(self, path):
        '''
        Encapsulate SDF files found at a particular path.

        Args:
        path - a file, directory or glob pattern for the SDF files
        '''
        print("Searching for {0}".format(path))
        self.glob(path)

    def __iter__(self):
        for name in self.files:
            print("Reading {0}".format(name))
            f = gzip.open(name)
            try:
                f.read(1) # Will raise if f is not GZip compressed
                f.rewind()
            except IOError:
                f = open(name)
            # Be as permissive as possible on read (sanitize, removeHs, strictParsing)
            # mol_to_dict can do cleanup.
            mols = rdkit.Chem.ForwardSDMolSupplier(f, sanitize=False, removeHs=False, strictParsing=False)
            molnum = 0
            for mol in mols:
                molnum += 1
                if mol is None:
                    print("Failed to read molecule {0} from {1}".format(molnum, name))
                    continue
                yield self.mol_to_dict(mol)
            f.close()

    def glob(self, path):
        path = os.path.abspath(path)
        if (os.path.isdir(path)):
            self.files = [d for d in [
                    os.path.join(path, f) for f in os.listdir(path)
                ] if os.path.isfile(d)]
        else:
            self.files = glob.glob(path)
        print("Found {0} files".format(len(self.files)))

    def canonicalize(self, mol):
        """
        Generate a canonical representation of a molecule.

        canonicalize( (rdkit.Chem.Mol)mol ) -> rdkit.Chem.Mol
        """
        # Attempt normalization; this may fail if there are ... oddities ... in the molecule.
        try:
            rdkit_mol = rdkit.Chem.RemoveHs(mol)
        except ValueError:
            rdkit_mol = None
        return rdkit_mol

    def cast(self, val):
        if val.lower() == str(True).lower():
            return True
        elif val.lower() == str(False).lower():
            return False

        try:
            return int(val)
        except ValueError:
            pass

        try:
            return float(val)
        except ValueError:
            pass

        return val
