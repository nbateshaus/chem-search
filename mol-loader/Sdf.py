"""
Encapsulate reading and formatting of an SDF collection
"""

import glob
import gzip
import os
import os.path

from rdkit import Chem


class Sdf:
    """
    Encapsulate SDF files found at a particular path.

    Subclasses must provide method mol_to_dict
    """

    def __init__(self, path, limit=None):
        """
        Encapsulate SDF files found at a particular path.

        Args:
        path - a file, directory or glob pattern for the SDF files
        """
        self.limit = limit
        print("Searching for {0}".format(path))
        self.files = self.glob(path)

    def __iter__(self):
        generated = 0
        for name in self.files:
            print("Reading {0}".format(name))
            f = gzip.open(name)
            try:
                f.read(1) # Will raise if f is not GZip compressed
                f.rewind()
            except IOError:
                f.close()
                f = open(name, 'rb')
            # Be as permissive as possible on read (sanitize, removeHs, strictParsing)
            # mol_to_dict can do cleanup.
            mols = Chem.ForwardSDMolSupplier(f, sanitize=False, removeHs=False, strictParsing=False)
            molnum = 0
            for mol in mols:
                molnum += 1
                if mol is None:
                    print("Failed to read molecule {0} from {1}".format(molnum, name))
                    continue
                yield self.mol_to_dict(mol)
                generated += 1
                if generated == self.limit:
                    break
            f.close()
            if generated == self.limit:
                break

    @staticmethod
    def glob(path):
        """
        Extend system glob to match all files in a directory.

        :param path: Directory or glob pattern
        :return: List of all files in or matching path
        """
        path = os.path.abspath(path)
        if os.path.isdir(path):
            files = [d for d in [
                    os.path.join(path, f) for f in os.listdir(path)
                ] if os.path.isfile(d)]
        else:
            files = glob.glob(path)
        print("Found {0} files".format(len(files)))
        return files

    @staticmethod
    def cast(val):
        """
        Given a string, guess the type of the value it represents, and return the string cast to that type.

        :param val: String representation of some value
        :return: bool, int, float or str derived from val
        """
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
