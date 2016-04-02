"""
Encapsulate reading and formatting of SureChEMBL
"""

import json

from Sdf import Sdf
from rdkit_utils import rdkit_standardize, rdkit_descriptors, rdkit_smiles, rdkit_fps_from_mol


class Surechembl(Sdf):
    def __init__(self, path, limit=None):
        Sdf.__init__(self, path, limit)

    def mol_to_dict(self, mol):
        """
        Capture all the information from a SureChEMBL molecule into a dictionary

        mol_to_dict( (Surechembl)self, (rdkit.Chem.Mol)mol) -> dict
        """
        props = set(mol.GetPropNames())
        d = {p.lower(): self.cast(mol.GetProp(p)) for p in props}
        d['id'] = 'https://www.surechembl.org/chemical/' + mol.GetProp('ID')

        smiles = []
        rdkit_mol = rdkit_standardize(mol)
        if rdkit_mol is not None:
            fps, fps_bits = rdkit_fps_from_mol(rdkit_mol)
            d['rdkit_fingerprint'] = fps
            d['rdkit_fingerprint_bits'] = fps_bits

            rs = rdkit_smiles(rdkit_mol)
            if rs is not None:
                d['rdkit_smiles'] = rs
                smiles.append(rs)

            descs = rdkit_descriptors(rdkit_mol)
            for name in descs:
                d['rdkit_' + name.lower()] = descs[name]

        if smiles:
            d['smiles'] = list(set(smiles))

        return d


def test_mol_to_dict():
    mols = Surechembl(path="resources/surechembl-test-data.sdf*")
    for mol in mols:
        print(json.dumps(mol))

if __name__ == "__main__":
    test_mol_to_dict()
