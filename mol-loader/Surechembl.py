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
        d = {'SureChEMBL_' + p: self.cast(mol.GetProp(p)) for p in props}
        d['id'] = d['SureChEMBL_ID']
        d['URL'] = 'https://www.surechembl.org/chemical/' + d['SureChEMBL_ID']

        smiles = []
        rdkit_mol = rdkit_standardize(mol)
        if rdkit_mol is not None:
            fps, fps_bits = rdkit_fps_from_mol(rdkit_mol)
            d['RDKit_Fingerprint'] = fps
            d['RDKit_Fingerprint_bits'] = fps_bits

            rs = rdkit_smiles(rdkit_mol)
            if rs is not None:
                d['RDKit_SMILES'] = rs
                smiles.append(rs)

            descs = rdkit_descriptors(rdkit_mol)
            for name in descs:
                d['RDKit_' + name] = descs[name]

        if smiles:
            d['SMILES'] = list(set(smiles))

        return d


def test_mol_to_dict():
    mols = Surechembl(path="resources/surechembl-test-data.sdf*")
    for mol in mols:
        print(json.dumps(mol))

if __name__ == "__main__":
    test_mol_to_dict()
