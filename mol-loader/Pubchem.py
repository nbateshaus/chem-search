"""
Encapsulate reading and formatting of PubChem
"""

import json

from Sdf import Sdf
from rdkit_utils import rdkit_descriptors, rdkit_standardize, rdkit_smiles, rdkit_fps_from_mol


class Pubchem(Sdf):
    def __init__(self, path, limit=None):
        Sdf.__init__(self, path, limit)

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
        Capture all the information from a PubChem molecule into a dictionary

        mol_to_dict( (Pubchem)self, (rdkit.Chem.Mol)mol) -> dict
        """
        props = set(mol.GetPropNames())
        d = {p.lower(): self.cast(mol.GetProp(p)) for p in props}
        d['id'] = 'https://pubchem.ncbi.nlm.nih.gov/compound/' + mol.GetProp('PUBCHEM_COMPOUND_CID')

        smiles = [mol.GetProp(p) for p in self.SMILES_PROPS if p in props]
        smiles.append(rdkit_smiles(mol))  # Original SMILES

        if smiles:
            d['smiles'] = list(set(smiles))

        synonyms = [mol.GetProp(p) for p in self.SYNONYM_PROPS if p in props]
        if synonyms:
            d['synonyms'] = list(set(synonyms))

        if 'pubchem_compound_canonicalized' in d:
            canon = d['pubchem_compound_canonicalized']
            if canon == 1:
                d['pubchem_compound_canonicalized'] = True
            elif canon == 0:
                d['pubchem_compound_canonicalized'] = False
            else:
                del d['pubchem_compound_canonicalized']

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

        return d


def test_mol_to_dict():
    p = Pubchem(path="resources/pubchem-test-data.sdf*")
    for mol in p:
        print(json.dumps(mol))


if __name__ == "__main__":
    test_mol_to_dict()
