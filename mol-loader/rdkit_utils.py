import hashlib

from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolKey import InchiInfo

INTERESTING_DESCRIPTORS = dict(
    ExactMolWt=Descriptors.ExactMolWt,
    FractionCSP3=Descriptors.FractionCSP3,
    HeavyAtomCount=Descriptors.HeavyAtomCount,
    LabuteASA=Descriptors.LabuteASA,
    MolLogP=Descriptors.MolLogP,
    MolWt=Descriptors.MolWt,
    NHOHCount=Descriptors.NHOHCount,
    NOCount=Descriptors.NOCount,
    NumAliphaticCarbocycles=Descriptors.NumAliphaticCarbocycles,
    NumAliphaticHeterocycles=Descriptors.NumAliphaticHeterocycles,
    NumAliphaticRings=Descriptors.NumAliphaticRings,
    NumAromaticCarbocycles=Descriptors.NumAromaticCarbocycles,
    NumAromaticHeterocycles=Descriptors.NumAromaticHeterocycles,
    NumAromaticRings=Descriptors.NumAromaticRings,
    NumHAcceptors=Descriptors.NumHAcceptors,
    NumHDonors=Descriptors.NumHDonors,
    NumRotatableBonds=Descriptors.NumRotatableBonds,
    NumSaturatedCarbocycles=Descriptors.NumSaturatedCarbocycles,
    NumSaturatedHeterocycles=Descriptors.NumSaturatedHeterocycles,
    NumSaturatedRings=Descriptors.NumSaturatedRings,
    RingCount=Descriptors.RingCount,
    TPSA=Descriptors.TPSA,
    NumAmideBonds=rdMolDescriptors.CalcNumAmideBonds,
    NumBridgeheadAtoms=rdMolDescriptors.CalcNumBridgeheadAtoms,
    NumSpiroAtom=rdMolDescriptors.CalcNumSpiroAtoms
)


def rdkit_descriptors(mol):
    """
    Given a molecule, return a dict of lots of interesting descriptors.

    :param mol: the molecule to process
    """
    dd = {}
    if mol is not None:
        try:
            dd = {k: fn(mol) for k, fn in INTERESTING_DESCRIPTORS.items()}
            inchi = Chem.MolToInchi(mol, options='/SUU')
            inchi_info = InchiInfo.InchiInfo(inchi).get_sp3_stereo()
            (n_stereo, n_undef_stereo, is_meso, dummy) = inchi_info['main']['non-isotopic']
            dd['NumChiralCenters'] = n_stereo
            dd['NumDefinedChiralCenters'] = n_stereo - n_undef_stereo
            dd['NumUndefinedChiralCenters'] = n_undef_stereo
            dd['IsMesoStructure'] = is_meso
        except ValueError:
            pass

    return dd


def rdkit_standardize(mol):
    """
    Generate a canonical representation of a molecule.

    rdkit_standardize( (rdkit.Chem.Mol)mol ) -> rdkit.Chem.Mol

    On error, returns None
    """
    std = None
    if mol is not None:
        try:
            std = Chem.RemoveHs(mol)
        except ValueError:
            pass
    return std


def rdkit_smiles(mol):
    smiles = None
    if mol is not None:
        try:
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        except (RuntimeError, ValueError):
            pass
    return smiles


def rdkit_mol_from_smiles(smiles):
    mol = None
    if smiles is not None:
        try:
            mol = Chem.MolFromSmiles(smiles)
        except ValueError:
            pass
    return mol


def rdkit_morgan_fps_from_mol(mol):
    """
    Generate a string representation of a Morgan fingerprint with radius 2 from mol

    :param mol: RDKit Molecule
    :return: String representation of RDKit fingerprint of mol, number of bits in the fingerprint
    """
    fps = None
    num_bits = None
    if mol is not None:
        try:
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            num_bits = fp.GetNumOnBits()
            fps = DataStructs.BitVectToFPSText(fp)
        except ValueError:
            pass
    return fps, num_bits
