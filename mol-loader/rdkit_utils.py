from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
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
    for name, fn in INTERESTING_DESCRIPTORS.items():
        dd[name] = fn(mol)

    inchi = Chem.MolToInchi(mol, options='/SUU')
    inchi_info = InchiInfo.InchiInfo(inchi).get_sp3_stereo()
    (n_stereo, n_undef_stereo, is_meso, dummy) = inchi_info['main']['non-isotopic']
    dd['NumChiralCenters'] = n_stereo
    dd['NumDefinedChiralCenters'] = n_stereo - n_undef_stereo
    dd['NumUndefinedChiralCenters'] = n_undef_stereo
    dd['IsMesoStructure'] = is_meso

    return dd


def rdkit_standardize(mol):
    """
    Generate a canonical representation of a molecule.

    canonicalize( (rdkit.Chem.Mol)mol ) -> rdkit.Chem.Mol

    On error, returns None
    """
    canon = None
    try:
        canon = Chem.RemoveHs(mol)
    except ValueError:
        pass
    return canon


def rdkit_smiles(mol):
    smiles = None
    if mol is not None:
        try:
            smiles = Chem.MolToSmiles(mol, True)
        except RuntimeError:
            pass
    return smiles
