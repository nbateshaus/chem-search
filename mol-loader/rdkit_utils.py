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


def rdkit_fps_from_mol(mol):
    """
    Generate a string representation of the fingerprint of mol

    :param mol: RDKit Molecule
    :return: String representation of RDKit fingerprint of mol, number of bits in the fingerprint
    """
    fps = None
    num_bits = None
    if mol is not None:
        try:
            fp = Chem.RDKFingerprint(mol)
            num_bits = fp.GetNumOnBits()
            fps = DataStructs.BitVectToFPSText(fp)
        except ValueError:
            pass
    return fps, num_bits


def rdkit_minhash_signature(mol, width):
    """
    Generate a MinHash signature of the specified width from the given molecule

    :param mol: RDKit Molecule
    :param width: The number of items to include in the hash signature
    :return: List of up to the specified width of integers
    """
    sig = None
    if mol is not None:
        try:
            all_hashes = []
            fp = Chem.RDKFingerprint(mol)
            # For small molecules, this will result in too few hashes for a valid signature.
            # Suggestion: hash groups of features, instead of individual features.
            # See http://infolab.stanford.edu/~ullman/mmds/ch3.pdf
            #     Page 113, Section 3.8.4, "Matching Fingerprints"
            # Challenge: the following molecules have 0 bits set in their fingerprints:
            #   O    - Water
            #   [H+] - Proton
            #   C    - Methane
            #   Au   - Elemental Gold
            for i in fp.GetOnBits():
                hasher = hashlib.sha256()
                hasher.update(i.to_bytes(2, 'little', signed=False))
                digest = hasher.digest()
                # fold the hash down to a signed 64-bit integer, because Java
                val = 0
                for i in range(int(len(digest)/8)):
                    val ^= int.from_bytes(digest[8 * i:8 * (i + 1)], 'little', signed=True)
                all_hashes.append(val)
            all_hashes = list(set(all_hashes))
            all_hashes.sort()
            sig = all_hashes[:width]
        except ValueError:
            pass
    return sig


def minhash_similarity(sig1, sig2):
    return len(set(sig1).intersection(sig2)) / max(len(sig1), len(sig2))


# Tests


def test_minhash():
    mols = [rdkit_mol_from_smiles('CCOC'), rdkit_mol_from_smiles('CCO'), rdkit_mol_from_smiles('COC')]
    sigs = [rdkit_minhash_signature(mol, 10) for mol in mols]
    for i in range(len(sigs)):
        print("sig {0} len {1}".format(i, len(sigs[i])))
    print("0 <-> 1: {0}".format(minhash_similarity(sigs[0], sigs[1])))
    print("0 <-> 2: {0}".format(minhash_similarity(sigs[0], sigs[2])))
    print("1 <-> 2: {0}".format(minhash_similarity(sigs[1], sigs[2])))


if __name__ == '__main__':
    test_minhash()
