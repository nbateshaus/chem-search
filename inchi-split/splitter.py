
import re
from collections import namedtuple

# this is almost a validating expression, it could certainly be simpler by just using [^/]* inside the groups
chargeDef = r"(/q[\-\+0-9;\*mMnNi]*)?"
protonationDef = r"(/p[\-\+0-9,;]*)?"
isotopeDef = r"(/i[\-\+0-9,;HDT]*(?:/h[0-9HDT]+)*)?"
stereoBondDef=r"(/b[\-\+0-9,\?\*;mNnNi]*)?"
stereoTetDef=r"(/t[\-\+0-9,\?;\*mMnNi]*)?"
stereoMDef=r"(/m[\-\+0-9,;\.]*)?"
stereoSDef=r"(/s[\-\+0-9,;]*)?"
inchiLayers=(
r"(InChI=1S?)",
r"(/[a-zA-Z0-9\.]*)", # formula
r"(/c[0-9\(\)\-\,\*;]*)?", # skeleton
r"(/h[0-9,\-\Hh\*\(\);]*)?", # hydrogens
chargeDef, # charge
protonationDef, # protonation
stereoBondDef, # stereo_bond
stereoTetDef, #stereo_tet
stereoMDef, #stereo_m
stereoSDef, #stereo_s
isotopeDef, #isotope
stereoBondDef, #isotope_stereo_bond
stereoTetDef, #isotope_stereo_tet
stereoMDef, #isotope_stereo_m
stereoSDef, #isotope_stereo_s
r"(/f[a-zA-Z0-9\.]*(?:/h[0-9,\-\Hh\*\(\);]*)?)?", # fixed_h
chargeDef, # fixedh_charge
protonationDef, # fixedh_protonation
stereoBondDef, #fixedh_stereo_bond
stereoTetDef, #fixedh_stereo_tet
stereoMDef, #fixedh_stereo_m
stereoSDef, #fixedh_stereo_s
isotopeDef, #fixedh_isotope
stereoBondDef, #fixedh_isotope_stereo_bond
stereoTetDef, #fixedh_isotope_stereo_tet
stereoMDef, #fixedh_isotope_stereo_m
stereoSDef, #fixedh_isotope_stereo_s
r"(/o[\(\)0-9,]*)?", # transposition
r"(/r.*)?", # reconnected_main # <- FIX: we punt on this
)
coreExpr=re.compile(''.join(inchiLayers))
Layers=namedtuple("Layers",['start','formula','skeleton','hydrogens',
                            # pos 4
                            'charge','protonation',
                            # pos 6
                            'stereo_bond','stereo_tet','stereo_m','stereo_s',
                            # pos 10
                            'isotope','isotope_stereo_bond','isotope_stereo_tet','isotope_stereo_m','isotope_stereo_s',
                            # pos 15
                            'fixedh','fixedh_charge','fixedh_protonation',
                            # pos 18
                            'fixedh_stereo_bond','fixedh_stereo_tet','fixedh_stereo_m','fixedh_stereo_s',
                            # pos 22
                            'fixedh_isotope','fixedh_isotope_stereo_bond','fixedh_isotope_stereo_tet','fixedh_isotope_stereo_m','fixedh_isotope_stereo_s',
                            # pos 27
                            'transposition',
                            'reconnected_main'
                            ])
layerGroups = {
'main':tuple(range(4)),
'charge':tuple(range(4,6)),
'stereo':tuple(range(6,10)),
'isotope':tuple(range(10,15)),
'fixedh':tuple(range(15,27)),
}

def formulaGrouping(tpl):
    return (tpl[0],tpl[1],)
def mainGrouping(tpl):
    return (tpl[x] for x in layerGroups['main'])
def chargeGrouping(tpl):
    return (tpl[x] for x in layerGroups['main']+layerGroups['charge'])
def stereoGrouping(tpl):
    return (tpl[x] for x in layerGroups['main']+layerGroups['charge']+layerGroups['stereo'])
def isotopeGrouping(tpl):
    return (tpl[x] for x in layerGroups['main']+layerGroups['charge']+layerGroups['isotope'][0:1])
def isotopestereoGrouping(tpl):
    return (tpl[x] for x in layerGroups['main']+layerGroups['charge']+layerGroups['isotope'])
def stereo_isotopeGrouping(tpl):
    return (tpl[x] for x in layerGroups['main']+layerGroups['charge']+layerGroups['stereo']+layerGroups['isotope'][0:1])
def stereo_isotopestereoGrouping(tpl):
    return (tpl[x] for x in layerGroups['main']+layerGroups['charge']+layerGroups['stereo']+layerGroups['isotope'])



def extractLayers(inchi):
    """

    >>> tpl=extractLayers('InChI=1S/C16H20N4O3/c1-9(21)19-15(18-4)20-13-11-7-10(8-17)5-6-12(11)23-16(2,3)14(13)22/h5-7,13-14,22H,1-4H3,(H2,18,19,20,21)/t13?,14-/m0/s1')
    >>> tpl.start
    'InChI=1S'
    >>> tpl.formula
    'C16H20N4O3'
    >>> tpl.skeleton
    'c1-9(21)19-15(18-4)20-13-11-7-10(8-17)5-6-12(11)23-16(2,3)14(13)22'
    >>> tpl.hydrogens
    'h5-7,13-14,22H,1-4H3,(H2,18,19,20,21)'
    >>> tpl.charge
    ''
    >>> tpl.protonation
    ''
    >>> tpl.stereo_bond
    ''
    >>> tpl.stereo_tet
    't13?,14-'
    >>> tpl.stereo_m
    'm0'
    >>> tpl.stereo_s
    's1'
    >>> tpl.isotope
    ''
    >>> tpl.fixedh
    ''

    Charge layers:
        From [O-]CCCC[NH3+]
    >>> tpl = extractLayers('InChI=1S/C4H10NO/c5-3-1-2-4-6/h1-5H2/q-1/p+1')
    >>> tpl.charge
    'q-1'
    >>> tpl.protonation
    'p+1'

    Stereochemistry:
        From [O-][C@H](Cl)/C=C/C=C(/CC(O)=N)CC(=O)N
    >>> tpl = extractLayers('InChI=1S/C9H12ClN2O3/c10-7(13)3-1-2-6(4-8(11)14)5-9(12)15/h1-3,7H,4-5H2,(H2,11,14)(H2,12,15)/q-1/b3-1+/t7-/m0/s1')
    >>> tpl.stereo_bond
    'b3-1+'
    >>> tpl.stereo_tet
    't7-'
    >>> tpl.stereo_m
    'm0'
    >>> tpl.stereo_s
    's1'

    Isotopes:
       From: [13CH3]O
    >>> tpl = extractLayers('InChI=1S/CH4O/c1-2/h2H,1H3/i1+1')
    >>> tpl.isotope
    'i1+1'
    >>> tpl.isotope_stereo_tet
    ''

    Isotope + stereo
       From: [13CH3]O[C@H](C)O
    >>> tpl = extractLayers('InChI=1S/C3H7ClO/c1-3(4)5-2/h3H,1-2H3/t3-/m1/s1/i2+1')
    >>> tpl.isotope
    'i2+1'
    >>> tpl.stereo_tet
    't3-'
    >>> tpl.isotope_stereo_tet
    ''

    Isotope causes stereo
       From: [13CH3][C@H](C)O
    >>> tpl = extractLayers('InChI=1S/C3H8O/c1-3(2)4/h3-4H,1-2H3/i1+1/t3-/m1/s1')
    >>> tpl.isotope
    'i1+1'
    >>> tpl.stereo_tet
    ''
    >>> tpl.isotope_stereo_tet
    't3-'

    Isotope causes stereo + standard stereo
        From: [13CH3][C@H](C)O[C@H](C)O
    >>> tpl = extractLayers('InChI=1S/C5H12O2/c1-4(2)7-5(3)6/h4-6H,1-3H3/t5-/m1/s1/i1+1/t4-,5-')
    >>> tpl.isotope
    'i1+1'
    >>> tpl.stereo_tet
    't5-'
    >>> tpl.isotope_stereo_tet
    't4-,5-'

    Fixed Hs and Isotopes
        From: O=C([18O])/C=C/C(=[18O])O
    >>> tpl = extractLayers('InChI=1/C4H3O4/c5-3(6)1-2-4(7)8/h1-2H,(H,5,6)/b2-1+/i5+2,7+2/f/h5H/i6+2,7+2')
    >>> tpl.isotope
    'i5+2,7+2'
    >>> tpl.fixedh_isotope
    'i6+2,7+2'

    Fixed Hs causes stereo_bond
        From: F[C@H](Cl)/C=C/C=C(/CC(O)=N)CC(=O)N
    >>> tpl = extractLayers('InChI=1/C9H12ClFN2O2/c10-7(11)3-1-2-6(4-8(12)14)5-9(13)15/h1-3,7H,4-5H2,(H2,12,14)(H2,13,15)/b3-1+/t7-/m0/s1/f/h12,14H,13H2/b3-1+,6-2-,12-8?')
    >>> tpl.fixedh
    'f/h12,14H,13H2'
    >>> tpl.fixedh_stereo_bond
    'b3-1+,6-2-,12-8?'

    Fixed Hs causes stereo
        From: C[C@H](Cl)[C@H](/CC(O)=N)CC(=O)N
    >>> tpl = extractLayers('InChI=1/C7H13ClN2O2/c1-4(8)5(2-6(9)11)3-7(10)12/h4-5H,2-3H2,1H3,(H2,9,11)(H2,10,12)/t4-/m0/s1/f/h9,11H,10H2/t4-,5+')
    >>> tpl.fixedh
    'f/h9,11H,10H2'
    >>> tpl.fixedh_stereo_tet
    't4-,5+'

    Fixed Hs cause a new formula
        From: C[C@H](CCC[C@@H](SCCC(C)(C)O)c1cccc(\C=C\c2ccc3ccc(Cl)cc3n2)c1)C(=O)[O-]  # from ChEMBL
    >>> tpl = extractLayers('InChI=1/C29H34ClNO3S/c1-20(28(32)33)6-4-9-27(35-17-16-29(2,3)34)23-8-5-7-21(18-23)10-14-25-15-12-22-11-13-24(30)19-26(22)31-25/h5,7-8,10-15,18-20,27,34H,4,6,9,16-17H2,1-3H3,(H,32,33)/p-1/b14-10+/t20-,27-/m1/s1/fC29H33ClNO3S/q-1')
    >>> tpl.formula
    'C29H34ClNO3S'
    >>> tpl.fixedh
    'fC29H33ClNO3S'
    >>> tpl.fixedh_charge
    'q-1'


    Disconnected parts + Fixed Hs causes stereo_bond + isotopes cause stereo
        From: [13CH3][C@H](C)O[C@H](C)O.F[C@H](Cl)/C=C/C=C(/CC(O)=N)CC(=O)N
    >>> tpl = extractLayers('InChI=1/C9H12ClFN2O2.C5H12O2/c10-7(11)3-1-2-6(4-8(12)14)5-9(13)15;1-4(2)7-5(3)6/h1-3,7H,4-5H2,(H2,12,14)(H2,13,15);4-6H,1-3H3/b3-1+;/t7-;5-/m01/s1/i;1+1/t;4-,5-/f/h12,14H,13H2;/b3-1+,6-2-,12-8?;')
    >>> tpl.stereo_bond
    'b3-1+;'
    >>> tpl.isotope
    'i;1+1'
    >>> tpl.isotope_stereo_tet
    't;4-,5-'
    >>> tpl.fixedh_stereo_bond
    'b3-1+,6-2-,12-8?;'

    Fixed Hs causes stereo + (FixedHs + isotopes) causes stereo (this is the most dependent example I can think of)
        From: N=C(NC)C(/C(=NC)N)=C/CC/C=C(/C1=NC=C[15NH]1)C1NC=C[15N]=1
    >>> tpl = extractLayers('InChI=1/C16H22N8/c1-19-13(17)11(14(18)20-2)5-3-4-6-12(15-21-7-8-22-15)16-23-9-10-24-16/h5-10H,3-4H2,1-2H3,(H2,17,19)(H2,18,20)(H,21,22)(H,23,24)/i21+1,23+1/f/h17,19,21,23H,18H2/b11-5-,17-13?,20-14?/i21+1,24+1/b11-5-,12-6-,17-13?,20-14?')
    >>> tpl.isotope
    'i21+1,23+1'
    >>> tpl.isotope_stereo_bond
    ''
    >>> tpl.fixedh
    'f/h17,19,21,23H,18H2'
    >>> tpl.fixedh_stereo_bond
    'b11-5-,17-13?,20-14?'
    >>> tpl.fixedh_isotope_stereo_bond
    'b11-5-,12-6-,17-13?,20-14?'

    Transposition:
        From the InChI tech manual Fig A3-3
    >>> tpl = extractLayers('InChI=1/2CH2O2/c2*2-1-3/h2*1H,(H,2,3)/i2+1;2-1/f/h2*2H/i3-1;2+1/o(1,2)')
    >>> tpl.transposition
    'o(1,2)'


    Edge cases:
    >>> tpl=extractLayers('InChI=1S/H2/h1H')
    >>> tpl.start
    'InChI=1S'
    >>> tpl.formula
    'H2'
    >>> tpl.skeleton
    ''
    >>> tpl.hydrogens
    'h1H'
    >>> tpl=extractLayers('InChI=1S/H')
    >>> tpl.start
    'InChI=1S'
    >>> tpl.formula
    'H'
    >>> tpl.skeleton
    ''
    >>> tpl.hydrogens
    ''

    """
    match = coreExpr.match(inchi)
    if not match:
        return None
    gps = list(match.groups())
    res = []
    for e in gps:
        if not e:
            res.append('')
        elif e[0]=='/':
            res.append(e[1:])
        else:
            res.append(e)
    res = Layers(*res)
    return res

if __name__=='__main__':
    import doctest
    doctest.testmod()
