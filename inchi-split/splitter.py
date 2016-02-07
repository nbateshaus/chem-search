
import re
from collections import namedtuple

# this is almost a validating expression, it could certainly be simpler by just using [^/]* inside the groups
inchiLayers=(
r"(InChI=1S?)",
r"(/[a-zA-Z0-9\.]*)", # formula
r"(/c[0-9\(\)\-\,;]*)?", # skeleton
r"(/h[0-9,\-\H\(\);]*)?", # hydrogens
r"(/q[\-\+0-9;]*)?", # charge
r"(/p[\-\+0-9,;]*)?", # protonation
r"(/b[\-\+0-9,\?;]*)?", # stereo_bond
r"(/t[\-\+0-9,\?;]*)?", #stereo_tet  FIX: probably could be tightened up
r"(/m[\-\+0-9,;]*)?", #stereo_m    FIX: probably could be tightened up
r"(/s[\-\+0-9,;]*)?", #stereo_s    FIX: probably could be tightened up
r"(/f/h[0-9,\-\H\(\);]*)?", # fixed_h
r"(/b[\-\+0-9,\?;]*)?", #fixedh_stereo_bond
r"(/t[\-\+0-9,\?;]*)?", #fixedh_stereo_tet  FIX: probably could be tightened up
r"(/m[\-\+0-9,;]*)?", #fixedh_stereo_m    FIX: probably could be tightened up
r"(/s[\-\+0-9,;]*)?", #fixedh_stereo_s    FIX: probably could be tightened up

)
coreExpr=re.compile(''.join(inchiLayers))
Layers=namedtuple("Layers",['start','formula','skeleton','hydrogens','charge','protonation','stereo_bond','stereo_tet','stereo_m','stereo_s',
                            'fixedh','fixedh_stereo_bond','fixedh_stereo_tet','fixedh_stereo_m','fixedh_stereo_s'])
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
    >>> tpl.fixedh
    ''

    Charge layers:
    # from [O-]CCCC[NH3+]
    >>> tpl = extractLayers('InChI=1S/C4H10NO/c5-3-1-2-4-6/h1-5H2/q-1/p+1')
    >>> tpl.charge
    'q-1'
    >>> tpl.protonation
    'p+1'

    Stereochemistry:
    # from [O-][C@H](Cl)/C=C/C=C(/CC(O)=N)CC(=O)N
    >>> tpl = extractLayers('InChI=1S/C9H12ClN2O3/c10-7(13)3-1-2-6(4-8(11)14)5-9(12)15/h1-3,7H,4-5H2,(H2,11,14)(H2,12,15)/q-1/b3-1+/t7-/m0/s1')
    >>> tpl.stereo_bond
    'b3-1+'
    >>> tpl.stereo_tet
    't7-'
    >>> tpl.stereo_m
    'm0'
    >>> tpl.stereo_s
    's1'

    Fixed Hs:
    # From: F[C@H](Cl)/C=C/C=C(/CC(O)=N)CC(=O)N
    >>> tpl = extractLayers('InChI=1/C9H12ClFN2O2/c10-7(11)3-1-2-6(4-8(12)14)5-9(13)15/h1-3,7H,4-5H2,(H2,12,14)(H2,13,15)/b3-1+/t7-/m0/s1/f/h12,14H,13H2/b3-1+,6-2-,12-8?')
    >>> tpl.fixedh
    'f/h12,14H,13H2'
    >>> tpl.fixedh_stereo_bond
    'b3-1+,6-2-,12-8?'

    # From: C[C@H](Cl)[C@H](/CC(O)=N)CC(=O)N
    >>> tpl = extractLayers('InChI=1/C7H13ClN2O2/c1-4(8)5(2-6(9)11)3-7(10)12/h4-5H,2-3H2,1H3,(H2,9,11)(H2,10,12)/t4-/m0/s1/f/h9,11H,10H2/t4-,5+')
    >>> tpl.fixedh
    'f/h9,11H,10H2'
    >>> tpl.fixedh_stereo_tet
    't4-,5+'



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

    return Layers(*res)

if __name__=='__main__':
    import doctest
    doctest.testmod()
