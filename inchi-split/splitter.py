
import re
from collections import namedtuple

# this is almost a validating expression, it could certainly be simpler by just using [^/]* inside the groups
coreExpr=re.compile(r"(InChI=1S?)(/[a-z,A-Z,0-9,\.]*)(/c[0-9,\(,\),\-,\,]*)?(/h[0-9,\,\-,\H,\(,\)]*)?")

Core=namedtuple("Core",['start','formula','skeleton','hydrogens'])
def extractCore(inchi):
    """

    >>> tpl=extractCore('InChI=1S/C16H20N4O3/c1-9(21)19-15(18-4)20-13-11-7-10(8-17)5-6-12(11)23-16(2,3)14(13)22/h5-7,13-14,22H,1-4H3,(H2,18,19,20,21)/t13?,14-/m0/s1')
    >>> tpl.start
    'InChI=1S'
    >>> tpl.formula
    'C16H20N4O3'
    >>> tpl.skeleton
    'c1-9(21)19-15(18-4)20-13-11-7-10(8-17)5-6-12(11)23-16(2,3)14(13)22'
    >>> tpl.hydrogens
    'h5-7,13-14,22H,1-4H3,(H2,18,19,20,21)'

    Edge cases:
    >>> tpl=extractCore('InChI=1S/H2/h1H')
    >>> tpl.start
    'InChI=1S'
    >>> tpl.formula
    'H2'
    >>> tpl.skeleton
    ''
    >>> tpl.hydrogens
    'h1H'
    >>> tpl=extractCore('InChI=1S/H')
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

    return Core(*res)

if __name__=='__main__':
    import doctest
    doctest.testmod()
