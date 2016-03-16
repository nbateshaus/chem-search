import argparse

import Chembl
import Pubchem
import Solr
import Surechembl

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--solr', required=True, help="URL of Solr instance")
    parser.add_argument('--solr_core', default='molecules', help="Name of the Solr core to which to post (default='molecules')")
    parser.add_argument('--chembl', help="DSN connection string for ChEMBLdb")
    parser.add_argument('--pubchem', help="Directory or file glob for PubChem SDF files")
    parser.add_argument('--surechembl', help="Directory or file glob for SureChEMBL SDF files")
    parser.add_argument('--limit', type=int, help="Maximum number of records to process")
    args = parser.parse_args()

    s = Solr.Solr(url=args.solr, core=args.solr_core)
    if args.chembl is not None:
        chembl = Chembl.Chembl(dsn=args.chembl, limit=args.limit)
        s.postall(chembl)

    if args.pubchem is not None:
        pubchem = Pubchem.Pubchem(path=args.pubchem, limit=args.limit)
        s.postall(pubchem)

    if args.surechembl is not None:
        surechembl = Surechembl.Surechembl(path=args.surechembl, limit=args.limit)
        s.postall(surechembl)
