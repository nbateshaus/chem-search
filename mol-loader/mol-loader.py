import Chembl
import Pubchem
import Solr
import Surechembl

if __name__=='__main__':
    s = Solr.Solr()
    chembl = Chembl.Chembl(dbname="chembl_20")
    s.postall(chembl)
    pubchem = Pubchem.Pubchem(path="/Users/nik/Data/PubChem/*.sdf.gz")
    s.postall(pubchem)
    surechembl = Surechembl.Surechembl(path="/Users/nik/Data/PubChem/*.sdf.gz")
    s.postall(surechembl)
