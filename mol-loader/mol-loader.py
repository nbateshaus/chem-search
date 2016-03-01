import Chembl
import Solr

if __name__=='__main__':
    c = Chembl.Chembl(dbname="chembl_20")
    s = Solr.Solr()
    s.postall(c)
