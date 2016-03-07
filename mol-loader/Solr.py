"""
Encapsulate interaction with Solr
"""

import decimal
import json
import requests

class Solr:
    def __init__(self, chunk_size=1000, **kwargs):
        self.chunk_size = chunk_size
        self.kwargs = kwargs

    # Built-in JSON can't serialize Decimals
    # Use like this: json.dumps(o, cls=__JsonDecimalEncoder)
    class JsonDecimalEncoder(json.JSONEncoder):
        def default(self, o):
            if isinstance(o, decimal.Decimal):
                return float(o)
            return super(DecimalEncoder, self).default(o)

    def postall(self, mols):
        '''Post all molecules in a sequence'''
        chunk = []
        for mol in mols:
            chunk.append(mol)
            if len(chunk) == self.chunk_size:
                self.__post_chunk(chunk)
                chunk = []
        if len(chunk) != 0:
                self.__post_chunk(chunk)

    def __post_chunk(self, chunk):
        url = 'http://localhost:8983/solr/chem-search/update?commit=true'
        headers = {'Content-Type' : 'application/json'}
        data = json.dumps(chunk, cls=Solr.JsonDecimalEncoder)
        requests.post(url, headers=headers, data=data)

def test_postall():
    mols = [{"id" : decimal.Decimal(i)} for i in range(1500)]
    solr = Solr()
    solr.postall(mols)

if __name__=='__main__':
    test_postall()
