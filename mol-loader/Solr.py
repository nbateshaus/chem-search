"""
Encapsulate interaction with Solr
"""

import decimal
import json
import requests
import os.path

class Solr:
    def __init__(self, url, core, chunk_size=1000):
        self.url = url
        self.core = core
        self.chunk_size = chunk_size

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
        url = os.path.join(self.url, self.core, 'update') + '?commit=true'
        headers = {'Content-Type' : 'application/json'}
        data = json.dumps(chunk, cls=Solr.JsonDecimalEncoder)
        print(data)
        print("Posting {0} molecules to Solr at {1}".format(len(chunk), url))
        requests.post(url, headers=headers, data=data)

def test_postall():
    mols = [{"id" : decimal.Decimal(i)} for i in range(1500)]
    solr = Solr()
    solr.postall(mols)

if __name__=='__main__':
    test_postall()
