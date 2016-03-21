"""
Encapsulate interaction with Solr
"""

# TODO: use pysolr where appropriate

import os.path
from decimal import Decimal
from json import dumps, loads, JSONEncoder

import requests


class Solr:
    def __init__(self, url, core, chunk_size=1000):
        self.url = url
        self.core = core
        self.chunk_size = chunk_size
        self.fields = { }

    # Built-in JSON can't serialize Decimals
    # Use like this: json.dumps(o, cls=__JsonDecimalEncoder)
    class JsonDecimalEncoder(JSONEncoder):
        def default(self, o):
            if isinstance(o, Decimal):
                return float(o)
            return super(Solr.JsonDecimalEncoder, self).default(o)

    def postall(self, mols):
        """Post all molecules in a sequence"""
        chunk = []
        self._get_fields()
        new_fields = []
        for mol in mols:
            chunk.append(mol)
            new_fields += self._new_fields(mol)
            if len(chunk) == self.chunk_size:
                self._post_fields(new_fields)
                new_fields = []
                self._post_chunk(chunk)
                chunk = []
        if len(chunk) != 0:
                self._post_chunk(chunk)

    def _post_chunk(self, chunk):
        url = os.path.join(self.url, self.core, 'update') + '?commit=true'
        headers = {'Content-Type' : 'application/json'}
        data = dumps(chunk, cls=Solr.JsonDecimalEncoder)
        print("Posting {0} molecules to Solr at {1}".format(len(chunk), url))
        requests.post(url, headers=headers, data=data)

    def _get_fields(self):
        url = os.path.join(self.url, self.core, 'schema/fields')
        response = requests.get(url, params={'wt': 'json'})
        fields = loads(response.text)['fields']
        self.fields = {f['name']: f for f in fields}
        print("Solr has {0} fields defined.".format(len(self.fields)))

    def _new_fields(self, mol):
        new_fields = []
        for key, val in mol.items():
            if key not in self.fields.keys():
                type = self._guess_type(val)
                if type is None:
                    # Skip this one; hopefully, another record will have a value for this field.
                    continue
                field = {
                    'name': key,
                    'indexed': True,
                    'stored': True,
                    'type': type
                }
                self.fields[key] = field
                new_fields.append(field)
        return new_fields

    def _guess_type(self, o):
        if o is None:
            return None
        if isinstance(o, list):
            subtype = self._guess_type(o[0])
            if subtype != 'text_general':
                return subtype + 's'
            return subtype
        if isinstance(o, bool):
            return 'boolean'
        if isinstance(o, int):
            # Use trie-encoded integers for range filtering
            return 'tlong'
        if isinstance(o, float):
            # Use trie-encoded floats for range filtering
            return 'tdouble'
        if isinstance(o, Decimal):
            if (o - int(o)) == 0:
                return 'tlong'
            return 'tdouble'
        # For short strings, use "string", because they will be used mainly as labels.
        # For long strings, use "text_general", because they will be used mainly for search.
        # '64' is a somewhat arbitrary cut-off between 'short' and 'long'
        if len(o) < 64:
            return 'string'
        return 'text_general'

    def _post_fields(self, new_fields):
        url = os.path.join(self.url, self.core, 'schema/fields') + '?commit=true'
        headers = {'Content-Type' : 'application/json'}
        # As of Solr version 5.4.1:
        # The documentation at https://cwiki.apache.org/confluence/display/solr/Schema+API is wrong.
        # It says we need to send something like:
        #    { "add-field": {"name": ..., "type": ..., ...} }
        # In reality, we need to send:
        #    [ {"name": ..., "type": ..., ...}, {"name": ..., "type": ..., ...}, ...]
        # Apparently, "add-field" is implied.
        data = dumps(new_fields, cls=Solr.JsonDecimalEncoder)
        requests.post(url, headers=headers, data=data)
        print("Posting {0} fields to Solr at {1}".format(len(new_fields), url))
        # Update our local field cache
        self._get_fields()


def test_guess_type():
    solr = Solr("#", 'core')
    assert(solr._guess_type(None) is None)
    assert(solr._guess_type(True) == 'boolean')
    assert(solr._guess_type([True, False]) == 'booleans')
    assert(solr._guess_type(1) == 'tlong')
    assert(solr._guess_type([1, 2, 3]) == 'tlongs')
    assert(solr._guess_type(1.0) == 'tdouble')
    assert(solr._guess_type([1.0, 2.0, 3.0]) == 'tdoubles')
    assert(solr._guess_type(Decimal('1.0')) == 'tlong')
    assert(solr._guess_type(Decimal('1.1')) == 'tdouble')
    assert(solr._guess_type('helios') == 'string')
    assert(solr._guess_type(['helios', 'centauri', 'betelgeuse']) == 'strings')
    assert(solr._guess_type('A ford is a fine kind of truck, whether it be a small, pickup style or a large semi-trailer. Either will serve you reliably for many years to come.')
           == 'text_general')
    assert(solr._guess_type(['A ford is a fine kind of truck, whether it be a small, pickup style or a large semi-trailer. Either will serve you reliably for many years to come.'])
           == 'text_general')

def test_postall():
    mols = [{
                "id": Decimal(i),
                'beauty':True,
                'beauties':[True, False],
                'ingot': 11,
                'ingots': [12, 13, 14],
                'fish': 3.14,
                'fishes': [1.1, 2.2, 3.3],
                'sun': 'helios',
                'suns': ['centauri', 'betelgeuse', 'sirius'],
                'truck': 'A ford is a fine kind of truck, whether it be a small, pickup style or a large semi-trailer. Either will serve you reliably for many years to come.'
            } for i in range(1500)]
    solr = Solr('http://localhost:8983/solr/', 'test')
    solr.postall(mols)

if __name__ == '__main__':
    test_guess_type()
    test_postall()
