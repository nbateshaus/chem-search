"""
Encapsulate interaction with Solr
"""

# TODO: use pysolr where appropriate

import os.path
from decimal import Decimal
from json import dump, dumps, load, loads, JSONEncoder

import requests


class Solr:
    SCHEMA_FILE = 'schema.json'
    SOLR_SCHEMA_PROPS = {
        "name",
        "type",
        "default",
        "indexed",
        "stored",
        "docValues",
        "sortMissingFirst",
        "sortMissingLast",
        "multiValued",
        "omitNorms",
        "omitTermFreqAndPositions",
        "omitPositions",
        "termVectors",
        "termPositions",
        "termOffsets",
        "termPayloads",
        "required",
        "uniqueKey",
        "useDocValuesAsStored"
    }

    def __init__(self, url, core, chunk_size=1000):
        self.url = url
        self.core = core
        self.chunk_size = chunk_size
        self.fields = { }
        self._load_fields()

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
        new_fields = {}
        for mol in mols:
            chunk.append(mol)
            new_fields.update(self._new_fields_in_mol(mol))
            if len(chunk) == self.chunk_size:
                self._reconcile_schema(new_fields)
                self._post_chunk(chunk)
                chunk = []
        if len(chunk) != 0:
                self._post_chunk(chunk)

    def _load_fields(self):
        if os.access(self.SCHEMA_FILE, os.R_OK):
            with open(self.SCHEMA_FILE, 'r') as f:
                tmp_fields = load(f)
            self.fields = {field['name']: field for field in tmp_fields}

    def _save_fields(self):
        keys = sorted([key for key in self.fields.keys()])
        tmp_fields = [self.fields[key] for key in keys]
        with open(self.SCHEMA_FILE, 'w') as f:
            dump(tmp_fields, f, indent=4, sort_keys=True)

    def _post_chunk(self, chunk):
        url = os.path.join(self.url, self.core, 'update') + '?commit=true'
        headers = {'Content-Type' : 'application/json'}
        data = dumps(chunk, cls=Solr.JsonDecimalEncoder)
        print("Posting {0} molecules to Solr at {1}".format(len(chunk), url))
        requests.post(url, headers=headers, data=data)

    def _reconcile_schema(self, missing_from_file):
        for f in missing_from_file:
            if f not in self.fields:
                self.fields[f] = missing_from_file[f]

        # Find any fields not in our save file, add and save them
        url = os.path.join(self.url, self.core, 'schema/fields')
        response = requests.get(url, params={'wt': 'json'})
        solr_fields = loads(response.text)['fields']
        solr_fields = {f['name']: f for f in solr_fields}
        for f in solr_fields.keys():
            if f not in self.fields:
                missing_from_file[f] = solr_fields[f]
                self.fields[f] = solr_fields[f]
        if missing_from_file:
            print('Saving fields {0}'.format(missing_from_file.keys()))
            self._save_fields()

        # Find any fields not in Solr, add and post them
        missing_from_solr = []
        for f in self.fields.keys():
            if f not in solr_fields:
                missing_from_solr.append(self.fields[f])
        if missing_from_solr:
            self._post_fields(missing_from_solr)

    def _new_fields_in_mol(self, mol):
        new_fields = {}
        deferred = {}
        for key, val in mol.items():
            if key not in self.fields:
                t, m = self._guess_type(val)
                auth = self._guess_authority(key)
                field = {
                    'name': key,
                    'indexed': True,
                    'stored': True,
                    'type': t,
                    'multiValued': m,
                    'list': False,
                    'facet': True,
                    'details': True,
                    'authority': auth
                }
                if t is not None:
                    new_fields[key] = field
                    if key in deferred:
                        deferred.pop(key)
                else:
                    deferred[key] = field
        for field in deferred:
            field['type'] = 'text_general'
            field['multiValued'] = True
            new_fields[field['name']] = field
        return new_fields

    def _guess_type(self, o):
        if o is None:
            return None, None
        if isinstance(o, list):
            subtype, ignored = self._guess_type(o[0])
            return subtype, True
        if isinstance(o, bool):
            return 'boolean', False
        if isinstance(o, int):
            # Use trie-encoded integers for range filtering
            return 'tlong', False
        if isinstance(o, float):
            # Use trie-encoded floats for range filtering
            return 'tdouble', False
        if isinstance(o, Decimal):
            if (o - int(o)) == 0:
                return 'tlong', False
            return 'tdouble', False
        # For short strings, use "string", because they will be used mainly as labels.
        # For long strings, use "text_general", because they will be used mainly for search.
        # '64' is a somewhat arbitrary cut-off between 'short' and 'long'
        if len(o) < 64:
            return 'string', False
        return 'text_general', True  # text_general is always multiValued

    @staticmethod
    def _guess_authority(name):
        parts = name.split('_')
        if parts:
            auth = parts[0]
            if auth == 'CHEMBL':
                auth = 'ChEMBL'
            elif auth == 'PUBCHEM':
                auth = 'PubChem'
            elif auth == 'SURECHEMBL':
                auth = 'SureChEMBL'
            elif auth == 'RDKIT':
                auth = 'RDKit'
            return auth
        else:
            return None

    def _post_fields(self, new_fields):
        # We store more than Solr consumes. Strip out stuff not supported by Solr.
        new_fields = [{k: f[k] for k in f.keys() if k in self.SOLR_SCHEMA_PROPS} for f in new_fields]
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
