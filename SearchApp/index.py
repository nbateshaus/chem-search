import pysolr
from flask import render_template, request, g
from rdkit.Chem import MolFromSmiles, MolToSmiles


from SearchApp import app

PAGE_SIZE = 10


class Pagination:
    def __init__(self, page, pages, page_size, hits):
        self.page = page
        self.pages = pages
        self.page_size = page_size
        self.hits = hits


def include_field(field):
    return field.startswith('RDKit_') or field == 'id'


@app.route('/')
@app.route('/index')
def index():
    args = request.args.to_dict()
    page = 1
    if 'page' in args:
        page = int(args['page'])
    q='*:*'
    if 'q' in args and args['q'] != '':
        q = args['q']
        # See if q is SMILES
        try:
            mol = MolFromSmiles(q)
            if mol is not None:
                # Use RawQueryParser, so all the special riff-raff in SMILES don't confuzzle Solr
                # https://cwiki.apache.org/confluence/display/solr/Other+Parsers#OtherParsers-RawQueryParser
                q = '{{!raw f=SMILES}}{0}'.format(q)
        except ValueError:
            pass
    start = int((page - 1) * PAGE_SIZE)
    solr = pysolr.Solr('http://localhost:8983/solr/molecules/')
    results = solr.search(q=q, start=start, rows=PAGE_SIZE)
    headers = list(
        {header for doc in results.docs for header in doc.keys() if include_field(header)}
    )
    headers.sort()
    rows = [{header: doc[header] if header in doc else '' for header in headers} for doc in results.docs]
    pagination = Pagination(page, int(results.hits / PAGE_SIZE) + 1, PAGE_SIZE, results.hits)
    g.qtime = results.qtime

    return render_template(
        'index.html',
        title='ChemSearch',
        headers=headers,
        rows=rows,
        pagination=pagination
    )
