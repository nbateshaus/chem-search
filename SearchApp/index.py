import pysolr
from flask import render_template, request

from SearchApp import app

PAGE_SIZE = 10

@app.route('/')
@app.route('/index')
def index():
    args = request.args.to_dict()
    page=1
    if 'page' in args:
        page=int(args['page'])
    start=int((page - 1) * PAGE_SIZE)
    solr = pysolr.Solr('http://localhost:8983/solr/molecules/')
    results = solr.search(q="*:*", start=start, rows=PAGE_SIZE)
    headers = list({header for doc in results.docs for header in doc.keys()})
    headers.sort()
    rows = [[doc[header] if header in doc else None for header in headers] for doc in results.docs]
    pagination = lambda: None
    pagination.page = page
    pagination.pages = int(results.hits / PAGE_SIZE)

    return render_template(
        'index.html',
        title='ChemSearch',
        headers=headers,
        rows=rows,
        pagination=pagination
    )
