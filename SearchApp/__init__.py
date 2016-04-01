import json

from flask import Flask

app = Flask(__name__)


def _read_schema():
    with open('schema.json', 'r') as f:
        tmp_fields = json.load(f)
    return {field['name']: field for field in tmp_fields}

app.config['schema'] = _read_schema()

from .index import index
from .render import render
