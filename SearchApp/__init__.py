from flask import Flask

app = Flask(__name__)

from .index import index
from .render import render
