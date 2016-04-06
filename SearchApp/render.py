from flask import abort, request, make_response
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from SearchApp import app


@app.route('/render')
def render():
    args = request.args.to_dict()
    if 'smiles' in args:
        try:
            if 'width' in args:
                width = int(args['width'])
            else:
                width = 300
            if 'height' in args:
                height = int(args['height'])
            else:
                height = 300
            if 'kekulize' in args:
                kekulize = bool(args['kekulize'])
            else:
                kekulize = True
            mol = Chem.MolFromSmiles(args['smiles'])
            mc = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=kekulize)
            drawer = rdMolDraw2D.MolDraw2DSVG(width,height)
            drawer.DrawMolecule(mc)
            drawer.FinishDrawing()
            response = make_response(drawer.GetDrawingText())
            response.headers['Content-Type'] = 'image/svg+xml'
            return response
        except:
                abort(500)
    abort(400)
