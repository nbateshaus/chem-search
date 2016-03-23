from io import BytesIO

from flask import abort, request, send_file
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

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
            mol = Chem.MolFromSmiles(args['smiles'])
            AllChem.Compute2DCoords(mol)
            img=Draw.MolToImage(mol, size=(width, height))
            img_io = BytesIO()
            img.save(img_io, 'PNG')
            img_io.seek(0)
            return send_file(img_io, mimetype='image/png')
        except:
            abort(500)
    abort(400)
