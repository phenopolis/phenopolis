from views import *
from lookups import *
import requests
import re
from utils import *
import itertools
from flask import request
import orm


@app.route('/exomiser_prioritise/',methods=['GET'])
def exomiser_prioritise():
    #phenotypes=HP:0001156,HP:0001363,HP:0011304,HP:0010055
    #prioritiser=hiphive
    #genes=341640,2263,4920,3909,10743
    #prioritiser-params=human,mouse,fish
    print(request.args)
    r=requests.get('http://localhost:8085/exomiser/api/prioritise/',params=request.args)
    #?phenotypes=HP:0001156,HP:0001363,HP:0011304,HP:0010055&prioritiser=hiphive&genes=341640,2263,4920,3909,10743&prioritiser-params=human,mouse,fish')
    return jsonify(result=r.json())

