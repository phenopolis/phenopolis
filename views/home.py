from views import *
from lookups import *
import requests
import re
from utils import *
import itertools
from config import config
if config.IMPORT_PYSAM_PRIMER3:
    import pysam
import csv
#hpo lookup
import orm
import subprocess



@app.route('/', methods=['GET'])
def homepage():
    cache_key = 't-home'
    return render_template('home.html', title='Phenopolis - Home Page')