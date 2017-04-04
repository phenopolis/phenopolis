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

@app.route('/register',methods=['POST'])
def register():
    name=request.form.get('name').replace(' ','')
    affiliation=request.form.get('affiliation')
    email=request.form.get('email')
    groups=request.form.getlist('group[]')
    user=orm.User(user_db=get_db(app.config['DB_NAME_USERS']),user=name,groups=groups,email=email,affiliation=affiliation)
    print(user.json())
    print(user.status)
    return jsonify(message=user.status['message']), user.status['http_code']


@app.route('/', methods=['GET'])
def homepage():
    cache_key = 't-home'
    db=get_db()
    patients_db=get_db(app.config['DB_NAME_PATIENTS']) 
    total_variants=db.variants.count()
    total_variants=db.variants.count()
    print('total_variants',total_variants,)
    total_patients=patients_db.patients.count()
    print('total_patients',total_patients,)
    male_patients=patients_db.patients.find( {'sex':'M'}).count()
    print('male_patients',male_patients,)
    female_patients=patients_db.patients.find( {'sex':'F'}).count()
    print('female_patients',female_patients,)
    unknown_patients=patients_db.patients.find( {'sex':'U'}).count()
    if config.LOCAL:
        hpo_json={}
    else:
        hpo_file='uclex_stats/overall_hpo_2016_Aug_2.json'
        hpo_json = json.load(open(hpo_file,'r'))
    exac_variants=0
    print('exac_variants',exac_variants,)
    pass_variants=db.variants.find({'FILTER':'PASS'}).count()
    print('pass_variants',pass_variants,)
    #pass_exac_variants=db.variants.find({'in_exac':True,'filter':'PASS'}).count()
    #pass_exac_variants=db.variants.find({'in_exac':True,'filter':'PASS'}).count()
    pass_exac_variants=0
    print('pass_exac_variants',pass_exac_variants,)
    #pass_exac_variants=db.variants.find({'in_exac':True,'filter':'PASS'}).count()
    pass_exac_variants=0
    #nonexac_variants=db.variants.find({'in_exac':False}).count()
    nonexac_variants=0
    #pass_nonexac_variants=db.variants.find({'in_exac':False,'filter':'PASS'}).count()
    pass_nonexac_variants=0
    nonpass_variants=(total_variants-pass_variants)
    nonpass_nonexac_variants=nonexac_variants-pass_nonexac_variants
    try:
        version_number = subprocess.check_output(['git', 'describe', '--exact-match'])
    except:
        version_number = None
    print('Version number is:-')
    print(version_number)
    labnames= ['black','brogan','elliott','gosgene','hardcastle','humphries','kelsell',
               'lachmann','marks','mead','moosajee','nejentsev','rahman','segal',
               'sisodiya','arvc','syrris','ukirdc','vulliamy','webster']
    username = ''
    if session and 'user' in session:
        username = session['user']
    return render_template('home.html', title='Phenopolis - Home Page',
        total_patients=total_patients,
        male_patients=male_patients,
        female_patients=female_patients,
        unknown_patients=unknown_patients,
        hpo_json=json.dumps(hpo_json),
        total_variants=total_variants,
        exac_variants=exac_variants,
        pass_variants=pass_variants,
        nonpass_variants=nonpass_variants,
        pass_exac_variants=pass_exac_variants,
        pass_nonexac_variants=pass_nonexac_variants,
        version_number=version_number,
        labnames=labnames,
        username=username)

