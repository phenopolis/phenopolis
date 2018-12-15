
#flask import
from flask import Flask
from flask import session
from flask.ext.session import Session
from flask import Response
from flask import stream_with_context
from flask import request
from flask import make_response
from flask import request
from flask import send_file
from flask import g
from flask import redirect
from flask import url_for
from flask import abort
from flask import render_template
from flask import flash
from flask import jsonify
from flask import send_from_directory
from flask.ext.compress import Compress
from flask.ext.runner import Runner
from flask_errormail import mail_on_500
from flask_debugtoolbar import DebugToolbarExtension 
from flask.ext.cache import Cache
import sys
import StringIO
import urllib, base64 
import numpy as np
from jinja2_extensions import *
import md5 
from scipy.stats import chisquare
import math 
from Bio import Entrez
from bson.json_util import loads
from mongodb import *
# fizz: hpo lookup
import phizz
import itertools
import json
import os
import pymongo
from config import config
if config.IMPORT_PYSAM_PRIMER3:
    import pysam
import gzip
import logging
import lookups
import random
from utils import * 
from collections import defaultdict, Counter
from collections import OrderedDict
from werkzeug.contrib.cache import SimpleCache 
from multiprocessing import Process
import glob
import sqlite3
import traceback
import time 
from functools import wraps 
from werkzeug.exceptions import default_exceptions, HTTPException 
import pandas
import csv
import time
import StringIO 
from urlparse import urlparse
import pickle 
#import pdb 
# handles live plotting if necessary
import plotly
print plotly.__version__  # version >1.9.4 required
from plotly.graph_objs import Scatter, Layout 
# connect to R session
#import pyRserve 
import numpy
import subprocess
import datetime
from Crypto.Cipher import DES
import base64
from binascii import b2a_base64, a2b_base64
from werkzeug.security import generate_password_hash, check_password_hash
from passlib.hash import argon2
import orm
from lookups import *
from config import config
import regex
import requests

logging.getLogger().addHandler(logging.StreamHandler())
logging.getLogger().setLevel(logging.INFO)

# Load default config and override config from an environment variable
if config.LOCAL:
    print 'LOCAL'
    app = Flask(__name__,static_url_path='/static', static_folder='../static', template_folder='../templates')
    app.config.from_pyfile('../local.cfg')
else:
    print 'SERVER'
    app = Flask(__name__, template_folder='../templates')
    app.config.from_pyfile('../phenopolis.cfg')

ADMINISTRATORS = ( 'n.pontikos@ucl.ac.uk',)
mail_on_500(app, ADMINISTRATORS)
Compress(app)
#app.config['COMPRESS_DEBUG'] = True
#cache = SimpleCache(default_timeout=70*60*24)
cache = Cache(app,config={'CACHE_TYPE': 'simple'})

#from flask_httpauth import HTTPBasicAuth
#auth = HTTPBasicAuth()

REGION_LIMIT = 1E5
EXON_PADDING = 50

# Check Configuration section for more details
SESSION_TYPE = 'mongodb'
app.config.from_object(__name__)
sess=Session()
sess.init_app(app)

print app.root_path

# Minifys the HTML when app.config['HTML_COMPRESS'] is True; otherwise prettifies the HTML
if app.config['HTML_COMPRESS']:
    from minify_output import uglify
    render_template = uglify(render_template)
else:
    from minify_output import prettify
    render_template = prettify(render_template)

def check_auth(username, password):
    """
    This function is called to check if a username / password combination is valid.
    """
    q={'statements':[{'statement': "MATCH (u:User {user:'%s'}) RETURN u" % username}]}
    print(q)
    resp=requests.post('http://localhost:57474/db/data/transaction/commit',auth=('neo4j', '1'),json=q)
    if not resp: return False
    r=resp.json()['results'][0]['data'][0]['row'][0]
    print(r)
    session['user']=username
    return argon2.verify(password, r['argon_password'])


def authenticate():
    """Sends a 401 response that enables basic auth"""
    return Response( 'Could not verify your access level for that URL.\n' 'You have to login with proper credentials', 401, {'WWW-Authenticate': 'Basic realm="Login Required"'})


def requires_auth(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        if session:
          if 'user' in session: 
             return f(*args, **kwargs)
        if request.method == 'POST':
          username=request.form['user']
          password=request.form['password']
          if check_auth(username,password):
             return f(*args, **kwargs)
        print 'Not Logged In - Redirect to home to login'
        if config.LOCAL:
           return redirect('/')
        else:
           return redirect('https://uclex.cs.ucl.ac.uk/')
    return decorated


@app.before_request
def make_session_timeout():
    session.permanent = True
    app.permanent_session_lifetime = datetime.timedelta(hours=2)

# 
@app.route('/login', methods=['POST'])
def login():
    username=request.form['name']
    password=request.form['password']
    if not check_auth(username,password):
       print 'Login Failed'
       return jsonify(error='Invalid Credentials. Please try again.'), 401
    else:
        print 'LOGIN SUCCESS'
        return jsonify(success="Authenticated"), 200

# 
@app.route('/logout')
def logout():
    print('DELETE SESSION')
    session.pop('user',None)
    if config.LOCAL:
        return redirect('/')
    else:
        return redirect('https://uclex.cs.ucl.ac.uk/')


# 
@app.route('/change_password', methods=['POST'])
def change_password():
    username = request.form['change_pwd_name']
    password = request.form['current_password']
    new_password_1 = request.form['new_password_1']
    if username == 'demo': 
        return jsonify(error='You do not have permission to change the password for username \'demo\'.'), 401
    elif not check_auth(username,password):
        print 'Change password:- Login Failed'
        return jsonify(error='Username and current password incorrect. Please try again.'), 401
    else:
        print 'LOGIN SUCCESS, CHANGING PASSWORD'
        argon_password = argon2.hash(new_password_1)
        db_users = get_db(app.config['DB_NAME_USERS'])
        #db_users.users.update_one({'user':username},{'$set':{'password':hash}})
        q={'query':'MATCH (u:User {user: $user}) SET u.password=$password','parameters':{'user':username,'password':password}}
        resp=requests.post('http://localhost:57474/db/data/cypher',auth=('neo4j', '1'),json=q)
        print(resp.json())
        #db_users.users.update_one({'user':username},{'$set':{'argon_password':hash}})
        q={'query':'MATCH (u:User {user: $user}) SET u.argon_password=$password','parameters':{'user':username,'password':argon_password}}
        resp=requests.post('http://localhost:57474/db/data/cypher',auth=('neo4j', '1'),json=q)
        msg = 'Password for username \''+username+'\' changed. You are logged in as \''+username+'\'.' 
        return jsonify(success=msg), 200


@app.route('/set/<query>')
def set(query):
    value = query
    session['key'] = value
    return value

@app.route('/get/')
def get():
    return session.get('key', 'not set')


def get_db(dbname=None):
    """
    Opens a new database connection if there is none yet for the
    current application context.
    """
    if dbname is None: dbname=app.config['DB_NAME']
    if not hasattr(g, 'db_conn'):
        g.db_conn=dict()
        g.db_conn[dbname] = connect_db(dbname)
    elif dbname not in g.db_conn:
        g.db_conn[dbname] = connect_db(dbname)
    return g.db_conn[dbname]

def get_hpo_graph():
    """
    """
    if not hasattr(g, 'hpo_graph'):
        from hpo_similarity.ontology import Ontology
        from hpo_similarity.similarity import CalculateSimilarity
        ontology=Ontology(app.config['HPO_OBO'])
        g.hpo_graph=ontology.get_graph()
    return g.hpo_graph


def connect_db(dbname=None):
    """
    Connects to the specific database.
    """
    if dbname=='neo4j':
        from neo4j.v1 import GraphDatabase, basic_auth
        neo4j=GraphDatabase.driver("bolt://localhost:7687", auth=basic_auth("neo4j", "1"))
        return neo4j.session()
    print(app.config['DB_HOST'], app.config['DB_PORT'])
    client = pymongo.MongoClient(host=app.config['DB_HOST'], port=app.config['DB_PORT'])
    print(client)
    if not dbname: dbname=app.config['DB_NAME']
    print(dbname)
    return client[dbname]



def parse_tabix_file_subset(tabix_filenames, subset_i, subset_n, record_parser):
    """
    Returns a generator of parsed record objects (as returned by record_parser) for the i'th out n subset of records
    across all the given tabix_file(s). The records are split by files and contigs within files, with 1/n of all contigs
    from all files being assigned to this the i'th subset.

    Args:
        tabix_filenames: a list of one or more tabix-indexed files. These will be opened using pysam.Tabixfile
        subset_i: zero-based number
        subset_n: total number of subsets
        record_parser: a function that takes a file-like object and returns a generator of parsed records
    """
    start_time = time.time()
    print(tabix_filenames)
    open_tabix_files = [pysam.Tabixfile(tabix_filename) for tabix_filename in tabix_filenames]
    tabix_file_contig_pairs = [(tabix_file, contig) for tabix_file in open_tabix_files for contig in tabix_file.contigs]
    # get every n'th tabix_file/contig pair
    tabix_file_contig_subset = tabix_file_contig_pairs[subset_i : : subset_n]
    short_filenames = ", ".join(map(os.path.basename, tabix_filenames))
    print(short_filenames)
    num_file_contig_pairs = len(tabix_file_contig_subset)
    print(("Loading subset %(subset_i)s of %(subset_n)s total: %(num_file_contig_pairs)s contigs from %(short_filenames)s") % locals())
    counter = 0
    for tabix_file, contig in tabix_file_contig_subset:
        header_iterator = tabix_file.header
        records_iterator = tabix_file.fetch(contig, 0, 10**9, multiple_iterators=True)
        for parsed_record in record_parser(itertools.chain(header_iterator, records_iterator)):
            counter += 1
            yield parsed_record
            if counter % 100000 == 0:
                seconds_elapsed = int(time.time()-start_time)
                print(("Loaded %(counter)s records from subset %(subset_i)s of %(subset_n)s from %(short_filenames)s " "(%(seconds_elapsed)s seconds)") % locals())
    print("Finished loading subset %(subset_i)s from  %(short_filenames)s (%(counter)s records)" % locals())



def create_cache():
    """
    This is essentially a compile step that generates all cached resources.
    Creates files like autocomplete_entries.txt
    Should be run on every redeploy.
    """
    # create autocomplete_entries.txt
    autocomplete_strings = []
    for gene in get_db().genes.find():
        autocomplete_strings.append(gene['gene_name'])
        if 'other_names' in gene:
            autocomplete_strings.extend(gene['other_names'])
    f = open(os.path.join(app.config['UCLEX_FILES_DIRECTORY'],'autocomplete_strings.txt'), 'w')
    for s in sorted(autocomplete_strings):
        f.write(s+'\n')
    f.close()
    # create static gene pages for genes in
    if not os.path.exists(app.config['GENE_CACHE_DIR']): os.makedirs(app.config['GENE_CACHE_DIR'])
    # get list of genes ordered by num_variants
    for gene_id in app.config['GENES_TO_CACHE']:
        try:
            page_content = get_gene_page_content(gene_id)
        except Exception as e:
            print e
            continue
        f = open(os.path.join(app.config['GENE_CACHE_DIR'], '{}.html'.format(gene_id)), 'w')
        f.write(page_content)
        f.close()


def precalculate_metrics():
    db = get_db()
    print 'Reading %s variants...' % db.variants.count()
    metrics = defaultdict(list)
    binned_metrics = defaultdict(list)
    progress = 0
    start_time = time.time()
    for variant in db.variants.find(fields=['quality_metrics', 'site_quality', 'allele_num', 'allele_count']):
        for metric, value in variant['quality_metrics'].iteritems():
            metrics[metric].append(float(value))
        qual = float(variant['site_quality'])
        metrics['site_quality'].append(qual)
        if variant['allele_num'] == 0: continue
        if variant['allele_count'] == 1:
            binned_metrics['singleton'].append(qual)
        elif variant['allele_count'] == 2:
            binned_metrics['doubleton'].append(qual)
        else:
            for af in AF_BUCKETS:
                if float(variant['allele_count'])/variant['allele_num'] < af:
                    binned_metrics[af].append(qual)
                    break
        progress += 1
        if not progress % 100000:
            print 'Read %s variants. Took %s seconds' % (progress, int(time.time() - start_time))
    print 'Done reading variants. Dropping metrics database... '
    db.metrics.drop()
    print 'Dropped metrics database. Calculating metrics...'
    for metric in metrics:
        bin_range = None
        data = map(numpy.log, metrics[metric]) if metric == 'DP' else metrics[metric]
        if metric == 'FS':
            bin_range = (0, 20)
        elif metric == 'VQSLOD':
            bin_range = (-20, 20)
        elif metric == 'InbreedingCoeff':
            bin_range = (0, 1)
        if bin_range is not None:
            data = [x if (x > bin_range[0]) else bin_range[0] for x in data]
            data = [x if (x < bin_range[1]) else bin_range[1] for x in data]
        hist = numpy.histogram(data, bins=40, range=bin_range)
        edges = hist[1]
        # mids = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
        lefts = [edges[i] for i in range(len(edges)-1)]
        db.metrics.insert({
            'metric': metric,
            'mids': lefts,
            'hist': list(hist[0])
        })
    for metric in binned_metrics:
        hist = numpy.histogram(map(numpy.log, binned_metrics[metric]), bins=40)
        edges = hist[1]
        mids = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
        db.metrics.insert({
            'metric': 'binned_%s' % metric,
            'mids': mids,
            'hist': list(hist[0])
        })
    db.metrics.ensure_index('metric')
    print 'Done pre-calculating metrics!'

def response(POS, REF, ALT, index, geno, chrom, pos):
    homozygous_genotype='/'.join([str(index),str(index)])
    heterozygous_genotype='/'.join(['0',str(index)])
    variant=dict()
    variant['pos']=POS
    variant['ref']=REF
    variant['alt']=ALT
    variant['hom_samples']=[h for h in geno if geno[h].split(':')[0]==homozygous_genotype][0:100]
    variant['HOM_COUNT']=len(variant['hom_samples'])
    variant['het_samples']=[h for h in geno if geno[h].split(':')[0]==heterozygous_genotype][0:100]
    variant['HET_COUNT']=len(variant['het_samples'])
    variant['wt_samples']=[h for h in geno if geno[h].split(':')[0]=='0/0'][1:100]
    variant['WT_COUNT']=len([h for h in geno if geno[h].split(':')[0]=='0/0'])
    variant['MISS_COUNT']=len([h for h in geno if geno[h].split(':')[0]=='./.'])
    variant['allele_num']= 2*(variant['HOM_COUNT'] + variant['HET_COUNT']+variant['WT_COUNT'])
    variant['allele_count']=2*variant['HOM_COUNT'] + variant['HET_COUNT']
    #variant['site_quality'] = variant['QUAL']
    #variant['filter'] = variant['FILTER']
    if variant['WT_COUNT']==0:
        variant['allele_freq'] = None
    else:
        variant['allele_freq'] = float(variant['HET_COUNT']+2*variant['HOM_COUNT']) / float(2*variant['WT_COUNT'])
    var2='-'.join([str(chrom),str(pos),variant['ref'],variant['alt']])
    variant['variant_id']=var2
    samples=variant['het_samples']+variant['hom_samples']
    print(samples)
    variant['hpo']=[p for p in get_db(app.config['DB_NAME_PATIENTS']).patients.find({'external_id':{'$in':samples}},{'_id':0,'features':1,'external_id':1})]
    return(jsonify(result=variant))



def get_rathergood_suggestions(query):
    """
    This generates autocomplete suggestions when user
    query is the string that user types
    If it is the prefix for a gene, return list of gene names
    """
    regex = re.compile(re.escape(query), re.IGNORECASE)
    patient_results = [x['external_id'] for x in get_db(app.config['DB_NAME_PATIENTS']).patients.find(
      {'external_id':regex}, {'score': {'$meta': 'textScore'}}
      ).sort([('score', {'$meta': 'textScore'})])
    ]
    gene_results = [x['gene_name'] for x in get_db().genes.find(
      {'gene_name':regex}, {'score': {'$meta': 'textScore'}}
      ).sort([('score', {'$meta': 'textScore'})])
    ]
    hpo_results = [x['name'][0] for x in get_db(app.config['DB_NAME_HPO']).hpo.find(
      {'name':regex}, {'score': {'$meta': 'textScore'}}
      ).sort([('score', {'$meta': 'textScore'})])
    ]
    results = patient_results+gene_results+hpo_results
    results = itertools.islice(results, 0, 20)
    return list(results)

@app.route('/phenotype_suggestions/<hpo>')
def get_phenotype_suggestions(hpo):
    # pattern = '.*?'.join(re.escape(hpo))   # Converts 'kid' to 'k.*?i.*?d'
    regex = re.compile(re.escape(hpo), re.IGNORECASE)  # Compiles a regex.
    suggestions = [x['name'][0] for x in get_db(app.config['DB_NAME_HPO']).hpo.find(
            {['name'][0]:regex}, {"score": {"$meta": "textScore"}} 
      ).sort([('score', {'$meta': 'textScore'})])]
    return json.dumps(suggestions[0:20])

@app.route('/gene_suggestions/<gene>')
def get_gene_suggestions(gene):
    # pattern = '.*?'.join(re.escape(gene))   # Converts 'kid' to 'k.*?i.*?d'
    regex = re.compile(re.escape(gene), re.IGNORECASE)  # Compiles a regex.
    suggestions = [x['gene_name'] for x in get_db().genes.find(
        {'gene_name':regex}, {"score": {"$meta": "textScore"}} 
      ).sort([('score', {'$meta': 'textScore'})]
    )]
    return json.dumps(suggestions[0:20])

def get_rathergood_result(db, query):
    """
    Similar to the above, but this is after a user types enter
    We need to figure out what they meant - could be gene, variant, region
    Return tuple of (datatype, identifier)
    Where datatype is one of 'gene', 'variant', or 'region'
    And identifier is one of:
    - ensembl ID for gene
    - variant ID string for variant (eg. 1-1000-A-T)
    - region ID string for region (eg. 1-1000-2000)
    Follow these steps:
    - if query is an ensembl ID, return it
    - if a gene symbol, return that gene's ensembl ID
    - if an RSID, return that variant's string
    Finally, note that we don't return the whole object here - only it's identifier.
    This could be important for performance later
    """
    query = query.strip()
    print 'Query: %s' % query
    # phenotype
    if query.startswith('HP:'):
        description=phizz.query_hpo([query])
        #description=hpo_db.hpo.find_one({'hpo_id':query})
        return 'hpo', query
    hpo=get_db(app.config['DB_NAME_HPO']).hpo.find_one({'name':query})
    if hpo:
        hpo_id=hpo['id'][0]
        return 'hpo', hpo_id
    if query.startswith('MIM'):
        disease=phizz.query_disease([query])
        return 'mim', query
    # patient
    patient=get_db(app.config['DB_NAME_PATIENTS']).patients.find_one({'external_id':query})
    if patient:
        return 'individual', patient['external_id']
    # Variant
    variant = orm.get_variants_by_rsid(db, query.lower())
    if variant:
        if len(variant) == 1:
            return 'variant', variant[0]['variant_id']
        else:
            return 'dbsnp_variant_set', variant[0]['rsid']
    variant = get_variants_from_dbsnp(db, query.lower())
    if variant:
        return 'variant', variant[0]['variant_id']
    # variant = get_variant(db, )
    # TODO - https://github.com/brettpthomas/exac_browser/issues/14
    gene = get_gene_by_name(db, query)
    if gene:
        return 'gene', gene['gene_id']
    # From here out, all should be uppercase (gene, tx, region, variant_id)
    query = query.upper()
    gene = get_gene_by_name(db, query)
    if gene:
        return 'gene', gene['gene_id']
    # Ensembl formatted queries
    if query.startswith('ENS'):
        # Gene
        gene = get_gene(db, query)
        if gene:
            return 'gene', gene['gene_id']
        # Transcript
        transcript = get_transcript(db, query)
        if transcript:
            return 'transcript', transcript['transcript_id']
    # From here on out, only region queries
    if query.startswith('CHR'):
        query = query.lstrip('CHR')
    # Region
    m = R1.match(query)
    if m:
        if int(m.group(3)) < int(m.group(2)):
            return 'region', 'invalid'
        return 'region', '{}-{}-{}'.format(m.group(1), m.group(2), m.group(3))
    m = R2.match(query)
    if m:
        return 'region', '{}-{}-{}'.format(m.group(1), m.group(2), m.group(2))
    m = R3.match(query)
    if m:
        return 'region', '{}'.format(m.group(1))
    m = R4.match(query)
    if m:
        return 'variant', '{}-{}-{}-{}'.format(m.group(1), m.group(2), m.group(3), m.group(4))
    return 'not_found', query



@app.route('/autocomplete/<query>')
def rathergood_autocomplete(query):
    suggestions = get_rathergood_suggestions(query)
    return Response(json.dumps(suggestions),  mimetype='application/json')


@app.route('/awesome')
@requires_auth
def awesome():
    db = get_db()
    query = str(request.args.get('query'))
    #for n in dir(request): print(n, getattr(request,n))
    #print(request.HTTP_REFERER)
    print(request.referrer)
    if request.referrer:
        referrer=request.referrer
        u = urlparse(referrer)
        referrer='%s://%s' % (u.scheme,u.hostname,)
        if u.port: referrer='%s:%s' % (referrer,u.port,)
    else:
        referrer=''
    #u.netloc
    print(referrer)
    datatype, identifier = get_rathergood_result(db, query)
    print "Searched for %s: %s" % (datatype, identifier)
    if datatype == 'gene':
        return redirect('{}/gene/{}'.format(referrer,identifier))
    elif datatype == 'transcript':
        return redirect('{}/transcript/{}'.format(referrer,identifier))
    elif datatype == 'variant':
        return redirect('{}/variant/{}'.format(referrer,identifier))
    elif datatype == 'region':
        return redirect('{}/region/{}'.format(referrer,identifier))
    elif datatype == 'dbsnp_variant_set':
        return redirect('{}/dbsnp/{}'.format(referrer,identifier))
    elif datatype == 'hpo':
        return redirect('{}/hpo/{}'.format(referrer,identifier))
    elif datatype == 'mim':
        return redirect('{}/mim/{}'.format(referrer,identifier))
    elif datatype == 'individual':
        return redirect('{}/individual/{}'.format(referrer,identifier))
    elif datatype == 'error':
        return redirect('{}/error/{}'.format(referrer,identifier))
    elif datatype == 'not_found':
        return redirect('{}/not_found/{}'.format(referrer,identifier))
    else:
        raise Exception


@app.route('/patient/<patient_str>')
def get_patient(patient_str):
    return patient_str

# AJAX
# Not finished
@app.route('/chisqu/<variant_str>',methods=['GET','POST'])
@requires_auth
def chisq(variant_str):
    if request.method=='POST':
        hpo_patients=request.form['patients'].strip().split(',')
    else:
        hpo_patients=request.args.get('patients').strip().split(',')
    print('hpo_patients',hpo_patients,)
    variant_str=str(variant_str).strip().replace('_','-')
    chrom, pos, ref, alt = variant_str.split('-')
    tb=pysam.TabixFile('chr%s.vcf.gz' % chrom,)
    region=str('%s:%s-%s'%(chrom, pos, int(pos),))
    headers=[h for h in tb.header]
    headers=(headers[len(headers)-1]).strip().split('\t')
    print(region)
    records=tb.fetch(region=region)
    geno=dict(zip(headers, [r.split('\t') for r in records][0]))
    samples=[h for h in geno if geno[h].split(':')[0]=='0/1' or geno[h].split(':')[0]=='1/1']
    res=jsonify(result=hpo_patients)
    return res


def stream_template(template_name, **context):
    app.update_template_context(context)
    t = app.jinja_env.get_template(template_name)
    rv = t.stream(context)
    rv.enable_buffering(5)
    return rv

@app.route('/my-large-page.html')
def render_large_template():
    rows = iter_all_rows()
    return Response(stream_template('the_template.html', rows=rows))



@app.route('/stream')
def streamed_response():
    def generate():
         yield 'Hello '
         yield request.args['name']
         yield '!'
    return Response(stream_with_context(generate()))

def generate_patient_table():
    def get_variants(variant_ids):
        for vid in db.variants.find({'variant_id':{'$in':variant_ids}}):
            yield 

'''
serve the Vincent annotated csv files
'''
@app.route('/download/send_csv', methods=['GET','POST'])
@requires_auth
def download_csv():
    p_id = request.args.get('p_id')
    if not lookup_patient(session['user'],p_id): return 'Sorry you are not permitted to see this patient, please get in touch with us to access this information.'
    folder = request.args.get('folder')
    path = DROPBOX
    csv_file = os.path.join(path,folder, p_id + '.csv')
    filename = folder+'_'+p_id+'.csv'
    if not os.path.isfile(csv_file):
        return 'Oops, file not found!'
    return send_file(csv_file,
                     mimetype='text/csv',
                     attachment_filename=filename,
                     as_attachment=True)

@app.route('/download', methods=['GET','POST'])
@requires_auth
def download():
    p_id = request.args.get('p_id')
    if not lookup_patient(session['user'],p_id): return 'Sorry you are not permitted to see this patient, please get in touch with us to access this information.'
    filetype = request.args.get('filetype')
    index = request.args.get('index')
    path=app.config[str(filetype)]
    print(path)
    if p_id:
        filename=os.path.join(path, p_id)
    else:
        filename=os.path.join(path)
    if filetype=='VCF_DIR':
        if index=='true':
            filename=os.path.join(path, p_id,'all.vcf.gz.tbi')
            attachment_filename=p_id+'.vcf.gz.tbi'
        else:
            filename=os.path.join(path, p_id,'all.vcf.gz')
            attachment_filename=p_id+'.vcf.gz'
    elif filetype=='BAM_DIR':
        if index=='true':
            filename=os.path.join(path, p_id+'_sorted_unique.bam.bai')
            attachment_filename=p_id+'.bam.bai'
        else:
            filename=os.path.join(path, p_id+'_sorted_unique.bam')
            attachment_filename=p_id+'.bam'
    elif filetype=='IRDC_VARIANT_FILES':
        filename=os.path.join(path)
        attachment_filename='IRDC_VARIANTS.zip'
    elif filetype=='IRDC_CNV_FILES':
        filename=os.path.join(path)
        attachment_filename='IRDC_CNV.zip'
    return send_file(filename, mimetype='*/*', attachment_filename=attachment_filename, as_attachment=True)


def encrypt(s):
    obj=DES.new(session['password'][:8], DES.MODE_ECB)
    s=s+(8-(len(s) % 8))*' '
    s=obj.encrypt(s)
    s=base64.urlsafe_b64encode(s)
    return s

def decrypt(s):
    obj=DES.new(session['password'][:8], DES.MODE_ECB)
    s=base64.urlsafe_b64decode(str(s))
    s=obj.decrypt(s)
    s=s.replace(' ','')
    return s

@app.route('/research_pubmed', methods=['POST'])
def research_pubmed():
    # use new search terms to update the individual-pubmedbatch table
    patient_id = request.form['p_id']
    search_term = request.form['OR']
    # update patient pubmed result status as running (1)
    db=get_db()
    db.patients.update({'external_id':patient_id},{'$set': {'pubmedbatch_status': 1}})
    # do the actual update
    exit_status = subprocess.call(['python','offline_analysis/pubmedScore/pubmedScore.py', '-p', patient_id, '--keywords', search_term])
    #exit_status=0
    # reset update status to 0
    db.patients.update({'external_id':patient_id},{'$set': {'pubmedbatch_status': 0}})
    return str(exit_status)

# AJAX
# fetch patients iwth hpo term
@app.route('/fetch_hpo',methods=['GET','POST'])
def fetch_hpo():
    if request.method=='POST':
        hpo_ids=request.form['hpo_ids'].strip().split(',')
    else:
        hpo_ids=request.args.get('hpo_ids').strip().split(',')
    hpo_id=hpo_ids[0]
    print('HPO',hpo_id)
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    patients_db=get_db(app.config['DB_NAME_PATIENTS'])
    hpo_patients=[p['external_id'] for p in lookups.get_hpo_patients(hpo_db,patients_db,hpo_id)]
    print('num patients',len(hpo_patients))
    res=jsonify(result=hpo_patients)
    return res

# AJAX
# fetch variants private to patients
# That is variants which are only seen in these patients and no one else.
@app.route('/fetch_private_variants',methods=['GET','POST'])
def fetch_private_variants():
    if request.method=='POST':
        hpo_patients=request.form['patients'].strip().split(',')
    else:
        hpo_patients=request.args.get('patients').strip().split(',')
    print('hpo_patients',hpo_patients,)
    db=get_db()
    if len(hpo_patients)==1:
        variants=db.variants.find({'PRIVATE_MUT':hpo_patients})
    else:
        #rsession=get_R_session()
        variants=rsession.r.private_variants(hpo_patients)
        #variants=[]
        print('private variants', variants)
        if type(variants) is str:
            variants=[variants]
        else:
            variants=variants.tolist()
    print('num of private variants',len(variants),)
    res=jsonify(result=variants)
    return res

# AJAX
# fetch common variants to patients
# That is variants which are seen in all these patients.
@app.route('/fetch_common_variants',methods=['GET','POST'])
def fetch_common_variants():
    if request.method=='POST':
        hpo_patients=request.form['patients'].strip().split(',')
    else:
        hpo_patients=request.args.get('patients').strip().split(',')
    print('hpo_patients',hpo_patients,)
    #rsession=get_R_session()
    #variants=rsession.r.common_variants(hpo_patients)
    variants=[]
    print('common variants', variants)
    if type(variants) is str:
        variants=[variants]
    else:
        variants=variants.tolist()
    print('num of common variants',len(variants),)
    res=jsonify(result=variants)
    return res


# AJAX
# fetches variant record from db
@app.route('/fetch_variant',methods=['GET','POST'])
def fetch_variant():
    if request.method=='POST':
        variants=request.form['variants'].strip().split(',')
    else:
        variants=request.args.get('variants').strip().split(',')
    db=get_db()
    req_len=len(variants)
    variant_ids=map(lambda x: x.replace('_','-'),variants)
    variants=[v for v in db.variants.find({'variant_id':{'$in':variant_ids}}, projection={'_id': False})]
    ans_len=len(variants)
    print(req_len==ans_len)
    res=jsonify(result=variants)
    return res


# AJAX
# fetches information from db
@app.route('/variant_count',methods=['GET','POST'])
def variant_count():
    if request.method=='POST':
        external_id=request.form['external_id'].strip()
    else:
        external_id=request.args.get('external_id').strip()
    #rsession=get_R_session()
    #res=jsonify(result={'variant_count':rsession.eval('sum(as.logical(variants[["%s"]]))' % external_id) , 'external_id':external_id})
    #return res

# AJAX
# fetches information from db
@app.route('/private_variant_count',methods=['GET','POST'])
def private_variant_count():
    if request.method=='POST':
        external_id=request.form['external_id'].strip()
    else:
        external_id=request.args.get('external_id').strip()
    db=get_db(app.config['DB_NAME_PATIENTS'])
    p=db.patients.find_one({'external_id':external_id})
    if 'PRIVATE_MUT' not in p: private_variant_count=0
    else: private_variant_count=len(p['PRIVATE_MUT'])
    res=jsonify(result={'variant_count': private_variant_count, 'external_id':external_id})
    return res


@app.route('/mim/<mim_id>')
def mim_page(mim_id):
    db=get_db(app.config['DB_NAME_PATIENTS'])
    print(str(mim_id))
    patients=[p for p in db.patients.find( { 'features': {'$elemMatch':{'id':str(hpo_id)}} } )]
    patient_ids=[p['external_id'] for p in patients]
    print(phizz.query_disease([hpo_id]))
    print(len([v['VARIANT_ID'] for v in db.variants.find({'HET' : { '$in': patient_ids }})]))
    print(len([v['VARIANT_ID'] for v in db.variants.find({'HOM' : { '$in': patient_ids }})]))
    return render_template('test.html')

@app.route('/patient/<patient_id>')
def patient_page(patient_id):
    db=get_db()
    patients=[p for p in db.patients.find({'external_id': str(patient_id)})]
    print(patients)
    return None

@app.route('/Exomiser/<path:path>')
@requires_auth
def exomiser_page(path):
    #is this user authorized to see this patient?
    return send_from_directory('Exomiser', path)

@app.route('/example/')
@requires_auth
def example():
    return send_from_directory('templates', 'temp-plot.html')



@app.route('/region/<region_id>')
def region_page(region_id):
    db = get_db()
    try:
        region = region_id.split('-')
        cache_key = 't-region-{}'.format(region_id)
        t = cache.get(cache_key)
        print 'Rendering %sregion: %s' % ('' if t is None else 'cached ', region_id)
        if t is None:
            chrom = region[0]
            start = None
            stop = None
            if len(region) == 3:
                chrom, start, stop = region
                start = int(start)
                stop = int(stop)
            if start is None or stop - start > REGION_LIMIT or stop < start:
                return render_template(
                    'region.html',
                    genes_in_region=None,
                    variants_in_region=None,
                    chrom=chrom,
                    start=start,
                    stop=stop,
                    coverage=None,
                    csq_order=csq_order,
                )
            if start == stop:
                start -= 20
                stop += 20
            genes_in_region = lookups.get_genes_in_region(db, chrom, start, stop)
            variants_in_region = lookups.get_variants_in_region(db, chrom, start, stop)
            xstart = get_xpos(chrom, start)
            xstop = get_xpos(chrom, stop)
            coverage_array = lookups.get_coverage_for_bases(db, xstart, xstop)
            t = render_template(
                'region.html',
                genes_in_region=genes_in_region,
                variants_in_region=variants_in_region,
                chrom=chrom,
                start=start,
                stop=stop,
                coverage=coverage_array,
                csq_order=csq_order,
            )
            cache.set(cache_key, t)
        return t
    except Exception, e:
        print 'Failed on region:', region_id, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/dbsnp/<rsid>')
def dbsnp_page(rsid):
    db = get_db()
    try:
        variants = lookups.get_variants_by_rsid(db, rsid)
        chrom = None
        start = None
        stop = None
        print 'Rendering rsid: %s' % rsid
        return render_template(
            'region.html',
            rsid=rsid,
            variants_in_region=variants,
            chrom=chrom,
            start=start,
            stop=stop,
            coverage=None,
            genes_in_region=None,
            csq_order=csq_order,
        )
    except Exception, e:
        print 'Failed on rsid:', rsid, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/not_found/<query>')
def not_found_page(query):
    return render_template(
        'not_found.html',
        query=query
    )


@app.route('/error/<query>')
@app.errorhandler(404)
def error_page(query):
    return render_template(
        'error.html',
        query=query
    )



@app.route('/about')
def about_page():
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
    return render_template('about.html',total_patients=total_patients)


@app.route('/participants')
def participants_page():
    return render_template('about.html')


@app.route('/terms')
def terms_page():
    return render_template('terms.html')


@app.route('/contact')
def contact_page():
    return render_template('contact.html')


@app.route('/faq')
def faq_page():
    patients_db=get_db(app.config['DB_NAME_PATIENTS']) 
    total_patients=patients_db.patients.count()
    return render_template('faq.html',total_patients=total_patients)

@app.route('/samples')
def samples_page():
    samples=pandas.read_csv('HPO/hpo.txt')
    return render_template('samples.html',samples=samples.to_html(escape=False))


@app.route('/text')
def text_page():
    db = get_db()
    query = request.args.get('text')
    datatype, identifier = get_awesomebar_result(db, query)
    if datatype in ['gene', 'transcript']:
        gene = lookups.get_gene(db, identifier)
        link = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr%(chrom)s%%3A%(start)s-%(stop)s" % gene
        output = '''Searched for %s. Found %s.
%s; Canonical: %s.
%s''' % (query, identifier, gene['full_gene_name'], gene['canonical_transcript'], link)
        output += '' if 'omim_accession' not in gene else '''
In OMIM: %(omim_description)s
http://omim.org/entry/%(omim_accession)s''' % gene
        return output
    elif datatype == 'error' or datatype == 'not_found':
        return "Gene/transcript %s not found" % query
    else:
        return "Search types other than gene transcript not yet supported"


@app.route('/read_viz/<path:path>')
def read_viz_files(path):
    full_path = os.path.abspath(os.path.join(app.config["READ_VIZ_DIR"], path))
    # security check - only files under READ_VIZ_DIR should be accsessible
    if not full_path.startswith(app.config["READ_VIZ_DIR"]):
        return "Invalid path: %s" % path
    logging.info("path: " + full_path)
    # handle igv.js Range header which it uses to request a subset of a .bam
    range_header = request.headers.get('Range', None)
    if not range_header:
        return send_from_directory(app.config["READ_VIZ_DIR"], path)
    m = re.search('(\d+)-(\d*)', range_header)
    if not m:
        error_msg = "ERROR: unexpected range header syntax: %s" % range_header
        logging.error(error_msg)
        return error_msg
    size = os.path.getsize(full_path)
    offset = int(m.group(1))
    length = int(m.group(2) or size) - offset
    data = None
    with open(full_path, 'rb') as f:
        f.seek(offset)
        data = f.read(length)
    rv = Response(data, 206, mimetype="application/octet-stream", direct_passthrough=True)
    rv.headers.add('Content-Range', 'bytes {0}-{1}/{2}'.format(offset, offset + length - 1, size))
    logging.info("GET range request: %s-%s %s" % (m.group(1), m.group(2), full_path))
    return rv


@app.after_request
def apply_caching(response):
    response.headers['Cache-Control'] = 'no-cache'
    response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
    # prevent click-jacking vulnerability identified by BITs
    response.headers["X-Frame-Options"] = "SAMEORIGIN"
    return response

### all the mongodb reading/writing code


def load_db():
    """
    Load the database
    """
    # Initialize database
    # Don't need to explicitly create tables with mongo, just indices
    confirm = raw_input('This will drop the database and reload. Are you sure you want to continue? [no] ')
    if not confirm.startswith('y'):
        print('Exiting...')
        sys.exit(1)
    all_procs = []
    for load_function in [load_variants_file, load_dbsnp_file, load_base_coverage, load_gene_models, load_constraint_information]:
        procs = load_function()
        all_procs.extend(procs)
        print("Started %s processes to run %s" % (len(procs), load_function.__name__))
    [p.join() for p in all_procs]
    print('Done! Loading MNPs...')
    load_mnps()
    print('Done! Creating cache...')
    #create_cache()
    print('Done!')


def load_base_coverage():
    """ """
    def load_coverage(coverage_files, i, n, db):
        coverage_generator = parse_tabix_file_subset(coverage_files, i, n, get_base_coverage_from_file)
        try:
            db.base_coverage.insert(coverage_generator, w=0)
        except pymongo.errors.InvalidOperation, e:
            print(e)
            # handle error when coverage_generator is empty
            pass  
    db = get_db()
    db.base_coverage.drop()
    print("Dropped db.base_coverage")
    # load coverage first; variant info will depend on coverage
    db.base_coverage.ensure_index('xpos')
    procs = []
    coverage_files = app.config['BASE_COVERAGE_FILES']
    num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    random.shuffle(app.config['BASE_COVERAGE_FILES'])
    for i in range(num_procs):
        p = Process(target=load_coverage, args=(coverage_files, i, num_procs, db))
        p.start()
        procs.append(p)
    return procs
    #print 'Done loading coverage. Took %s seconds' % int(time.time() - start_time)


def load_variants_file():
    def load_variants(sites_file, i, n, db):
        for f in sites_file:
            print(f)
            variants_generator = parse_tabix_file_subset([f], i, n, get_variants_from_sites_vcf)
            try:
                db.variants.insert(variants_generator, w=0)
            except pymongo.errors.InvalidOperation:
                pass  # handle error when variant_generator is empty
    db = get_db('exac')
    db.variants.drop()
    print("Dropped db.variants")
    # grab variants from sites VCF
    db.variants.ensure_index('xpos')
    db.variants.ensure_index('xstart')
    db.variants.ensure_index('xstop')
    db.variants.ensure_index('rsid')
    db.variants.ensure_index('genes')
    db.variants.ensure_index('transcripts')
    db.variants.ensure_index('variant_id')
    #sites_vcfs = app.config['SITES_VCFS']
    sites_vcfs=['ExAC.r0.3.1.sites.vep.vcf.gz']
    print(sites_vcfs)
    #if len(sites_vcfs) > 1: raise Exception("More than one sites vcf file found: %s" % sites_vcfs)
    procs = []
    num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    #pdb.set_trace()
    for i in range(num_procs):
        p = Process(target=load_variants, args=(sites_vcfs, i, num_procs, db))
        p.start()
        procs.append(p)
    return procs

    #print 'Done loading variants. Took %s seconds' % int(time.time() - start_time)


def load_constraint_information():
    db = get_db()
    db.constraint.drop()
    print 'Dropped db.constraint.'
    start_time = time.time()
    with gzip.open(app.config['CONSTRAINT_FILE']) as constraint_file:
        for transcript in get_constraint_information(constraint_file):
            db.constraint.insert(transcript, w=0)
    db.constraint.ensure_index('transcript')
    print 'Done loading constraint info. Took %s seconds' % int(time.time() - start_time)


def load_mnps():
    db = get_db()
    start_time = time.time()
    db.variants.ensure_index('has_mnp')
    print 'Done indexing.'
    while db.variants.find_and_modify({'has_mnp' : True}, {'$unset': {'has_mnp': '', 'mnps': ''}}): pass
    print 'Deleted MNP data.'
    with gzip.open(app.config['MNP_FILE']) as mnp_file:
        for mnp in get_mnp_data(mnp_file):
            variant = lookups.get_raw_variant(db, mnp['xpos'], mnp['ref'], mnp['alt'], True)
            db.variants.find_and_modify({'_id': variant['_id']}, {'$set': {'has_mnp': True}, '$push': {'mnps': mnp}}, w=0)
    db.variants.ensure_index('has_mnp')
    print 'Done loading MNP info. Took %s seconds' % int(time.time() - start_time)


@app.route('/load_gene_models/')
def load_gene_models():
    db = get_db()
    db.genes.drop()
    db.transcripts.drop()
    db.exons.drop()
    print 'Dropped db.genes, db.transcripts, and db.exons.'
    start_time = time.time()
    canonical_transcripts = {}
    with gzip.open(app.config['CANONICAL_TRANSCRIPT_FILE']) as canonical_transcript_file:
        for gene, transcript in get_canonical_transcripts(canonical_transcript_file):
            canonical_transcripts[gene] = transcript
    omim_annotations = {}
    with gzip.open(app.config['OMIM_FILE']) as omim_file:
        for fields in get_omim_associations(omim_file):
            if fields is None:
                continue
            gene, transcript, accession, description = fields
            omim_annotations[gene] = (accession, description)
    dbnsfp_info = {}
    with gzip.open(app.config['DBNSFP_FILE']) as dbnsfp_file:
        for dbnsfp_gene in get_dbnsfp_info(dbnsfp_file):
            other_names = [other_name.upper() for other_name in dbnsfp_gene['gene_other_names']]
            dbnsfp_info[dbnsfp_gene['ensembl_gene']] = (dbnsfp_gene['gene_full_name'], other_names)
    print 'Done loading metadata. Took %s seconds' % int(time.time() - start_time)
    # grab genes from GTF
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        for gene in get_genes_from_gencode_gtf(gtf_file):
            gene_id = gene['gene_id']
            if gene_id in canonical_transcripts:
                gene['canonical_transcript'] = canonical_transcripts[gene_id]
            if gene_id in omim_annotations:
                gene['omim_accession'] = omim_annotations[gene_id][0]
                gene['omim_description'] = omim_annotations[gene_id][1]
            if gene_id in dbnsfp_info:
                gene['full_gene_name'] = dbnsfp_info[gene_id][0]
                gene['other_names'] = dbnsfp_info[gene_id][1]
            db.genes.insert(gene, w=0)
    print 'Done loading genes. Took %s seconds' % int(time.time() - start_time)
    start_time = time.time()
    db.genes.ensure_index('gene_id')
    db.genes.ensure_index('gene_name_upper')
    db.genes.ensure_index('gene_name')
    db.genes.ensure_index('other_names')
    db.genes.ensure_index('xstart')
    db.genes.ensure_index('xstop')
    print 'Done indexing gene table. Took %s seconds' % int(time.time() - start_time)
    # and now transcripts
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        db.transcripts.insert((transcript for transcript in get_transcripts_from_gencode_gtf(gtf_file)), w=0)
    print 'Done loading transcripts. Took %s seconds' % int(time.time() - start_time)
    start_time = time.time()
    db.transcripts.ensure_index('transcript_id')
    db.transcripts.ensure_index('gene_id')
    print 'Done indexing transcript table. Took %s seconds' % int(time.time() - start_time)
    # Building up gene definitions
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        db.exons.insert((exon for exon in get_exons_from_gencode_gtf(gtf_file)), w=0)
    print 'Done loading exons. Took %s seconds' % int(time.time() - start_time)
    start_time = time.time()
    db.exons.ensure_index('exon_id')
    db.exons.ensure_index('transcript_id')
    db.exons.ensure_index('gene_id')
    print 'Done indexing exon table. Took %s seconds' % int(time.time() - start_time)
    return []


def load_dbsnp_file():
    db = get_db()
    def load_dbsnp(dbsnp_file, i, n, db):
        if os.path.isfile(dbsnp_file + ".tbi"):
            dbsnp_record_generator = parse_tabix_file_subset([dbsnp_file], i, n, get_snp_from_dbsnp_file)
            try:
                db.dbsnp.insert(dbsnp_record_generator, w=0)
            except pymongo.errors.InvalidOperation:
                pass  # handle error when coverage_generator is empty
        else:
            with gzip.open(dbsnp_file) as f:
                db.dbsnp.insert((snp for snp in get_snp_from_dbsnp_file(f)), w=0)
    db.dbsnp.drop()
    db.dbsnp.ensure_index('rsid')
    db.dbsnp.ensure_index('xpos')
    start_time = time.time()
    dbsnp_file = app.config['DBSNP_FILE']
    print "Loading dbsnp from %s" % dbsnp_file
    if os.path.isfile(dbsnp_file + ".tbi"): num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    else:
        # see if non-tabixed .gz version exists
        if os.path.isfile(dbsnp_file):
            print(("WARNING: %(dbsnp_file)s.tbi index file not found. Will use single thread to load dbsnp."
                "To create a tabix-indexed dbsnp file based on UCSC dbsnp, do: \n"
                "   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp141.txt.gz \n"
                "   gzcat snp141.txt.gz | cut -f 1-5 | bgzip -c > snp141.txt.bgz \n"
                "   tabix -0 -s 2 -b 3 -e 4 snp141.txt.bgz") % locals())
            num_procs = 1
        else:
            raise Exception("dbsnp file %s(dbsnp_file)s not found." % locals())
    procs = []
    for i in range(num_procs):
        p = Process(target=load_dbsnp, args=(dbsnp_file, i, num_procs, db))
        p.start()
        procs.append(p)
    return procs
    #print 'Done loading dbSNP. Took %s seconds' % int(time.time() - start_time)
    #start_time = time.time()
    #db.dbsnp.ensure_index('rsid')
    #print 'Done indexing dbSNP table. Took %s seconds' % int(time.time() - start_time)


"""
Get the most recent common ancestor between two sets of hpo terms.
"""
def mrc_hpo():
    hpo_graph=get_hpo_graph()
    db=get_db()
    for var in db.variants.find():
        hpo_anc=[]
        for eid in list(set(var['HET']+var['HOM'])):
            patient=db.patients.find_one({'external_id':eid})
            if not patient: continue
            if 'features' not in patient: continue
            for f in patient['features']:
                fid=f['id']
                if not fid.startswith('HP'): continue
                hpo_anc.append(set(hpo_graph.get_ancestors(fid)))
        if not hpo_anc: continue
        if 'SYMBOL' not in var: continue
        var['ALL_HPO']=list(set(set.union(*hpo_anc)))
        var['SHARED_HPO']=list(set.intersection(*hpo_anc))
        print(var['VARIANT_ID'],var['SYMBOL'],len(var['HET']+var['HOM']),var['SHARED_HPO'],var['ALL_HPO'])
        db.variants.update({'VARIANT_ID':var['VARIANT_ID']},var,upsert=True)


#progressbar
'''
{
    'random_p_id':{
        'total':456,
        'count':123,
        'status':['running','done']
    },
    ...
}
'''
PROGRESS_BAR = {}

'''
initiate a progress instance
arg: total length of genes
return: progress_id
'''

def init_progress_bar(id,length):
    # check id
    if id in PROGRESS_BAR:
        if PROGRESS_BAR[id]['status'] != 'done':
            return 'the id already exists in PROGRESS_BAR'

    # initialise progress_bar
    PROGRESS_BAR[id] = {
            'total': length,
            'count':0,
            'message': '',
            'status':'running'
    }
    return 0

'''
update progress
arg: {
    id: id, 
    message: message,
    step: 1
    }
default step 1
'''

def update_progress_bar(obj):
    # check if id in PROGRESS_BAR
    if not obj['id'] in PROGRESS_BAR:
        return 'ID does not exist in PROGRESS_BAR'

    # update progress
    if not 'step' in obj:
        obj['step'] = 1
    PROGRESS_BAR[obj['id']]['count'] += obj['step']

    PROGRESS_BAR[obj['id']]['message'] = obj['message']
    # done?
    if PROGRESS_BAR[obj['id']]['count'] == PROGRESS_BAR[obj['id']]['total']:
        PROGRESS_BAR[obj['id']]['status'] = 'done'

'''
kill a progress
'''

def kill_progress_bar(key):
    if key in PROGRESS_BAR:
        del PROGRESS_BAR[key]


'''
to check if an iterable is empty
'''
def peek(iterable):
    try:
        first = next(iterable)
    except RuntimeError:
        return None
    except StopIteration:
        return None
    return first, itertools.chain([first], iterable)


'''
find the freaking PID, Title or Abstract no matter what!
'''
def find_item(obj, key):
    if key in obj:
        return obj[key]
    if isinstance(obj, dict):
        for k in obj:
            if isinstance(obj[k], dict):
                item = find_item(obj[k], key)
                if item is not None:
                    return item
            elif isinstance(obj[k], list):
                for i in obj[k]:
                    if isinstance(i, str):
                        continue
                    item = find_item(i, key)
                    if item is not None:
                        return item
    elif isinstance(obj, list):
        for k in obj:
            if isinstance(k, dict):
                item = find_item(k, key)
                if item is not None:
                    return item
            elif isinstance(k, list):
                for i in k:
                    if isinstance(i, str):
                        continue
                    item = find_item(i, key)
                    if item is not None:
                        return item


"""
for pubmedBatch
check title and abstract is truely relevant. Assign to both this gene and each ref
"""
def scrutinise(obj):
    print obj['smashed_all']
    if obj['lag']:
        obj['lag'] = obj['lag']/3600/24 # convert it to days
            # need to update
        search_results = Entrez.read(Entrez.esearch(db='pubmed', term=obj['smashed_all'], reldate=obj['lag'], datetype='pdat', usehistory='y'))
    else:
        # just search
        search_results = Entrez.read(Entrez.esearch(db='pubmed',retmax=50, term=obj['smashed_all'], usehistory='y'))
    # now done the search. let's get results
    count = int(search_results["Count"])
    print count
    results = {'results':[], 'total_score':0}
    # get search content
    attempt = 1
    while attempt <= 10:
        try:
            handle = Entrez.efetch("pubmed",
                                   restart=0,
                                   retmax=50,
                                   retmode="xml",
                                   webenv=search_results['WebEnv'],
                                   query_key=search_results['QueryKey']
                                   )
            break
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print('Received error from server %s' % err)
            else:
                print('Something is wrong while efetch..')
            print('Attempt %i of 10' % attempt)
            attempt += 1
            time.sleep(5)
    record = Entrez.parse(handle)
    if peek(record):
        # got something. let's do some calculation
        for r in record:
            # calculate score
            score = 0
            pid = str(find_item(r, 'PMID'))
            abstract_list = find_item(r, 'AbstractText')
            # parse abstract
            abstract = ''
            if abstract_list:
                for a in abstract_list:
                    if hasattr(a, 'attributes') and 'Label' in a.attributes:
                        abstract = abstract + '<b>' + a.attributes['Label'] + ': </b>'
                        abstract = abstract + a + '<br/>'
                    else:
                        abstract = abstract + a
        
            title = find_item(r, 'ArticleTitle')
            if title:
                score = score + len(obj['reg'].findall(title))
            if abstract:
                score = score + len(obj['reg'].findall(abstract))
        
            # add result to genes[gene_name]
            if score:
                results['results'].append({
                    'id': pid,
                    'title': title,
                    'abstract': abstract,
                    'score': score
                })
                results['total_score'] = results['total_score'] + score
    results['results'] = sorted(results['results'], key=lambda k: k['score'], reverse=True)
    return results

def get_pred_score(obj):
    # for the batch_pubmed route.
    # calculate the pred score
    # [D/A].each = 10, [P].each = 5, [C].each = 6, [T/B/N].each = -1. If there is a splicing/insertion/deletion event, the score is set as 1000. Not given is set as 0
    # ref: https://github.com/plagnollab/DNASeq_pipeline/blob/master/GATK_v2/filtering.md
    pred = 0
    if ('Func' in obj and re.search('splic', obj['Func'])) or ('ExonicFunc' in obj and re.search(r'stop|frame|del|insert', obj['ExonicFunc'])):
        pred = 1000;
    else:
        for k in obj:
            if re.search('Pred', k):
                if obj[k] == 'D' or obj[k] == 'A':
                    pred = pred + 10
                elif obj[k] == 'P':
                    pred = pred + 5
                elif obj[k] == 'C':
                    pred = pred + 6
                elif obj[k] == 'T' or obj[k] == 'B' or obj[k] == 'N':
                    pred = pred - 1
                else:
                    pass
    return pred;


@app.route('/plot/<gene>')
def plot(gene):
    #db = get_db()
    #var=db.variants.find_one({'VARIANT_ID':'3_8775295_C_T'})
    d=csv.DictReader(file('CARDIO/assoc_3.csv','r'),delimiter=',')
    x=[i for i, r, in enumerate(d)]
    d=csv.DictReader(file('CARDIO/assoc_3.csv','r'),delimiter=',')
    y=[-math.log10(float(r['HCM.chisq.p'])) for r in d]
    print(x)
    print(y)
    d=csv.DictReader(file('CARDIO/assoc_3.csv','r'),delimiter=',')
    #layout = dict( yaxis = dict( type = 'log', tickvals = [ 1.5, 2.53, 5.99999 ]), xaxis = dict( ticktext = [ "green eggs", "& ham", "H2O", "Gorgonzola" ], tickvals = [ 0, 1, 2, 3, 4, 5 ]))
    labels=[r['VARIANT_ID'] for r in d]
    layout = Layout( xaxis = dict( ticktext=labels, tickvals=x ), title="p-value plot" )
    #Layout( title="p-value plot")
    plotly.offline.plot({
        "data": [
                Scatter(
                    x=x,
                    y=y
                    )
                ],
        "layout": layout
        }, filename='genes/%s-pvalues.html' % (gene,), auto_open=False)
    return send_from_directory('genes', '%s-pvalues.html' % gene,)

""" JINJA2 filer """
def highlight(text, list, myclass):
    # wrap list element in text (case insensitive) with <span>
    # note that gene description has to be split by ','
    #  with class to do highlighting
    for l in list:
        # remove (.*), escape +?.*
        l = re.sub(r'\(.*\)', '', l)
        l = re.sub(r'\+','\\+',l)
        l = re.sub(r'\?','\\?',l)
        l = re.sub(r'\.','\\.',l)
        l = re.sub(r'\*','\\*',l)
        l = re.sub(r'\[.*\]','',l)
        l = re.sub(r'\\', '\\\\',l)
        words = l.split(',')
        for w in words:
            # wrap w with brackets to be a catch group
            text = re.sub(r'(\b%s\b)' % w, r'<span class="%s">\1</span>' % myclass, text, flags=re.I)
    return text
jinja2.filters.FILTERS['highlight'] = highlight

def highlight2(text, kw, myclass):
    # wrap list element in text (case insensitive) with <span>
    # note that gene description has to be split by ','
    #  with class to do highlighting
    # remove (.*), escape +?.*
    for w in kw:
        # wrap w with brackets to be a catch group
        text = re.sub(r'\b(%s)\b'%w, r'<span class="%s">\1</span>' % myclass, text, flags=re.I)
    return text
jinja2.filters.FILTERS['highlight2'] = highlight2



import views.my_patients
import views.gene
import views.variant
import views.individual
import views.igv
import views.hpo
import views.search
import views.home
import views.exomiser
# work in progress, comment out if not needed
import views.pheno4j



