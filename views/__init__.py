
#flask import
from flask import Flask, session
from flask.ext.session import Session
from flask import Response
from flask import stream_with_context, request, Response, make_response
from flask import Flask
from flask import request
from flask import send_file
from flask import session
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
import sys
import StringIO
import urllib, base64 
import numpy as np
from jinja2_extensions import *
import md5 
from scipy.stats import chisquare
import math 
from Bio import Entrez
from phenotips_python_client import PhenotipsClient
from phenotips_python_client import browser
from bson.json_util import loads
from mongodb import *
# fizz: hpo lookup
import phizz
import itertools
import json
import os
import pymongo
import pysam
import gzip
from parsing import *
import logging
import lookups
import random
import sys
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
import math
import plotly
print plotly.__version__  # version >1.9.4 required
from plotly.graph_objs import Scatter, Layout 
# connect to R session
#import pyRserve 
import numpy
import subprocess
from flask import Flask, render_template, redirect, url_for, request


logging.getLogger().addHandler(logging.StreamHandler())
logging.getLogger().setLevel(logging.INFO)

#ADMINISTRATORS = ( 'n.pontikos@ucl.ac.uk',)
#app = Flask(__name__)
#mail_on_500(app, ADMINISTRATORS)
#Compress(app)
#app.config['COMPRESS_DEBUG'] = True
##cache = SimpleCache(default_timeout=60*60*24)



app = Flask(__name__)
ADMINISTRATORS = ( 'n.pontikos@ucl.ac.uk',)
app = Flask(__name__)
mail_on_500(app, ADMINISTRATORS)
Compress(app)
app.config['COMPRESS_DEBUG'] = True
cache = SimpleCache(default_timeout=60*60*24)

REGION_LIMIT = 1E5
EXON_PADDING = 50
# Load default config and override config from an environment variable
#app.config.from_pyfile('uclex.cfg')
app.config.from_pyfile('uclex-old.cfg')

# Check Configuration section for more details
#SESSION_TYPE = 'null'
SESSION_TYPE = 'mongodb'
#SESSION_USE_SIGNER=True
app.config.from_object(__name__)
sess=Session()
sess.init_app(app)



def check_auth(username, password):
    """
    This function is called to check if a username / password combination is valid.
    Will try to connect to phenotips instance.
    """
    conn=PhenotipsClient()
    response=conn.get_patient(auth='%s:%s' % (username, password,),number=1)
    if response:
        session['password2'] = password
        password=md5.new(password).hexdigest()
        session['user'] = username
        session['password'] = password
        return True
    else: return False
    # check that user name and hash of password exist in database
    db_users=get_db('users')
    # setting a session key for pubmedBatch to save result
    session['password2'] = password
    password=md5.new(password).hexdigest()
    session['user'] = username
    session['password'] = password
    r=db_users.users.find_one({'user':username})
    if r is None:
        return False
    elif md5.new(r['password']).hexdigest() == md5.new(password).hexdigest():
        print('LOGIN', session['user'])
        return True
    else:
        return False


def authenticate():
    """Sends a 401 response that enables basic auth"""
    return Response( 'Could not verify your access level for that URL.\n' 'You have to login with proper credentials', 401, {'WWW-Authenticate': 'Basic realm="Login Required"'})


def requires_auth(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        if session:
          if 'user' in session and 'password2' in session and check_auth(session['user'],session['password2']):
             return f(*args, **kwargs)
          else:
             return redirect('https://uclex.cs.ucl.ac.uk/login')
             #return render_template('login.html', error='Invalid Credentials. Please try again.')
        print 'method', request.method
        error=None
        if request.method == 'POST':
          username=request.form['username']
          password=request.form['password']
          if check_auth(username,password):
             return f(*args, **kwargs)
          else:
             # doesn't redirect
             #return render_template('login.html', error='Invalid Credentials. Please try again.')
             #return login()
             return redirect('https://uclex.cs.ucl.ac.uk/login')
    return decorated


# 
@app.route('/login', methods=['GET','POST'])
def login():
    print request.method
    error = None
    print 'login', request.method
    print request.form
    if request.method == 'POST':
       username=request.form['username']
       password=request.form['password']
       if not check_auth(username,password):
          error = 'Invalid Credentials. Please try again.'
       else:
           return redirect('https://uclex.cs.ucl.ac.uk')
    return render_template('login.html', error=error)

# 
@app.route('/logout')
def logout():
    try:
        print session
        del session['user']
        del session['password']
        del session['password2']
        del session
    except NameError:
        return redirect('https://uclex.cs.ucl.ac.uk/login')
    return render_template('login.html', error="You have been logged out")



@app.route('/')
@requires_auth
def homepage():
    cache_key = 't-homepage'
    #t = cache.get(cache_key)
    #if t: return t
    db=get_db()
    total_variants=db.variants.count()
    print('total_variants',total_variants,)
    total_patients=db.patients.count()
    print('total_patients',total_patients,)
    male_patients=db.patients.find( {'sex':'M'}).count()
    print('male_patients',male_patients,)
    female_patients=db.patients.find( {'sex':'F'}).count()
    print('female_patients',female_patients,)
    unknown_patients=db.patients.find( {'sex':'U'}).count()
    dotfile='static/dot/patients.dot'
    DOT=file(dotfile,'r').read().replace('\n','\\n')
    # replace single quote
    DOT=re.sub("'", '&#39;', DOT)
    #fontsize=7
    # change fontsize to 7
    #DOT=re.sub(r'fontsize="\d+"', 'fontsize="%d"' % fontsize, DOT)
    exac_variants=db.variants.find({'in_exac':True}).count()
    print('exac_variants',exac_variants,)
    pass_variants=db.variants.find({'filter':'PASS'}).count()
    print('pass_variants',pass_variants,)
    pass_exac_variants=db.variants.find({'in_exac':True,'filter':'PASS'}).count()
    print('pass_exac_variants',pass_exac_variants,)
    pass_exac_variants=db.variants.find({'in_exac':True,'filter':'PASS'}).count()
    nonexac_variants=db.variants.find({'in_exac':False}).count()
    pass_nonexac_variants=db.variants.find({'in_exac':False,'filter':'PASS'}).count()
    nonpass_variants=(total_variants-pass_variants)
    nonpass_nonexac_variants=nonexac_variants-pass_nonexac_variants
    #labels = 'PASS', 'non-PASS',
    #sizes =[100*pass_variants/float(total_variants),100*(nonpass_variants)/float(total_variants)]
    #print(sizes)
    #colors = ['yellowgreen', 'red']
    #explode = (0.1, 0)
    #plt.figure(figsize=(5,5))
    #plt.margins(1, 1)
    #plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=True, startangle=90)
     ## Set aspect ratio to be equal so that pie is drawn as a circle.
    #plt.axis('equal')
    #plt.axis('off')
    #plt.show()
    # word cloud
    #from os import path
    #from wordcloud import WordCloud
    #text = 'HPO HPO HPO HPO all day'
    ## Read the whole text.
    ## take relative word frequencies into account, lower max_font_size
    #wordcloud = WordCloud().generate(text)
    #plt.figure()
    #plt.imshow(wordcloud)
    #plt.axis("off")
    #plt.show()
    #imgdata = StringIO.StringIO()
    #plt.savefig(imgdata, format='svg')
    #imgdata.seek(0)  # rewind the data
    #import urllib
    #image=urllib.quote(base64.b64encode(imgdata.buf))
    #image=imgdata.buf
    #image = '<svg' + image.split('<svg')[1]
    t = render_template('homepage.html',
        total_patients=total_patients,
        male_patients=male_patients,
        female_patients=female_patients,
        unknown_patients=unknown_patients,
        DOT=DOT,
        total_variants=total_variants,
        exac_variants=exac_variants,
        pass_variants=pass_variants,
        pass_exac_variants=pass_exac_variants,
        pass_nonexac_variants=pass_nonexac_variants,
        #image=image.decode('utf8'))
        image="")
    #cache.set(cache_key, t)
    return t



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

def get_R_session():
    if not hasattr(g, 'R_session'): g.R_session=pyRserve.connect()
    return g.R_session

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
    client = pymongo.MongoClient(host=app.config['DB_HOST'], port=app.config['DB_PORT'])
    print(client)
    if not dbname: dbname=app.config['DB_NAME']
    print(dbname)
    return client[dbname]


# @app.teardown_appcontext
# def close_db(error):
#     """Closes the database again at the end of the request."""
#     if hasattr(g, 'db_conn'):
#         g.db_conn.close()
def get_gene_page_content(gene_id,hpo_id=None,hpo_string=""):
    db = get_db()
    patients_db=get_db('patients')
    hpo_db=get_db('hpo')
    # scrape exac
    b=browser.Browser('exac.broadinstitute.org')
    p=b.get_page('/gene/%s'%gene_id)
    m = re.compile('window.table_variants\s*=\s*(.*)\s*;')
    exac_table_variants=dict()
    if m: exac_table_variants=loads(m.search(p).group(1))
    #gene = lookups.get_gene(db, gene_id)
    gene=db.genes.find_one({'gene_id':gene_id},{'_id':0})
    gene_name=gene['gene_name']
    if gene is None: abort(404)
    # gene and hpo term ideally
    cache_key = 't-gene-{}'.format(gene_id)
    t = cache.get(cache_key)
    print 'Rendering %sgene: %s' % ('' if t is None else 'cached ', gene_id)
    if t is not None:
        print 'Rendering gene: %s' % gene_id
        #print(phizz.query_gene(ensembl_id=str(gene_id)))
        return t
    # 1. get variants in gene from db
    #variants_in_gene = [v for v in db.variants.find({'variant_id':{'$in':gene['variant_ids']}},{'_id':0})]
    variants_in_gene=lookups.get_variants_in_gene(db,gene_id)
    print('variants_in_gene',len(variants_in_gene))
    # which transcript contains most of the variants
    transcript_counter=Counter([t for v in variants_in_gene for t in v['transcripts'] ])
    print(transcript_counter)
    transcript_with_most_variants=transcript_counter.most_common(1)[0][0]
    print('transcript with most variants',transcript_with_most_variants)
    transcripts_in_gene = lookups.get_transcripts_in_gene(db, gene_id)
    print('transcripts_in_gene',len(transcripts_in_gene))
    # Get some canonical transcript and corresponding info
    transcript_id = gene['canonical_transcript']
    #transcript_id = transcript_with_most_variants
    # if none of the variants are on the canonical transcript use the transcript with the most variants on it
    transcript = lookups.get_transcript(db, transcript_id)
    variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)
    #print('variants_in_transcript',len(variants_in_transcript))
    coverage_stats = lookups.get_coverage_for_transcript(db, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING)
    print('coverage_stats',len(coverage_stats))
    add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)
    constraint_info = lookups.get_constraint_for_transcript(db, transcript_id)
    #print('constraint_info',constraint_info)
    # 2. get patients with hpo term and patients without
    #patients=[p for p in patients_db.patients.find( { 'features': {'$elemMatch':{'id':str(hpo_id)}} } )]
    #patient_ids=[p['external_id'] for p in patients]
    #hpo=phizz.query_hpo([hpo_id])[0]
    # samples
    everyone=frozenset(file('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/headers.txt','r').read().strip().split('\t'))
    #everyone=frozenset([p['external_id'] for p in patients_db.patients.find()])
    #for p in patients_db.patients.find(): if 'external_id' not in p: print(p)
    if hpo_id is None:
        hpo="HP:0000001"
        hpo_name='All'
        cases=frozenset()
        print('num cases',len(cases))
    else:
        #cases
        hpo_name=hpo_db.hpo.find_one({'id':hpo_id})['name'][0]
        print(hpo_name)
        cases=[p['external_id'] for p in lookups.get_hpo_patients(hpo_db,patients_db,hpo_id)]
        cases=frozenset(cases)
        print('num cases',len(cases))
    #controls
    #everyone=frozenset(everyone.tolist())
    controls=everyone-cases
    #everyone=everyone & headers
    #controls=frozenset(everyone) - frozenset(cases)
    print('num controls',len(controls))
    # 3. for each variant tabix,  get counts and chisq
    chrom, pos, ref, alt = variants_in_gene[0]['variant_id'].split('-')
    tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % chrom,)
    region ='%s:%s-%s' % (str(gene['chrom']), str(gene['start']), str(gene['stop']),)
    headers=[h for h in tb.header]
    headers=(headers[len(headers)-1]).strip('#').strip().split('\t')
    records=[dict(zip(headers,r.strip().split('\t'))) for r in tb.fetch(region)]
    print(len(records))
    records=dict([('%s-%s-%s-%s' % (r['CHROM'], r['POS'], r['REF'], r['ALT'],),r,) for r in records])
    for i,_, in enumerate(variants_in_gene):
        v=variants_in_gene[i]
        variant_str=v['variant_id']
        print(variant_str)
        variant_str=str(variant_str).strip().replace('_','-')
        chrom, pos, ref, alt = variant_str.split('-')
        v['pos_coding_noutr']= get_xpos(chrom, pos)
        #region=str('%s:%s-%s'%(chrom, pos, int(pos),))
        #records=tb.fetch(region=region)
        #geno=dict(zip(headers, [r.split('\t') for r in records][0]))
        if variant_str not in records:
            v["-log10pvalue"]=0
            continue
        geno=records[variant_str]
        #samples=[h for h in geno if geno[h].split(':')[0]=='0/1' or geno[h].split(':')[0]=='1/1']
        geno_cases=[geno[ca].split(':')[0] for ca in cases]
        caco=Counter(geno_cases)
        geno_controls=[geno[co].split(':')[0] for co in controls]
        coco=Counter(geno_controls)
        ca_mut=caco.get('1/1',0)+caco.get('0/1',0)
        ca_wt=caco.get('0/0',0)
        co_mut=coco.get('1/1',0)+coco.get('0/1',0)
        co_wt=coco.get('0/0',0)
        v['co_wt']=co_wt
        v['ca_wt']=ca_wt
        v['co_mut']=co_mut
        v['ca_mut']=ca_mut
        if ca_mut==0:
            v["-log10pvalue"]=0
        else:
            counts=numpy.array([[ca_mut,ca_wt],[co_mut,co_wt]])
            print(counts)
            stat=chisquare(counts)
            print(stat)
            v['-log10pvalue']=-math.log10(stat['p.value'])
            #print(chisquare([ca_mut,ca_wt,co_mut,co_wt]))
            #d=csv.DictReader(file('/data/uclex_files/UCLexInfo/uclex-samples.csv','r'),delimiter=',')
            #headers=file('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/headers.txt','r').read().strip().replace('#','').split('\t')
            #get EXAC info
        for j,_, in enumerate(exac_table_variants):
            exac_v=exac_table_variants[j]
            if exac_v['variant_id']!=variant_str: continue
            v['EXAC']=exac_v
            v['HGVS']=exac_v['HGVS']
            v['HGVSp']=exac_v['HGVSp']
            v['HGVSc']=exac_v['HGVSc']
            v['major_consequence']=exac_v['major_consequence']
            v['exac_allele_freq']=exac_v['allele_freq']
            break
        variants_in_gene[i]=v
    # gene hpo
    # dotfile='/slms/UGI/vm_exports/vyp/phenotips/HPO/dot/%s.dot' % gene_name
    # temporary dotfile for test
    dotfile='static/dot/literature_phenotype/%s.dot' % gene_name
    if os.path.isfile(dotfile):
        literature_DOT=file(dotfile,'r').read().replace('\n','\\n')
        # replace single quote
        literature_DOT=re.sub("'", '&#39;', literature_DOT)
        #fontsize=7
        # change fontsize to 7
        #DOT=re.sub(r'fontsize="\d+"', 'fontsize="%d"' % fontsize, DOT)
    else:
        literature_DOT=''
    simreg_DOT=''
    print('all good')
    # this is to set the pos
    for i,_, in enumerate(variants_in_gene):
        variant_id=variants_in_gene[i]['variant_id']
        for j,_, in enumerate(variants_in_transcript):
            variant_id2=variants_in_transcript[j]['variant_id']
            if variant_id==variant_id2:
                variants_in_transcript[j]['-log10pvalue']=variants_in_gene[i]['-log10pvalue']
                break
    print variants_in_transcript
    if not variants_in_transcript: variants_in_transcript=variants_in_gene
    t=render_template( 'gene.html',
            gene=gene,
            transcript=transcript,
            variants_in_gene=variants_in_gene,
            variants_in_transcript=variants_in_transcript,
            transcripts_in_gene=transcripts_in_gene,
            constraint=constraint_info,
            csq_order=csq_order,
            literature_DOT=hpo_string.replace('\n','\\n'),
            simreg_DOT=simreg_DOT,
            hpo_name=hpo_name, coverage_stats=coverage_stats)
    #cache.set(cache_key, t, timeout=1000*60)
    return t


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
    variant['hpo']=[p for p in get_db('patients').patients.find({'external_id':{'$in':samples}},{'_id':0,'features':1,'external_id':1})]
    return(jsonify(result=variant))

@app.route('/variant_json/<variant_str>')
def variant_json(variant_str):
    variant_str=str(variant_str).strip().replace('_','-')
    chrom, pos, ref, alt = variant_str.split('-')
    tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % chrom,)
    #mainset_February2016_chrX_filtered.vcf.gz
    region=str('%s:%s-%s'%(chrom, pos, int(pos),))
    headers=[h for h in tb.header]
    headers=(headers[len(headers)-1]).strip().split('\t')
    records=tb.fetch(region=region)
    records=[r.split('\t') for r in records]
    for r in records:
        geno=dict(zip(headers, r))
        POS=geno['POS']
        REF=geno['REF']
        print 'POS', POS
        print 'REF', REF
        for i, ALT, in enumerate(geno['ALT'].split(',')):
            print 'ALT', ALT
            # insertion
            if ref=='-' and REF+alt==ALT:
                return response(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            # deletion
            # replace leftmost
            elif alt=='-' and ALT==REF.replace(ref,''):
                return reponse(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            # replace rightmost
            elif alt=='-' and ALT==REF[::-1].replace(ref[::-1], "", 1)[::-1]:
                return response(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            # 
            elif alt=='-' and ref==REF and ALT=='*':
                return response(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            elif alt==ALT and ref==REF:
                return response(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            continue


@app.route('/set_variant_causal/<individual>/<variant_str>')
def set_variant_causal(individual, variant_str):
    print individual, variant_str
    db=get_db()
    #get_db().patients.update({'patient_id':individual},{'$addToSet':{'causal_variants':variant_str}})
    var=db.variants.find_one({'variant_id':variant_str})
    gene_id=var['genes'][0]
    gene_name=db.genes.find_one({'gene_id':gene_id})['gene_name_upper']
    print 'GENE_NAME', gene_name
    # update Gene in phenotips
    conn=PhenotipsClient()
    auth='%s:%s' % (session['user'],session['password2'],)
    p=conn.get_patient(eid=individual,auth=auth)
    p['genes']=p.get('genes',[])+[{'gene':gene_name}]
    print conn.update_patient( eid=p['external_id'], auth=auth, patient=p )
    print get_db('patients').patients.update({'external_id':individual},{'$set':p},w=0)
    p=db.patients.find_one({'external_id':individual})
    p['causal_variants']=list(frozenset(p.get('causal_variants',[])+[variant_str]))
    db.patients.update({'external_id':individual},{'$set':{'causal_variants':p['causal_variants']}},w=0)
    if request.referrer:
        referrer=request.referrer
        u = urlparse(referrer)
        referrer='%s://%s' % (u.scheme,u.hostname,)
        if u.port: referrer='%s:%s' % (referrer,u.port,)
    return redirect(referrer+'/individual/'+individual)

@app.route('/unset_variant_causal/<individual>/<variant_str>')
def unset_variant_causal(individual, variant_str):
    print individual, variant_str
    db=get_db()
    p=db.patients.find_one({'external_id':individual})
    if 'causal_variants' in p and not p['causal_variants']: p['causal_variants']=[]
    if variant_str in p.get('causal_variants',[]):
        p['causal_variants']=p['causal_variants'].remove(variant_str)
    db.patients.update({'external_id':individual},{'$set':{'causal_variants':p['causal_variants']}},w=0)
    conn=PhenotipsClient()
    auth='%s:%s' % (session['user'],session['password2'],)
    p2=conn.get_patient(eid=individual,auth=auth)
    p2['genes']=[]
    for var in p['causal_variants']:
        var=db.variants.find_one({'variant_id':var})
        gene_id=var['genes'][0]
        gene_name=db.genes.find_one({'gene_id':gene_id})['gene_name_upper']
        print 'GENE_NAME', gene_name
        p2['genes']=list(frozenset(p2.get('genes',[])+[{'gene':gene_name}]))
    # update Gene in phenotips
    print conn.update_patient( eid=p2['external_id'], auth=auth, patient=p2 )
    print get_db('patients').patients.update({'external_id':individual},{'$set':p2},w=0)
    if request.referrer:
        referrer=request.referrer
        u = urlparse(referrer)
        referrer='%s://%s' % (u.scheme,u.hostname,)
        if u.port: referrer='%s:%s' % (referrer,u.port,)
    return redirect(referrer+'/individual/'+individual)

@app.route('/set_variant_status/<individual>/<variant_str>/<status>')
def set_variant_status(individual, variant_str, status):
    print individual, variant_str, status
    db=get_db()
    #print get_db().patients.update({'patient_id':individual},{'$addToSet':{'variant_status':{variant_str:status}}})
    rare_variants=db.patients.find_one({'external_id':individual},{'rare_variants':1})['rare_variants']
    for rv in rare_variants:
        if rv['variant_id']==variant_str:
            rv['status']=status
    print db.patients.update({'external_id':individual},{'$set':{'rare_variants':rare_variants}})
    return status


import views.uclex
import views.uclex_gene
import views.pubmedbatch



