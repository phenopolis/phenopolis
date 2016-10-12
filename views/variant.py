import flask
from views import *
from lookups import *
import rest as annotation
import requests
import primer3
import myvariant
import re
from utils import *
import itertools
import pysam
import csv
#hpo lookup
import phizz
import random
import orm
import vcf
import subprocess
import os


@app.route('/variant/<variant_str>')
@requires_auth
def variant_page(variant_str):
    db = get_db()
    variant=orm.Variant(db=db,variant_id=variant_str)
    # pos, ref, alt = get_minimal_representation(pos, ref, alt)
    #v=load_variant(db,variant_id)
    #xpos = get_xpos(chrom, pos)
    return jsonify(result=variant.__dict__)
    if variant is None:
        variant = {
            'chrom': chrom,
            'pos': pos,
            'xpos': xpos,
            'ref': ref,
            'alt': alt
        }
    consequences = []
    ordered_csqs = []
    # Adds major_consequence
    #base_coverage = lookups.get_coverage_for_bases(db, xpos, xpos + len(ref) - 1)
    base_coverage = []
    #any_covered = any([x['has_coverage'] for x in base_coverage])
    any_covered = any([x['has_coverage'] for x in base_coverage])
    # check the appropriate sqlite db to get the *expected* number of
    # available bams and *actual* number of available bams for this variant
    print 'Rendering variant: %s' % variant_str
    return render_template(
        'variant.html',
        variant=variant,
        base_coverage=base_coverage,
        consequences=consequences,
        any_covered=any_covered,
        ordered_csqs=ordered_csqs,
        metrics=[]
    )


@app.route('/variant_json/<variant_str>')
def variant_json(variant_str):
    return jsonify(result=vcf.vcf_query(variant_str=variant_str))

@app.route('/variant_json_db/<variant_str>')
def variant_json_db(variant_str):
    variant=orm.Variant(variant_id=variant_str,db=get_db())
    return jsonify(result=variant.__dict__)

@app.route('/variant_json_db_new/<variant_str>')
def variant_json_db_new(variant_str):
    variant=orm.Variant(variant_id=variant_str,db=get_db('uclex'))
    return jsonify(result=variant.__dict__)

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


@app.route('/private_variants/<individual>')
def private_variants(individual):
    pv=[]
    cmd="bgt view -s,"+individual+" -s 'name!=\""+individual+"\"' -f 'AC1>0&&AC2==0' -G "+ "/slms/gee/research/vyplab/UCLex/mainset_July2016/bgt/mainset_July2016.bgt"
    print(cmd)
    s=subprocess.check_output([cmd],shell=True)
    for l in s.split('\n'):
        if len(l)<5: continue
        if l.startswith('##'): continue
        if l.startswith('#'):
            headers=l.split('\t')
            continue
        d=dict(zip(headers,l.split('\t')))
        d.update(dict([x.split('=') for x in d['INFO'].split(';')]))
        del d['INFO']
        d['variant_id']='-'.join([d['#CHROM'],d['POS'],d['REF'],d['ALT']])
        pv.append(d)
    return jsonify(result=pv)

@app.route('/rare_variants/<individual>/<AC>')
def rare_variants(individual,AC=10):
    pv=[]
    cmd="bgt view -s,"+individual+" -s 'name!=\""+individual+"\"' -f 'AC1>0&&AC2<%s' "%str(AC)+ "-G /slms/gee/research/vyplab/UCLex/mainset_July2016/bgt/mainset_July2016.bgt" 
    print(cmd)
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True)
    def generate():
        for l in iter(proc.stdout.readline,''):
            l=l.strip()
            print(l)
            if len(l)<5: continue
            if l.startswith('##'): continue
            if l.startswith('#'):
                headers=l.split('\t')
                continue
            d=dict(zip(headers,l.split('\t')))
            d.update(dict([x.split('=') for x in d['INFO'].split(';')]))
            del d['INFO']
            if ',' in d['ALT']: d['ALT']=d['ALT'].split(',')[0]
            d['variant_id']='-'.join([d['#CHROM'],d['POS'],d['REF'],d['ALT']])
            try:
                var=orm.Variant(variant_id=d['variant_id'],db=get_db('uclex'))
            except Exception, e:
                print(e)
                print(d)
                continue
            yield flask.json.dumps(var.__dict__)+'\n'
            #yield l+'\n'
    #return Response(stream_with_context(generate()),mimetype='application/json')
    return Response(stream_with_context(generate()),mimetype='text/plain')

@app.route('/common_private_variants/<individual>/<individual2>')
def common_private_variants(individual,individual2):
    pv=[]
    s=subprocess.check_output(["bgt view -s,"+individual+" -s,"+individual2+" -s 'name!=\""+individual+"\"&&name!=\""+individual2+"\"' -f 'AC1>0&&AC2>0&&AC3==0' -G /slms/gee/research/vyplab/UCLex/mainset_July2016/bgt/mainset_July2016.bgt" ],shell=True)
    #bgt view -s,IRDC_batch6_LON_2055 -s,WebsterURMD_Sample_06G02870 -s 'name!="IRDC_batch6_LON_2055"&&name!="WebsterURMD_Sample_06G02870"' -f 'AC1>0&&AC2>0&&AC3==0' -G mainset_July2016_chr1.bgt
    for l in s.split('\n'):
        if len(l)<5: continue
        if l.startswith('##'): continue
        if l.startswith('#'):
            headers=l.split('\t')
            continue
        d=dict(zip(headers,l.split('\t')))
        d.update(dict([x.split('=') for x in d['INFO'].split(';')]))
        del d['INFO']
        d['variant_id']='-'.join([d['#CHROM'],d['POS'],d['REF'],d['ALT']])
        pv.append(d)
    return jsonify(result=pv)

@app.route('/common_rare_variants/<individual>/<individual2>/<AC>')
def common_rare_variants(individual,individual2,AC=1):
    pv=[]
    s=subprocess.check_output(["bgt view -s,"+individual+" -s,"+individual2+" -s 'name!=\""+individual+"\"&&name!=\""+individual2+"\"' -f 'AC1>0&&AC2>0&&AC3<%s' "%AC+ "-G /slms/gee/research/vyplab/UCLex/mainset_July2016/bgt/mainset_July2016.bgt" ],shell=True)
    #bgt view -s,IRDC_batch6_LON_2055 -s,WebsterURMD_Sample_06G02870 -s 'name!="IRDC_batch6_LON_2055"&&name!="WebsterURMD_Sample_06G02870"' -f 'AC1>0&&AC2>0&&AC3==0' -G mainset_July2016_chr1.bgt
    for l in s.split('\n'):
        if len(l)<5: continue
        if l.startswith('##'): continue
        if l.startswith('#'):
            headers=l.split('\t')
            continue
        d=dict(zip(headers,l.split('\t')))
        d.update(dict([x.split('=') for x in d['INFO'].split(';')]))
        del d['INFO']
        #d['variant_id']='-'.join([d['#CHROM'],d['POS'],d['REF'],d['ALT']])
        #pv.append(d)
        d['variant_id']='-'.join([d['#CHROM'],d['POS'],d['REF'],d['ALT']])
        try:
            var=orm.Variant(variant_id=d['variant_id'],db=get_db('uclex'))
        except Exception, e:
            print(e)
            print(d)
            continue
        pv.append(var.__dict__)
    return jsonify(result=pv)



