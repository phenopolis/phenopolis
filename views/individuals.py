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

@app.route('/individuals_json')
@requires_auth
def individuals_json():
    page=int(request.args.get('page',0))
    number=int(request.args.get('number',200))
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    def f(p):
        print p['external_id']
        p['features']=[f for f in p.get('features',[]) if f['observed']=='yes']
        if 'solved' in p:
            if 'gene' in p['solved']:
                p['solved']=[p['solved']['gene']]
            else:
                p['solved']=[]
        else: p['solved']=[]
        if 'genes' in p: p['genes']=[x['gene'] for x in p['genes'] if 'gene' in x]
        else: p['genes']=[]
        p['genes']=list(frozenset(p['genes']+p['solved']))
        p2=get_db().patients.find_one({'external_id':p['external_id']},{'rare_homozygous_variants_count':1,'rare_compound_hets_count':1, 'rare_variants_count':1,'total_variant_count':1})
        if not p2: return p
        p['rare_homozygous_variants_count']=p2.get('rare_homozygous_variants_count','')
        p['rare_compound_hets_count']=p2.get('rare_compound_hets_count','')
        p['rare_variants_count']=p2.get('rare_variants_count','')
        p['total_variant_count']=p2.get('total_variant_count','')
        #p['all_variants_count']=get_db().patients.find_one({'external_id':p['external_id']},{'_id':0,'all_variants_count':1})['all_variants_count']
        #db.cache.find_one({"key" : "%s_blindness,macula,macular,retina,retinal,retinitis,stargardt_" % })
        return p
    #conn=PhenotipsClient()
    patients_db=get_db(app.config['DB_NAME_PATIENTS'])
    #all_patients=conn.get_patient(session=session).get('patientSummaries',[]) 
    all_patients=[p for p in patients_db.patients.find()]
    all_eids=[p['external_id'] for p in patients_db.patients.find() if 'external_id' in p]
    #all_eids=[p['eid'] for p in all_patients if p['eid']]
    total=len(all_eids)
    print('TOTAL NUMBER OF PATIENTS',total)
    #patients=conn.get_patient(session=session,start=page*number,number=number).get('patientSummaries',[])
    patients=all_patients
    eids=[p['external_id'] for p in patients if 'external_id' in p]
    print(eids)
    patients=patients_db.patients.find({'external_id':{'$in':eids}},projection={'_id':False})
    #patients=get_db(app.config['DB_NAME_PATIENTS']).patients.find({'external_id':re.compile('^IRDC')},{'pubmedBatch':0})
    individuals=[f(p) for p in patients if 'external_id' in p]
    # family_history":{"consanguinity":true}
    #if session['user']=='demo': for ind in individuals: ind['external_id']=encrypt(ind['external_id'])
    return(jsonify(result=individuals))



# shows each individual, 
# all_individuals
@app.route('/individuals')
@requires_auth
def individuals_page():
    page=int(request.args.get('page',0))
    number=int(request.args.get('number',200))
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    def f(p):
        print p['external_id']
        p['features']=[f for f in p.get('features',[]) if f['observed']=='yes']
        if 'solved' in p:
            if 'gene' in p['solved']:
                p['solved']=[p['solved']['gene']]
            else:
                p['solved']=[]
        else: p['solved']=[]
        if 'genes' in p: p['genes']=[x['gene'] for x in p['genes'] if 'gene' in x]
        else: p['genes']=[]
        p['genes']=list(frozenset(p['genes']+p['solved']))
        p2=get_db().patients.find_one({'external_id':p['external_id']},{'rare_homozygous_variants_count':1,'rare_compound_hets_count':1, 'rare_variants_count':1,'total_variant_count':1})
        if not p2: return p
        p['rare_homozygous_variants_count']=p2.get('rare_homozygous_variants_count','')
        p['rare_compound_hets_count']=p2.get('rare_compound_hets_count','')
        p['rare_variants_count']=p2.get('rare_variants_count','')
        p['total_variant_count']=p2.get('total_variant_count','')
        #p['all_variants_count']=get_db().patients.find_one({'external_id':p['external_id']},{'_id':0,'all_variants_count':1})['all_variants_count']
        #db.cache.find_one({"key" : "%s_blindness,macula,macular,retina,retinal,retinitis,stargardt_" % })
        return p
    users_db=get_db(app.config['DB_NAME_USERS'])
    patients_db=get_db(app.config['DB_NAME_PATIENTS'])
    user=users_db.users.find_one({'user':session['user']})
    if 'individuals' in user:
        individuals=user['individuals']
        return render_template('individuals_page.html',individuals=individuals,page=page,number=number,total=len(individuals))
    eids=user['external_ids']
    patients=patients_db.patients.find({'external_id':{'$in':eids}})
    individuals=[f(p) for p in patients if 'external_id' in p]
    users_db.users.update_one({'user':session['user']},{'$set':{'individuals':individuals}})
    return render_template('individuals_page.html',individuals=individuals,page=page,number=number,total=len(individuals))



# shows each individual, 
# all_individuals
@app.route('/individuals_csv')
@requires_auth
def individuals_csv():
    page=int(request.args.get('page',0))
    number=int(request.args.get('number',200))
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    def f(p):
        print p['external_id']
        p['features']=[f for f in p.get('features',[]) if f['observed']=='yes']
        if 'solved' in p:
            if 'gene' in p['solved']:
                p['solved']=[p['solved']['gene']]
            else:
                p['solved']=[]
        else: p['solved']=[]
        if 'genes' in p: p['genes']=[x['gene'] for x in p['genes'] if 'gene' in x]
        else: p['genes']=[]
        p['genes']=list(frozenset(p['genes']+p['solved']))
        p2=get_db().patients.find_one({'external_id':p['external_id']},{'rare_homozygous_variants_count':1,'rare_compound_hets_count':1, 'rare_variants_count':1,'total_variant_count':1})
        if not p2: return p
        p['rare_homozygous_variants_count']=p2.get('rare_homozygous_variants_count','')
        p['rare_compound_hets_count']=p2.get('rare_compound_hets_count','')
        p['rare_variants_count']=p2.get('rare_variants_count','')
        p['total_variant_count']=p2.get('total_variant_count','')
        #p['all_variants_count']=get_db().patients.find_one({'external_id':p['external_id']},{'_id':0,'all_variants_count':1})['all_variants_count']
        #db.cache.find_one({"key" : "%s_blindness,macula,macular,retina,retinal,retinitis,stargardt_" % })
        return p
    conn=PhenotipsClient()
    all_patients=conn.get_patient(session=session).get('patientSummaries',[]) 
    all_eids=[p['eid'] for p in all_patients if p['eid']]
    total=len(all_eids)
    print('TOTAL NUMBER OF PATIENTS',total)
    patients=conn.get_patient(session=session,start=page*number,number=number).get('patientSummaries',[])
    eids=[p['eid'] for p in patients if p['eid']]
    print(eids)
    patients=get_db(app.config['DB_NAME_PATIENTS']).patients.find({'external_id':{'$in':eids}})
    #patients=get_db(app.config['DB_NAME_PATIENTS']).patients.find({'external_id':re.compile('^IRDC')},{'pubmedBatch':0})
    individuals=[f(p) for p in patients if 'external_id' in p]
    # family_history":{"consanguinity":true}
    #if session['user']=='demo': for ind in individuals: ind['external_id']=encrypt(ind['external_id'])
    #return render_template('individuals_page.html',individuals=individuals,page=page,number=number,total=total)
    return '\n'.join([','.join([ind['external_id'],ind['total_variant_count'],ind['rare_variants_count']]) for ind in individuals])



