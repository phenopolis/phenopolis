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


def get_individuals():
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
        if '_id' in p: del p['_id']
        return p
    users_db=get_db(app.config['DB_NAME_USERS'])
    patients_db=get_db(app.config['DB_NAME_PATIENTS'])
    user=users_db.users.find_one({'user':session['user']})
    if 'individuals' in user:
        individuals=user['individuals']
        individuals=[{k:ind.get(k,'') for k in ['external_id','sex','specificity','features','solved','genes','rare_homozygous_variants_count','rare_compound_hets_count','rare_variants_count','total_variant_count']} for ind in individuals]
        return individuals
    eids=user['external_ids']
    patients=patients_db.patients.find({'external_id':{'$in':eids}},{'_id':False})
    individuals=[f(p) for p in patients if 'external_id' in p]
    users_db.users.update_one({'user':session['user']},{'$set':{'individuals':individuals}})
    return individuals


@app.route('/individuals_json')
@requires_auth
def individuals_json():
    individuals=get_individuals()
    return(jsonify(result=individuals))


# shows each individual, 
# all_individuals
@app.route('/individuals')
@requires_auth
def individuals_page():
    individuals=get_individuals()
    return render_template('individuals_page.html',individuals=individuals)


