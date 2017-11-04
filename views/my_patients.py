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


def individuals_update(external_ids):
    patients_db=get_db(app.config['DB_NAME_PATIENTS'])
    users_db=get_db(app.config['DB_NAME_USERS'])
    def f(eid):
        p=patients_db.patients.find_one({'external_id':eid},{'_id':False})
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
    new_individuals=[f(eid) for eid in external_ids]
    old_individuals=users_db.users.find_one({'user':session['user']}).get('individuals',[])
    old_individuals=[ind for ind in old_individuals if ind['external_id'] not in external_ids]
    individuals=new_individuals+old_individuals
    users_db.users.update_one({'user':session['user']},{'$set':{'individuals':individuals}})
    return individuals

def merge_dicts(*dict_args): # TODO LMTW use common function
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

def old_get_individuals(user):
    s="""
    MATCH (u:User {user:'%s'})--(p:Person)-[:PersonToObservedTerm]->(t:Term),
    (p)-[:CandidateGene]-(g:Gene)
    RETURN p.personId as individual,
    p.gender as gender,
    collect(DISTINCT t) as phenotypes,
    p.score as phenotypeScore,
    size((p)<-[:HomVariantToPerson]-()) as hom_count,
    size((p)<-[:HetVariantToPerson]-()) as het_count,
    collect(DISTINCT g.gene_name) as genes;
    """ % user
    print(s)
    uri='http://'+app.config['NEO4J_HOST']+':'+str(app.config['NEO4J_PORT'])+'/db/data/cypher'
    data=requests.post(uri,auth=('neo4j',app.config['NEO4J_PWD']),json={'query':s})
    temp = data.json() # TODO LMTW remove
    print('######################### data.json() ###############')
    print(temp)


    return data.json()

def get_individuals(user):
    s="""
    MATCH (u:User {user:'%s'})--(p:Person)-[:PersonToObservedTerm]->(t:Term),
    (p)-[:CandidateGene]-(g:Gene)
    RETURN p.personId as individual,
    p.gender as gender,
    collect(DISTINCT t) as phenotypes,
    p.score as phenotypeScore,
    size((p)<-[:HomVariantToPerson]-()) as hom_count,
    size((p)<-[:HetVariantToPerson]-()) as het_count,
    collect(DISTINCT g.gene_name) as genes;
    """ % user

    db_session = neo4j_driver.session()
    result=db_session.run(s)
    print('######################### start res $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
    records = []
    for r in result:
        records.append({
            'individual': r['individual'],
            'gender': r['gender'],
            'phenotypes': [dict(x) for x in r['phenotypes']],
            'phenotypeScore': r['phenotypeScore'],
            'hom_count': r['hom_count'],
            'het_count': r['het_count'],
            'genes': [y for y in r['genes']]
        })

    print(records)
    d=[]
    d.append( {
        'first_name': 'Jodi',
        'second_name': 'Hunt',
        'titles': ['Dr', 'Developer'],
    })
    d.append( {
        'first_name': 'Lou',
        'second_name': 'Sams',
        'titles': ['Dr', 'Developer'],
    })

    j=(jsonify(result=d))
    return jsonify(result=records)


@app.route('/my_patients_json')
@requires_auth
def my_patients_json():
    users_db=get_db(app.config['DB_NAME_USERS'])
    user=users_db.users.find_one({'user':session['user']})
    individuals=get_individuals(user['user'])
    return(jsonify(result=individuals))


# shows each patients, 
# all_individuals
@app.route('/my_patients')
@requires_auth
def my_patients():
    return render_template('my_patients.html')

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
