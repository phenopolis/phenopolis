import json
from views import *
from lookups import *
from orm import *
import rest as annotation
import requests
from config import config
import myvariant
from vcf import vcf_query
import hashlib
from bson.json_util import dumps
from neo4j.v1 import GraphDatabase, basic_auth


@app.route('/pheno4j/',methods=['GET'])
@requires_auth
def pheno4j():
    neo=get_neo4j()
    result = neo.run("MATCH (a:Person) return a.personId as personId ")
    s='\n'.join([ "%s" % (record["personId"]) for record in result ])
    return s


@app.route('/rv_sharing/<individual_id>/<thresh>/<allele_freq>/<limit>')
@requires_auth
def rv_sharing(individual_id,thresh,allele_freq,limit):
    #thresh=0.05
    #allele_freq=0.001
    print individual_id
    print float(thresh)
    print float(allele_freq)
    neo=get_db('neo4j')
    q= """ MATCH (k:Person)
    WITH count(k) as numberOfPeople
    MATCH (p:Person {{personId:"{personId}"}})<-[:PRESENT_IN]-(gv:GeneticVariant)
    WHERE (gv.allele_freq < {allele_freq} or gv.hasExac = false)
    WITH size(()<-[:PRESENT_IN]-(gv)) as count , gv, p, numberOfPeople
    WHERE count > 1 
    AND ((count / toFloat(numberOfPeople))  <= {thresh})
    MATCH (gv)-[:PRESENT_IN]->(q:Person)
    WHERE p <> q
    WITH p,q,count(gv) as intersection, numberOfPeople
    ORDER BY intersection DESC limit {limit}
    MATCH (x:Person)<-[:PRESENT_IN]-(v:GeneticVariant)
    WHERE (x.personId = p.personId or x.personId = q.personId)
    AND (v.allele_freq < {allele_freq} or v.hasExac = false)
    AND ((size(()<-[:PRESENT_IN]-(v)) / toFloat(numberOfPeople))  <= {thresh})
    WITH p, q, v, intersection
    RETURN p.personId, q.personId, intersection, size(collect(distinct v)) as unionSum, (round((intersection/toFloat(size(collect(distinct v))))*100.0*10)/10) as PercentShared
    ORDER BY PercentShared DESC;
    """.format(thresh=float(thresh), personId=individual_id,allele_freq=float(allele_freq),limit=int(limit))
    result = neo.run(q)
    get_db('neo4j').close()
    return json.dumps([r.__dict__ for r in result], indent=4)


