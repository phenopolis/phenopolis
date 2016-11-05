
import json
import urllib 
import pymongo
from itertools import chain

conn=pymongo.MongoClient(host='localhost',port=27017)
db=conn['hpo']

for r in db.tiger_eye_genes.find():
    refseq=r['RefSeq']
    url = 'http://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&format=row&input=refseqmrnaaccession&inputValues={}&outputs=ensembltranscriptid,ensemblgeneid&,ensemblgeneidtaxonId=9606'.format( refseq )
    u = urllib.urlopen(url)
    s=json.loads(u.read())
    d=dict()
    d['ensembl_transcript_id'] = list(chain(*[x['Ensembl Transcript ID'].split('//') for x in s]))
    d['ensembl_gene_id'] = list(chain(*[x['Ensembl Gene ID'].split('//') for x in s]))
    print(refseq,d)
    print(db.tiger_eye_genes.update({'RefSeq':refseq},{'$set':d},w=0))

