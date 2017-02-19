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



@app.route('/search', methods=['GET','POST'])
@requires_auth
def search():
    cache_key = 't-homepage'
    #t = cache.get(cache_key)
    #if t: return t
    db=get_db()
    patients_db=get_db(app.config['DB_NAME_PATIENTS']) 
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

    try:
        version_number = subprocess.check_output(['git', 'describe', '--exact-match'])
    except:
        version_number = None
    print('Version number is:-')
    print(version_number)

    t = render_template('search.html',
        title='home',
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
        #image=image.decode('utf8'))
        image="",
        version_number=version_number)
    #cache.set(cache_key, t)
    return t


