from views import *
from lookups import *
from orm import *



@app.route('/transcript_json/<transcript_id>')
@requires_auth
def transcript_json(transcript_id):
    db = get_db()
    def f(v):
        del v['_id']
        if session['user']=='demo':
            del v['het_samples']
            del v['hom_samples']
            del v['wt_samples']
        return v
    variants=[f(v) for v in db.variants.find({'canonical_transcript':transcript_id})]
    #cache_key = 't-transcript-{}'.format(transcript_id)
    #t = cache.get(cache_key)
    #print 'Rendering %stranscript: %s' % ('' if t is None else 'cached ', transcript_id)
    #if t: return t
    #cache.set(cache_key, t)
    return jsonify(result={'variants':variants,'count':len(variants)})


@app.route('/transcript/<transcript_id>')
@requires_auth
def transcript(transcript_id):
    db = get_db()
    transcript=db.transcripts.find_one({'transcript_id':transcript_id})
    transcript['variants']=[Variant(variant_id=v['variant_id'],db=db) for v in db.variants.find({'canonical_transcript':transcript_id})]
    #cache_key = 't-transcript-{}'.format(transcript_id)
    #t = cache.get(cache_key)
    #print 'Rendering %stranscript: %s' % ('' if t is None else 'cached ', transcript_id)
    #if t: return t
    #cache.set(cache_key, t)
    individuals=dict()
    for v in transcript['variants']:
        if v.canonical_hgvsc[0]:
            v.cdna_pos=v.canonical_hgvsc[0].split(':')[1].split('.')[1]
        else:
            v.cdna_pos=''
        v.canonical_hgvs=dict(zip( v.canonical_hgvsp, v.canonical_hgvsc))
        v.__dict__['protein_mutations']=dict([(p,p.split(':')[1],) for p in v.canonical_hgvsp if ':' in p])
        for csq in v.transcript_consequences:
            if csq['transcript_id']!=transcript_id: continue
            v.distance=csq.get('distance','')
        for s in v.het_samples:
            if v.HET_COUNT < 10: individuals[s]=individuals.get(s,[])+[v]
    table_headers=re.findall("<td class='?\"?(.*)-cell'?\"?>",file('templates/transcript.html','r').read())
    return render_template('transcript.html',transcript=transcript,individuals=individuals,table_headers=table_headers)



