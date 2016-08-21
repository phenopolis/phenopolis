import rest
from variant import Variant


class Transcript(object):
    db=None
    def __init__(self, transcript_id, db=None):
        self.transcript_id=transcript_id
        if db is None: return
        Transcript.db=db
        self.data=Transcript.db.transcripts.find_one({'transcript_id':self.transcript_id},fields={'_id':False})
        if not self.data: raise Exception(transcript_id, 'not in db')
        #self.xpos = get_xpos(self.chrom, self.pos)
        #self.variant = lookups.get_variant(self.db, self.xpos, self.ref, self.alt)
        self.__dict__.update(self.data) 
    def get_variants_in_transcript(self):
        """
        """
        variants = []
        print('transcripts')
        for variant in Transcript.db.variants.find({'transcripts': self.transcript_id}, fields={'_id': False}):
            variant['vep_annotations'] = [x for x in variant['vep_annotations'] if x['Feature'] == self.transcript_id]
            #add_consequence_to_variant(variant)
            #remove_extraneous_information(variant)
            variants.append(variant)
        return variants


