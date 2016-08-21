import rest
import json
import orm


class Gene(object):
    def __init__(self, gene_id, db=None):
        self.gene_id=gene_id
        Gene.db=db
        self.data=Gene.db.genes.find_one({'gene_id':gene_id},{'_id':False})
        self.__dict__.update(self.data)
    @property
    def variants(self):
        #if 'variant_ids' not in self.__dict__:
        variants=[v for v in Gene.db.variants.find({'genes': self.gene_id}, fields={'_id': False})]
        print('number of variants', len(variants))
        self.__dict__['variant_ids']=[v['variant_id'] for v in variants]
        #self.save()
        orm.Variant.db=Gene.db
        for variant in variants:
            try:
                v=orm.Variant(variant_id=variant['variant_id'],data=variant)
            except Exception, e:
                print e
                continue
            yield v
    @property
    def variant_ids(self):
        if 'variant_ids' in self.__dict__: return self.__dict__['variant_ids']
        variants=[v for v in Gene.db.variants.find({'genes': self.gene_id}, fields={'_id': False})]
        self.__dict__['variant_ids']=[v['variant_id'] for v in variants]
        self.save()
        return self.__dict__['variant_ids']
    @property
    def transcripts(self):
        if 'transcripts' in self.__dict__: return self.__dict__['transcripts']
        self.__dict__['transcripts']=rest.mg.getgene(self.gene_id,'ensembl.transcripts')
        print(self.save())
        return self.__dict__['transcripts']
    @property
    def canonical_transcript(self):
        if 'transcripts' in self.__dict__: return self.__dict__['transcripts'][0]
        self.__dict__['transcripts']=rest.mg.getgene(self.gene_id,'ensembl.transcripts')
        print(self.save())
        return self.__dict__['transcripts'][0]
    @property
    def summary(self):
        if 'summary' in self.__dict__: return self.__dict__['summary']
        self.__dict__['summary']=rest.mg.getgene(self.gene_id,'reporter.summary')
        #self.save()
        return self.__dict__['summary']
    @property
    def exons(self):
        if 'exons' in self.__dict__: return self.__dict__['exons']
        self.__dict__['exons']=rest.mg.getgene(self.gene_id,'exons')
        #self.save()
        return self.__dict__['exons']
    def save(self):
        print('writing', self.gene_id, 'to database')
        return Gene.db.genes.update({'gene_id':self.gene_id},self.__dict__,upsert=True)







