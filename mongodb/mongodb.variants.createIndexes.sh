
mongo
use uclex
db.variants.createIndex({'VARIANT_ID':1})
db.variants.createIndex({'Transript':1})
db.variants.createIndex({'Gene':1})
db.variants.createIndex({'genes':1})

db.variants.createIndex({'SYMBOL':1})
db.variants.createIndex({'HET':1})
db.variants.createIndex({'HOM':1})


[
{
"v" : 1,
"key" : {
"_id" : 1
},
"name" : "_id_",
"ns" : "uclex.variants"
},

{
"v" : 1,
"unique" : true,
"key" : {
"variant_id" : 1
},
"name" : "variant_id_1",
"ns" : "uclex.variants"
},

{
"v" : 1,
"key" : {
"CHROM" : 1
},
"name" : "CHROM_1",
"ns" : "uclex.variants"
},

{
"v" : 1,
"key" : {
"canonical_gene_name_upper" : 1
},
"name" : "canonical_gene_name_upper_1",
"ns" : "uclex.variants"
},

{
"v" : 1,
"key" : {
"filter" : 1
},
"name" : "filter_1",
"ns" : "uclex.variants"
},
{
"v" : 1,
"key" : {
"het_samples" : 1
},
"name" : "het_samples_1",
"ns" : "uclex.variants"
},
{
"v" : 1,
"key" : {
"hom_samples" : 1
},
"name" : "hom_samples_1",
 "ns" : "uclex.variants"
  },
   {
         "v" : 1,
           "key" : {
              "canonical_cadd" : 1
                },
                  "name" : "canonical_cadd_1",
                    "ns" : "uclex.variants"
                     },
                      {
                            "v" : 1,
                              "key" : {
                                 "FILTER" : 1
                                   },
                                     "name" : "FILTER_1",
                                       "ns" : "uclex.variants"
                                        },
                                         {
                                               "v" : 1,
                                                 "key" : {
                                                    "EXAC" : "hashed"
                                                      },
                                                        "name" : "EXAC_hashed",
                                                          "ns" : "uclex.variants"
                                                           },
                                                            {
                                                                  "v" : 1,
                                                                    "key" : {
                                                                       "most_severe_consequence" : 1
                                                                         },
                                                                           "name" : "most_severe_consequence_1",
                                                                             "ns" : "uclex.variants"
                                                                              }
                                                                              ]
