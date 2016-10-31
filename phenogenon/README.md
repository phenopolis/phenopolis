# Phenogenon

Unrelated patients with observed HPO terms are included in the analysis, the number of which is denoted as Pata. The number of patients affected by a given HPO term h is denoted as Pat<sub>h</sub>. The number of patients who have a specific genotype in a gene g is denoted as Pat<sub>g</sub> (this may be a single variant or more than one for compound hets). The variants are filtered by ExAC allele frequency and CADD phred score  (M Kircher 2014) according to user specified thresholds so that only rare, predicted, damaging variants are included (frameshift variants which do not have a CADD score are assigned a default score of 50). The number of patients having both HPO term h and filtered variants in g is denoted as Pat<sub>gh</sub>.
Therefore one can construct a 2 × 2 table as follows, per gene g and HPO term h:



|                                         | Number of patients not affected by h  | Number of patients affected by h |
|-----------------------------------------|---------------------------------------|----------------------------------|
| Number of patients without variant in g |   Pat<sub>a</sub> - Pat<sub>h</sub> - Pat<sub>g</sub> + Pat<sub>gh</sub> |  Pat<sub>h</sub> - Pat<sub>gh</sub>               |
| Number of patients with variant in g    |    Pat<sub>g</sub> - Pat<sub>gh</sub>        |  Pat<sub>gh</sub>            |



Fisher’s Exact test is used to test the non-independence between g and h, and the phi correlation coefficient is used to quantify the correlation between g and h. Fisher’s Exact test produces two p values: right-tail and left-tail p values. A small right-tail p value indicates a positive correlation between g and h. The sign and magnitude of phi also determines the type of correlation. A strong positive phi indicates that rare damaging variants in the gene appear more often than expected by chance in patients with that HPO term.  This suggests the disrupted gene may cause the phenotype. A negative phi would indicate a gene has less rare damaging variants than expected by chance, which might suggest a phenotype is driven by a more functional gene.


Obviously this approach is currently susceptible to sequencing bias.  For instance, if all individuals which share a particular HPO term have improved coverage for a given gene then there will be an enrichment of rare variants in that gene.  We are currently working on extending the method to account for this.


This method is applied to the gene with two possible inheritance modes, dominant where patient has to have at least one qualified variant, and recessive where the patient has to have at least two qualified variants. The significance of the top p-values can also help infer whether a gene is more likely to cause dominant or recessive mode of the phenotype.


The HPO tree view also nicely illustrates syndromic genes which cause a range of phenotypes, for example <i>USH2A</i>:



