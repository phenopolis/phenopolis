# PubmedScore

Genes are prioritised by PubmedScore calculated using literature support. PubMed publications are mined for interesting phenotypes given relevant keywords. Users may also supply a list of genes of interest that are known to be associated with the phenotype. PubmedScore will assign at least a score of 1 to the genes on the given list. This will help to keep those genes on top of the result table if the user choose to sort it by the Pubmed Relevance Score.


Genes are prioritised by PubmedScore calculated using literature support. PubMed publications are mined for interesting phenotypes given relevant keywords. Users may also supply a list of genes of interest that are known to be associated with the phenotype. PubmedScore will assign at least a score of 1 to the genes on the given list. This will help to keep those genes on top of the result table if the user choose to sort it by the Pubmed Relevance Score.


Define an indicator function B<sub>known</sub>, given a gene (g):

B<sub>known</sub>(g) = { if the gene is on a given list of genes known to be associated with the phenotype : 1 ,  else : 0 }


The pubmedscore(g) is calculated by counting the appearance of each keyword (user defined) in every returned pubmed paper’s title and abstract, and taking the sum.


The overall Pubmed Relevance (Sp) score of a given gene (g) is calculated as:

Sp(g) = max(B<sub>known</sub>(g), pubmedscore(g))


One can also insert an additional function adapted to problems at hand. Taking patients with retinal dystrophy as an example: [Retnet] (https://sph.uth.edu/Retnet/) has detailed information for retinal dystrophy associated genes, such as if a gene may cause dominant or recessive inheritance mode of a disease. This can then be used to help calculate Sp. Given:


B<sub>ret</sub>(g) = { if the gene is reported on RetNet : 1, else : 0 }

B<sub>retMode</sub>(g) = { if the gene is registered on RetNet and changes on the gene is described therein as possibly causing Dominant or X-linked or M-linked diseases : 1, else : 0 }


B<sub>searchMode</sub>(g) = { if it is searched in the context of recessive inheritance mode : 1, else : 0 }

Then 

Sp(g) = max(B<sub>known</sub>(g), pubmedscore(g), fc(g))

where

fc(g) = B<sub>ret</sub>(g) ✕ 100 ✕ [B<sub>retMode</sub>(g) or B<sub>searchMode</sub>(g)]

The function is then guaranteed to rank the known retinal dystrophy associated genes on top of the table that matches the inheritance mode.


Take gene DRAM2 as an example (can be found as candidate gene for two of the demo patients): Pubmed search with keywords ‘blindness macula macular pigmentosa retina retinal retinitis stargardt’ returns three publications (pubmed IDs: 27518550, 26720460, 25983245). The app then counted the appearance of the keywords in the titles and abstracts of the three publications, which equaled to 16. The pubmedscore(DRAM2) was then assigned as 16. However, since DRAM2 was registered on RETNET as a recessive gene, and this patient has a relevant homozygous change on DRAM2, therefore

* B<sub>ret</sub>(DRAM2)=1
* B<sub>retMode</sub>(DRAM2)=0
* B<sub>searchMode</sub>(DRAM2)=1 

One can derive that

fc(DRAM2)=1 ✕ 100 ✕ [0 or 1]=100
	
And hence the Pubmed Relevance score of <i>DRAM2</i> is

Sp(g)=max(0,16,100)=100




