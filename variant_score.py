
# given a pubmed result for a gene/search term, calculate
#  each paper's pubmed score, and sum them up to calculate 
# a total score for the gene
# search_mode = 1 if search hom / comp_het. else 0
# reg are the search terms, refer to pubmed code  since this modifies gene_pubmed in place, , it won't return anything
'''
gene_pubmed = {
    results:[{
        id: pubmed_id,
        title: title,
        abstract: abstract,
        bad: 0 or 1 #
        },
        ...]
    }
'''
def pubmed_scoring(symbol, gene_pubmed, reg, known_gene=None, retnet=None, search_mode=1):
    # get known_genes and retnet
    known_genes = known_genes or open('ret_known_genes.txt', 'r').readline().strip().split()
    retnet  = retnet or json.load(open('retnet.json', 'r'))
    total_score = 0
    # calculate pubmed scores for each record
    for r,v in gene_pubmed['result'].iteritems():
        if v['bad']:
            v['score'] = 0
            continue
        score = 0
        if 'title' in v:
            score += len(reg.findall(v['title']))
        if 'abstract' in v:
            score += len(reg.findall(v['abstract']))
        v['score'] = score
        total_score += score
    # sort
    gene_pubmed['results'] = sorted(gene_pubmed['results'], key=lambda k: k['score'], reverse=True)
    # calculate the total score
    known = 1 if symbol in known_genes else 0
    in_retnet = 1 if symbol in retnet else 0
    ret_mode = 1 if symbol in retnet and retnet[symbol][mode] in ['d', 'x', 'm'] else 0
    gene_pubmed['total_score'] = max(known, total_score, in_retnet * 100**(ret_mode or search_mode))


