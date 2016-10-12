from operator import itemgetter
import lookups

AF_BUCKETS = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]
METRICS = [
    'BaseQRankSum',
    'ClippingRankSum',
    'DP',
    'FS',
    'InbreedingCoeff',
    'MQ',
    'MQRankSum',
    'QD',
    'ReadPosRankSum',
    'VQSLOD'
]


def add_transcript_coordinate_to_variants(db, variant_list, transcript_id):
    """
    Each variant has a 'xpos' and 'pos' positional attributes.
    This method takes a list of variants and adds a third position: the "transcript coordinates".
    This is defined as the distance from the start of the transcript, in coding bases.
    So a variant in the 7th base of the 6th exon of a transcript will have a transcript coordinate of
    the sum of the size of the first 5 exons) + 7
    This is 0-based, so a variant in the first base of the first exon has a transcript coordinate of 0.

    You may want to add transcript coordinates for multiple transcripts, so this is stored in a variant as
    variant['transcript_coordinates'][transcript_id]

    If a variant in variant_list does not have a `transcript_coordinates` dictionary, we create one

    If a variant start position for some reason does not fall in any exons in this transcript, its coordinate is 0.
    This is perhaps logically inconsistent,
    but it allows you to spot errors quickly if there's a pileup at the first base.
    `None` would just break things.

    Consider the behavior if a 20 base deletion deletes parts of two exons.
    I think the behavior in this method is consistent, but beware that it might break things downstream.

    Edits variant_list in place; no return val
    """
    # make sure exons is sorted by (start, end)
    exons = sorted(lookups.get_exons_in_transcript(db, transcript_id), key=itemgetter('start', 'stop'))
    # offset from start of base for exon in ith position (so first item in this list is always 0)
    exon_offsets = [0 for i in range(len(exons))]
    for i, exon in enumerate(exons):
        for j in range(i+1, len(exons)):
            exon_offsets[j] += exon['stop'] - exon['start']

    for variant in variant_list:
        if 'transcript_coordinates' not in variant:
            variant['transcript_coordinates'] = {}
        variant['transcript_coordinates'][transcript_id] = 0
        for i, exon in enumerate(exons):
            if exon['start'] <= variant['pos'] <= exon['stop']:
                variant['transcript_coordinates'][transcript_id] = exon_offsets[i] + variant['pos'] - exon['start']




protein_letters_1to3 = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu',
    'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
    'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
    'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
    'X': 'Ter', '*': 'Ter'
}


def get_proper_hgvs(csq):
    # Needs major_consequence
    if csq['major_consequence'] in ('splice_donor_variant', 'splice_acceptor_variant', 'splice_region_variant'):
        return get_transcript_hgvs(csq)
    else:
        return get_protein_hgvs(csq)


def get_transcript_hgvs(csq):
    return csq['HGVSc'].split(':')[-1]

CHROMOSOMES = ['chr%s' % x for x in range(1, 23)]
CHROMOSOMES.extend(['chrX', 'chrY', 'chrM'])
CHROMOSOME_TO_CODE = { item: i+1 for i, item in enumerate(CHROMOSOMES) }

def get_single_location(chrom, pos):
    """
    Gets a single location from chromosome and position
    chr must be actual chromosme code (chrY) and pos must be integer

    Borrowed from xbrowse
    """
    return CHROMOSOME_TO_CODE[chrom] * int(1e9) + pos


def get_xpos(chrom, pos):
    """
    Borrowed from xbrowse
    """
    if not chrom.startswith('chr'):
        chrom = 'chr{}'.format(chrom)
    return get_single_location(chrom, int(pos))



