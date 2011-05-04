'''
It takes as input Tyler generates gff and converts it into a sequence ontology
recomendations formated gff
'''
from franklin.gff import GffFile, FEATURE, COMMENT, create_id_to_name_mapper
from optparse import OptionParser
from operator import itemgetter

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--ingff3', dest='ingff',
                      help='Input GFF3 file')
    parser.add_option('-o', '--outgff3', dest='outgff',
                      help='Output GFF3 file')
    return parser

def get_parameters():
    'It reads and fixes the parsed command line options'
    parser  = parse_options()
    options = parser.parse_args()[0]

    if not options.ingff:
        parser.error('An input GFF3 file is required')
    else:
        ingff3_fpath = options.ingff
    if not options.outgff:
        parser.error('An output GFF3 file is required')
    else:
        outgff3_fpath = options.outgff

    return ingff3_fpath, outgff3_fpath

def main():
    'It runs the script'
    ingff3_fpath, outgff3_fpath = get_parameters()
    gff_out = GffFile(fpath=outgff3_fpath, mode='w')

    for gene in genes_in_gff(ingff3_fpath):
        gene = reorder_gene(gene)
        join_same_exons(gene)
        join_same_cds(gene)
        add_polypeptide_feature(gene)
        write_gene(gene, gff_out)

def _get_feat_loc(feature):
    'it gets the feature loc'
    start = feature['start']
    end = feature['end']
    strand = feature['strand']
    phase = feature['phase']
    return (start, end, strand, phase)

def _change_feature_name(feature, name):
    'It changes name and feature with the given name'
    feature['name'] = name
    feature['id'] = name
    feature['attributes']['ID'] = name
    feature['attributes']['Name'] = name

def get_feature_fingerprint(feature):
    'it extracts the feature fingerprint'
    if isinstance(feature, list):
        loc = []
        for cds_part in feature:
            loc.append(_get_feat_loc(cds_part))
    else:
        loc = _get_feat_loc(feature)
    return tuple(loc)

def join_same_exons(gene):
    'It joins same exons in exons'

    exons = gene['exon']
    exon_dict = {}
    finger_exons = [(get_feature_fingerprint(exon), exon) for exon in exons]
    for finger,exon in finger_exons:
        if finger not in exon_dict:
            exon_dict[finger] = exon
        else:
            parent = exon['attributes']['Parent']
            exon_dict[finger]['attributes']['Parent'] += ',' + parent
    exons_sorted = sorted(exon_dict.values(), key=itemgetter('start'))
    gene_name = gene['gene']['name']
    exons = []
    for index, exon in enumerate(exons_sorted):
        name = '%s.exon%d' % (gene_name, index + 1)
        _change_feature_name(exon, name)
        exons.append(exon)
    gene['exon'] = exons

def join_same_cds(gene):
    'it joins same cds'
    cdss = gene['CDS']
    cds_dict = {}
    finger_cds = [(get_feature_fingerprint(cds), cds) for cds in cdss.values()]
    for finger, cds in finger_cds:
        if finger not in cds_dict:
            cds_dict[finger] = cds
        else:

            parent = cds[0]['attributes']['Parent']
            for cds_ in cds_dict[finger]:
                cds_['attributes']['Parent'] += ',' + parent
    gene_name = gene['gene']['name']

    cds_sorted = cds_dict.values()
    cdss = []
    for index, cds in enumerate(cds_sorted):
        name = '%sC%d' % (gene_name, index + 1)
        for cds_part in cds:
            _change_feature_name(cds_part, name)
            # remove Target for protein
            protein = cds_part['attributes']['Target'].split()[0]
            del(cds_part['attributes']['Target'])
            cds_part['attributes']['protein_id'] = protein
        cdss.append(cds)
    gene['CDS'] = cdss
def add_polypeptide_feature(gene):
    'Taking into account the cds create a polypeptide feature'
    for cds in gene['CDS']:
        start, end = None, None
        for cds_part in cds:
            part_start = cds_part['start']
            part_end  = cds_part['end']
            if start is None or start > part_start:
                start = part_start
            if end is None or end < part_end:
                end = part_end
        name   = cds[0]['attributes']['protein_id']
        seqid  = cds[0]['seqid']
        strand = cds[0]['strand']
        source = cds[0]['source']
        parent = cds[0]['attributes']['Parent']
        polypeptide_feat = {'start':start, 'end':end, 'name':name, 'id':name,
                            'seqid':seqid, 'phase':'.', 'strand':strand,
                            'source':source, 'score':'.', 'type':'polypeptide',
                            'attributes':{'Derives_from':parent}}
        if 'polypeptide' not in gene:
            gene['polypeptide'] = []
        gene['polypeptide'].append(polypeptide_feat)

def write_gene(gene, out_gff):
    'It writes the gene feature'
    out_gff.write(FEATURE, gene['gene'])

    for transcript in gene['transcript']:
        out_gff.write(FEATURE, transcript)
    for exon in gene['exon']:
        out_gff.write(FEATURE, exon)
    for cds in gene['CDS']:
        for cds_part in cds:
            out_gff.write(FEATURE, cds_part)
    for polypeptide in gene['polypeptide']:
        out_gff.write(FEATURE, polypeptide)

    out_gff.write(COMMENT, '##')

def reorder_gene(features):
    'It organices the features inside a list'
    gene = {'gene':None, 'CDS':{}, 'exon':[], 'transcript':[]}
    for feature in features:
        kind = feature['type']
        if kind == 'gene':
            gene[kind] = feature
        elif kind == 'CDS':
            id_ = feature['id']
            if id_ not in gene['CDS']:
                gene['CDS'][id_] = []
            gene['CDS'][id_].append(feature)
        else:
            gene[kind].append(feature)
    return gene

def genes_in_gff(fpath):
    '''It yields the features containing genes. It suposes to have a sorted gff
    file: All exons, cdss and popypeptides are after the gene feature'''
    mapper = create_id_to_name_mapper()
    gff = GffFile(fpath=fpath, feature_mappers=[mapper])
    gene_feats = []
    for feature in gff.features:
        kind = feature['type']
        if kind == 'gene':
            if gene_feats:
                yield gene_feats
            gene_feats = []
        gene_feats.append(feature)
    else:
        if gene_feats[0]['type'] == 'gene':
            yield gene_feats

def test():
    from tempfile import NamedTemporaryFile
    from franklin.utils.misc_utils import TEST_DATA_DIR
    import os
    gff_fpath = os.path.join(TEST_DATA_DIR, 'previous.gff3')
    gff_in = GffFile(fpath=gff_fpath)

    gff_out_fhand = NamedTemporaryFile()
    gff_out = GffFile(fpath=gff_out_fhand.name, mode='w')

    for feature in gff_in.features:
        gff_out.write(FEATURE, feature)
    gff_out_fhand.flush()

    # test genes_in_gff
    genes = list(genes_in_gff( gff_out_fhand.name))

    assert len(genes)    == 3
    assert len(genes[0]) == 11
    assert len(genes[1]) == 11

    # test_reorder features inside genes
    gene_ordered = reorder_gene(genes[0])
    assert len(gene_ordered['transcript']) == 2
    assert len(gene_ordered['exon']) == 4
    assert len(gene_ordered['CDS'].keys()) == 2
    assert len(gene_ordered['CDS']['MELO3A000483C1']) == 2


    #test join same exones
    join_same_exons(gene_ordered)
    assert len(gene_ordered['exon']) == 3

    # test join same cds
    join_same_cds(gene_ordered)
    assert len(gene_ordered['CDS']) == 1

    # create polypeptide feature
    add_polypeptide_feature(gene_ordered)
    assert len(gene_ordered['polypeptide']) == 1

    # write gene
    fhand = NamedTemporaryFile()
    write_gene(gene_ordered, out_gff=GffFile(fpath=fhand.name, mode='w'))
    #print open(fhand.name).read()

if __name__ == '__main__':
    test()
    main()



