'''
It creates the relation file with the id and the dbxref accesion of a gbrowse
database. Taking into account the seqid and the start and end
'''
from optparse import OptionParser
from franklin.gff import GffFile

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--ingff3', dest='ingff',
                      help='Input GFF3 file')
    parser.add_option('-o', '--output', dest='output',
                      help='Output relations file')


    return parser

def get_parameters():
    'It reads and fixes the parsed command line options'
    parser  = parse_options()
    options = parser.parse_args()[0]

    if not options.ingff:
        parser.error('An input GFF3 file is required')
    else:
        ingff3_fpath = options.ingff
    if not options.output:
        parser.error('An output GFF3 file is required')
    else:
        output_fhand = open(options.output, 'w')

    return ingff3_fpath, output_fhand

def main():
    'The main part'
    ingff3_fpath, output_fhand = get_parameters()
    gff = GffFile(ingff3_fpath)

    for feature in gff.features:
        id_ = feature['id']
        gbrowse_acc = '%s:%s..%s' % (feature['seqid'], str(feature['start']),
                                     str(feature['end']))
        output_fhand.write('%s\t%s\n' % (id_, gbrowse_acc))

if __name__ == '__main__':
    main()






