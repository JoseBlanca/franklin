'''
Given an input GFF3 file, the database name and a file with the relationship
between the gff3 IDs and the database accessions for each feature the script
adds the dbxref to the GFF3 file
'''

from optparse import OptionParser
from franklin.gff import modify_gff3, create_go_annot_mapper

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--ingff3', dest='ingff',
                      help='Input GFF3 file')
    parser.add_option('-o', '--outgff3', dest='outgff',
                      help='Output GFF3 file')
    parser.add_option('-a', '--annots', dest='annots',
                      help='b2g annot file')

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
    if not options.annots:
        parser.error('b2g annot file is required')
    else:
        annots_fhand = open(options.annots)
    return ingff3_fpath, outgff3_fpath, annots_fhand


def main():
    'It runs the script'
    ingff3_fpath, outgff3_fpath, annots_fhand = get_parameters()

    mappers = []
    mappers.append(create_go_annot_mapper(annots_fhand))

    modify_gff3(ingff3_fpath, outgff3_fpath, mappers)

if __name__ == '__main__':
    main()
