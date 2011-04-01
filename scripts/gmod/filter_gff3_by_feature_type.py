'''
Created on 2011 mar 29

@author: peio
'''
from optparse import OptionParser
from franklin.gff import modify_gff3, create_feature_type_filter

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--ingff3', dest='ingff',
                      help='Input GFF3 file')
    parser.add_option('-o', '--outgff3', dest='outgff',
                      help='Output GFF3 file')
    parser.add_option('-t', '--types', dest='types',
                    help='types to filter by. separated by commas')
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
    if not options.types:
        parser.error('needs typed to filter by')
    else:
        types = options.types.split(',')
    return ingff3_fpath, outgff3_fpath, types

def main():
    'It runs the script'
    ingff3_fpath, outgff3_fpath, types = get_parameters()

    filters = []
    filters.append(create_feature_type_filter(types))

    modify_gff3(ingff3_fpath, outgff3_fpath, filters=filters)

if __name__ == '__main__':
    main()
