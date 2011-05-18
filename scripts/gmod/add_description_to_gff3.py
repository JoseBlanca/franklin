'''
Given an input GFF3 file, the database name and a file with the relationship
between the gff3 IDs and the database accessions for each feature the script
adds the dbxref to the GFF3 file
'''

from optparse import OptionParser
from franklin.gff import modify_gff3, create_add_description_mapper

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--ingff3', dest='ingff',
                      help='Input GFF3 file')
    parser.add_option('-o', '--outgff3', dest='outgff',
                      help='Output GFF3 file')
    parser.add_option('-d', '--descriptions', dest='desc',
                      help='Descriptions file')

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
    if not options.desc:
        parser.error('Description file is required')
    else:
        desc_fhand = open(options.desc)
    return ingff3_fpath, outgff3_fpath, desc_fhand


def main():
    'It runs the script'
    ingff3_fpath, outgff3_fpath, desc_fhand = get_parameters()

    mappers = []
    mappers.append(create_add_description_mapper(desc_fhand))

    modify_gff3(ingff3_fpath, outgff3_fpath, mappers)

if __name__ == '__main__':
    main()
