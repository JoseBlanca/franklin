#!/usr/bin/env python
'''It merges the genetic markers in the gff to chado. This gff is the gff used
populate the gbrowse's phisical map
Inputs: cmap_gff and physical_map
outputs: 3 gff3 to chado. One with the physical map and other with the markers
         not located inthe physical map

'''
from optparse import OptionParser

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-p', '--physical', dest='physical',
                      help='physical map file')
    parser.add_option('-g', '--genetic', dest='genetic',
                      help='genetic_map file')
    parser.add_option('-o', '--output', dest='output',
                      help='output file')
    parser.add_option('-m', '--markers', dest='orphan',
                      help='markers output file')
    return parser

def set_parameters():
    'Set parameters'
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.physical is None:
        parse.error('In file required')
    else:
        physical_fhand = open(options.physical)

    if options.genetic is None:
        parse.error('In file required')
    else:
        genetic_fhand = open(options.genetic)

    if options.output is None:
        outfhand = open('output.gff3', 'w')
    else:
        outfhand = open(options.output, 'w')

    if options.orphan is None:
        orphan_fhand = open('orphaned_markers.gff3', 'w')
    else:
        orphan_fhand = open(options.orphan, 'w')

    return physical_fhand, genetic_fhand, outfhand, orphan_fhand

def main():
    'The main part'
    # get parameters
    physical_fhand, genetic_fhand, outfhand, orphan_fhand = set_parameters()

    # get genetic markers
    genetic_markers = index_genetic_markers(genetic_fhand)

    # insert genetic_markers in physical markers
    orphan_markers = merge_markers(physical_fhand, outfhand, genetic_markers)


def merge_markers(physical_fhand, outfhand,genetic_markers):
    '''It merges the markers in the given outfhand. The orpahned ones are
    returned as a list'''
    genetic_marker_names = genetic_markers.keys()

    for line in physical_fhand:
        line = line.split()
        if not line:
            continue
        if line[0] != '#':
            line.split()



        outfhand.write(line + '\n')



def parse_gff3_feat(line):
    'It parses a gff3'

defindex_genetic_markers(fhand):
    'It indexes the genetic markers'
    markers = {}
    for line in fhand:
        line = line.strip()
        if not line or line[0] == '#':
            continue
        name = get_value_from_gff3_attrs(line.split('\t')[8], 'ID')
        if name not in markers:
            markers[name] = []
        markers[name].append()
    print markers

def get_value_from_gff3_attrs(attrs, field='ID'):
    'It gets the the dbxref from a field of the attrs field. Fallbacks to ID'
    for items in attrs.split(';'):
        key, value = items.split('=')
        if key == field:
            return value
    raise ValueError("the field  %s is not in this gff3 attr field" % field)




if __name__ == '__main__':
    main()




