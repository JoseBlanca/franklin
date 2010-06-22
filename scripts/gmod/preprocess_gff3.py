#!/usr/bin/env python
'''This scripts preprocesses a gff3 file and coverts some feature_types to sofa
compatible ones'''
from optparse import OptionParser
from franklin.gff import gff_parser, write_gff
import sys

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--infile', dest='infile', help='infile')
    parser.add_option('-o', '--outfile', dest='outfile', help='outfile')
    return parser

def set_parameters():
    'Set parameters'
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.infile is None:
        parser.error('In file required')
    else:
        infhand = open(options.infile)

    if options.outfile is None:
        outfhand = sys.stdout
    else:
        outfhand =  open(options.outfile, 'w')

    return infhand, outfhand

def main():
    'The main part'
    # set parameters
    infhand, outfhand = set_parameters()

    #get features
    features = gff_parser(infhand, version=3)

    # correct_feature_type
    features = correct_feature_type(features)

    # write_features
    write_gff(features, outfhand)

def correct_feature_type(features):
    'it corrects  the feature_type'
    for feature in features:
        feature_type = feature['type']
        if feature_type == 'Chromosome':
            feature['type'] = 'chromosome'
        elif feature_type in('frameworkmarker', 'marker', 'placementmarker'):
            feature['type'] = 'genetic_marker'
        yield feature


def correct_feature_type2(line):
    'It changes the feature names to sofa compatible ones'
    items = line.split()
    feature_type = items[2]
    if feature_type == 'Chromosome':
        items[2] = 'chromosome'
    elif feature_type in('frameworkmarker', 'marker', 'placementmarker'):
        items[2] = 'genetic_marker'
    return '\t'.join(items)




if __name__ == '__main__':
    main()

