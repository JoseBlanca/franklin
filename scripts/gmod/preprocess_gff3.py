#!/usr/bin/env python
'''This scripts preprocesses a gff3 file and coverts some feature_types to sofa
compatible ones'''

from optparse import OptionParser
from franklin.gff import features_in_gff, write_gff, get_gff_header
import sys
from itertools import imap

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

    # get header
    header = get_gff_header(infhand)

    #get features
    features = features_in_gff(infhand, version=3)

    # corrections
    corrections = [correct_feature_type, remove_attributes]

    for correction in corrections:
        features = imap(correction, features)

    # write_features
    write_gff(features, outfhand, header)

def correct_feature_type(feature):
    'it corrects  the feature_type'
    feature_type = feature['type']
    if feature_type == 'Chromosome':
        feature['type'] = 'chromosome'
    elif feature_type in('frameworkmarker', 'marker', 'placementmarker'):
        feature['type'] = 'genetic_marker'
    elif feature_type == 'rflp':
        feature['type'] = 'RFLP_fragment'
    elif feature_type in('snp', 'SNP-CAPS', 'SNP-Snapshot'):
        feature['type'] = 'SNP'
    elif feature_type == 'SNP-INDEL':
        feature['type'] = 'indel'
    elif feature_type in ('issr', 'rapd', 'aflp', 'morphological',
                          'isozyme', 'spelling error', 'Isozyme', 'RFLP'):
        feature['type'] = 'genetic_marker'
        feature['attributes']['marker_type'] = feature_type
    elif feature_type == 'sequence feature':
        feature['type'] = 'sequence_feature'
    elif feature_type in ('SSR', 'EST-SSR'):
        feature['type'] = 'microsatellite'
    return feature

def remove_attributes(feature):
    'It removes the feature tag from attributes'
    tags = ['Contig_hit', 'BAC', 'Sequence', 'Marker_hit']
    for tag in tags:
        if tag in feature['attributes']:
            del feature['attributes'][tag]
    return feature

if __name__ == '__main__':
    main()

