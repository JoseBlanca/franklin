#!/usr/bin/env python

'''
It creates a marker centralized file to be used for the rest of the scripts to
normalize all marker name and type fields

@author: peio
'''
from optparse import OptionParser
import os
from franklin.gff import features_in_gff

SOFA_TRADUCTOR = {'rflp': 'RFLP_fragment',
                  'snp' : 'SNP',
                  'SNP-CAPS': 'SNP',
                  'SNP-Snapshot': 'SNP',
                  'SNP-INDEL':'indel',
                  'SSR':'microsatellite',
                  'EST-SSR':'microsatellite',
                  'aflp':'genetic_marker',
                  'isozyme':'genetic_marker',
                  'issr':'genetic_marker',
                  'morphological':'genetic_marker',
                  'spelling error':'genetic_marker',
                  'rapd':'genetic_marker',
                  'placementmarker': 'biological_region' ,
                  'frameworkmarker': 'biological_region',
                  'marker':'biological_region',
                  }
ACCEPTED_MARKERS = ['genetic_marker', 'frameworkmarker', 'indel', 'marker',
                        'microsatellite', 'placementmarker', 'rflp', 'snp',
                        'isozyme', 'aflp', 'issr', 'morphological', 'rapd',
                        'spelling error']

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--infiles', dest='infiles', help='infiles')
    parser.add_option('-m', '--markersfile', dest='markersfile',
                      help='markers file')
    parser.add_option('-n', '--name_correlations', dest='name_cor',
                      help='Evaluable file with name correlations')
    parser.add_option('-t', '--type_correlations', dest='type_cor',
                      help='Evaluable file with type correlations')
    parser.add_option('-s', '--sequence_correlations', dest='sequence_cor',
                      help='Evaluable file with sequence correlations')
    return parser

def set_parameters():
    'Set parameters'
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.infiles is None:
        parser.error('In file required')
    else:
        infhands = [open(infile)for infile in options.infiles.split(',')]

    if options.markersfile is None:
        parser.error('We need a markers file to read or write markers to')
    else:
        markersfname =  options.markersfile

    if options.name_cor is None:
        name_cor = None
    else:
        name_cor = correlations_in_file(open(options.name_cor))

    if options.type_cor is None:
        type_cor = None
    else:
        type_cor = correlations_in_file(open(options.type_cor))
    if options.sequence_cor is None:
        sequence_cor = None
    else:
        sequence_cor = correlations_in_file(open(options.sequence_cor))

    return infhands, markersfname, name_cor, type_cor, sequence_cor

def main():
    'The main part'
    # set parameters
    infhands, markersfname, name_cor, type_cor, sequence_cor = set_parameters()

    #read markes file if exists
    if os.path.exists(markersfname):
        markersfhand = open(markersfname)
        markers = parse_markersfile(markersfhand)
        markersfhand.close()
    else:
        markers = {}

    # add files to marker structure
    for infhand in infhands:
        add_markers(infhand, markers, name_cor, type_cor, sequence_cor)

    # write markers to file:
    markersfhand = open(markersfname, 'w')
    write_markers(markersfhand, markers)
    markersfhand.close()

def correlations_in_file(fhand):
    'It indexes a correlation file'
    correlations = {}
    for line in fhand:
        line = line.strip()
        if not line:
            continue
        name, correlation = line.split()
        correlations[name] = correlation
    return correlations

def add_markers(infhand, markers, name_cor=None, type_cor=None,
                sequence_cor=None):
    'It adds markers from fhand correction the names and adding them as alias'
    # guess file type
    file_type = _guess_file_type(infhand)
    if file_type == 'gff3':
        _add_markers_gff3(infhand, markers, name_cor, type_cor, sequence_cor)

def _add_markers_gff3(infhand, markers, name_cor=None, type_cor=None,
                      sequence_cor = None):
    'It adds markers using a gff3 file'
    for feature in features_in_gff(infhand, 3):
        if feature['type'] not in ACCEPTED_MARKERS:
            continue
        original_name = feature['name'].lower()
        original_type = feature['type']
        if name_cor is not None and original_name in name_cor:
            name  = name_cor[original_name].lower()
            alias = original_name

        else:
            name  = original_name
            alias = None

        # Use the types found in correlation file
        if type_cor is not None and name in type_cor:
            type_ = type_cor[name]
        else:
            type_ = original_type
        # traduce from type_ to sofa
        if type_ in SOFA_TRADUCTOR:
            sofa_type = SOFA_TRADUCTOR[type_]
        else:
            sofa_type = type_

        if (('sequence' not in  feature or not feature['sequence'])
            and sequence_cor is not None and name in sequence_cor):
            sequence = sequence_cor[name]
        else:
            sequence = None

        if 'Publications' in feature['attributes']:
            pub = feature['attributes']['Publications']
        else:
            pub = None
        # name, alias , type_, sofa_type, sequence, pub
        add_marker(markers, name, alias , original_type, sofa_type, sequence, pub)

def add_marker(markers, name, alias , type_, sofa_type, sequence, pub):
    'It adds a marker to markers'
    if name not in markers:
        markers[name] = {}

    if 'alias' not in markers[name]:
        markers[name]['alias'] = set()
    if alias is not None:
        markers[name]['alias'].add(alias)

    for field, value in (('type', type_), ('sofa', sofa_type),
                  ('sequence', sequence), ('publication', pub), ('name',name)):
        if field not in markers[name] or markers[name][field] is None:

            markers[name][field] = value

def _guess_file_type(infhand):
    'It guesses file type'
    first_line = infhand.readline()
    if first_line.strip() == '##gff-version 3':
        file_type = 'gff3'

    infhand.seek(0)
    return file_type

def parse_markersfile(markersfhand):
    'It parses the markers file'
    markers = {}
    for line in markersfhand:
        line = line.strip()
        if not line or line[0] == '#':
            continue

        marker_name, sofa, type_, alias, sequence, pub = line.split('\t')
        if sofa == '.':
            sofa = None
        if type_ == '.':
            type_ = None
        if alias == '.':
            alias = set()
        else:
            alias = set(alias.split(','))
        if sequence == '.':
            sequence = None
        if pub == '.':
            pub = None
        markers[marker_name] = {'sofa':sofa, 'type':type_,
                                'alias':alias, 'sequence':sequence,
                                'publication':pub, 'name':marker_name}
    return markers

def write_markers(outfhand, markers):
    'It writes the markers to a file'
    outfhand.write('#markers file\n')
    outfhand.write('#name\tsofa\toriginal_type\talias\tsequence\tpublications\n')
    for marker_name in markers:
        if 'sofa' in markers[marker_name]:
            sofa = markers[marker_name]['sofa']
        else:
            sofa = '.'
        if 'type' in markers[marker_name]:
            type_ = markers[marker_name]['type']
        else:
            type_ = '.'
        if ('alias' in markers[marker_name] and markers[marker_name]['alias']):
            alias = ','.join(markers[marker_name]['alias'])
        else:
            alias = '.'
        if ('sequence' in markers[marker_name] and
                                markers[marker_name]['sequence'] is not None):
            sequence = markers[marker_name]['sequence']
        else:
            sequence = '.'
        if ('publication' in markers[marker_name] and
                                markers[marker_name]['publication'] is not None):
            publication = markers[marker_name]['publication']
        else:
            publication = '.'
        outfhand.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (marker_name, sofa, type_,
                                                  alias, sequence, publication))

if __name__ == '__main__':
    main()
