'''
This script provides the tools to deal with markers. It provides tools to read,
write and search the markers data structure

Created on 08/07/2010

@author: peio
'''
from franklin.gff import features_in_gff
import csv

SOFA_TRADUCTOR = {'rflp': 'RFLP_fragment',
                  'snp' : 'SNP',
                  'RFLP': 'RFLP_fragment',
                  'SNP-CAPS': 'SNP',
                  'SNP-Snapshot': 'SNP',
                  'SNP-INDEL':'indel',
                  'ssr':'microsatellite',
                  'SSR':'microsatellite',
                  'EST-SSR':'microsatellite',
                  'aflp':'genetic_marker',
                  'AFLP':'genetic_marker',
                  'isozyme':'genetic_marker',
                  'issr':'genetic_marker',
                  'morphological':'genetic_marker',
                  'Morphological':'genetic_marker',
                  'spelling error':'genetic_marker',
                  'rapd':'genetic_marker',
                  'RAPD':'genetic_marker',
                  'RAPDS':'genetic_marker',
                  'placementmarker': 'biological_region' ,
                  'frameworkmarker': 'biological_region',
                  'marker':'biological_region',
                  'Indel':'indel',
                  'IMA':'genetic_marker',
                  'trait':'genetic_marker',
                  'Trait':'genetic_marker',
                  'sequence':'sequence_feature'


                  }
ACCEPTED_MARKERS = ['genetic_marker', 'frameworkmarker', 'indel', 'marker',
                        'microsatellite', 'placementmarker', 'rflp', 'snp',
                        'isozyme', 'aflp', 'issr', 'morphological', 'rapd',
                        'spelling error']

def get_marker_from_correlations(name, correlations=None):
    'It gets the marker from correlations. It look in the name and alias field'
    if correlations is None:
        return None
    for marker in correlations.values():
        if name == marker['name'] or name in marker['alias']:
            return marker
    return None

def parse_markersfile(markersfhand):
    'It parses the markers file'
    markers = {}
    for line in markersfhand:
        line = line.strip()
        if not line or line[0] == '#':
            continue
        try:
            marker_name, sofa, type_, alias, sequence, pub = line.split('\t')
        except ValueError:
            #print line
            raise ValueError(line + ': Malformed')
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
        if 'sofa' in markers[marker_name] and markers[marker_name]['sofa'] is not None:
            sofa = markers[marker_name]['sofa']
        else:
            sofa = '.'
        if 'type' in markers[marker_name] and markers[marker_name]['type'] is not None:
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

def add_markers_gff3(infhand, markers, name_cor=None, type_cor=None,
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


def add_markers_custom_csv(infhand, markers):
    '''It adds the markes from a custom csv file. It is not a general function
    for all csvs.
    The format of the csv is:
    name       type      pub
    '''

    for marker in csv.reader(infhand, delimiter='\t'):
        alias = marker[0]
        name  = marker[0].lower()
        type_ = marker[1]
        if type_ in SOFA_TRADUCTOR:
            sofa_type = SOFA_TRADUCTOR[type_]
            original_type = type_
        else:
            sofa_type = type_
            original_type = None

        pub = marker[2]
        if not pub:
            pub = None
        sequence = None
        add_marker(markers, name, alias, original_type, sofa_type, sequence,
                   pub)
