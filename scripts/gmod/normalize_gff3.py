#!/usr/bin/env python
'''
    This script correct names and types. It uses a previously created feature
correlation file. It looks in the correlation file if there is information to
modify in the file.

Created on 05/07/2010
@author: peio
'''

from optparse import OptionParser
from centralize_markers import parse_markersfile
from franklin.gff import create_feature, write_gff, get_gff_header
import sys

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--infile', dest='infile', help='infile')
    parser.add_option('-o', '--outfile', dest='outfile', help='outfile')
    parser.add_option('-c', '--correlations', dest='correlations',
                      help='Marker correlations file')

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

    if options.correlations is None:
        parser.error('Correlations file is required')
    else:
        correlations_fhand = open(options.correlations)

    return infhand, outfhand, correlations_fhand

def main():
    'The main part'
    # set parameters
    infhand, outfhand, correlations_fhand = set_parameters()
    run(infhand, outfhand, correlations_fhand)

def run(infhand, outfhand, correlations_fhand):
    'It runs the script'
    # correlations =
    correlations = parse_markersfile(correlations_fhand)

    # get header
    header = get_gff_header(infhand)

    # fix features
    features = fix_features(infhand, correlations)

    # write_features
    write_gff(features, outfhand, header)

def fix_features(fhand, correlations):
    'It fix features'
    #get features
    for line in fhand:
        line = line.strip()
        if not line:
            continue
        if line.startswith('##'):
            yield line
        else:
            feature = create_feature(line, 3)
            feature = remove_attributes(feature)
            yield correct_feature(feature, correlations)

def _get_marker_from_correlations(name, correlations=None):
    'It gets the marker from correlations. It look in the name and alias field'
    if correlations is None:
        return None
    for marker in correlations.values():
        if name == marker['name'] or name in marker['alias']:
            return marker
    return None

def correct_feature(feature, correlations=None):
    'It correct the feature giving a correlation file'
    marker =  _get_marker_from_correlations(feature['name'], correlations)
    if marker is not None:
        feature['name'] = marker['name']
        feature['id']   = marker['name']
        feature['type'] = marker['sofa']
        feature['attributes']['original_type'] = marker['type']
        if 'publication' not in feature['attributes'] and marker['publication'] is not None:
            feature['attributes']['publication'] = marker['publication']

    else:
        feature['name'] = feature['name'].lower()
        feature['id'] = feature['id'].lower()

    return feature

def remove_attributes(feature):
    'It removes the feature tag from attributes'
    tags = ['Contig_hit', 'BAC', 'Sequence', 'Marker_hit']
    for tag in tags:
        if tag in feature['attributes']:
            del feature['attributes'][tag]
    return feature

def test():
    'It tests the script'
    from StringIO import StringIO
    gff3	=	'''deleu_IX	CMap	snp	0	0	.	.	.	ID=mc088_0;Name=mc088;est=MRGH21;publication=van Leeuwen et al 2005
deleu_IX	CMap	rflp	81	81	.	.	.	ID=mc092_0;Name=mc092;est=MC092;publication=Oliver et al 2001
##cmap lalalal
deleu_IX	CMap	microsatellite	134	134	.	.	.	ID=cmtc47_0;Name=cmxth4;publication=Danin-Poleg et al 2001
'''

    correlations = '''#markers	file
#name	sofa	original_type	alias	sequence	publications
mc88	RFLP_fragment	placementmarker	mc088	.	.
cmxth4	SNP	frameworkmarker	.	.	.
'''
    infhand = StringIO(gff3)
    outfhand = StringIO()
    correlations_fhand = StringIO(correlations)
    run(infhand, outfhand, correlations_fhand)
    res = outfhand.getvalue()
    assert 'publication=van%20Leeuwen%20et%20al%202005;' in res
    assert 'original_type=placementmarker;est=MRGH21;ID=mc88;Name=mc88' in res
    assert '##cmap lalalal' in  res


if	__name__	==	'__main__':
    test()
    main()
