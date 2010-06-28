#!/usr/bin/env python

'''
It collapses the features by its name.
'''
from optparse import OptionParser
import sys
from franklin.gff import features_in_gff, write_gff, get_gff_header
from franklin.utils.misc_utils import OrderedDict

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
    run(infhand, outfhand)

def run(infhand, outfhand):
    'It runs the script'
    # get header
    header = get_gff_header(infhand)

    # get features
    features = features_in_gff(infhand, version=3)

    # join features
    joined_features_ = join_features(features)

    # write_gff
    write_gff(joined_features_, outfhand, header)

def join_features(features):
    'It joins features with the same'
    # join features by its name
    feature_dict = OrderedDict()
    for feature in features:
        name = feature['name']
        if name not in feature_dict:
            feature_dict[name] = []
        feature_dict[name].append(feature)

    for feature_list in feature_dict.values():
        yield joined_features(feature_list)

def joined_features(feature_list):
    'It joins various features in a list'
    if len(feature_list) == 1:
        new_feature = feature_list[0]
    else:
        new_attributes = {}
        for feature in feature_list:
            for attr_key, attr_value in feature['attributes'].items():
                if attr_key not in new_attributes:
                    new_attributes[attr_key] = attr_value
                elif attr_value != new_attributes[attr_key]:
                    new_attributes[attr_key] +=',%s' % str(attr_value)

        new_feature = feature_list[0]
        new_feature['id'] = new_feature['name']
        new_feature['attributes'] = new_attributes

    return new_feature
def test():
    'It tests the script'
    from StringIO import StringIO
    gff_in	=	'''##gff-version	3
Chrctg0	assembly	chromosome	1	140722177	.	.	.	Dbxref=gbrowse_melon:Chrctg0%253A1..140722177;Name=Chrctg0;ID=Chrctg0
Chrctg0	FPC	contig	1	140722177	.	.	.	contig=ctg0;Dbxref=gbrowse_melon:Chrctg0%253A1..140722177;Name=ctg0;ID=ctg0
Chrctg0	FPC	genetic_marker	52303873	52303873	.	.	.	marker=mc133;Dbxref=gbrowse_melon:Chrctg0%253A52303873..52303873;Name=mc133;ID=mc133
Chr0	FPC	genetic_marker	225281	225281	.	.	.	marker=mc133;Dbxref=gbrowse_melon:Chr0%253A225281..225281;Name=mc133;ID=mc133_2
'''
    gff_out	=	'''##gff-version	3
Chrctg0	assembly	chromosome	1	140722177	.	.	.	Dbxref=gbrowse_melon:Chrctg0%253A1..140722177;Name=Chrctg0;ID=Chrctg0
Chrctg0	FPC	contig	1	140722177	.	.	.	contig=ctg0;Dbxref=gbrowse_melon:Chrctg0%253A1..140722177;Name=ctg0;ID=ctg0
Chrctg0	FPC	genetic_marker	52303873	52303873	.	.	.	marker=mc133;Dbxref=gbrowse_melon:Chrctg0%253A52303873..52303873,gbrowse_melon:Chr0%253A225281..225281;Name=mc133;ID=mc133
'''
    infhand  = StringIO(gff_in)
    outfhand = StringIO()
    run(infhand, outfhand)
    result = outfhand.getvalue()
    assert  result == gff_out

if __name__ == '__main__':
    test()
    main()
