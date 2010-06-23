#!/usr/bin/env python
'''It merges the genetic markers in the gff to chado. This gff is the gff used
populate the gbrowse's phisical map
Inputs: cmap_gff and physical_map
outputs: 3 gff3 to chado. One with the physical map and other with the markers
         not located inthe physical map

'''
from optparse import OptionParser
from franklin.gff import (features_in_gff, get_gff_header,
                          add_dbxref_to_feature, write_gff)
from itertools import tee

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
    parser.add_option('-d', '--dbxref', dest='dbxref',
                      help='Name of the database to make the links')
    return parser

def set_parameters():
    'Set parameters'
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.physical is None:
        parser.error('In file required')
    else:
        physical_fhand = open(options.physical)

    if options.genetic is None:
        parser.error('In file required')
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

    dbxref_db = 'melon_cmap' if options.dbxref is None else options.dbxref

    return physical_fhand, genetic_fhand, outfhand, orphan_fhand, dbxref_db

def main():
    'The main part'
    # get parameters
    physical_fhand, genetic_fhand, outfhand, orphan_fhand, dbxref_db = \
                                                                set_parameters()
    run(physical_fhand, genetic_fhand, outfhand, orphan_fhand, dbxref_db)

def run(physical_fhand, genetic_fhand, outfhand, orphan_fhand, dbxref_db):
    'It runs the script'
    # get gff3 header
    physical_header = get_gff_header(physical_fhand)

    genetic_features = features_in_gff(genetic_fhand, 3)
    # get genetic markers
    genetic_markers = get_indexed_features(genetic_features)

    # insert genetic_markers in physical markers
    merged_features, orphan_features = merge_markers(physical_fhand,
                                                     genetic_markers, dbxref_db)

    # write merged gff
    write_gff(merged_features, outfhand, header=physical_header)

    # add rdbxref to orphan features
    orphan_features = add_dbxref_to_orphans(orphan_features, dbxref_db)

    # prepare orphan features to load in chado. Add a region were all features
    # are located
    orphan_features = prepare_to_chado(orphan_features)


    write_gff(orphan_features, orphan_fhand)

def prepare_to_chado(features):
    '''It prepares the features to load in chado. It adds a fake region were all
     features are located'''
    features, features2 = tee(features, 2)

    max_length = _get_longest_feature_lenght(features2)
    fake_region_name = 'fake_region'

    yield {'seqid':fake_region_name, 'source':'.', 'type':'region',
           'start':1, 'end': max_length + 10, 'score':'.',
           'strand':'.', 'phase':'.', 'attributes':{'ID':fake_region_name,
                                                    'Name':fake_region_name}}
    for feature in features:
        feature['seqid'] = fake_region_name

        feature['end'] = feature['end'] - feature['start'] + 1
        feature['start'] = 1
        yield feature

def _get_longest_feature_lenght(features):
    'It get the longest features length'
    length = 0
    for feature in features:
        new_length = int(feature['end']) - int(feature['start'])
        if new_length > length:
            length = new_length
    return length

def merge_markers(physical_fhand, genetic_markers, dbxref_db):
    '''It merges the markers in the given outfhand. The orpahned ones are
    returned as a list'''
    physical_features = features_in_gff(physical_fhand, version=3)

    merged_features = add_dbxref_to_common_features(physical_features,
                                                    genetic_markers, dbxref_db)
    merged_features, merged_features2 = tee(merged_features, 2)

    orphan_features = get_orphaned_features(merged_features2, genetic_markers)

    return merged_features, orphan_features

def get_orphaned_features(physical_features, genetic_markers):
    'It get orphaned features'
    physycal_feat_names = set([feat['name'] for feat in physical_features])
    genetic_feat_names = set(genetic_markers.keys())
    orphan_feat_names = genetic_feat_names.difference(physycal_feat_names)
    for orphan_feat_name in orphan_feat_names:
        for orphan_feature in genetic_markers[orphan_feat_name]:
            yield orphan_feature


def add_dbxref_to_orphans(features, dbxref_db):
    'It builds the dbxref within the genetic feature information'
    indexed_orphans = get_indexed_features(features)
    for feature_list in indexed_orphans.values():
        dbxref_ids = get_genetic_dbxrefs(feature_list)
        for index, feature in enumerate(feature_list):
            add_dbxref_to_feature(feature, dbxref_db, dbxref_ids[index])
            yield feature

def add_dbxref_to_common_features(physical_features, genetic_markers,
                                  dbxref_db):
    'It yields common genetic_markers'

    genetic_marker_names = genetic_markers.keys()
    for feature in physical_features:
        name = feature['name']
        if name in genetic_marker_names:
            genet_features = genetic_markers[name]

            genetic_dbxref_ids =  get_genetic_dbxrefs(genet_features)
            for genetic_dbxref_id in genetic_dbxref_ids:
                add_dbxref_to_feature(feature, dbxref_db, genetic_dbxref_id)
        yield feature

def get_genetic_dbxrefs(features):
    'It builds the dbxrefs taking into account the genetic_feature_data'
    feature_dbxrefs = []
    if len(features) > 1:
        for index, feature in enumerate(features):
            letter = chr(65 + index)
            feature_dbxref = '%s_%s_%s' % (feature['seqid'], feature['name'],
                                           letter)
            feature_dbxrefs.append(feature_dbxref)
    else:
        feature = features[0]
        feature_dbxref = '%s_%s' % (feature['seqid'], feature['name'])
        feature_dbxrefs.append(feature_dbxref)
    return feature_dbxrefs

def get_indexed_features(features):
    'It indexes the genetic markers'
    markers = {}
    for feature in features:
        name = feature['name']
        if name not in markers:
            markers[name] = []
        markers[name].append(feature)
    return markers

def test():
    'It test the script'
    from StringIO import StringIO
    import os
    from franklin.utils.misc_utils import DATA_DIR

    physical_fhand = open(os.path.join(DATA_DIR, 'map_fis.gff3'))
    genetic_fhand  = open(os.path.join(DATA_DIR, 'map_fis.gff3'))
    outfhand       = StringIO()
    orphan_fhand   = StringIO()
    dbxref = 'cmap_melon'

    run(physical_fhand, genetic_fhand, outfhand, orphan_fhand, dbxref)
    result = outfhand.getvalue()
    assert 'Dbxref=cmap_melon:Chrctg0_Cm13_B04;' in result
    print 'Test OK'


if __name__ == '__main__':
    #test()
    main()




