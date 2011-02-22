'''
Given an input GFF3 file, the database name and a file with the relationship
between the gff3 IDs and the database accessions for each feature the script
adds the dbxref to the GFF3 file
'''

import sys
from optparse import OptionParser
import unittest
from StringIO import StringIO
from franklin.gff import (features_in_gff, get_gff_header, write_gff,
                          add_dbxref_to_feature)

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--ingff3', dest='ingff',
                      help='Input GFF3 file')
    parser.add_option('-o', '--outgff3', dest='outgff',
                      help='Output GFF3 file')
    parser.add_option('-r', '--relations', dest='relations',
                      help='Relationships between accessions file')
    parser.add_option('-d', '--database', dest='database',
                      help='Database name')

    return parser

def get_parameters():
    'It reads and fixes the parsed command line options'
    parser  = parse_options()
    options = parser.parse_args()[0]
    params = {}
    if not options.ingff:
        parser.error('An input GFF3 file is required')
    else:
        params['ingff3_fhand'] = open(options.ingff)
    if not options.outgff:
        parser.error('An output GFF3 file is required')
    else:
        params['outgff3_fhand'] = open(options.outgff, 'w')
    if not options.relations:
        parser.error('A file with the relations between accessions is required')
    else:
        params['rels_fhand'] = open(options.ingff)
    if not options.database:
        parser.error('A database name for the dbxref is required')
    else:
        params['database'] = options.database
    return params

def _add_dbxref_to_features(features, dbxref_db, acc_relations):
    'It adds dbxref to feature taking into account the kind of dbxref to add'

    for feature in  features:
        add_dbxref_to_feature(feature, dbxref_db, acc_relations[feature['id']])
        yield feature

def _get_relations(rels_fhand):
    'It returns a dict with the relations between accessions'

    rels_fhand.seek(0)

    acc_relations = {}
    for line_index, line in enumerate(rels_fhand):
        line = line.strip()
        if not line:
            continue
        try:
            acc1, acc2 = line.split()
        except ValueError:
            msg = 'Malformed relations file in line number %i: %s' % \
                                                          (line_index + 1, line)
            raise ValueError(msg)
        acc_relations[acc1] = acc2
    return acc_relations

def add_dbxref_to_gff3(ingff3_fhand, outgff3_fhand, database, rels_fhand):
    'It adds a new dbxref to a GFF3 file'

    header = get_gff_header(ingff3_fhand)

    features = features_in_gff(ingff3_fhand, 3)

    acc_relations = _get_relations(rels_fhand)

    features = _add_dbxref_to_features(features, database, acc_relations)

    write_gff(features, outgff3_fhand, header)

def main():
    'It runs the script'
    params = get_parameters()
    add_dbxref_to_gff3(**params)


class AddDbxrefToGff3Test(unittest.TestCase):
    'It tests the addition of dbxrefs to gff3'
    def test_dbxref_addition(self):
        'It adds a dbxref to a gff3 without'

        #no previous dbxref
        in_gff3 = '''##gff-version 3
##sequence-region ctg123 1 1497228
ctg123\t.\tgene\t1000\t9000\t.\t.\t.\tID=gene00001;Name=EDEN\n'''
        ingff3_fhand = StringIO(in_gff3)
        outgff3_fhand = StringIO()
        database = 'database'
        relations = 'gene00001\tacc1\n'
        rels_fhand = StringIO(relations)
        add_dbxref_to_gff3(ingff3_fhand=ingff3_fhand,
                           outgff3_fhand=outgff3_fhand,
                           database=database,
                           rels_fhand=rels_fhand)
        assert 'Dbxref=database:acc1' in outgff3_fhand.getvalue()

        #we add a new dbxref
        ingff3_fhand = StringIO(outgff3_fhand.getvalue())
        outgff3_fhand = StringIO()
        add_dbxref_to_gff3(ingff3_fhand=ingff3_fhand,
                           outgff3_fhand=outgff3_fhand,
                           database='database2',
                           rels_fhand=rels_fhand)
        assert 'database:acc1' in outgff3_fhand.getvalue()
        assert 'database2:acc1' in outgff3_fhand.getvalue()

        #we add a dbxref that was already there
        ingff3_fhand = StringIO(outgff3_fhand.getvalue())
        outgff3_fhand = StringIO()
        add_dbxref_to_gff3(ingff3_fhand=ingff3_fhand,
                           outgff3_fhand=outgff3_fhand,
                           database='database2',
                           rels_fhand=rels_fhand)
        result = outgff3_fhand.getvalue()
        assert 'database:acc1' in result
        assert 'database2:acc1' in result
        assert result.count('database2:acc1') == 1

        #test error malformed line
        relations = 'gene00001\tacc1\n1\t2\t3\t4\n'
        rels_fhand = StringIO(relations)
        ingff3_fhand = StringIO(in_gff3)
        outgff3_fhand = StringIO()
        try:
            add_dbxref_to_gff3(ingff3_fhand=ingff3_fhand,
                       outgff3_fhand=outgff3_fhand,
                       database='database2',
                       rels_fhand=rels_fhand)
            self.fail('ValueError expected')
        except ValueError:
            pass

def _test():
    'It tests the script'
    unittest.main()

if __name__ == '__main__':
    if len(sys.argv) == 2 and sys.argv[1].lower() == 'test':
        sys.argv.pop(1)
        _test()
    else:
        main()
