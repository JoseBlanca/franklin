#!/usr/bin/env python
'''
Created on 2011 uzt 7

@author: peio
'''

from optparse import OptionParser
import sys
from franklin.seq.alignment_result import (filter_alignments,
                                           get_alignment_parser)

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-o', '--outfile', dest='outfile',
                      help='output file')
    parser.add_option('-b', '--blast_result', dest='blast',
                      help='blast results file')
    parser.add_option('-a', '--array_filter', dest='array_filters',
                      action='store_true', help='Use array filters')
    return parser

def set_parameters():
    'Set parameters'
    # Set parameters
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.outfile is None:
        out_fhand = sys.stdout
    else:
        out_fhand = open(options.outfile)

    if options.blast is None:
        parser.error('blast output file required')
    else:
        blast_fhand = open(options.blast)

    array_filters = options.array_filters
    return out_fhand, blast_fhand, array_filters

def similar_sequences_for_blast(blast_fhand, filters):
    'It look fro similar sequences ina blast result'
    #now we parse the blast
    blast_parser = get_alignment_parser('blast+')
    blast_result = blast_parser(blast_fhand)

    alignments = filter_alignments(blast_result, config=filters)
    for alignment in alignments:
        query_name = alignment['query'].name
        for match in alignment['matches']:
            #to which sequence our query is similar?
            name = match['subject'].name
            subj_desc = match['subject'].description

            if 'expect' in match['scores']:
                evalue = str(match['scores']['expect'])
            else:
                evalue = None
            if 'identity'in match['scores']:
                identity = str(match['scores']['identity'])
            else:
                identity = None
            if 'similarity' in match['scores']:
                similarity = str(match['scores']['similarity'])
            else:
                similarity = None

            yield{'name':name,
                  'subject_description':subj_desc,
                  'query_name':query_name,
                  'subject_start': match['subject_start'],
                  'subject_end':   match['subject_end'],
                  'query_start':   match['start'],
                  'query_end':     match['end'],
                  'evalue':        evalue,
                  'identity':      identity,
                  'similarity':    similarity
                  }

def main():
    'The main part'
    out_fhand, blast_fhand, array_filters = set_parameters()

    header  = 'query\tsubject\tdescription\tquery_start\tquery_end\tsubject_start\t'
    header += 'subject_end\tevalue\n'
    out_fhand.write(header)
    if array_filters:
        filters = [{'kind'     : 'score_threshold',
                    'score_key': 'similarity',
                    'min_score': 100}, #95
                   {'kind'            : 'min_length',
                    'min_percentage': 80,
                    'length_in_query' : True }]
    else:
        filters =  [{'kind'     : 'score_threshold',
                     'score_key': 'similarity',
                     'min_score': 90},
                    {'kind'            : 'min_length',
                     'min_num_residues': 200,
                     'length_in_query' : True }]

    for similar_seq in similar_sequences_for_blast(blast_fhand, filters=filters):
        query_name    = similar_seq['query_name']
        subject_name  = similar_seq['name']
        subject_description = similar_seq['subject_description']
        query_start   = similar_seq['query_start']
        query_end     = similar_seq['query_end']
        subject_start = similar_seq['subject_start']
        subject_end   = similar_seq['subject_end']
        evalue = similar_seq['evalue'] if similar_seq['evalue'] else ''
        #similarity = similar_seq['similarity'] if similar_seq['similarity'] else ''
        #identity = similar_seq['identity'] if similar_seq['identity'] else ''
        out_fhand.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (query_name, subject_name,
                                                     subject_description,
                                                     query_start, query_end,
                                                     subject_start,
                                                     subject_end,
                                                     evalue))
    out_fhand.close()


if __name__ == '__main__':
    main()
