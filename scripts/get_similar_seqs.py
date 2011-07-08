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

    return out_fhand, blast_fhand
def similar_sequences_for_blast(blast_fhand, filters=None):
    'It look fro similar sequences ina blast result'
    #now we parse the blast
    blast_parser = get_alignment_parser('blast+')
    blast_result = blast_parser(blast_fhand)

    # We filter the results with appropiate  filters
    if filters is None:
        filters = [{'kind'     : 'score_threshold',
                    'score_key': 'similarity',
                    'min_score': 90,
                   },
                   {'kind'            : 'min_length',
                    'min_num_residues': 100,
                    'length_in_query' : True
                   }
                  ]
    alignments = filter_alignments(blast_result, config=filters)
    for alignment in alignments:
        query_name = alignment['query'].name
        for match in alignment['matches']:
            #to which sequence our query is similar?
            name = match['subject'].name
            subj_desc = match['subject'].description
            yield{'name':name,
                  'subject_description':subj_desc,
                  'query_name':query_name,
                  'subject_start': match['subject_start'],
                  'subject_end':   match['subject_end'],
                  'query_start':   match['start'],
                  'query_end':     match['end']}

def main():
    'The main part'
    out_fhand, blast_fhand = set_parameters()

    header  = 'query\tsubject\tdescription\tquery_start\tquery_end\tsubject_start\t'
    header += 'subject_end\n'
    out_fhand.write(header)


    filters = None

    for similar_seq in similar_sequences_for_blast(blast_fhand, filters=filters):
        query_name    = similar_seq['query_name']
        subject_name  = similar_seq['name']
        subject_description = similar_seq['subject_description']
        query_start   = similar_seq['query_start']
        query_end     = similar_seq['query_end']
        subject_start = similar_seq['subject_start']
        subject_end   = similar_seq['subject_end']
        out_fhand.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (query_name, subject_name,
                                                     subject_description,
                                                     query_start, query_end,
                                                     subject_start,
                                                     subject_end))
    out_fhand.close()


if __name__ == '__main__':
    main()
