'''
This script clean(mask and strip) and filter sequence depending in a lot of
variables in sanger sequences.
'''
from optparse import OptionParser
from itertools import imap, ifilter
import re
from biolib.biolib_utils import checkpoint, seqs_in_file
from biolib.seq_cleaner import (create_vector_striper_by_alignment,
                                create_striper_by_quality_trimpoly,
                                create_striper_by_quality_lucy,
                                create_masker_repeats_by_repeatmasker )
from biolib.seq_filters import create_length_filter

def main():
    'The main function'
    parser = OptionParser('usage: %prog -f fastafile [-q quality_file]')
    parser.add_option('-f', '--fastafile', dest='fastafile',
                      help='Sanger input fasta')
    parser.add_option('-f', '--qualfile', dest='qualfile',
                      help='Quality fasta file')
    options = parser.parse_args()[0]

    # Input checks and file preparations
    if options.fastafile is None:
        parser.error('Script at least needs an input file (fasta)')
    else:
        fhand_seqs = open(options.fastafile, 'r')
    if options.qualfile is None:
        fhand_qual = None
    else:
        fhand_qual = open(options.qualfile, 'r')


    analisis_steps = [{'function':create_vector_striper_by_alignment,
                       'arguments':{'vectors':'UniVec', 'aligner':'blast'},
                       'type': 'cleaner', 'statistics': None,
                       'name': 'Remove sanger vectors'},

                      {'function':create_vector_striper_by_alignment,
                       'arguments':{'vectors':'XXXXX', 'aligner':'exonerate'},
                       'type': 'cleaner','statistics': None,
                       'name': 'Remove our adaptors'},

                      {'function': create_striper_by_quality_lucy,
                        'type':'cleaner', 'statistics': None,
                       'name':'Strip low quality with lucy'},

                      {'function': create_striper_by_quality_trimpoly,
                       'type':'cleaner', 'statistics': None,
                       'name':'Strip low quality with trimpoly'},

                      {'function':create_masker_repeats_by_repeatmasker ,
                       'arguments':{'species':'eudicotyledons'},
                       'type': 'cleaner', 'statistics':None ,
                       'name': 'Mask repeats'},


                      {'function': create_length_filter,
                       'arguments':{'lenght':100, 'count_masked': False},
                       'type':'filter' ,
                       'statistics': None ,
                       'name':'Remove seqs shorter than 100 nt' }
                    ]

    seq_iter = seqs_in_file(fhand_seqs, fhand_qual)

    for analisis_step in analisis_steps:
        function  = analisis_step['function']
        if 'arguments' in analisis_step:
            arguments = analisis_step['arguments']
        else:
            analisis_step = None
        type_     = analisis_step['type']
        name      = analisis_step['name']
        ungapped_name = re.sub(' ', '_', name)

        if type_ == 'cleaner':
            # here we prepare the cleaner with iterator that are going to be
            # executed when the checkpoint is arrived
            if arguments is None:
                cleaner_function = function()
            else:
                cleaner_function = function(**arguments)
            filtered_seqs = imap(seq_iter, cleaner_function)

        if type_ == 'filter':
            filter_function = function(**arguments)
            filtered_seqs = ifilter(seq_iter, filter_function)
        # Now we need to create again the iterator. And this is going to be
        # useful to actually run the previous filter. Unitil now everithing is
        # an iterator mega structure and noting is executed

        # First we create the new seq and quality files
        fhand_seqs_new = open('%s_seq_file.fasta' % ungapped_name, 'wt')
        if fhand_qual is not None:
            fhand_qual_new = open('%s_qual_file.fasta' % ungapped_name, 'wt')
        else:
            fhand_qual_new = None

        seq_iter = checkpoint(filtered_seqs, fhand_seqs_new, fhand_qual_new)
        fhand_seqs = fhand_seqs_new
        fhand_qual = fhand_qual_new






if __name__ == '__main__':
    main()
