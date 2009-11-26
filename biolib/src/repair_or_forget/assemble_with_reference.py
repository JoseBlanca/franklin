#!/usr/bin/env python
'''
With this script you can assemble a bunch of different reads and types  into
contigs.

The script is modular. You can run each of the steps independently.
It needs a configuration file.There is a configuration file in data directory
    Configuration file format:
        1)reference    = (reference sequences, fasta format)
        2) working dir = (working dir)
        3)For each sequence file you need to create a dictionary with this
         information: file-name = (file_name)
                      format    = (file_format)
                      seq_type  = (454|sanger|solexa)
                      aligner   = (aligner program to use)
                      align parameters = (parameters)
'''

from optparse import OptionParser
from biolib.libassemble import (check_and_fix_config, load_seqs_in_bank)

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -c config_file [-s steps]')
    parser.add_option('-c', '--configfile', dest='cfile',
                    help='COnfiguration file')
    parser.add_option('-s', '--steps', dest='steps',
                     help='Steps that you want to run. Comma separated [1,2,3]')
    return parser
def set_parameters():
    '''It sets the parameters for the script.'''
    parser = parse_options()
    options = parser.parse_args()[0]

    if options.cfile is None:
        parser.error('Script needs the configuration file')
    else:
        cfile = options.cfile

    if options.steps is None:
        steps = [1, 2, 3]
    else:
        steps = [int(step) for step in options.steps.split(',')]
    return cfile, steps
def prepare_env():
    from biolib.utils.misc_utils import NamedTemporaryDir
    import os
    import biolib
    data_dir = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

    temp_dir = NamedTemporaryDir()
    ref_fhand  = open(os.path.join(data_dir, 'seq.fasta'), 'r')
    read_fhand = open(os.path.join(data_dir, 'solexa.fastq'), 'r')
    config =  {'work_dir' : temp_dir.get_name(),
               'reference': [{'name'     : 'reference',
                              'seq_fpath': ref_fhand.name}],
                'reads'   : [{'name'     : 'solexa',
                              'seq_fpath': read_fhand.name,
                              'format'   :'fastq'}]}
    return config, temp_dir

def main():
    'The main function'
    #configuration
    cfile, steps  = set_parameters()
    configuration = check_and_fix_config(eval(open(cfile).read()))


    #configuration, temp_dir  = prepare_env()

    # load seqs in AMOS bank
    if 1 in steps:
        load_seqs_in_bank(configuration)



if __name__ == '__main__':
    main()