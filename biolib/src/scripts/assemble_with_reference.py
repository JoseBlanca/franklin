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
from biolib.libassemble import check_cfile

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

def main():
    'The main function'

    #configuration
    cfile, steps  = set_parameters()
    configuration = check_cfile(open(cfile))

    # load seqs in AMOS bank
    if 1 in steps:
        load_seqs_in_bank(configuration)





if __name__ == '__main__':
    main()