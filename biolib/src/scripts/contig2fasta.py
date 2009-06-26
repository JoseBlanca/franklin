'''
Created on 2009 mai 13

@author: peio
'''
from optparse import OptionParser
from biolib.contig_io import get_parser_by_name
from biolib.biolib_utils import fasta_str
import sys

def main():
    '''The main function '''
    parser = OptionParser('usage: %prog -i [caffile|acefile] -o fastafile')
    parser.add_option('-o', '--outfile', dest='outfile', help='Output file')
    parser.add_option('-i', '--infile', dest='infile', help='Maximun limit')
    parser.add_option('-c', '--contigs', dest='contig_list', \
                      help='Contig list to print. Comma separated')
    options = parser.parse_args()[0]

    if options.infile is None:
        parser.error('Script at least needs an input file (caf|ace)')
    else:
        infile = options.infile

    if options.outfile is None:
        outfile = sys.stdout
    else:
        outfile = options.outfile
        outfile = open(outfile, 'w')

    parser = get_parser_by_name( infile)

    if options.contig_list is None:
        contig_list = None
    else:
        contig_list = options.contig_list.split(',')

    for contig in parser.contigs():
        name       = contig.consensus.sequence.name
        print_this = False
        if contig_list is  None:
            print_this = True
        elif name in contig_list:
            print_this = True

        if print_this:
            sequence = contig.consensus.sequence
            toprint  = fasta_str(sequence, name)
            outfile.write(toprint)
    outfile.close()

if __name__ == '__main__':
    main()
