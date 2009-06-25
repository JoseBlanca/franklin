'''
Created on 2009 eka 25

@author: peio
'''
from biolib.contig_io import CafParser, AceParser
from optparse import OptionParser
from biolib.biolib_utils import fasta_str

def main():
    ''' Main function where we find snps'''

    parser = OptionParser('usage: %prog -i infile -o outfile')
    parser.add_option('-i', '--infile', dest='infile',
                      help='Input file')
    parser.add_option('-o', '--outfile', dest='outfile',
                      help='Output file')
    (options, args) = parser.parse_args()
    if options.infile is None:
        parser.error('Script at least needs an input file (caf|ace)')
    else:
        infile = options.infile
    if options.outfile is None:
        parser.error('The script needs the output file')
    else:
        outfile = options.outfile

    # File Parsing
    if infile[-3:].lower() == 'ace':
        print "File type: ace"
        parser = AceParser(infile)
    elif infile[-3:].lower() == 'caf':
        print "File type: caf"
        parser = CafParser(infile)

    for contig in parser.contigs():
        seq  = contig.consensus.sequence
        name = contig.consensus.sequence.name
        outfile.write(fasta_str(seq, name))
    outfile.close()


if __name__ == '__main__':
    main()
