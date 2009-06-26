'''
Created on 2009 eka 26

@author: peio
'''
from optparse import OptionParser
import os
from biolib.contig_io import contig_to_fasta
from biolib.biolib_utils import call

def main():
    '''main section '''
    parser = OptionParser('usage: %prog -i infile [-t]...', version='%prog 1.0')
    parser.add_option('-i', '--infile', dest='infile', help='Input file')
    parser.add_option('-d', '--directory', dest='directory', help='directory')
    (options, args) = parser.parse_args()

    if options.infile is None:
        parser.error('Script at least needs an input file (caf|ace)')
    else:
        infile = options.infile
    if options.directory is None:
        parser.error('Script needs the working directory')
    else:
        directory = options.directory

    if not os.path.exists(directory):
        os.mkdir(directory)

    fasta   = contig_to_fasta(infile)
    fasta_name = os.path.basename(infile)[:-3] + 'fasta'
    fasta_fhand = open(os.path.join(directory, fasta_name))
    fasta_fhand.write(fasta)
    fasta_fhand.close()

    os.chdir(directory)
    cmd = ['formatdb', '-i', fasta_fhand.name, '-V', '-p', 'F', '-o']
    stdout, stderr, retcode = call(cmd)
    print stdout

    cmd = ['blastall', '-p', 'blastn', '-d', fasta_name, '-e' , '0.0001', '-n',
           '-b', '30', '-v', '30']
    stdout, stderr, retcode = call(cmd)
    print stdout



if __name__ == '__main__':
    main()
