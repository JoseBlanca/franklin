'''
Created on 2009 mai 4

@author: peio
'''


from biolib.SeqVariation import seqvariations_in_alignment, calculate_pic, \
                                cap_enzime,  _allele_count
from biolib.contig_parser import CafParser, AceParser
from biolib.contig_cleaner import contig_strip, water_alignment_strip
from biolib.biolib_utils import get_start_end
from optparse import OptionParser
from os.path import basename


def main():
    ''' Main function where we find snps'''
    
    parser = OptionParser('usage: %prog -i infile [-t]...', version='%prog 1.0')
    parser.add_option('-s', '--soutfile',  dest='soutfile', 
                      help='snp information output file')
    parser.add_option('-c', '--coutfile',  dest='coutfile', 
                      help='contig information output file')
    parser.add_option('-i', '--infile', dest='infile',
                      help='Input file')
    parser.add_option('-t', '--contig_strip', dest = 'c_strip', 
                     action='store_true', help = 'Strip contig?'  )
    (options, args) = parser.parse_args()
  

    if options.infile is None:
        parser.error('Script at least needs an input file (caf|ace)')
    else:
        infile = options.infile
        
    if options.soutfile is None:
        soutfile = "".join(infile.split('.')[:-1]) + '.snp_out'
        soutfile = basename(soutfile)
    else:
        soutfile = options.soutfile
    
    if options.coutfile is None:
        coutfile = "".join(infile.split('.')[:-1]) + '.contig_out'
        coutfile = basename(coutfile)
    else:
        coutfile = options.coutfile
      
    print "Starting file indexing"  
    if infile[-3:].lower() == 'ace':
        print "File type: ace"
        parser = AceParser(infile)
    elif infile[-3:].lower() == 'caf':
        print "File type: caf"
        parser = CafParser(infile)
    print "File indexing finished"
        
    
    soutfileh = open(soutfile, 'wt')
    soutfileh.write('#Id\tLoc_start\tLoc_end\tType\tpic\tcap\n')
    
    coutfileh = open(coutfile, 'wt')
    coutfileh.write('#Id\tsnps in contig\t%snps in contig\n')
    for contig in parser.contigs():
        contig_name = contig.consensus.sequence.name
        if options.c_strip is not None:
            contig = contig_strip(contig, 13) # este numero es aproximado
        contig = water_alignment_strip(contig)
        print ".-Searching in contig: %s" % contig_name
        var_count = 0
        for seqvar in seqvariations_in_alignment(contig):
            sorted_alleles = seqvar.sorted_alleles()
            second_alleles = _allele_count(sorted_alleles[1][1]) 
            
            kind      = seqvar.kind()
            loc_start, loc_end = get_start_end(seqvar.location)
            id_snp      = contig_name + '_' +  str(loc_start)
            pik         = calculate_pic(seqvar)
            enzymes     = cap_enzime(seqvar)
            if enzymes is not None:
                enzyme_list = ",".join(enzymes)
            else:
                enzyme_list = ''
            toprint =  "%s\t%d\t%d\t%s\t%d\t%f\t%s\n" % (id_snp, loc_start,
                                                            loc_end, kind, 
                                                            second_alleles,
                                                            pik, 
                                                            enzyme_list)
            print toprint.strip()
            soutfileh.write(toprint)
            
            var_count += 1
        if var_count == 0:
            percentaje_snp = 0
        else:
            percentaje_snp = (var_count/float(len(contig.consensus))) * 100
            
        toprint = '%s\t%d\t%f\n' % (contig_name, var_count, percentaje_snp)
        coutfileh.write(toprint)
        
    print "Finished"    


if __name__ == '__main__':
    main()