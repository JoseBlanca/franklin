'''This module holds utilities to filter sequences.

The filtering can be divided in two kinds of functions, the ones that return
a bool meaning if the sequence should be removed or not and the ones that
return a masked sequence or trimmed sequence. In this latter case the filter
can also return None if no sequence is left after the filtering process.
'''

from itertools import imap
from biolib.biolib_utils import call, temp_fasta_file
from biolib.blast_summary import BlastSummary
from Bio.Blast import NCBIXML
from StringIO import StringIO

def ifiltering_map(func, *iterators):
    """Version of imap that only yield the True items."""
    for result in imap(func, *iterators):
        if result:
            yield result

class BlastRunner(object):
    'A blast runner for sequence objects'
    def __init__(self, database, program, expect, nhits=20,
                 use_megablast=False):
        '''The init.
        
        keyword arguments
        database -- The blast database name
        program  -- The blast program (e.g blastn or blastp)
        expect   -- The expect to filter the result
        nhits    -- number of results to keep (default 20)
        use_megablast -- default(False)
        '''
        #pylint: disable-msg=R0913
        #blast requires this arguments, and maybe some more

        #why don't we just use a run_blast function?
        #because this interface would also work with if we already have a 
        #multiblast result file and we want to access it in a random way.
        self._database = database
        self._program  = program
        self._expect   = expect
        self._nhits    = nhits
        self._use_megablast = use_megablast

    def get_blast(self, sequence):
        'Given a sequence it returns the xml blast result as a string'
        #we create the fasta file
        fastah = temp_fasta_file(sequence)
        #we run the blast
        nhits = str(self._nhits)
        cmd = ['blastall', '-i', fastah.name, '-p', self._program,
               '-e', str(self._expect), '-m', '7', '-v', nhits, '-b', nhits,
               '-d', self._database]
        if self._use_megablast:
            cmd.append('-n')
            cmd.append('T')
        stdout, stderr, retcode = call(cmd)
        if retcode:
            raise RuntimeError('Problem running blastall: ' + stderr)
        fastah.close()
        return stdout

def create_blast_filter(expect, database, program, keep_better_hits=True,
                        use_megablast=True):
    '''A function factory factory that creates blast filters.
    
    It returns a function that will accept a sequence and it will return
    True or False depending on the blast outcome.
    database is the blast database name.
    program is the blast program (e.g. blastn or blastp)
    '''
    xml_blast_source = BlastRunner(database=database, program=program,
                                    expect=expect, use_megablast=use_megablast)
    def blast_filter(sequence):
        'Given a sequence it returns True or False depending on the blast'
        #first we need the xml blast result
        xml_blast = StringIO(xml_blast_source.get_blast(sequence))
        #now we want a summary
        biopython_blast = NCBIXML.parse(xml_blast)
        summary = BlastSummary(biopython_blast.next())
        #we filter it
        summary.filter_expect_threshold(expect)
        #is the filter positive or not?
        if len(summary.hits):
            result = True
        else:
            result = False
        if not keep_better_hits:
            result = not result
        return result
    return blast_filter
