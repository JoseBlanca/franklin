'''
Created on 26/04/2011

@author: jose
'''

import unittest
import os.path, subprocess, sys, random
from tempfile import NamedTemporaryFile

import franklin
from franklin.utils.misc_utils import NamedTemporaryDir
from franklin.seq.writers import create_temp_seq_file
from franklin.seq.readers import seqs_in_file
from franklin.seq.seqs import SeqWithQuality, Seq

CLEAN_READS = os.path.join(os.path.split(franklin.__path__[0])[0],
                                         'scripts', 'clean_reads')

def _call(cmd):
    'It runs the command and it returns stdout, stderr and retcode'
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    retcode = process.returncode
    return stdout, stderr, retcode

def _call_python(cmd):
    'It runs a python executable'
    cmd.insert(0, sys.executable)
    return _call(cmd)

def create_random_seq(length, gc=50, qual_range=40):
    'It returns a random sequence'
    nucl_choice = {'at': {0:'A', 1:'T'}, 'gc': {0:'G', 1:'C'}}
    gc = gc / 100.0
    if isinstance(qual_range, int):
        qual_range = [qual_range, qual_range]
    nucls = []
    quals = []
    for index in range(length):
        gc_choice = random.uniform(0, 1)
        gc_choice = 'gc' if gc_choice <= gc else 'at'
        at_choice = round(random.uniform(0, 1))
        nucl = nucl_choice[gc_choice][at_choice]
        qual = random.randint(qual_range[0], qual_range[1])
        nucls.append(nucl)
        quals.append(qual)
    assert len(nucls) == length
    return ''.join(nucls), quals

def create_random_seqwithquality(length, gc=50, qual_range=40):
    'It returns a random seqwithquality'
    seq, qual = create_random_seq(length, gc, qual_range)
    name = list('holacaracola')
    random.shuffle(name)
    name = ''.join(name)
    return SeqWithQuality(Seq(seq), qual=qual, name=name)

class CleanReadsTest(unittest.TestCase):
    'It tests the clean_reads script'

    def test_script_exists(self):
        'clean_reads exists'
        if not os.path.exists(CLEAN_READS):
            self.fail('clean_reads does not exists')

    def test_help(self):
        'The help text is generated'
        #print _call(['ls'])
        stdout = _call_python([CLEAN_READS])[0]
        assert 'SEQ_IN' in stdout

        stdout = _call_python([CLEAN_READS , '-h'])[0]
        assert 'SEQ_IN' in stdout

    def test_sanger(self):
        'It tests the basic sanger cleaning'
        seq1 = create_random_seqwithquality(500, qual_range=55)
        seq2 = create_random_seqwithquality(50, qual_range=15)
        seqs = [seq1 + seq2]
        inseq_fhand, inqual_fhand = create_temp_seq_file(seqs, format='qual')
        outseq_fhand = NamedTemporaryFile()
        outqual_fhand = NamedTemporaryFile()

        #platform is required
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name]
        stderr = _call_python(cmd)[1]
        assert 'required' in stderr

        #a correct platform is required
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', 'hola']
        stderr = _call_python(cmd)[1]
        assert 'choice' in stderr

        #we can clean a sanger sequence with quality
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-q', inqual_fhand.name,
               '-o', outseq_fhand.name, '-u', outqual_fhand.name,
               '-p', 'sanger']
        retcode = _call_python(cmd)[2]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     qual_fhand=open(outqual_fhand.name)))
        assert out_seqs[0].qual[-1] == 55

        #we can clean a sanger sequence without quality
        seq1 = create_random_seqwithquality(500, qual_range=55)
        seqs = [SeqWithQuality(seq1.seq + Seq('NNNNNNNNNNNNNN'), name='Ns')]
        inseq_fhand = create_temp_seq_file(seqs, format='fasta')[0]
        outseq_fhand = NamedTemporaryFile()
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', 'sanger']
        retcode = _call_python(cmd)[2]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name)))
        assert not str(out_seqs[0].seq).lower().endswith('nnnnn')

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'RunnerFactorytest.test_create_lucy_runner']
    unittest.main()
