'''
Created on 26/04/2011

@author: jose
'''

import unittest
import os.path, subprocess, sys
from tempfile import NamedTemporaryFile

import franklin
from franklin.utils.misc_utils import TEST_DATA_DIR
from franklin.seq.writers import create_temp_seq_file
from franklin.seq.readers import seqs_in_file
from franklin.seq.seqs import SeqWithQuality, Seq
from franklin.utils.test_utils import create_random_seqwithquality

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

class CleanReadsTest(unittest.TestCase):
    'It tests the clean_reads script'

    def test_script_exists(self):
        'clean_reads exists'
        if not os.path.exists(CLEAN_READS):
            self.fail('clean_reads does not exists')

    def test_help(self):
        'The help text is generated'

        stdout = _call_python([CLEAN_READS])[0]
        assert 'SEQ_IN' in stdout

        stdout = _call_python([CLEAN_READS , '-h'])[0]
        assert 'SEQ_IN' in stdout

        stdout = _call_python([CLEAN_READS , '--version'])[0]
        assert 'clean_reads' in stdout


    def test_error(self):
        'It tests that we can capture an unexpected error'
        error_fhand = NamedTemporaryFile()
        cmd = [CLEAN_READS, 'testerror', '--error_log', error_fhand.name]
        stdout, stderr = _call_python(cmd)[:-1]
        error_log = open(error_fhand.name).read()
        assert 'An unexpected error happened' in stderr
        assert 'function calls leading' in  error_log

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

        #disable quality trimming and lucy_splice are incompatible
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', 'sanger', '-x', '--lucy_splice', 'splice.fasta']
        stderr = _call_python(cmd)[1]
        assert 'incompatible' in stderr

        #we can clean a sanger sequence with quality
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-q', inqual_fhand.name,
               '-o', outseq_fhand.name, '-u', outqual_fhand.name,
               '-p', 'sanger']
        retcode = _call_python(cmd)[2]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     qual_fhand=open(outqual_fhand.name)))
        assert out_seqs[0].qual[-1] == 55

        #disable quality trimming
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-q', inqual_fhand.name,
               '-o', outseq_fhand.name, '-u', outqual_fhand.name,
               '-p', 'sanger', '-x']
        retcode = _call_python(cmd)[2]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     qual_fhand=open(outqual_fhand.name)))
        assert seqs[0].seq == out_seqs[0].seq

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

    def test_illumina(self):
        'It tests the Illumina cleaning'
        seq1 = create_random_seqwithquality(50, qual_range=35)
        seq2 = create_random_seqwithquality(10, qual_range=15)
        seqs = [seq1 + seq2]
        inseq_fhand = create_temp_seq_file(seqs, format='fastq')[0]
        outseq_fhand = NamedTemporaryFile()
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', 'illumina', '-f', 'fastq']
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     format='fastq'))
        assert out_seqs[0].qual[-2] == 35

        #disable quality trimming
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', 'illumina', '-f', 'fastq', '-x']
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     format='fastq'))
        assert seqs[0].seq == out_seqs[0].seq


        #illumina format
        inseq_fhand = create_temp_seq_file(seqs, format='fastq-illumina')[0]
        outseq_fhand = NamedTemporaryFile()
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', 'illumina', '-f', 'fastq-illumina']
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     format='fastq-illumina'))
        assert out_seqs[0].qual[-2] == 35

    def test_solid(self):
        'It tests the solid cleaning'
        #csfasta
        cs_seq  = '''# Thu Jul 15 13:25:41 2010 /share/apps/corona/bin/filter_
# Cwd: /home/pipeline
# Title: solid0065_20100630_FRAG
>2_21_490_F3
T3.23121101332.0133.2221.23.2.2103.330320302..32320
>2_22_386_F3
T3.00222003211.1011.2122.30.0.3210.013012201..20222
>2_22_431_F3
T0.03020122002.2022.2122.21.2.2122.222102322..12221
>8_25_1748_F3
T0..11031202101103031103110303212300122113032213202
'''
        cs_qual = '''# Thu Jul 15 13:25:41 2010 /share/apps/corona/bin/filter_
# Cwd: /home/pipeline
# Title: solid0065_20100630_FRAG
>2_21_490_F3
31 -1 12 24 17 20 29 21 16 18 30 22 24 -1 24 10 26 22 -1 19 26 23 14 -1 27 26 -1 13 -1 20 6 10 11 -1 15 30 19 15 22 4 18 31 4 -1 -1 33 14 9 8 5
>2_22_386_F3
33 -1 23 27 30 24 31 30 32 30 33 14 33 -1 27 18 28 27 -1 31 27 27 30 -1 26 27 -1 17 -1 27 28 26 28 -1 4 17 21 33 14 28 14 17 26 -1 -1 30 30 30 15 7
>2_22_431_F3
29 -1 18 29 4 16 14 19 26 24 16 4 22 -1 21 26 4 30 -1 22 17 6 24 -1 24 32 -1 27 -1 23 17 21 27 -1 5 19 19 6 4 16 6 18 12 -1 -1 14 25 13 12 5
>8_25_1748_F3
31 -1 -1 29 30 30 28 27 28 24 29 24 24 31 20 32 31 18 28 15 32 28 31 29 31 29 32 27 30 29 27 24 31 32 23 27 28 14 30 17 31 20 7 30 29 23 30 8 29 29
'''
        cs_seq_fhand  = NamedTemporaryFile()
        cs_qual_fhand = NamedTemporaryFile()
        cs_seq_fhand.write(cs_seq)
        cs_qual_fhand.write(cs_qual)
        cs_seq_fhand.flush()
        cs_qual_fhand.flush()
        out_fhand = NamedTemporaryFile()
        cmd = [CLEAN_READS, '-i', cs_seq_fhand.name, '-q', cs_qual_fhand.name,
               '-o', out_fhand.name, '-p', 'solid', '-f', 'csfasta',
               '-g', 'fastq']
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(out_fhand.name),
                                     format='fastq'))
        assert not out_seqs

        #we allow more than one missing calls
        cmd = [CLEAN_READS, '-i', cs_seq_fhand.name, '-q', cs_qual_fhand.name,
               '-o', out_fhand.name, '-p', 'solid', '-f', 'csfasta',
               '-g', 'fastq', '--solid_allow_missing_call']
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(out_fhand.name),
                                     format='fastq'))
        assert out_seqs[0].seq.startswith('..1103120210110')
        assert out_seqs[0].qual[2] == 29

        #no quality trimming
        cmd = [CLEAN_READS, '-i', cs_seq_fhand.name, '-q', cs_qual_fhand.name,
               '-o', out_fhand.name, '-p', 'solid', '-f', 'csfasta',
               '-g', 'fastq', '--solid_allow_missing_call', '-x']
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(out_fhand.name),
                                     format='fastq'))
        assert len(out_seqs) == 4
        assert len(out_seqs[0]) == 49

        #double encoding
        #we allow more than one missing calls
        cmd = [CLEAN_READS, '-i', cs_seq_fhand.name, '-q', cs_qual_fhand.name,
               '-o', out_fhand.name, '-p', 'solid', '-f', 'csfasta',
               '-g', 'fastq', '--solid_allow_missing_call', '--double_encoding']
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(out_fhand.name),
                                     format='fastq'))
        assert out_seqs[0].seq.startswith('NNCCATCGAGCACC')

    def test_adaptors(self):
        'It removes adaptors'
        seq1 = create_random_seqwithquality(5, qual_range=35)
        adaptor = create_random_seqwithquality(15, qual_range=35)
        seq2 = create_random_seqwithquality(50, qual_range=35)
        seqs = [seq1 + adaptor + seq2]
        inseq_fhand = create_temp_seq_file(seqs, format='fastq')[0]
        outseq_fhand = NamedTemporaryFile()
        adaptor_fhand = create_temp_seq_file([adaptor], format='fasta')[0]
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', 'illumina', '-f', 'fastq', '-a', adaptor_fhand.name]
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     format='fastq'))
        assert seq2.seq == out_seqs[0].seq


        seq1 = create_random_seqwithquality(5, qual_range=35)
        adaptor = create_random_seqwithquality(15, qual_range=35)
        seq2 = create_random_seqwithquality(50, qual_range=35)
        seqs = [seq1 + adaptor + seq2]
        inseq_fhand = create_temp_seq_file(seqs, format='fastq')[0]
        outseq_fhand = NamedTemporaryFile()
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', '454', '-f', 'fastq', '-a']
        stdout, stderr, retcode = _call_python(cmd)
        assert retcode == 0
        assert  "--adaptors_file: {'454': '" in stdout

        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', 'illumina', '-f', 'fastq', '-a']
        stdout, stderr, retcode = _call_python(cmd)
        assert 'clean_reads does not have default adaptors file' in stderr
        assert  retcode == 14


    def test_vector(self):
        'It removes the vector'
        seq1 = create_random_seqwithquality(5, qual_range=35)
        vector = create_random_seqwithquality(3000, qual_range=35)
        seq2 = create_random_seqwithquality(250, qual_range=35)
        seqs = [seq1 + vector[30:60] + seq2]
        inseq_fhand = create_temp_seq_file(seqs, format='fastq')[0]
        outseq_fhand = NamedTemporaryFile()
        vector_fhand = create_temp_seq_file([vector], format='fasta')[0]
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', '454', '-f', 'fastq', '-v', vector_fhand.name]
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     format='fastq'))
        assert (len(seq2.seq) - len(out_seqs[0].seq)) < 5

    def test_vectordb(self):
        'It removes the vector from a vector database'
        seq1 = create_random_seqwithquality(5, qual_range=35)
        vector = 'CACTATCTCCGACGACGGCGATTTCACCGTTGACCTGATTTCCAGTTGCTACGTCAAGTTC'
        vector = SeqWithQuality(Seq(vector), name='vect', qual=[30]*len(vector))
        seq2 = create_random_seqwithquality(250, qual_range=35)
        seqs = [seq1 + vector + seq2]
        inseq_fhand = create_temp_seq_file(seqs, format='fastq')[0]
        outseq_fhand = NamedTemporaryFile()
        vector_db = os.path.join(TEST_DATA_DIR, 'blast', 'arabidopsis_genes+')
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', '454', '-f', 'fastq', '-d', vector_db]
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     format='fastq'))
        assert (len(seq2.seq) - len(out_seqs[0].seq)) < 5

        seq1 = create_random_seqwithquality(5, qual_range=35)
        vector = 'GGTGCCTCCGGCGGGCCACTCAATGCTTGAGTATACTCACTAGACTTTGCTTCGCAAAG'
        vector = SeqWithQuality(Seq(vector), name='vect', qual=[30]*len(vector))
        seq2 = create_random_seqwithquality(250, qual_range=35)
        seqs = [seq1 + vector + seq2]
        inseq_fhand = create_temp_seq_file(seqs, format='fastq')[0]
        outseq_fhand = NamedTemporaryFile()
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', '454', '-f', 'fastq', '-d']
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     format='fastq'))
        assert (len(seq2.seq) - len(out_seqs[0].seq)) < 5

    def test_words(self):
        'It trims re words'
        vector = 'ACTG'
        vector = SeqWithQuality(Seq(vector), name='vect', qual=[30]*len(vector))
        seq2 = create_random_seqwithquality(250, qual_range=35)
        seqs = [vector + seq2]
        inseq_fhand = create_temp_seq_file(seqs, format='fastq')[0]
        outseq_fhand = NamedTemporaryFile()
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', '454', '-f', 'fastq', '-r', '"^ACTG","TTTTTTTTTTTTTT"']
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     format='fastq'))
        assert seq2.seq == out_seqs[0].seq

    def test_edge_trim(self):
        'It trims the sequence edges'
        seq2 = create_random_seqwithquality(250, qual_range=35)
        seqs = [seq2]
        inseq_fhand = create_temp_seq_file(seqs, format='fastq')[0]
        outseq_fhand = NamedTemporaryFile()
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', '454', '-f', 'fastq', '-e', '10,10']
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     format='fastq'))
        assert len(seq2.seq) - len(out_seqs[0].seq) == 20

    def test_min_length(self):
        'Filtering by length'
        seq1 = create_random_seqwithquality(250, qual_range=35)
        seq2 = create_random_seqwithquality(50, qual_range=35)
        seq3 = create_random_seqwithquality(250, qual_range=35)
        seqs = [seq1, seq2, seq3]
        inseq_fhand = create_temp_seq_file(seqs, format='fastq')[0]
        outseq_fhand = NamedTemporaryFile()
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', '454', '-f', 'fastq', '-m', '51']
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     format='fastq'))
        assert len(out_seqs) == 2
        assert len(out_seqs[0]) == 250
        assert len(out_seqs[1]) == 250

    def test_filter(self):
        'Filtering by blast similarity'
        seq1 = create_random_seqwithquality(150, qual_range=35)
        seq2 = 'CACTATCTCCGACGACGGCGATTTCACCGTTGACCTGATTTCCAGTTGCTACGTCAAGTTCTC'
        seq2 += 'TACGGCAAGAATATCGCCGGAAAACTCAGTTACGGATCTGTTAAAGACGTCCGTGGAATCCA'
        seq2 += 'AGCTAAAGAAGCTTTCCTTTGGCTACCAATCACCGCCATGGAATCGGATCCAAGCTCTGCCA'
        seq2 = SeqWithQuality(Seq(seq2), name='ara', qual=[30]*len(seq2))
        seq3 = create_random_seqwithquality(150, qual_range=35)
        seqs = [seq1, seq2, seq3]
        inseq_fhand = create_temp_seq_file(seqs, format='fastq')[0]
        outseq_fhand = NamedTemporaryFile()
        ara_db = os.path.join(TEST_DATA_DIR, 'blast', 'arabidopsis_genes+')
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', '454', '-f', 'fastq',
               '--filter_dbs', ','.join((ara_db, ara_db))]
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     format='fastq'))
        assert len(out_seqs) == 2

class ParallelTest(unittest.TestCase):
    'It tests the clean_reads script parallel operation'

    def test_fasta_qual(self):
        'Cleaning fasta and qual seqs in parallel'
        seq1 = create_random_seqwithquality(500, qual_range=55)
        seq2 = create_random_seqwithquality(50, qual_range=15)
        seq3 = create_random_seqwithquality(500, qual_range=55)
        seq4 = create_random_seqwithquality(50, qual_range=15)
        seq5 = create_random_seqwithquality(500, qual_range=55)
        seq6 = create_random_seqwithquality(50, qual_range=15)
        seqs = [seq1 + seq2, seq3 + seq4, seq5 + seq6]
        inseq_fhand, inqual_fhand = create_temp_seq_file(seqs, format='qual')
        outseq_fhand = NamedTemporaryFile()
        outqual_fhand = NamedTemporaryFile()
        #we can clean a sanger sequence with quality
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-q', inqual_fhand.name,
               '-o', outseq_fhand.name, '-u', outqual_fhand.name,
               '-p', 'sanger', '-t', '2']
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     qual_fhand=open(outqual_fhand.name)))
        assert out_seqs[0].qual[-1] == 55

    def test_fastq(self):
        'Cleaning fastq seqs in parallel'
        seq1 = create_random_seqwithquality(500, qual_range=55)
        seq2 = create_random_seqwithquality(50, qual_range=15)
        seq3 = create_random_seqwithquality(500, qual_range=55)
        seq4 = create_random_seqwithquality(50, qual_range=15)
        seq5 = create_random_seqwithquality(500, qual_range=55)
        seq6 = create_random_seqwithquality(50, qual_range=15)
        seqs = [seq1 + seq2, seq3 + seq4, seq5 + seq6]
        inseq_fhand = create_temp_seq_file(seqs, format='fastq')[0]
        outseq_fhand = NamedTemporaryFile()
        #we can clean a sanger sequence with quality
        cmd = [CLEAN_READS, '-i', inseq_fhand.name, '-o', outseq_fhand.name,
               '-p', 'sanger', '-t', '4', '-f', 'fastq']
        retcode = _call_python(cmd)[-1]
        assert retcode == 0
        out_seqs = list(seqs_in_file(seq_fhand=open(outseq_fhand.name),
                                     format='fastq'))
        assert out_seqs[0].qual[-1] == 55

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'ParallelTest.test_fasta_qual']
    unittest.main()
