'''
Command utilities for franklin -- tests

This module provides utilities to run external commands into franklin
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of franklin.
# franklin is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# franklin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with franklin. If not, see <http://www.gnu.org/licenses/>.

import unittest

from franklin.utils.cmd_utils import (_process_parameters, create_runner,
                                      _which_binary, b2gpipe_runner, call)
from franklin.seq.seqs import SeqWithQuality, Seq
from franklin.utils.misc_utils import DATA_DIR
import os, tempfile

class ProcessParametersTest(unittest.TestCase):
    'tests the parameter processing'

    def test_basic_parameter_processing(self):
        'test the most basic parameter processing options'
        #default value
        param_def = {'hola': {'default':'caracola', 'option':'-h'}}
        params = {}
        cmd_params = _process_parameters(params, param_def)
        assert cmd_params == ['-h', 'caracola']

        #an error if a required parameter is not given
        param_def = {'hola': {'required':True}}
        params = {}
        try:
            _process_parameters(params, param_def)
            self.fail()
            #pylint: disable-msg=W0704
        except ValueError:
            pass

        #we can add custom cmd options
        param_def = {'hola': {'default':'caracola', 'option':'-h'}}
        params = {'bin':['-j', 'k']}
        cmd_params = _process_parameters(params, param_def)
        assert cmd_params == ['-j', 'k', '-h', 'caracola']

        #the parameter is a number
        param_def = {'hola': {'default':1, 'option':'-h'}}
        params = {}
        cmd_params = _process_parameters(params, param_def)
        assert cmd_params == ['-h', '1']

        param_def = {'caracola': {'default':1.0, 'option':'-c'}}
        params = {}
        cmd_params = _process_parameters(params, param_def)
        assert cmd_params == ['-c', '1.0']

        param_def = {'caracola': {'default':[1.0, 2.0], 'option':'-c'}}
        params = {}
        cmd_params = _process_parameters(params, param_def)
        assert cmd_params == ['-c', '1.0', '2.0']

        assert _which_binary('est2genome') == '/usr/bin/est2genome'


class RunnerFactorytest(unittest.TestCase):
    'test the creation of external binary runners'
    @staticmethod
    def test_create_blast_runner():
        'We can create a runner class for blast'
        blastpath = os.path.join(DATA_DIR, 'blast')
        run_blast_for_seq = create_runner(tool='blast',
                                    parameters={'database':'arabidopsis_genes',
                                                'program':'blastn'},
                                          environment={'BLASTDB':blastpath})
        seq = 'AACTACGTAGCTATGCTGATGCTAGTCTAGCTAGTCGTAGTCTGATCGTAGTCAGTT'
        seq = Seq(seq)
        seq1 = SeqWithQuality(seq)
        result = run_blast_for_seq(seq1)['blast']
        assert result.read()[0] == '<'
    @staticmethod
    def test_create_mdust_runner():
        'We can create a runner class for mdust'
        run_mdust__for_seq = create_runner(tool='mdust')
        seq  = 'AACTACGTAGCTATGCTGATGCTAGTCTAGAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        seq = Seq(seq)
        seq1 = SeqWithQuality(seq)
        result = run_mdust__for_seq(seq1)['sequence']
        assert result.read()[-10:-1] == 'aaaaaaaaa'
    @staticmethod
    def test_create_lucy_runner():
        'We can create a runner class for lucy'
        fastafile = os.path.join(DATA_DIR, 'seq.fasta')
        run_lucy_for_seq = create_runner(tool='lucy',
                                    parameters={'vector':(fastafile,fastafile)})
        seq  = 'AACTACGTAGCTATGCTGATGCTAGTCTAGAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        qual = [30] * len(seq)
        seq = Seq(seq)
        seq1 = SeqWithQuality(seq, qual=qual)
        seqs = [seq1, seq1]
        result = run_lucy_for_seq(seqs)['sequence']
        assert result[0].read() == ''

    @staticmethod
    def test_run_b2g4pipe():
        'It test the runner of b2g4pipe'

        blast = open(os.path.join(DATA_DIR, 'blast2.xml'))
        fhand, annot_fpath = tempfile.mkstemp()
        os.close(fhand)
        fhand, dat_fpath = tempfile.mkstemp()
        os.close(fhand)
        b2gpipe_runner(blast, annot_fpath, dat_fpath)

        assert os.path.exists(annot_fpath)
        assert os.path.exists(dat_fpath)
        os.remove(annot_fpath)
        os.remove(dat_fpath)

class TestCall(unittest.TestCase):
    'It test call'

    def test_call(self):
        'Test call'
        cmd = ['ls', '/@#@#@#@']

        ## When fails
        #without stdout file. Raise False
        stdout, stderr, retcode = call(cmd)
        assert '/@#@#@#@' in stderr

        try:
            call(cmd, raise_on_error=True)
            self.fail()
        except RuntimeError as error:
            assert '/@#@#@#@' in str(error)


        #with stdout file
        stdout = tempfile.NamedTemporaryFile()
        stderr = tempfile.NamedTemporaryFile()
        stdout_str, stderr_str, retcode = call(cmd, stdout=stdout,
                                               stderr=stderr)
        assert not stdout_str
        assert not stderr_str

        assert '/@#@#@#@' in stderr.read()
        try:
            call(cmd, stdout=stdout, stderr=stderr, raise_on_error=True)
            self.fail()
        except RuntimeError as error:
            assert '/@#@#@#@' in stderr.read()


        ## when it do right
        cmd = ['ls', '/']
        #without stdout file.
        stdout, stderr, retcode = call(cmd)
        assert 'root'in stdout

        #without stdout file.
        stdout = tempfile.NamedTemporaryFile()
        stderr = tempfile.NamedTemporaryFile()
        stdout_str, stderr_str, retcode = call(cmd, stdout=stdout, stderr=stderr)
        assert  stdout_str is None

        assert 'root' in stdout.read()









if __name__ == "__main__":
    import sys;sys.argv = ['', 'TestCall.test_call']
    unittest.main()
