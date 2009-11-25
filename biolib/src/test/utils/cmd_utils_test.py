'''
Command utilities for biolib -- tests

This module provides utilities to run external commands into biolib
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of biolib.
# biolib is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# biolib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with biolib. If not, see <http://www.gnu.org/licenses/>.

import unittest

from biolib.utils.cmd_utils import _process_parameters, create_runner
from biolib.seq.seqs import SeqWithQuality
import biolib
import os

DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

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


class RunnerFactorytest(unittest.TestCase):
    'test the creation of external binary runners'
    @staticmethod
    def test_create_blast_runner():
        'We can create a runner class for blast'
        run_blast_for_seq = create_runner(bin_='blast2', kind='blast',
                                          parameters={'database':'tair7_cdna',
                                                      'program':'blastn'})
        seq = 'AACTACGTAGCTATGCTGATGCTAGTCTAGCTAGTCGTAGTCTGATCGTAGTCAGTT'
        seq1 = SeqWithQuality(seq)
        result = run_blast_for_seq(seq1)[0]
        assert result.read()[0] == '<'
    @staticmethod
    def test_create_mdust_runner():
        'We can create a runner class for mdust'
        run_mdust__for_seq = create_runner(kind='mdust')
        seq  = 'AACTACGTAGCTATGCTGATGCTAGTCTAGAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        seq1 = SeqWithQuality(seq)
        result = run_mdust__for_seq(seq1)[0]
        assert result.read()[-10:-1] == 'aaaaaaaaa'
    @staticmethod
    def test_create_lucy_runner():
        'We can create a runner class for lucy'
        fastafile = os.path.join(DATA_DIR, 'seq.fasta')
        run_lucy_for_seq = create_runner(kind='lucy',
                                    parameters={'vector':(fastafile,fastafile)},
                                    multiseq=True)
        seq  = 'AACTACGTAGCTATGCTGATGCTAGTCTAGAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        qual = [30] * len(seq)
        seq1 = SeqWithQuality(seq, qual=qual)
        seqs = [seq1, seq1]
        result = run_lucy_for_seq(seqs)[0]
        assert result.read() == ''

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()


