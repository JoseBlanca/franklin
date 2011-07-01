'''
This module is part the ngs_backbone. It performas statistics related to snvs

Created on 27/06/2011

@author: dani
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

from franklin.backbone.analysis import Analyzer

from franklin.snv.snv_statistics import do_snv_stats

class SnvStatsAnalyzer(Analyzer):
    'It performs the statistic analysis for ngs_backbone'
    def run(self):
        'It runs the analysis.'
        output_dir = self._create_output_dirs()['result']
        inputs = self._get_input_fpaths()
        pickle_paths = inputs['pickle']

        # do analysis
        for seq_path in pickle_paths:
            do_snv_stats(seq_path, output_dir)

DEFINITIONS = {
    'snv_stats':
        {'inputs':{
            'pickle':
                {'directory': 'annotation_dbs',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'result':{'directory': 'snv_stats'}},
         'analyzer': SnvStatsAnalyzer}
    }
