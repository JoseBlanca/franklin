'''
This module is part of ngs_backbone. This module provide mapping related
analyses

Created on 15/03/2010

@author: peio
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

import os, shutil
from tempfile import NamedTemporaryFile

from franklin.backbone.analysis import (Analyzer, scrape_info_from_fname,
                                        _LastAnalysisAnalyzer)
from franklin.mapping import map_reads
from franklin.utils.misc_utils import NamedTemporaryDir, VersionedPath
from franklin.backbone.specifications import BACKBONE_BASENAMES
from franklin.sam import (bam2sam, add_header_and_tags_to_sam, merge_sam,
                          sam2bam, sort_bam_sam, standardize_sam, realign_bam,
                          bam_distribs)

class SetAssemblyAsReferenceAnalyzer(Analyzer):
    'It sets the reference assembly as mapping reference'
    def run(self):
        '''It runs the analysis.'''
        contigs_path = self._get_input_fpaths()['contigs']
        contigs_ext = contigs_path.extension
        reference_dir = self._create_output_dirs()['result']
        reference_fpath = os.path.join(reference_dir,
                          BACKBONE_BASENAMES['mapping_reference'] + '.' + \
                                                                    contigs_ext)
        if os.path.exists(reference_fpath):
            os.remove(reference_fpath)
        os.symlink(contigs_path.last_version, reference_fpath)

def _get_basename(fpath):
    'It returns the base name without path and extension'
    return os.path.splitext(os.path.basename(fpath))[0]

class MappingAnalyzer(Analyzer):
    'It performs the mapping of the sequences to the reference'
    def run(self):
        '''It runs the analysis.'''
        self._log({'analysis_started':True})
        settings = self._project_settings['Mappers']
        inputs = self._get_input_fpaths()
        reads_fpaths = inputs['reads']
        reference_fpath = inputs['reference']
        output_dir = self._create_output_dirs(timestamped=True)['result']

        #memory for the java programs
        java_mem = self._project_settings['Other_settings']['java_memory']
        picard_path = self._project_settings['Other_settings']['picard_path']

        for read_fpath in reads_fpaths:
            read_info = scrape_info_from_fname(read_fpath)
            platform = read_info['pl']
            #which maper are we using for this platform
            mapper = settings['mapper_for_%s' % platform]
            out_bam_fpath = os.path.join(output_dir,
                                         read_fpath.basename + '.bam')
            mapping_parameters = {}
            if platform in ('454', 'sanger'):
                mapping_parameters['reads_length'] = 'long'
            else:
                mapping_parameters['reads_length'] = 'short'
            if not os.path.exists(out_bam_fpath):
                map_reads(mapper,
                          reads_fpath=read_fpath.last_version,
                          reference_fpath=reference_fpath.last_version,
                          out_bam_fpath=out_bam_fpath,
                          parameters=mapping_parameters,
                          threads=self.threads,
                          java_conf={'java_memory':java_mem,
                                     'picard_path':picard_path})

        # Now we run the select _last mapping
        self._spawn_analysis(PRIVATE_DEFINITIONS['select_last_mapping'])

        self._log({'analysis_finished':True})

class MergeBamAnalyzer(Analyzer):
    'It performs the merge of various bams into only one'
    def run(self):
        '''It runs the analysis.'''
        self._log({'analysis_started':True})
        settings = self._project_settings
        project_path = settings['General_settings']['project_path']
        os.chdir(project_path)
        inputs = self._get_input_fpaths()
        bam_paths = inputs['bams']
        reference_path = inputs['reference']

        output_dir = self._create_output_dirs()['result']
        merged_bam_path = VersionedPath(os.path.join(output_dir,
                                        BACKBONE_BASENAMES['merged_bam']))

        merged_bam_fpath = merged_bam_path.next_version

        #Do we have to add the default qualities to the sam file?
        #do we have characters different from ACTGN?
        add_qualities = settings['Sam_processing']['add_default_qualities']
        #memory for the java programs
        java_mem = settings['Other_settings']['java_memory']
        picard_path = settings['Other_settings']['picard_path']

        if add_qualities:
            default_sanger_quality = settings['Other_settings']['default_sanger_quality']
            default_sanger_quality = int(default_sanger_quality)
        else:
            default_sanger_quality = None

        temp_dir = NamedTemporaryDir()
        for bam_path in bam_paths:
            bam_basename = bam_path.basename
            temp_sam = NamedTemporaryFile(prefix='%s.' % bam_basename,
                                          suffix='.sam')
            sam_fpath = os.path.join(temp_dir.name, bam_basename + '.sam')
            bam2sam(bam_path.last_version, temp_sam.name)
            sam_fhand = open(sam_fpath, 'w')
            # First we need to create the sam with added tags and headers
            add_header_and_tags_to_sam(temp_sam, sam_fhand)
            temp_sam.close()
            sam_fhand.close()
            #the standardization
            temp_sam2 = NamedTemporaryFile(prefix='%s.' % bam_basename,
                                           suffix='.sam', delete=False)
            standardize_sam(open(sam_fhand.name), temp_sam2,
                            default_sanger_quality,
                            add_def_qual=add_qualities,
                            only_std_char=True)
            temp_sam2.flush()
            shutil.move(temp_sam2.name, sam_fhand.name)

            temp_sam2.close()

        get_sam_fpaths = lambda dir_: [os.path.join(dir_, fname) for fname in os.listdir(dir_) if fname.endswith('.sam')]

        # Once the headers are ready we are going to merge
        sams = get_sam_fpaths(temp_dir.name)
        sams = [open(sam) for sam in sams]

        temp_sam = NamedTemporaryFile(suffix='.sam')
        reference_fhand = open(reference_path.last_version)
        merge_sam(sams, temp_sam, reference_fhand)
        reference_fhand.close()

        # close files
        for sam in sams:
            sam.close()
        # Convert sam into a bam,(Temporary)
        temp_bam = NamedTemporaryFile(suffix='.bam')
        sam2bam(temp_sam.name, temp_bam.name)

        # finally we need to order the bam
        sort_bam_sam(temp_bam.name, merged_bam_fpath,
                     java_conf={'java_memory':java_mem,
                                'picard_path':picard_path})
        temp_bam.close()
        temp_sam.close()
        self._log({'analysis_finished':True})

class RealignBamAnalyzer(Analyzer):
    'It realigns the bam using GATK'
    def run(self):
        '''It runs the analysis.'''
        self._log({'analysis_started':True})
        settings = self._project_settings
        project_path = settings['General_settings']['project_path']
        os.chdir(project_path)
        inputs = self._get_input_fpaths()
        bam_path = inputs['bam']
        bam_fpath = bam_path.last_version
        reference_path = inputs['reference']

        #memory for the java programs
        osettings   = settings['Other_settings']
        java_mem    = osettings['java_memory']
        picard_path = osettings['picard_path']
        gatk_path = osettings['gatk_path']

        #we need a temporary path
        temp_bam = NamedTemporaryFile(suffix='.bam')
        temp_bam_fpath = temp_bam.name
        temp_bam.close()

        #do the realigment
        realign_bam(bam_fpath=bam_fpath,
                    reference_fpath=reference_path.last_version,
                    out_bam_fpath=temp_bam_fpath,
                    java_conf={'java_memory':java_mem,
                               'picard_path':picard_path,
                               'gatk_path':gatk_path},
                    threads=self.threads)
        #a new version for the original bam
        out_bam_fpath = bam_path.next_version
        shutil.move(temp_bam_fpath, out_bam_fpath)
        self._log({'analysis_finished':True})

class BamStatsAnalyzer(Analyzer):
    'It makes the stats of the mapping'

    def run(self):
        '''It runs the analysis.'''
        self._log({'analysis_started':True})
        settings = self._project_settings
        self._create_output_dirs()['result']
        project_path = settings['General_settings']['project_path']
        project_name = settings['General_settings']['project_name']
        sample_size = settings['Sam_stats']['sampling_size']
        os.chdir(project_path)
        inputs = self._get_input_fpaths()
        bam_path = inputs['bam']
        bam_fpath = bam_path.last_version
        bam_fhand = open(bam_fpath)

        out_dir = os.path.abspath(self._get_output_dirs()['result'])
        for kind in ('coverage', 'mapq'):
            basename = os.path.join(out_dir, "%s" % (project_name))
            bam_fhand.seek(0)
            bam_distribs(bam_fhand, kind, basename=basename,
                         sample_size=sample_size)
        bam_fhand.close()

DEFINITIONS = {
    'set_assembly_as_reference':
        {'inputs':{
                   'contigs':
                            {'directory': 'assembly_result',
                             'file': 'contigs'},
                   },
         'outputs':{'result':{'directory': 'mapping_reference'}},
         'analyzer': SetAssemblyAsReferenceAnalyzer,
        },
    'mapping':
        {'inputs':{
            'reads':
                {'directory': 'cleaned_reads',
                 'file_kinds': 'sequence_files'},
            'reference':
                {'directory': 'mapping_reference',
                'file': 'mapping_reference'},
            },
         'outputs':{'result':{'directory': 'mappings_by_readgroup'}},
         'analyzer': MappingAnalyzer,
        },
#    'select_last_mapping':
#        {'inputs':{'analyses_dir':{'directory': 'mappings'}},
#         'outputs':{'result':{'directory': 'mapping_result',
#                              'create':False}},
#         'analyzer': _LastAnalysisAnalyzer,
#        },
    'merge_bams':
        {'inputs':{
            'bams':
                {'directory': 'mappings_by_readgroup',
                 'file_kinds': 'bam'},
            'reference':
                {'directory': 'mapping_reference',
                'file': 'mapping_reference'},
            },
         'outputs':{'result':{'directory': 'mapping_result'}},
         'analyzer': MergeBamAnalyzer,
        },
    'realign_bam':
        {'inputs':{
            'bam':
                {'directory': 'mapping_result',
                 'file': 'merged_bam'},
            'reference':
                {'directory': 'mapping_reference',
                'file': 'mapping_reference'},
            },
         'outputs':{'result':{'directory': 'mapping_result'}},
         'analyzer': RealignBamAnalyzer,
        },
    'mapping_stats':
        {'inputs':{
            'bam':
                {'directory': 'mapping_result',
                 'file': 'merged_bam'},
            },
         'outputs':{'result':{'directory': 'mapping_stats'}},
         'analyzer': BamStatsAnalyzer,
        },
}
PRIVATE_DEFINITIONS ={
    'select_last_mapping':
         {'inputs':{'analyses_dir':{'directory': 'mappings'}},
         'outputs':{'result':{'directory': 'mapping_result',
                              'create':False}},
         'analyzer': _LastAnalysisAnalyzer,
         }
}
