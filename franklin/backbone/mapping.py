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
from franklin.utils.cmd_utils import call
from franklin.utils.misc_utils import (NamedTemporaryDir, VersionedPath,
                                       rel_symlink)
from franklin.backbone.specifications import (BACKBONE_BASENAMES,
                                              PLOT_FILE_FORMAT)
from franklin.sam import (bam2sam, add_header_and_tags_to_sam, merge_sam,
                          sam2bam, sort_bam_sam, standardize_sam, realign_bam,
                          bam_distribs, create_bam_index, bam_general_stats)

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
        rel_symlink(contigs_path.last_version, reference_fpath)

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
        output_dir = self._create_output_dirs(timestamped=True)['result']

        # define color and sequence references
        reference_path    = inputs['reference']
        mapping_index_dir = inputs['mapping_index']
        #print reference_path, mapping_index_dir

        #memory for the java programs
        java_mem = self._project_settings['Other_settings']['java_memory']
        picard_path = self._project_settings['Other_settings']['picard_path']

        for read_fpath in reads_fpaths:
            mapping_parameters = {}
            read_info = scrape_info_from_fname(read_fpath)
            platform = read_info['pl']
            #which maper are we using for this platform
            mapper = settings['mapper_for_%s' % platform]

            (reference_fpath,
             color_space) = self._prepare_mapper_index(mapping_index_dir,
                                                      reference_path,
                                                      platform, mapper)
            mapping_parameters['colorspace'] = color_space
            out_bam_fpath = os.path.join(output_dir,
                                         read_fpath.basename + '.bam')

            if platform in ('454', 'sanger'):
                mapping_parameters['reads_length'] = 'long'
            else:
                mapping_parameters['reads_length'] = 'short'

            if not os.path.exists(out_bam_fpath):
                mapping_parameters['threads']   = self.threads
                mapping_parameters['java_conf'] = {'java_memory':java_mem,
                                                   'picard_path':picard_path}
                map_reads(mapper,
                          reads_fpath=read_fpath.last_version,
                          reference_fpath=reference_fpath,
                          out_bam_fpath=out_bam_fpath,
                          parameters=mapping_parameters)

        # Now we run the select _last mapping
        self._spawn_analysis(DEFINITIONS['_select_last_mapping'],
                             silent=self._silent)

        self._log({'analysis_finished':True})

    def _prepare_mapper_index(self, mapping_index_dir, reference_path, platform,
                              mapper):
        'It creates reference_fpath depending on the mapper and the platform'
        kind        = 'color' if platform == 'solid' else 'sequence'
        color_space = True if kind == 'color' else False
        mapping_index_dir = mapping_index_dir[0].original_path
        index_dir   = mapping_index_dir % (mapper, kind)
        if not os.path.exists(index_dir):
            os.mkdir(index_dir)

        reference_fpath = reference_path.last_version
        reference_fname = os.path.basename(reference_fpath)
        index_fpath     = os.path.join(index_dir, reference_fname)
        if not os.path.exists(index_fpath):
            rel_symlink(reference_fpath, index_fpath)
        return index_fpath, color_space

class MergeBamAnalyzer(Analyzer):
    'It performs the merge of various bams into only one'
    def run(self):
        '''It runs the analysis.'''
        self._log({'analysis_started':True})
        settings = self._project_settings
        project_path = settings['General_settings']['project_path']
        tmp_dir      = settings['General_settings']['tmpdir']
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
        try:
            merge_sam(sams, temp_sam, reference_fhand)
        except Exception:
            if os.path.exists(merged_bam_fpath):
                os.remove(merged_bam_fpath)
            raise
        reference_fhand.close()

        # close files
        for sam in sams:
            sam.close()
        # Convert sam into a bam,(Temporary)
        temp_bam = NamedTemporaryFile(suffix='.bam')
        sam2bam(temp_sam.name, temp_bam.name)

        # finally we need to order the bam
        #print 'unsorted.bam', temp_bam.name
        #raw_input()
        sort_bam_sam(temp_bam.name, merged_bam_fpath,
                     java_conf={'java_memory':java_mem,
                                'picard_path':picard_path}, tmp_dir=tmp_dir )
        temp_bam.close()
        temp_sam.close()
        create_bam_index(merged_bam_fpath)

        self._log({'analysis_finished':True})

class CalmdBamAnalyzer(Analyzer):
    'It runs samtools calmd '
    def run(self):
        '''It runs the analysis.'''
        self._log({'analysis_started':True})
        inputs = self._get_input_fpaths()
        bam_path = inputs['bam']
        bam_fpath = bam_path.last_version
        reference_fpath = inputs['reference'].last_version
        out_fhand = open(bam_path.next_version, 'w')
        cmd = ['samtools', 'calmd', '-Abr', bam_fpath, reference_fpath]
        call(cmd, raise_on_error=True, stdout=out_fhand)
        out_fhand.close()
        self._log({'analysis_finished':True})


class RealignBamAnalyzer(Analyzer):
    'It realigns the bam using GATK'
    def run(self):
        '''It runs the analysis.'''
        self._log({'analysis_started':True})
        settings = self._project_settings
        project_path = settings['General_settings']['project_path']
        tmp_dir      = settings['General_settings']['tmpdir']
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
                    threads=self.threads,
                    tmp_dir=tmp_dir)
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
        summary_fname = os.path.join(out_dir,
                                     BACKBONE_BASENAMES['statistics_file'])
        summary_fhand = open(summary_fname, 'w')

        #The general statistics
        bam_general_stats(bam_fhand, summary_fhand)

        for kind in ('coverage', 'mapq'):
            basename = os.path.join(out_dir, "%s" % (project_name))
            bam_fhand.seek(0)
            bam_distribs(bam_fhand, kind, basename=basename,
                         sample_size=sample_size, summary_fhand=summary_fhand,
                         plot_file_format=PLOT_FILE_FORMAT)
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
            'mapping_index':
                {'directory': 'mapping_index'},
            },
         'outputs':{'result':{'directory': 'mappings_by_readgroup'}},
         'analyzer': MappingAnalyzer,
        },
    '_select_last_mapping':
        {'inputs':{'analyses_dir':{'directory': 'mappings'}},
         'outputs':{'result':{'directory': 'mapping_result',
                              'create':False}},
         'analyzer': _LastAnalysisAnalyzer,
        },
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
    'calmd_bam':
        {'inputs':{
            'bam':
                {'directory': 'mapping_result',
                 'file': 'merged_bam'},
            'reference':
                {'directory': 'mapping_reference',
                'file': 'mapping_reference'},
            },
         'outputs':{'result':{'directory': 'mapping_result'}},
         'analyzer': CalmdBamAnalyzer,
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
