'''
Created on 12/03/2010

@author: peio
'''
import shutil, os
from franklin.backbone.analysis import Analyzer
from tempfile import NamedTemporaryFile
from franklin.pipelines.pipelines import seq_pipeline_runner
from franklin.backbone.blast_runner import backbone_blast_runner
from franklin.pipelines.annotation_steps import (annotate_cdna_introns,
                                                 annotate_orthologs,
                                                 annotate_with_descriptions,
                                                 annotate_microsatellites,
                                                 annotate_orfs,
                                                 annotate_gos)
from franklin.sam import create_bam_index
from franklin.pipelines.snv_pipeline_steps import (
                                            unique_contiguous_region_filter,
                                            close_to_intron_filter,
                                            high_variable_region_filter,
                                            close_to_snv_filter,
                                            close_to_limit_filter,
                                            major_allele_freq_filter,
                                            kind_filter,
                                            cap_enzyme_filter,
                                            is_variable_filter)

def _get_basename(fpath):
    'It returns the base name without path and extension'
    return os.path.splitext(os.path.basename(fpath))[0]

class AnnotationAnalyzer(Analyzer):
    'It annotates the introns in cdna sequences'

    def _get_inputs_and_prepare_outputs(self):
        'It creates the output dir and it returns the inputs'
        output_dir   = self._create_output_dirs()['result']
        inputs       = self._get_input_fpaths()
        return inputs, output_dir

    def _run_annotation(self, pipeline, configuration, inputs, output_dir):
        'It runs the analysis.'
        repr_fpaths  = inputs['repr']
        try:
            seqs_fpaths  = inputs['input']
        except KeyError:
            seqs_fpaths = []
        seqs_fpaths  = self._get_seq_or_repr_fpath(seqs_fpaths, repr_fpaths)
        for seq_fpath in seqs_fpaths:
            temp_repr = NamedTemporaryFile(suffix='.repr', mode='a',
                                           delete=False)
            io_fhands = {'in_seq': open(seq_fpath),
                         'outputs':{'repr':temp_repr}}
            if seq_fpath in configuration:
                #there is a different configuration for every file to annotate
                config = configuration[seq_fpath]
            else:
                config = configuration
            seq_pipeline_runner(pipeline, configuration=config,
                                io_fhands=io_fhands)
            temp_repr.close()
            repr_fpath = os.path.join(output_dir,
                                      _get_basename(seq_fpath) + '.repr')
            shutil.move(temp_repr.name, repr_fpath)

    def _get_seq_or_repr_fpath(self, seqs_fpaths, repr_fpaths):
        'It returns for every file the repr or the seq file'
        repr_fpaths_ = dict([(self._get_basename(fpath),
                                               fpath) for fpath in repr_fpaths])
        if not seqs_fpaths:
            return repr_fpaths
        new_seq_fpaths = []
        for fpath in seqs_fpaths:
            basename = _get_basename(fpath)
            if basename in repr_fpaths_:
                new_seq_fpaths.append(repr_fpaths_[basename])
            else:
                new_seq_fpaths.append(fpath)
        return new_seq_fpaths

class AnnotateOrthologsAnalyzer(AnnotationAnalyzer):
    '''It annotates the orthologs for the working species against the species
    given in the configuration file'''

    def run(self):
        'It runs the analysis'
        inputs, output_dir = self._get_inputs_and_prepare_outputs()
        blast_settings = self._project_settings['blast']
        settings = self._project_settings['Annotation']
        ortholog_databases = settings['ortholog_databases']

        #first we need some blasts
        project_dir = self._project_settings['General_settings']['project_path']
        blasts = {}
        for input_ in inputs['input']:
            for database in ortholog_databases:
                db_kind = blast_settings[database]['kind']
                if db_kind == 'nucl':
                    blast_program = 'tblastx'
                else:
                    blast_program = 'blastx'
                blastdb = blast_settings[database]['path']
                #this could be different adding something to the settings
                blastdb_seq_fpath = blastdb
                blast = backbone_blast_runner(query_fpath=input_,
                                              project_dir=project_dir,
                                              blast_program=blast_program,
                                              blast_db=blastdb,
                                              dbtype=db_kind)
                reverse_blast = backbone_blast_runner(
                                              query_fpath=blastdb_seq_fpath,
                                              project_dir=project_dir,
                                              blast_program=blast_program,
                                              blast_db_seq=input_,
                                              dbtype=db_kind)
                if input_ not in blasts:
                    blasts[input_] = {}
                blasts[input_][database] = {'blast':blast,
                                            'reverse_blast':reverse_blast}

        pipeline = []
        configuration = {}
        for database in ortholog_databases:
            step = annotate_orthologs
            step['name_in_config'] = database
            #an annotation step for every ortholog database
            pipeline.append(step)
            for input_ in inputs['input']:
                reverse_blast = ''
                step_config = {
                    'blast':{'blast': blasts[input_][database]['blast']},
                    'reverse_blast':{'blast':
                                     blasts[input_][database]['reverse_blast']},
                    'species': database}
                configuration[input_] = {}
                configuration[input_][database] = step_config

        return self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=output_dir)

class AnnotateIntronsAnalyzer(AnnotationAnalyzer):
    'It annotates the introns in cdna sequences'

    def run(self):
        'It runs the analysis.'
        inputs, output_dir = self._get_inputs_and_prepare_outputs()
        pipeline = [annotate_cdna_introns]

        settings = self._project_settings['Annotation']
        if 'Cdna_intron_annotation' not in settings:
            msg = 'You should set up genomic_db and genomic_seqs in settings'
            raise RuntimeError(msg)
        genomic_db = settings['Cdna_intron_annotation']['genomic_db']
        genomic_seqs = settings['Cdna_intron_annotation']['genomic_seqs']
        configuration = {'annotate_cdna_introns': {'genomic_db':genomic_db,
                                       'genomic_seqs_fhand':open(genomic_seqs)}}

        return self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=output_dir)

class SnvCallerAnalyzer(AnnotationAnalyzer):
    'It performs the calling of the snvs in a bam file'

    def run(self):
        'It runs the analysis.'
        inputs, output_dir = self._get_inputs_and_prepare_outputs()
        merged_bam   = inputs['merged_bam']
        create_bam_index(merged_bam)



        pipeline = 'snv_bam_annotator'
        bam_fhand = open(merged_bam)
        configuration = {'snv_bam_annotator': {'bam_fhand':bam_fhand}}
        settings = self._project_settings
        if 'Snvs' in settings and 'min_quality' in settings['Snvs']:
            min_quality = settings['Snvs']['min_quality']
            configuration['snv_bam_annotator']['min_quality'] = int(min_quality)

        return self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=output_dir)

class WriteAnnotationAnalyzer(Analyzer):
    'It writes all the output annotation files'
    def run(self):
        'It runs the analysis.'
        output_dir   = self._create_output_dirs()['result']
        inputs       = self._get_input_fpaths()
        repr_fpaths  = inputs['repr']

        output_files = ['vcf', 'gff']
        for seq_fpath in repr_fpaths:
            outputs = {}
            for output_kind in output_files:
                output_fpath = os.path.join(output_dir,
                                   _get_basename(seq_fpath) + '.' + output_kind)
                if os.path.exists(output_fpath):
                    os.remove(output_fpath)
                output_fhand = open(output_fpath, 'a')
                outputs[output_kind] = output_fhand
            io_fhands = {'in_seq': open(seq_fpath),
                         'outputs': outputs}
            seq_pipeline_runner(pipeline=None,
                                configuration=None,
                                io_fhands=io_fhands)

class SnvFilterAnalyzer(AnnotationAnalyzer):
    'It performs the filtering analysis of the snvs'

    def run(self):
        'It runs the analysis.'
        inputs, output_dir = self._get_inputs_and_prepare_outputs()

        pipeline, configuration = self._get_snv_filter_specs()

        return self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=output_dir)
    def _get_snv_filter_specs(self):
        'It gets the pipeline and configuration from settings'
        configuration = {}
        settings = self._project_settings['snv_filters']
        for filter_data in settings.values():
            if not filter_data['use']:
                continue
            name   = filter_data['name']
            if 'step_name' in filter_data:
                step_name = filter_data['step_name']
            else:
                step_name = name
            filter_config = {}
            for argument, value in filter_data.items():
                if argument == 'use' or argument == 'name':
                    continue
                filter_config[argument] = value
            filter_config['step_name'] = step_name
            configuration[name] =  filter_config
        pipeline = self._get_pipeline_from_step_name(configuration)

        return pipeline, configuration

    def _get_pipeline_from_step_name(self, configuration):
        'It gets the pipeline steps from filter_names'
        translator = {'uniq_contiguous' : unique_contiguous_region_filter,
                      'close_to_intron':close_to_intron_filter,
                      'high_variable_region': high_variable_region_filter,
                      'close_to_snv':close_to_snv_filter,
                      'close_to_limit': close_to_limit_filter,
                      'maf': major_allele_freq_filter,
                      'by_kind' : kind_filter,
                      'cap_enzyme': cap_enzyme_filter,
                      'is_variable': is_variable_filter,}
        pipeline = []
        for name in configuration:
            step_name = configuration[name]['step_name']
            #we have to remove the step name from the configuration because this
            #dict will be used as **kargs
            del configuration[name]['step_name']
            step = translator[step_name]
            step['name_in_config'] = name
            pipeline.append(step)
        return pipeline

class AnnotateDescriptionAnalyzer(AnnotationAnalyzer):
    '''This class is used to annotate sequences using the first hit on the given
    blast databases'''
    def run(self):
        'It runs the analysis'
        inputs, output_dir = self._get_inputs_and_prepare_outputs()
        blast_settings = self._project_settings['blast']

        settings = self._project_settings['Annotation']
        annot_settings = settings['description_annotation']
        description_databases = annot_settings['description_databases']

        #first we need some blasts
        project_dir = self._project_settings['General_settings']['project_path']
        blasts = {}
        for input_ in inputs['input']:
            for database in description_databases:
                db_kind = blast_settings[database]['kind']
                if db_kind == 'nucl':
                    blast_program = 'tblastx'
                else:
                    blast_program = 'blastx'

                blastdb = blast_settings[database]['path']
                blast   = backbone_blast_runner(query_fpath=input_,
                                                project_dir=project_dir,
                                                blast_program=blast_program,
                                                blast_db=blastdb,
                                                dbtype=db_kind)
                if input_ not in blasts:
                    blasts[input_] = []
                blasts[input_].append({'blast':blast, 'modifier':None})

        pipeline = []
        configuration = {}
        for database in description_databases:
            step = annotate_with_descriptions
            step['name_in_config'] = database
            pipeline.append(step)
            for input_ in inputs['input']:
                step_config = {'blasts': blasts[input_]}
                configuration[input_] = {}
                configuration[input_][database] = step_config

        return self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=output_dir)

class AnnotateMicrosatelliteAnalyzer(AnnotationAnalyzer):
    '''This class is used to annotate the microsatellite of a sequence as a
    feature'''

    def run(self):
        'It runs the analysis.'
        inputs, output_dir = self._get_inputs_and_prepare_outputs()

        pipeline = [annotate_microsatellites]
        configuration = {}
        return self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=output_dir)

class AnnotateOrfAnalyzer(AnnotationAnalyzer):
    '''This class is used to annotate the orf of a sequence as a feature'''

    def run(self):
        'It runs the analysis.'
        matrix = self._project_settings['Annotation']['orf_annotation']['estscan_matrix']
        inputs, output_dir = self._get_inputs_and_prepare_outputs()

        pipeline = [annotate_orfs]

        configuration = {'annotate_orfs': {'parameters':{'matrix':matrix}}}
        return self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=output_dir)

class AnnotateGoAnalyzer(AnnotationAnalyzer):
    'This class is used to annotate gos using blast2go'

    def run(self):
        'It runs the analysis.'
        inputs, output_dir = self._get_inputs_and_prepare_outputs()
        blast_settings = self._project_settings['blast']

        annot_settings = self._project_settings['Annotation']
        go_settings = annot_settings['go_annotation']
        go_database = go_settings['blast_database']

        #first we need some blasts
        project_dir = self._project_settings['General_settings']['project_path']
        blasts = {}
        for input_ in inputs['input']:
            db_kind = blast_settings[go_database]['kind']
            if db_kind == 'nucl':
                blast_program = 'blastn'
            else:
                blast_program = 'blastp'

            blastdb = blast_settings[go_database]['path']
            blast   = backbone_blast_runner(query_fpath=input_,
                                            project_dir=project_dir,
                                            blast_program=blast_program,
                                            blast_db=blastdb,
                                            dbtype=db_kind)
            blasts[input_] = blast
        # prepare pipeline
        pipeline      = [annotate_gos]
        configuration = {}
        for input_ in  inputs['input']:
            step_config = {'blast': blasts[input_]}
            configuration[input_] = {'annotate_gos': step_config}

        return self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=output_dir)

DEFINITIONS = {
    'filter_snvs':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'result':{'directory': 'annotation_repr'}},
         'analyzer': SnvFilterAnalyzer},
    'annotate_snv':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            'merged_bam':
                {'directory':'mapping_result',
                 'file': 'merged_bam'},
            'input':
                {'directory': 'annotation_input',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'result':{'directory': 'annotation_repr'}},
         'analyzer': SnvCallerAnalyzer},
    'annotate_introns':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            'input':
                {'directory': 'annotation_input',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'result':{'directory': 'annotation_repr'}},
         'analyzer': AnnotateIntronsAnalyzer},
    'write_annotation':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'result':{'directory': 'annotation_result'}},
         'analyzer': WriteAnnotationAnalyzer},
    'annotate_orthologs':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            'input':
                {'directory': 'annotation_input',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'result':{'directory': 'annotation_repr'}},
         'analyzer': AnnotateOrthologsAnalyzer},
    'annotate_description':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            'input':
                {'directory': 'annotation_input',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'result':{'directory': 'annotation_repr'}},
         'analyzer': AnnotateDescriptionAnalyzer},

    'annotate_microsatellite':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            'input':
                {'directory': 'annotation_input',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'result':{'directory': 'annotation_repr'}},
         'analyzer': AnnotateMicrosatelliteAnalyzer},
    'annotate_orf':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            'input':
                {'directory': 'annotation_input',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'result':{'directory': 'annotation_repr'}},
         'analyzer': AnnotateOrfAnalyzer},
    'annotate_go':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            'input':
                {'directory': 'annotation_input',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'result':{'directory': 'annotation_repr'},  },
                    #'b2g_file':{'directory':'go_files'}},
         'analyzer': AnnotateGoAnalyzer},
    }
