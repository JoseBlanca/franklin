'''
Created on 15/03/2010

@author: peio
'''
from tempfile import NamedTemporaryFile
import os, shutil
from franklin.utils.misc_utils import NamedTemporaryDir
from franklin.backbone.analysis import (Analyzer, scrape_info_from_fname,
                                        LastAnalysisAnalyzer)
from franklin.utils.cmd_utils import call
from franklin.seq.readers import guess_seq_file_format
from franklin.backbone.specifications import (BACKBONE_DIRECTORIES,
                                              BACKBONE_BASENAMES)
from franklin.utils.seqio_utils import (seqio, cat, seqs_in_file,
                                        write_seqs_in_file)

class PrepareMiraAssemblyAnalyzer(Analyzer):
    'It assembles the cleaned reads'

    def run(self):
        '''It runs the analysis. It checks if the analysis is already done per
        input file'''
        self._log({'analysis_started':True})
        files_454 = []
        files_sanger_with_qual = []
        files_sanger_without_qual = []
        for path in self._get_input_fpaths()['reads']:
            fpath = path.last_version
            fhand = open(fpath)
            fname = os.path.split(fpath)[-1]
            if 'pl_454' in fname.lower():
                files_454.append(fhand)
            elif 'pl_sanger' in fname.lower():
                format_ = guess_seq_file_format(fhand)
                if format_ == 'fasta':
                    files_sanger_without_qual.append(fhand)
                elif format_ == 'fastq':
                    files_sanger_with_qual.append(fhand)

        #fastq are processed before
        files_sanger = files_sanger_with_qual[:]
        files_sanger.extend(files_sanger_without_qual)

        #all files should be fasta and fasta.qual
        output_dir = self._create_output_dirs()['assembly_input']
        project_name = self._get_project_name()
        for ext, files in (('_in.454', files_454),
                           ('_in.sanger', files_sanger)):
            base_name = os.path.join(output_dir, project_name + ext)
            fasta_fpath = base_name + '.fasta'
            qual_fpath = base_name + '.fasta.qual'
            if os.path.exists(fasta_fpath) or not files:
                continue
            fasta_fhand = open(fasta_fpath, 'w')
            qual_fhand = open(qual_fpath, 'w')
            self._cat_to_fasta(files, fasta_fhand, qual_fhand)
            fasta_fhand.close()
            qual_fhand.close()

        # close all files
        for file_ in files_454 + files_sanger:
            file_.close()
        self._log({'analysis_finished':True})

    @staticmethod
    def _files_to_temp_fasta(files):
        'It converts the given files to a temporary fasta and qual'
        fastas, quals = [], []
        for file_ in files:
            #are we dealing with a fastq file (with qual)
            if 'fastq' in os.path.splitext(file_.name)[-1]:
                qual  = NamedTemporaryFile(suffix='.qual')
                fasta = NamedTemporaryFile(suffix='.fasta')
                seqio(in_seq_fhand=file_, out_seq_fhand=fasta,
                        out_qual_fhand=qual, out_format='fasta')
            else:
                #the file is already fasta
                fasta = file_
                qual = None
            fastas.append(fasta)
            quals.append(qual)
        return fastas, quals

    def _cat_to_fasta(self, files, fasta_fhand, qual_fhand):
        'It cats the given files together and it returns a fasta and qual'
        #all files should be in fasta and qual
        fastas, quals = self._files_to_temp_fasta(files)
        cat(fastas, fasta_fhand)
        cat(quals , qual_fhand)

class MiraAssemblyAnalyzer(Analyzer):
    'It assembles the cleaned reads'

    def run(self):
        '''It runs the analysis. It checks if the analysis is already done per
        input file'''
        self._log({'analysis_started':True})
        proj_name = self._get_project_name()
        #we need an assembly dir for this run
        #in this case every assembly is done in a new directory
        #the directory for this assembly
        output_dirs = self._create_output_dirs(timestamped=True)
        assembly_dir = output_dirs['analysis']
        mira_dir  = os.path.join(assembly_dir, '%s_assembly' % proj_name)
        original_dir = os.getcwd()
        os.chdir(assembly_dir)
        #with technologies have we used?
        #each input file should have a link in the assembly dir
        techs = set()
        for path in self._get_input_fpaths()['reads']:
            fpath = path.last_version
            fname = os.path.split(fpath)[-1]
            mira_input_fname = os.path.join(assembly_dir, fname)
            os.symlink(fpath, mira_input_fname)
            if '454' in fname:
                techs.add('454')
            elif 'sanger' in fname:
                techs.add('sanger')

        settings = self._project_settings['Mira']
        #job part of the command

        job = settings['job_options']
        job.extend(techs)
        job = ','.join(job)
        job = '-job=' + job

        cmd = ['mira', '-project=%s' % proj_name, '-fasta', job]
        if 'general_settings' in settings:
            general_settings = settings['general_settings']
            cmd.extend(general_settings)

        cmd.append('-OUT:rld=yes')

        for tech in techs:
            tech_str = '%s_settings' % tech
            if tech_str in settings:
                cmd.append(tech_str.upper())
                cmd.extend(settings[tech_str])

        stdout = open(os.path.join(assembly_dir, 'stdout.txt'), 'w')
        stderr = open(os.path.join(assembly_dir, 'stderr.txt'), 'w')
        retcode = call(cmd, stdout=stdout, stderr=stderr)[-1]
        if retcode:
            stdout.flush()
            stdout.seek(0)
            stderr.flush()
            stderr.seek(0)
            raise RuntimeError(open(stdout.name).read() + \
                               open(stderr.name).read())

        stdout.close()
        stderr.close()
        #remove the log directory
        mirachkpt_dir = os.path.join(mira_dir,  "%s_d_chkpt" % proj_name)
        if os.path.exists(mirachkpt_dir):
            shutil.rmtree(mirachkpt_dir)

        #results
        mira_results_dir = os.path.join(mira_dir, '%s_d_results' % proj_name)
        results_dir = output_dirs['result']
        mira_basename = '%s_out' % proj_name
        file_exts = (('.ace', '.ace'),
                     ('.caf', '.caf'),
                     ('.unpadded.fasta', '.fasta'),
                     ('.unpadded.fasta.qual', '.qual'))
        for mira_ext, result_ext in file_exts:
            mira_fpath = os.path.join(mira_results_dir,
                                      mira_basename + mira_ext)
            res_fpath = os.path.join(results_dir,
                                     BACKBONE_BASENAMES['contigs'] + result_ext)
            if os.path.exists(mira_fpath):
                os.symlink(mira_fpath, res_fpath)


        mira_info_dir = os.path.join(mira_dir, '%s_d_info' % proj_name)
        if os.path.exists(mira_info_dir):
            os.symlink(mira_info_dir, os.path.join(results_dir,
                                                  BACKBONE_DIRECTORIES['info']))
        self._log({'analysis_finished':True})
        os.chdir(original_dir)

DEFINITIONS ={
    'prepare_mira_assembly':
        {'inputs':{
            'reads':
                {'directory': 'cleaned_reads',
                 'file_kinds': 'sequence_files'}
            },
         'outputs':{'assembly_input':{'directory': 'assembly_input'}},
         'analyzer': PrepareMiraAssemblyAnalyzer,
        },
    'mira_assembly':
        {'inputs':{
            'reads':
                {'directory': 'assembly_input',
                 'file_kinds': 'sequence_or_qual_files'}
            },
         'outputs':{
                    'analysis':  {'directory': 'assemblies'},
                    'result':  {'directory': 'assembly_result'},
                    },
         'analyzer': MiraAssemblyAnalyzer,
        },
    'select_last_assembly':
        {'inputs':{'analyses_dir':{'directory': 'assemblies'}},
         'outputs':{'result':{'directory': 'assembly_result',
                              'create': False}},
         'analyzer': LastAnalysisAnalyzer,
        },
   }
