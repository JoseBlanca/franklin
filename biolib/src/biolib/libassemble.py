'''
Created on 15/09/2009

@author: peio
'''
import os
from biolib.utils.seqio_utils import (guess_seq_file_format,
                                      quess_seq_type, seqio)
from biolib.utils.cmd_utils import call

def check_and_fix_config(config):
    'It checks and prepares the configurations'
    errors = []
    # check_working_dir
    if 'work_dir' not in config:
        errors.append('Working dir not configured')
    elif not os.path.exists(config['work_dir']):
        os.mkdir(os.path.abspath(config['work_dir']))

    # reference
    _check_and_fix_seq_definitions(config, errors, 'reference')
    # reads
    _check_and_fix_seq_definitions(config, errors, 'reads')

    # Is the cfile OK?
    if errors:
        error_msg = ["Please correct this errors in your config file:\n"]
        for error in errors:
            error_msg.append('\t.-%s\n' % error)
        raise ValueError("".join(error_msg))
    return config

def _check_and_fix_seq_definitions(config, errors, kind):
    'It checks seq definition'
    if kind not in config:
        errors.append('%s not in config file' % kind)
    else:
        for seq_config in config[kind]:
            _check_and_fix_seq_definition(seq_config, errors, kind)

def _check_and_fix_seq_definition(seq_config, errors, kind):
    'It checks if the seq definitions seq_config is correct'
    #defaults aligners:
    aligners = {'short_seqs':{'binary':'nucmer' , 'parameters': '-l 20 -c 20'},
                'long_seqs' :{'binary':'nucmer' , 'parameters': ''}}

    #check required fields
    for field in ['name', 'seq_fpath']:
        if field not in seq_config:
            errors.append('%s : required %s field' % (kind, field))

    #qual file
    if 'qual_fpath' not in seq_config:
        seq_config['qual_fpath'] = None

    #check optional
    # check if we can guess file format
    if 'format' not in seq_config and not errors:
        try:
            format = guess_seq_file_format(open(seq_config['seq_fpath']))
            seq_config['format'] = format
        except ValueError:
            errors.append('File Format not defined and not recognizable')

    # check the seq type and use the apropiated alignments
    if kind == 'reads' and 'seq_type' not in seq_config and not errors:
        seq_type = quess_seq_type(open(seq_config['seq_fpath']),
                                  seq_config['format'], 40)
        seq_config['seq_type'] = seq_type
        if 'aliger' not in seq_config:
            seq_config['aligner'] = aligners[seq_type]['binary']
            seq_config['aligner_parameters'] = aligners[seq_type]['parameters']

def load_seqs_in_bank(configuration):
    '''This function loads the data in AMOS bank directory/format. It takes
     the input files from configuration'''

    seqs = prepare_files_to_load_in_bank(configuration)

    work_dir = configuration['work_dir']
    os.chdir(work_dir)

    # Create message_file from seq files
    cmd = ['tarchive2amos', '-o', '%s/message_file' % work_dir ]
    cmd.extend(seqs)
    retcode = call(cmd)[2]
    if retcode:
        msg = open(os.path.join(work_dir, 'tarchive2amos.log')).read()
        raise RuntimeError(msg)

    # Load message file in bank
    cmd = ['bank-transact', '-c', '-z', '-b', os.path.join(work_dir, 'bank'),
           '-m', os.path.join(work_dir, 'message_file.afg')]
    stdout, stderr, retcode = call(cmd)
    if retcode:
        raise RuntimeError(stdout, stderr)

def prepare_files_to_load_in_bank(configuration):
    '''This functions prepare files to load in bank from configuration. It also
    saves a list of the original seq files

    If it is needed it converts the files to fasta'''
    seqs_to_bank      = []
    work_dir          = configuration['work_dir']

    for seqs_file in configuration['reference'] + configuration['reads'] :
        seq_fpath  = seqs_file['seq_fpath']
        if 'qual_fpath' in seqs_file:
            qual_fpath = seqs_file['qual_fpath']
        else:
            qual_fpath = None

        if 'format' not in seqs_file:
            format = guess_seq_file_format(open(seq_fpath))
        else:
            format     = seqs_file['format']
        if format == 'fasta':
            dst = os.path.join(work_dir, os.path.basename(seq_fpath)) + '.seq'
            os.symlink(seq_fpath, dst)
            seqs_to_bank.append(dst)
        else:
            in_seq_fhand = open(seq_fpath)

            if qual_fpath == None:
                in_qual_fhand = None
            else:
                in_qual_fhand = open(qual_fpath)

            basename = '.'.join(os.path.basename(seq_fpath).split('.')[:-1])
            out_seq_fpath = os.path.join(work_dir, basename + '.seq')
            out_seq_fhand = open(out_seq_fpath, 'w')
            if format == 'fasta' and qual_fpath == None:
                out_qual_fhand = None
            else:
                out_qual_fpath = os.path.join(work_dir, basename+'.qual')
                out_qual_fhand = open(out_qual_fpath, 'w')

            seqio(in_seq_fhand, in_qual_fhand, format,
                  out_seq_fhand, out_qual_fhand, 'fasta')
            seqs_to_bank.append(out_seq_fhand.name)

    return seqs_to_bank




