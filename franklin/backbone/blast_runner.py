'''
A blast runner hardly tied to backbone folder structure
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

import logging
from os.path import join, exists, splitext, basename, split
from os import makedirs, remove, listdir
from franklin.backbone.specifications import (BACKBONE_DIRECTORIES,
                                              BACKBONE_BASENAMES)
from franklin.utils.cmd_utils import call
from franklin.utils.misc_utils import get_num_threads, rel_symlink
from franklin.seq.readers import guess_seq_file_format
from franklin.utils.seqio_utils import seqio
from tempfile import NamedTemporaryFile
from franklin.backbone.analysis import scrape_info_from_fname

LOGGER_NAME = 'franklin'

def blast_runner_plus(seq_fpath, blast_db, blast_type, result_fpath,
                       threads=False):
    'It runs a blast giving a file and a database path'
    cmd = [blast_type, '-db', blast_db, '-num_alignments', '25',
           '-num_descriptions', '25', '-evalue', '0.0001', '-outfmt', '5',
           '-query', seq_fpath, '-out', result_fpath]
    if threads:
        cmd.extend(['-num_threads', str(threads)])
    call(cmd, raise_on_error=True, log=True)

def makeblastdb_plus(seq_fpath, dbtype, outputdb=None):
    'It creates the blast db database'
    cmd = ['makeblastdb', '-in', seq_fpath, '-dbtype', dbtype]
    if outputdb is not None:
        cmd.extend(['-out', outputdb])
    call(cmd, raise_on_error=True)

def _get_basename(fpath):
    'It returns the base name without path and extension'
    return splitext(basename(fpath))[0]

def _create_temp_fasta_file(fpath):
    'It creates a fasta file format temfile'
    fasta_fhand = NamedTemporaryFile(suffix='.fasta', mode='a')
    seqio(in_seq_fhand=open(fpath),
          out_seq_fhand=fasta_fhand,
          out_format='fasta')
    return fasta_fhand

def guess_blastdb_kind(blastdb):
    'it infers the kind of the blastdb'
    blastdir, basename_ = split(blastdb)
    for file_ in listdir(blastdir):
        if file_.startswith(basename_ + '.') and file_ != basename_ :
            if splitext(file_)[1][1] == 'n':
                return 'nucl'
            elif splitext(file_)[1][1] == 'p':
                return 'prot'

def guess_blast_program(seq_type, db_type, prefer_tblastx=False):
    'It guess the blast program to use for the given seq and db type'
    if seq_type == 'nucl' and db_type == 'nucl':
        if prefer_tblastx:
            blast_program = 'tblastx'
        else:
            blast_program = 'blastn'
    elif seq_type == 'nucl' and db_type == 'prot':
        blast_program = 'blastx'
    elif seq_type == 'prot' and db_type == 'nucl':
        blast_program = 'tblastn'
    elif seq_type == 'prot' and db_type == 'prot':
        blast_program = 'blastp'

    return blast_program

def backbone_blast_runner(query_fpath, project_dir, blast_program,
                          blast_db=None, blast_db_seq=None, dbtype='nucl',
                          threads=False):
    '''It returns the blast if the results doesn't exist'''
    if blast_db is None and blast_db_seq is None:
        raise RuntimeError('It needs a blast database or seqfile')

    #create a logger
    logger = logging.getLogger(LOGGER_NAME)
    query_basename = _get_basename(query_fpath)
    blast_dir = join(project_dir, BACKBONE_DIRECTORIES['blast_dir'])

    if blast_db:
        result_dir = join(blast_dir, query_basename, _get_basename(blast_db))
    else:
        result_dir = join(blast_dir, query_basename,
                          _get_basename(blast_db_seq))
    if not exists(result_dir):
        makedirs(result_dir)
    result_fpath = join(result_dir,
                        '%s.%s.xml' % (BACKBONE_BASENAMES['blast_basename'],
                                       blast_program))
    if exists(result_fpath):
        logger.info('Using the stored blast result %s' % result_fpath)
        return result_fpath

    #the input file should be fasta
    fasta_query_fhand = None
    fasta_db_fhand = None
    if guess_seq_file_format(open(query_fpath)) != 'fasta':
        fasta_query_fhand = _create_temp_fasta_file(query_fpath)
        query_fpath = fasta_query_fhand.name

    #we have to create a database in BACKBONE_DIRECTORIES['blast_databases']
    if blast_db_seq:
        blast_db = make_backbone_blast_db(project_dir, blast_db_seq, dbtype)

    logger.info('Running the blast %s' % result_fpath)
    try:
        blast_runner_plus(query_fpath, blast_db, blast_program,
                                 result_fpath, threads=threads)
    except RuntimeError as error:
        if exists(result_fpath):
            remove(result_fpath)
        msg = '%s \n database: %s\n database type: %s' % (str(error),
                                                             blast_db, dbtype)
        raise RuntimeError(msg)

    if fasta_query_fhand:
        fasta_query_fhand.close()
    if fasta_db_fhand:
        fasta_db_fhand.close()

    return result_fpath

def make_backbone_blast_db(project_dir, blast_db_seq, dbtype):
    'It formats a blastdb when need it'
    logger = logging.getLogger(LOGGER_NAME)
    #the name should be the basename of the blast_db_seq
    db_dir = join(project_dir, BACKBONE_DIRECTORIES['blast_databases'])
    if not exists(db_dir):
        makedirs(db_dir)
    db_seq_fpath = join(db_dir, _get_basename(blast_db_seq))
    if not exists(db_seq_fpath):
        #which is the name of the new databae?
        blast_db_seq_format = guess_seq_file_format(open(blast_db_seq))
        if blast_db_seq_format == 'fasta':
            rel_symlink(blast_db_seq, db_seq_fpath)
        else:
            seqio(in_seq_fhand=open(blast_db_seq),
                  out_seq_fhand=open(db_seq_fpath, 'w'),
                  out_format='fasta')
        logger.info('Formatting the database %s' % db_seq_fpath)
        try:
            makeblastdb_plus(db_seq_fpath, dbtype=dbtype)
        except RuntimeError:
            msg = 'Error making blastdb. db:%s\n dbtype:%s\n' % \
                                               (db_seq_fpath, dbtype)
            remove(db_seq_fpath)
            raise RuntimeError(msg)
    return db_seq_fpath


