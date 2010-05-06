#!/usr/bin/env python
'This script prepares the franklin tarballs'

from optparse import OptionParser
import tempfile, shutil, os, subprocess, tarfile
from os.path import join, abspath

def remove_download_link(index_fpath):
    'It removes the download link from the index page'
    mod_fhand = open(join(index_fpath + '.mod'), 'w')
    in_toc = False
    for line in open(index_fpath, 'r'):
        if '.. toctree::' in line:
            in_toc = True
        if in_toc and 'download' in line:
            line = ''
        mod_fhand.write(line)
    shutil.move(mod_fhand.name, index_fpath)

def prepare_release(indir):
    'It creates the tar.gz for the release'

    #we need a tempdir
    work_dir = tempfile.mkdtemp()

    #the version
    version = None
    franklin_init = join(indir, 'franklin/__init__.py')
    for line in open(franklin_init):
        if '__version__' in line:
            version = line.split('=')[1].strip().strip("'")
    if not version:
        raise RuntimeError('version not found in ' + franklin_init)

    #the output name
    release_name = 'franklin-%s' % version
    release_dir = join(work_dir, release_name)
    #copy the git dir
    shutil.copytree(indir, release_dir, ignore=shutil.ignore_patterns('.git*',
                                                                    '*project',
                                                        '*franklin_release.py'))

    #build the documentation
    doc_dir = join(release_dir, 'doc')
    #remove the download page from the doc
    os.remove(join(doc_dir, 'source', 'download.rst'))
    remove_download_link(join(doc_dir, 'source', 'index.rst'))

    cwd = os.getcwd()
    os.chdir(doc_dir)
    subprocess.check_call(["make", "clean"])
    subprocess.check_call(["make", "html"])
    os.chdir(cwd)
    temp_doc_dir = join(work_dir, 'html_doc')
    shutil.move(join(doc_dir, 'build', 'html'), temp_doc_dir)
    shutil.rmtree(doc_dir)
    shutil.move(temp_doc_dir, doc_dir)


    #now we can create the tar.gz
    tar_fpath = release_name + '.tar.gz'
    tar = tarfile.open(tar_fpath, 'w:gz')
    tar.add(release_dir, arcname='')
    tar.close()

    shutil.rmtree(work_dir)

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-d', '--directory', dest='directory',
                    help='franklin git directory')
    opts = parser.parse_args()[0]
    opts = {'indir': abspath(opts.directory)}
    if not opts['indir']:
        raise ValueError('a franklin git directory is required')
    return opts

if __name__ == '__main__':
    prepare_release(**parse_options())

