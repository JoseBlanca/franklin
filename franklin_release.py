#!/usr/bin/env python

from optparse import OptionParser
import tempfile, shutil, os, subprocess, tarfile
from os.path import join

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
    cwd = os.getcwd()
    os.chdir(doc_dir)
    subprocess.check_call(["make", "clean"])
    subprocess.check_call(["make", "html"])
    os.chdir(cwd)
    temp_doc_dir = join(work_dir, 'html_doc')
    shutil.move(join(doc_dir, 'build', 'html'), temp_doc_dir)
    shutil.rmtree(doc_dir)
    shutil.move(temp_doc_dir, doc_dir)

    #remove the download page
    os.remove(join(doc_dir, 'download.html'))
                
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
    options = parser.parse_args()[0]
    options = {'indir': options.directory}
    if not options['indir']:
        raise ValueError('a franklin git directory is required')
    return options

if __name__ == '__main__':
    options = parse_options()
    prepare_release(**options)
