#!/usr/bin/env python
'This script prepares the franklin tarballs'

from optparse import OptionParser
import tempfile, shutil, os, subprocess, tarfile
from os.path import join, abspath, split

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

def prepare_release(indir, program):
    'It creates the tar.gz for the release'

    #the relased tool
    tool_name = split(indir)[-1]

    #we need a tempdir
    work_dir = tempfile.mkdtemp()

    #the version
    version = None
    init_fname = join(indir, '%s/__init__.py' % tool_name)
    for line in open(init_fname):
        if '__version__' in line:
            version = line.split('=')[1].strip().strip("'")
    if not version:
        raise RuntimeError('version not found in ' + init_fname)

    #the output name
    if program:
        tool_name = program
        #tool_name = 'ngs_backbone'
    release_name = '%s-%s' % (tool_name, version)
    tar_dir = join(work_dir, tool_name)
    release_dir = join(tar_dir, release_name)
    #copy the git dir
    shutil.copytree(indir, release_dir, ignore=shutil.ignore_patterns('.git*',
                                                                    '*project',
                                                            '*tool_release.py',
                                                                    '*.pyc'))

    #build the documentation
    doc_dir = join(release_dir, 'doc')
    #which kind of doc is this? do we have a source directory in the doc?
    if program in os.listdir(doc_dir):
        program_doc_dir = join(doc_dir, program)
    else:
        program_doc_dir = doc_dir

    #remove the download page from the doc
    download_page = join(program_doc_dir, 'download.rst')
    os.remove(download_page)

    index_page = join(program_doc_dir, 'index.rst')
    remove_download_link(index_page)

    cwd = os.getcwd()
    os.chdir(program_doc_dir)
    subprocess.check_call(["make", "clean"])
    subprocess.check_call(["make", "html"])
    os.chdir(cwd)
    temp_doc_dir = join(work_dir, 'html_doc')

    build_dir = join(program_doc_dir, '_build', 'html')

    shutil.move(build_dir, temp_doc_dir)
    shutil.rmtree(doc_dir)
    shutil.move(temp_doc_dir, doc_dir)

    #remove the build directory
    build_dir = join(release_dir, 'build')
    if os.path.exists(build_dir):
        shutil.rmtree(build_dir)

    #now we can create the tar.gz
    tar_fpath = release_name + '.tar.gz'
    tar = tarfile.open(tar_fpath, 'w:gz')
    tar.add(tar_dir, '')
    tar.close()

    shutil.rmtree(work_dir)

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-d', '--directory', dest='directory',
                      help='tool git directory')
    parser.add_option('-p', '--program', dest='program', default=None,
                      help='program to choose in the given repository')
    opts = parser.parse_args()[0]
    opts = {'indir': abspath(opts.directory), 'program':opts.program}
    if not opts['indir']:
        raise ValueError('a tool git directory is required')
    return opts

if __name__ == '__main__':
    prepare_release(**parse_options())

