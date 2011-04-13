'''
Created on 25/03/2009

@author: jose blanca
'''
#taken from django-tagging

import imp
import os, sys, glob, fnmatch
from distutils.core import setup

import franklin


def opj(*args):
    path = os.path.join(*args)
    return os.path.normpath(path)

def find_data_files(srcdir, *wildcards, **kw):
    # get a list of all files under the srcdir matching wildcards,
    # returned in a format to be used for install_data
    def walk_helper(arg, dirname, files):
        if '.svn' in dirname or '.git' in dirname:
            return
        names = []
        lst, wildcards = arg
        for wc in wildcards:
            wc_name = opj(dirname, wc)
            for f in files:
                filename = opj(dirname, f)
                if fnmatch.fnmatch(filename, wc_name) and not os.path.isdir(filename):
                    names.append(filename)
        if names:
            lst.append( (dirname, names ) )

    file_list = []
    recursive = kw.get('recursive', True)
    if recursive:
        os.path.walk(srcdir, walk_helper, (file_list, wildcards))
    else:
        walk_helper((file_list, wildcards),
                    srcdir,
                    [os.path.basename(f) for f in glob.glob(opj(srcdir, '*'))])
    return file_list

def check_modules(modules):
    'It check the modules needed'
    all_optionals = True
    for module_name, module_info in  modules.items():
        try:
            # we use this and not a simple import because matpliotlib writes
            # a directory in the first import. with posible wrong permissions
            imp.find_module(module_name)
        except ImportError:
            if 'required' in module_info and module_info['required']:
                print '%s is required and is not installed for %s' % \
                                                   (module_name, sys.executable)
                print "Instalation Aborted"
                sys.exit(1)
            print "%s is not installed, %s" % (module_name, module_info['msg'])
            all_optionals = False

    if not all_optionals:
        msg = 'Some optional requirements were not met, do you want to continue'
        msg += ' the installation (y/n)? '
        answer = raw_input(msg)
        if answer.lower() not in ('y', 'yes'):
            print 'Installation aborted'
            sys.exit()

def fullsplit(path, result=None):
    """
    Split a pathname into components (the opposite of os.path.join) in a
    platform-neutral way.
    """
    if result is None:
        result = []
    head, tail = os.path.split(path)
    if head == '':
        return [tail] + result
    if head == path:
        return result
    return fullsplit(head, [tail] + result)

def _guess_packages(directory):
    # Compile the list of packages available, because distutils doesn't have
    # an easy way to do this.
    packages, data_files, modules = [], [], []

    pieces = fullsplit(directory)
    if pieces[-1] == '':
        len_root_dir = len(pieces) - 1
    else:
        len_root_dir = len(pieces)

    for dirpath, dirnames, filenames in os.walk(os.path.join(directory,
                                                             PACKAGE_DIR)):
        if '__init__.py' in filenames:
            package = '.'.join(fullsplit(dirpath)[len_root_dir:])
            packages.append(package)
            for filename in os.listdir(dirpath):
                if (filename.startswith('.') or filename.startswith('_') or
                    not filename.endswith('.py')):
                    continue
                modules.append(package + '.' + filename)
        elif filenames:
            data_files.append([dirpath, [os.path.join(dirpath, f) for f in filenames]])
    return packages

def _guess_scripts(directory, wanted_scripts):
    scripts = []
    for dirpath, dirnames, filenames in os.walk(os.path.join(root_dir,
                                                             SCRIPTS_DIR)):
        for filename in filenames:
            if filename in wanted_scripts:
                scripts.append(os.path.join(dirpath, filename))
    return scripts

import distutils.command.install_data

#Taken from http://wiki.python.org/moin/Distutils/Tutorial
# Specializations of some distutils command classes
## Code borrowed from wxPython's setup and config files
## I am not 100% sure what's going on, but it works!
class wx_smart_install_data(distutils.command.install_data.install_data):
    """need to change self.install_dir to the actual library dir"""
    def run(self):
        install_cmd = self.get_finalized_command('install')
        self.install_dir = getattr(install_cmd, 'install_lib')
        return distutils.command.install_data.install_data.run(self)

PACKAGE_DIR = 'franklin'
SCRIPTS_DIR = 'scripts'

#dependencies
MODULES = {'Bio':        {'required': True},
           'configobj':  {'required': True},
           'psubprocess':{'msg':'no parallel processing will be possible'},
           'pysam':      {'msg':'SNP calling will fail'},
           'matplotlib': {'msg':'some statistics will fail'}
           }

# check module dependencies
check_modules(MODULES)

root_dir = os.path.dirname(__file__)
packages = _guess_packages(root_dir)
wanted_scripts = ['backbone_analysis.py', 'seqio.py',
                  'backbone_create_project.py']
scripts = _guess_scripts(root_dir, wanted_scripts)
data_files = find_data_files('ext', '*')

dist = setup(
    # basic package data
    name         = "franklin",
    version      = franklin.__version__,
    author       ='Jose Blanca, Peio Ziarsolo',
    author_email ='jblanca@upv.es',
    description  ='Some genomics related classes',
    # package structure
    include_package_data = True,
    packages=packages,
    package_dir={'':'.'},

    package_data={'': ['data/*.*']},
    data_files = data_files,
    cmdclass = { 'install_data':    wx_smart_install_data },
    requires=['BioPython', 'matplotlib', 'configobj',
              'pysam', 'psubprocess'],
    scripts=scripts,
)
