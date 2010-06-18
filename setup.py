'''
Created on 25/03/2009

@author: jose blanca
'''
#taken from django-tagging

import os, sys
#from setuptools import setup
from distutils.core import setup
import franklin

PACKAGE_DIR = 'franklin'
SCRIPTS_DIR = 'scripts'

#dependencies

python_exec = sys.executable
try:
    import Bio
except:
    print 'BioPython is required and is not installed for %s' % python_exec
    print 'Installation aborted'
    sys.exit(1)
try:
    import configobj
except:
    print'ConfigObj is required and is not installed for %s' % python_exec
    print 'Installation aborted'
    sys.exit(1)
all_optionals = True
try:
    import psubprocess
except:
    print 'psubprocess is not installed for %s, no parallel processing will be possible' % python_exec
    all_optionals = False
try:
    import pysam
except:
    print 'pysam is not installed for %s, SNP calling will fail' % python_exec
    all_optionals = False
try:
    import matplotlib
except:
    print 'matplotlib is not installed for %s, some statistics will fail' % python_exec
    all_optionals = False

answer = raw_input('Some optional requirements were not met, do you want to continue the installation (y/n)? ')
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

# Compile the list of packages available, because distutils doesn't have
# an easy way to do this.
packages, data_files, modules = [], [], []
root_dir = os.path.dirname(__file__)
pieces = fullsplit(root_dir)
if pieces[-1] == '':
    len_root_dir = len(pieces) - 1
else:
    len_root_dir = len(pieces)

for dirpath, dirnames, filenames in os.walk(os.path.join(root_dir,
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

wanted_scripts = ['backbone_analysis.py', 'seqio.py',
                  'backbone_create_project.py']
scripts = []
for dirpath, dirnames, filenames in os.walk(os.path.join(root_dir,
                                                         SCRIPTS_DIR)):
    for filename in filenames:
        if filename not in wanted_scripts:
            continue
        elif filename.endswith('.py') or filename == 'blast+':
            scripts.append(os.path.join(dirpath, filename))

setup(
    # basic package data
    name = "franklin",
    version = franklin.__version__,
    author='Jose Blanca, Peio Ziarsolo',
    author_email='jblanca@btc.upv.es',
    description='Some genomics related classes',
    # package structure
    include_package_data = True,
    packages=packages,
    package_dir={'':'.'},
    #py_modules = modules,

    package_data={'': ['data/samtools/*', 'data/blast/*', 'data/*.*']},
    requires=['BioPython', 'sqlalchemy', 'matplotlib', 'configobj',
              'pysam', 'psubprocess'],
    scripts=scripts,
)
