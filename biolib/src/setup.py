'''
Created on 25/03/2009

@author: jose blanca
'''

from setuptools import setup
setup(
    # basic package data
    name = "biolib",
    version = "0.0.1",
    author='Jose Blanca, Peio Ziarsolo',
    author_email='jblanca@btc.upv.es',
    description='Some genomics related classes',
    # package structure
    packages=['biolib'],
    package_dir={'':'.'},
    package_data={'mypkg': ['data/*']},
    requires=['Nose', 'BioPython'],
    test_suite='nose.collector'
)
