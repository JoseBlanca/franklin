#!/usr/bin/env python2.6
'''
Created on 04/05/2010

@author: jose

This script is meant for franklin internal use only, we use it for the
parallelization.
'''
import sys, tempfile, traceback
from franklin.pipelines.pipelines import process_sequences_for_script

def main():

    #stdout = open('/tmp/script.test', 'w')
    args = sys.argv[1:]
    if len(args) == 6:
        (in_fpath_seq, file_format,
                           pipeline, configuration, out_fpath, temp_dir) = args
        in_fpath_qual = None
    else:
        (in_fpath_seq, in_fpath_qual, file_format,
                           pipeline, configuration, out_fpath, temp_dir) = args
    tempfile.tempdir = temp_dir

    if in_fpath_qual:
        process_sequences_for_script(in_fpath_seq, file_format,
                                     pipeline, configuration, out_fpath,
                                     in_fpath_qual=in_fpath_qual)
    else:
        process_sequences_for_script(in_fpath_seq, file_format,
                                     pipeline, configuration, out_fpath)

    #stdout.close()

if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        sys.stderr.write(str(error) + '\n')
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)
