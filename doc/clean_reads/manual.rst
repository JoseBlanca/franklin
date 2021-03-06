=============
 clean_reads
=============

--------------------------
Clean Sanger and NGS reads
--------------------------

:Author: jblanca@upv.es, pziarsolo@upv.es
:Date: 28-4-2011
:Version: 0.1
:Manual section: 1

SYNOPSIS
========

  **clean_reads** **-i** *input_seqs.fastq* **-o** *ouput_seqs.fastq* **-p** **sanger|454|illumina|solid**
  [**-f** *input_file_format*] [**-g** *output_file_format*] [**--mask_no_trim**]
  [**--double_encoding** *double_encode_solid*] [**-q** *input_quals.fasta*] [**-u** *output_quals.fasta*]
  [**-a** *adaptors_file.fasta*] [**-v** *vector_file.fasta*] [**-d** *blast_vector_database*]
  [**-r** *regular_expression_list*] [**-e** *bases_to_trim_at_edges*] [**-x**]
  [**-n** *allowable_percentage_ of_n*] [**--lucy_splice** *lucy_splice_file*]
  [**--lucy_error** *lucy_error*] [**-lucy_window** *lucy_window*] [**--lucy_bracket** *lucy_bracket*]
  [**--qual_window** *qual_window*] [**--qual_threshold** *qual_threshold*] [**--only_3_end** *True*]
  [**--solid_qual_length** *window length*] [**--solid_qual_threshold** *quality_threshold*]
  [**--solid_allow_missing_call**] [**-m** *minimum_seq_length*] [**--filter_evalue** *evalue*]
  [**--filter_idendity** *identity*] [**--filter_num_residues** *num_residues*]
  [**--filter_length_percentage** *percentage_of_length*] [**--filter_dbs** *list_of_blast_dbs*]
  [**-t** *number_of_threads*] [**--error_log** *error_log_file*]


DESCRIPTION
===========

clean_reads cleans Sanger, 454, Illumina and Solid reads taking into account the read quality and the vector and adaptor contaminations.

clean_reads uses internally several algorithms and third party tools to do the quality, adaptor and vector cleaning.
The exact algorithm and tools used for each set of reads will depend on the sequencing platform and on the knowledge of the cloning process that we have.
The sequencing platforms supported are:

  * Sanger (with or without quality)
  * 454
  * Illumina
  * solid

For quality trimming clean_reads is capable of using three different algorithms.
For 454 and Sanger reads with quality it uses lucy_ to remove the bad quality regions of the sequences.
For Illumina and solid reads uses a quality thresholded sliding window algorithm.
When no qualities none of these approaches can be used, but we still are able to infer the quality of the sequence by noticing how many **Ns** are in the reads.
clean_reads uses the trimpoly_ algorithm to get rid of the stretches with too many undetermined nucleotides.

The vector is removed by using BLAST_ and lucy_ (for Sanger and 454) while the adaptors will be trimmed by the blast-short algorithm.

clean_reads also supports the removal of regular expressions patterns by using an exact match algorithm.

The solid read cleaning is more limited than for the other platforms due to the color encoding.
Only the quality trimming and quality filtering will be allowed.
For the solid reads there is an extra quality filtering step that will remove reads that fail to meet minimum quality criteria for the first colors and for the presence of missing calls.

Also filtering out the reads by blasting them against a database is supported.

OPTIONS
=======

**-i**, **--seq_in** *input_sequence_file*
        The path to the sequence file to be cleaned.
**-o**, **--seq_out** *output_sequence_file*
        The path to the resulting output file.
**-p**, **--platform** *sanger|454|illumina|solid*
        Sequencing platform used to generate the sequences.
        This is an important parameter because some of the algorithms to be used will depend on the sequencing platform.

**-f**, **--format** *input_file_format*
        The input file format.
        The supported formats are: fasta, sanger fastq (fastq), illumina fastq (ifastq) and csfasta.
        Also any other `sequence format supported by Biopython <http://www.biopython.org/wiki/SeqIO>`_ could be also used.
        This parameter is optional because in some instances clean_reads will be able to guess the input file format.

**-g**, **--output_format** *output_file format*
        The output file format.
        By default the input and output file will have the same format, but this parameter allows to change that behaviour.
        The supported output formats are fasta and sanger fastq, but any of the Biopython supported formats can be used.
        csfasta is not supported as an format of the output file, so when the input is in csfasta the default output will be sanger fastq.

**--mask_no_trim**
		Do not trim the reads, just mask them.
		Instead of trimming the bad quality and vector regions it would lower-case them when this parameter is enabled.
		The default behavior is trimming.

**--double_encoding**
        Use double encoding for solid.
        Use the letters (NACGT) to encode the solid colors (.0123).

**-q**, **--qual_in**
        The path to the qualities input file.

**-u**, **--qual_out**
        The path to the qualities output file.

**-a**, **-adaptors** *adaptors_fasta_file*
        The path to the adaptors file.
        The adaptors to be trimmed should be given as a fasta file.
        These adaptors will be look for by using the blast-short algorithm.
        For the 454 platform clean_reads holds an internal file with some standard adaptors, so in this case this option can be used without providing any file name and those standard adaptors will be used.

**-v**, **-vector** *vector_fasta_file*
        The path to the vector file.
        A fasta file with the vector sequence.

**-d**, **-vector_db** *vector_blast_db*
        Vector BLAST database.
        A BLAST formated database like Univec can be provided to detect and remove unknown vectors.
        This option can be used without giving a vector BLAST database, in that case an internal Univec_Core database will be used.

**-r**, **--re_words** *regular_expression_list*
        A list of regular expressions to be trimmed.
        The regular expressions should be surrounded by double quotes and separated by commas. Valid examples could be: "^ACGT" and "^ACGT","AAAAA$".

**-e**, **-edge_trim** *edge trimming values*
        The number of nucleotides to be removed from the sequence beginning and end.
        This would be an unconditional trimming.
        The number of nucleotides should be given as two integers separated by a comma (e.g. 0,10).

**-x**, **--disable_quality_trimming**
    The quality trimming algorithms will not be used.
    This parameter will not affect the quality filtering, for instance the solid quality filtering will continue to work.

**-n**, **-n_percent**
        Allowed percent of Ns.
        Trim the regions with a percent of unknown bases above the given threshold.

**--lucy_splice**
        The lucy splice site file.
        The splice site definition used for the exact trimming of the vector and adaptor.
        Refer to the lucy man page for the format of this file.

**--lucy_error**
        lucy error parameter.
        Refer to the lucy man page for an explanation of this parameter.

**--lucy_window**
        lucy window sizes.
        Refer to the lucy man page for an explanation of this parameter.

**--lucy_bracket**
        lucy bracket parameter.
        Refer to the lucy man page for an explanation of this parameter.

**-qual_window** *integer*
        Length of the window used for quality trimming.
        This parameters affects the quality cleaning of Illumina and Solid sequences.

**--qual_threshold** *integer*
        Quality threshold.
        The phred scaled quality threshold used to discriminate good from bad nucleotides.

**--only_3_end** *boolean*
        Quality trim only from the 3' end of the read.
        This parameter will be set to True for solid and False for the rest of the platforms.

**--solid_qual_length** *integer*
        Number of 5' colors to consider to quality filtering.

**--solid_qual_threshold** *integer*
        Minimum mean quality allowable for solid reads.
        The mean quality of a number of 5' colors will discriminate if the read is to be completely removed and not only trimmed.

**--solid_allow_missing_call**
        Disable filtering out solid reads with more than 1 missing calls.

**-m**, **-min_len** *integer*
        Minimum number of nucleotides after the trimming.
        All sequences shorted than the given length will be filtered out.

**--filter_evalue** *float*
        Sequences with a better evalue against any of the given databases will be filtered out.
        This filter is not used by default.

**--filter_idendity** *float*
        Minimun identity to consider a BLAST hsp (default 95%)

**--filter_num_residues** *int*
        Sequences with BLAST matches longer than this length will be filtered out.
        The default behaviour is based on the percentage, so this filter is not used.

**--filter_length_percentage** *float*
        Sequences with BLAST matches longer than this length will be filtered out (default 75)
        The percentage is calculated as the matched region divided by the total read length.

**--filter_dbs** *database list*
        List of BLAST databases used for similarity filtering.

**-t**, **-threads**
        Number of threads to use.
        The reads can be processed in parallel using several processes.

**--tmpdir** *temporary/directory*
		Directory to be used for the temporary files.

**--error_log** *error_log_file*
        Path to the error log file to use (default clean_reads.error)

.. include:: ../ngs_backbone/links.txt
