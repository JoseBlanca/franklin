
Architecture
============

This section is intended to be read by python developers interested in fixing bugs or extending the functionality of franklin. A standard user does not need to read it.

The franklin functionality depends heavily on Biopython and on some external tools, these are the foundations that franklin requires to do its work. The functions and classes provided by franklin itself are collected in two main packages: franklin and backbone.

.. image:: figures/franklin_architecture.png

The package franklin provides the required functions to do the different analysis. The backbone package (located in franklin/backbone) is used to build the user interface to the application. For instance in the cleaning analysis:
  * the backbone will look at the configuration file, it will look for the input and output paths,
  * and the franklin package has the functions that read the files, clean the sequence objects (using several external applications), and write the output sequences.

franklin package architecture
-----------------------------

In franklin all the analyses are done on sequences so the main class is SeqWithQuality (franklin/seq/seqs.py). SeqWithQuality is just a Biopython SeqRecord with some modifications. These SeqWithQuality sequences when a sequence file is read by the function seqs_in_file (franklin/seq/readers.py). Once the analysis is completed the resulting sequence is written to the output by a writer class named SequenceWriter (franklin/seq/writer.py).

.. image:: figures/pipeline.png

An analysis is represented by a pipeline composed by several steps. All the analyses are run by the same function: seq_pipeline_runner (franklin/pipelines/pipelines.py). This function creates a list of analysis steps and a configuration for them and runs them. Each step is composed by a function and a configuration. These analysis step functions are classified in three types: mappers, filters and bulk mappers. The mapper and filter take a SeqWithQuality as input and generate either a new SeqWithQuality (the mapper) or a boolean (the filter). The bulk mappers take a stream of SeqWithQuality sequences and return a new streams.

The analysis steps are entries in a python dictionary that define the analysis function, the kind of function and the description of the analysis. These analysis steps are defined in several modules in the fraklin/pipelines package like seq_pipeline_steps.py (sequence cleaning), annotation_steps.py (sequence annotation) and snv_pipeline_steps.py (snv calling and filtering).

franklin is designed to deal with millions of sequences. Since most of the analysis performed can be done for each sequence independently only once sequence is keep in memory at each time. This is accomplished by using python generators, the sequence streams are in fact a python generator that yield one SeqWithQuality at a time. The function seqs_in_file is the main of such generators.


