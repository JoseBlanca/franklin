
Snv filters
===========

Once we have annotate snvs with annotate_snv analysis, we can filter those snvs using some of the available filters.

This analysis is performed by the filter_snv analysis.

The configuration of the snv filtering is handler by the :ref:`snv-filter` section of the backbone.conf file. By default only kind filtering is activated. To use any of the filters change the Use Option to True in each of the filter sections::

  use = True

TAny of the filters remove the snv, instead it keeps the information of the filter result and writes it to the vcf file.

A detailed view and configuration hooks for each filter are described below.

Available filters:
------------------
Close to intron
_______________
To use this filter you need first to annotate introns using annotate_introns :ref:`intron-annotation` analysis.

Use this filter if you want to filter snvs that are closer than X nucleotides to an intron.
This filter only have and option and it is the distance to the intron::

  distance = 30

Close to snv
____________
Use this filter if you want to filter snvs that are closer than X nucleotides of another snv.

This filter only have and option and is the distance to the intron::

  distance = 60

Close to limit
______________
Use this filter if you want to filter snvs that are closer than X nucleotides to the edges of the sequence

This filter only have and option and is the distance to the intron::

  distance = 60

High variable region
____________________
It filters the snv if it is in a variable region. The variability of a region is the frequency that we can find a snv per 100 nucleotides.

This filter have two options to configure. The frequency and the length of the region. By default all the sequence is taken::

  max_variability = 0.4
  window          = None  (All the sequence)

More frequent allele
_____________________
This filters rejects snvs in which the more frequent allele's frequency is bigger than the given one.

This filter only have one option to configure::

  frequency = 0.8

Kind filter
___________
Snv caller in annotate_snv detects if the snv is a SNP, INDEL OR COMPLEX. With this filter you can filter by the kind you want.

These are the correspondent codes to kinds::

  SNP       = 0
  INDEL     = 4
  COMPLEX   = 5

This filter only has one option to configure::

  kind = 0  (Choose your snv type)

Cap enzymes
___________
It filters the snv looking if it is detectable by restriction enzymes. If alleles in the snv behave different with the restriction enzyme then it filters it.

This filter uses remap from EMBOOS. The unique configurable option is to use most used enzymes or all of them::

  all_enzymes = True

Sequence filter
_______________

It filters the snv looking if the snv is in one of the given sequences(name).

This option only needs a path to the file with one sequence name per line::

  list_path = 'path to file with seq names'

Variable in group
_________________
In case you have more than one read_group library or sample and it it correcly defined in the sequence file names, you can filter looking if the snv is variable in a subgroup.

The configuration file provides three basic filter combinations form each one of the group. For example: if you want to filter see if the snv are variables is just some of your libraries::

  name       = 'is_variable_in_lb'
  step_name  = 'is_variable'
  group_kind = 'libraries'
  groups     = ['lib1', 'lib2']


First you have to use a uniq name for this filter step.

Group kind: one of samples, libraries or read_groups.

groups : a list of libraries in which you want to look at.


Unique contiguous
_________________
With this filter you can filter snvs that are in regions that are similar to other regions in the seq pool.

This filter have 3 configurable options::

  distance           = 'distance from each side os the snv to select a region'
  genomic_db         = 'path to seq pool'
  genomic_seqs_fpath = 'path to seq pool blast db'

