
Annotation
==========

There are different annotation analyses for the sequences, but they all operate in a similar way, so it is worth to explain the general annotation process. The sequences to annotate should be placed in one or several files at annotations/input/. When an annotation is done the results are stored in a kind of database at annotations/db/. There is a versioned file in annotations/db/ for each original file set to annotate. The annotations are cumulative, so if we annotate the ORF and after that the SSRs both will be stored at annotations/db/. At every time after an annotation the real output files can be generated running the write_annotation analysis. The output files will be found at annotations/features/. The output files are:

* VCF_ file. It stores the SNP and indel information.
* DNA and pep ORF fasta files. The DNA to translate and the translation.
* a csv file for the SSR (microsatellite) information.
* GFF_ file. It stores all found sequence features.
* annot and dat files for Blast2Go.
* a text file with the orthologs.


.. _snp-calling:

SNP calling
-----------

SNP calling is a form of sequence annotation. To annotate some sequences you need the input sequences to annotate, a bam file in mapping/bams/merged.bam and a reference genome in mapping/reference/reference.fasta. SNP calling is based on the samtools pileup with some filters on top. The ngs_backbone analysis is: annotate_snvs.

Configuration parameters
________________________

In the section Snvs:

min_quality
  Minimum allele quality. The alleles with a lower quality will not be consider. It is a phred quality.

min_mapq
  Reads mapped with a lower quality will not be considered. This is the phred quality found on the bam file.

min_num_alleles
  Usually it will be 1 or 2. If 1 the positions with only one allele different from the reference will be considered if 2 only the positions with at least 2 alleles will be considered.

default_bam_platform
  Default is None. Backbone takes platform information from the bam file. If your bam doesn't have readgroup information. You can provide the platform information using this option. Backbone needs platform information to perform the snv calling.

Also a subsection edge_removal inside Snvs could be defined. The parameters in this section are used to when we want to ignore the nucleotides close to the edged of the reads. This parameters are dependent on the platform and are: 454_left, 454_right, sanger_left, sanger_right, illumina_left and illumina_right.


.. _snv-filter:

SNV filtering
-------------
Once you have the snv annotated, you may not want all of them. Only some of then may be of your interest. ngs_backbone provide an analysis to filter snvs.

Yo can find more information on how to use this analysis :doc:`here <snv_filters>`.


.. _blast-databases:

Blast databases
---------------

Several annotations make use of blast. ngs_backbone requires some data about the blast databases to be able to run. This information should be set up in the ngs_backbone.conf file on the blast section. In this section for each blast database a subsection for every database should be prepared. An example for the nr would be::

 [blast]
 [[nr]]
 path = /absolute/path/to/blast/database
 species = 'all'  #this database is not species specific
 kind = 'prot'    #or nucl for nucleotide databases

.. _description-annotation:

Description annotation
----------------------

A description for the sequences can be created blasting some databases. Several blast databases can be used in a sequential way. Once a sequence has a blast hit in one of the databases, the description will be build from the description of that hit. For instance let's imagine that we have to annotate 10 sequences with the databases swissprot and tair. Both blasts would be carried out. If we find significant hits for 5 sequences in swissprot those sequence would be already annotated and a description for the 5  remaining sequences will be look for at the second database (tair in this case). The corresponding ngs_backbone analysis is annotate_description.


Configuration parameters
________________________

In the ngs_backbone.conf section Annotation, subsection description_annotation the parameter description_databases should have a list of at least one blast database. The blast databases used should be defined in the corresponding section, see ::`blast-databases`.


.. _go-annotation:

GO annotation
-------------

ngs_backbone annotates the `Gene Ontology`_ ontology terms by using blast2go_. Blast2GO uses the results of a blast nr search to infer the relevant GO terms for every sequence. The corresponding blast analysis is annotate_go.

Output files
____________

The Blast2GO executable creates two files that can be loaded by the graphical Blast2GO interface, the annot and the dat file. Both will be present at the directory annotations/features. Be aware that generating the dat file will require quite memory.


Configuration parameters
________________________

In the ngs_backbone.conf annotation section a go_annotation subsection should be present. The parameters are:

create_dat_file
  It can be set to True or False

java_memory
  This optional parameter is especially important if the creation of the dat file is required.

b2g_properties_file
  Optionaly you can configure blast2go's properties tunning this file.

.. _orf-annotation:

ORF annotation
--------------

To annotate the ORFs found in your sequences just run the ngs_backbone analysis annotate_orf. ESTScan will be used to look for the ORFs. The output will be (after running write_annotation) a couple of files for each input file, one for the DNA and another for the proteins.


Configuration parameters
________________________

In order to run this analysis in the section orf_annotation at the ngs_backbone.conf file the estscan_matrix matrix file should be defined. This is a valid specific matrix file for ESTScan


.. _ssr-annotation:

Microsatellite annotation
-------------------------

The SSRs can be annotated just by running the ngs_backbone analysis annotate_microsatellite. The result of this analysis will be shown in the gff file and in a csv microsatellite file.


.. _ortholog-annotation:

Ortholog annotation
-------------------

ngs_backbone can annotate the orthologs doing a reciprocal blast search. It can be done on one or several blast databases. The ngs_backbone analysis is called ortholog_annotation. The list of orthologs will be found in annotations/features/


Configuration parameters
________________________

In the ngs_backbone.conf section Annotation, subsection ortholog_annotation the parameter ortholog_databases should have a list of at least one blast database. The blast databases used should be defined in the corresponding section, see ::`blast-databases`.


.. _intron-annotation:

cDNA intron annotation
----------------------

When the sequences to annotate are cDNA ngs_backbone can guess where the introns were by using the analysis annotate_introns. To do it it aligns the cDNA with a genomic sequence using the emboss program est2genome. As a shortcut ngs_backbone before running est2genome with the whole genomic sequence it does a blast search to look for the relevant genome region and only after that the est2genome alignment is done.

Configuration parameters
________________________

In the ngs_backbone.conf section Annotation, subsection Cdna_intron_annotation the parameter genomic_db should have one blast database. The blast database used should be defined in the corresponding section, see ::`blast-databases`. Also in the same section the parameter genomic_seqs should have the absolute path to the fasta file with the genomic sequences that make up the employed database.


Annotation statistics
---------------------

Once we have annotated some sequences we can get a summary of the process by running the annotation_stats analysis. It will create one text file for each input annotation file with statistics about: description, microsatellites, SNVs, ORF, GO terms and orhologs.


.. include:: links.txt

