
Annotation
==========

There are different annotation analyses for the sequences, but they all operate in a similar way, so it is worth to explain the general annotation process. The sequences to annotate should be placed in one or several files at annotations/input/. When an annotation is done the results are stored in a kind of database at annotations/repr/. There is a versioned file in annotations/repr/ for each original file set to annotate. The annotations are cumulative, so if we annotate the ORF and after that the SSRs both will be stored at annotations/repr/. At every time after an annotation the real output files can be generated running the write_annotation analysis. The output files will be found at annotations/result/. The output files are:

* `vcf <http://1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcf3.3>`_. It stores the SNP and indel information.
* DNA and pep ORF fasta files. The DNA to translate and the translation.
* a csv file for the SSR (microsatellite) information.
* `gff <http://www.sequenceontology.org/resources/gff3.html>`_. Stores all found sequence features.
* annot and dat files for Blast2Go.
* a text file with the orthologs.


.. _snp-calling:

SNP calling
-----------

SNP calling is a form of sequence annotation. To annotate some sequences you need the input sequences to annotate, a bam file in mapping/results/merged.bam and a reference genome in mapping/reference/reference.fasta. SNP calling is based on the samtools pileup with some filters on top. The ngs_backbone analysis is: annotate_snv.

Configuration parameters
________________________

In the section Snvs:

min_quality
  Minimum allele quality. The alleles with a lower quality will not be consider. It is a phred quality.

min_mapq
  Reads mapped with a lower quality will not be considered. This is the phred quality found on the bam file.

min_num_alleles
  Usually it will be 1 or 2. If 1 the positions with only one allele different from the reference will be considered if 2 only the positions with at least 2 alleles will be considered.

Also a subsection edge_removal inside Snvs could be defined. The parameters in this section are used to when we want to ignore the nucleotides close to the edged of the reads. This parameters are dependent on the platform and are: 454_left, 454_right, sanger_left, sanger_right, illumina_left and illumina_right.

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

ngs_backbone annotates the `Gene Ontology <http://www.geneontology.org/>` ontology terms by using `Blast2GO <http://www.blast2go.org/>`. Blast2GO uses the results of a blast nr search to infer the relevant GO terms for every sequence. The corresponding blast analysis is annotate_go.

Output files
____________

The Blast2GO executable creates two files that can be loaded by the graphical Blast2GO interface, the annot and the dat file. Both will be present at the directory annotations/result. Be aware that generating the dat file will require quite memory.


Configuration parameters
________________________

In the ngs_backbone.conf annotation section a go_annotation subsection should be present. The parameters are:

create_dat_file
  It can be set to True or False

java_memory
  This optional parameter is especially important if the creation of the dat file is required.


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

ngs_backbone can annotate the orthologs doing a reciprocal blast search. It can be done on one or several blast databases. The ngs_backbone analysis is called ortholog_annotation. The list of orthologs will be found in annotations/result/


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


