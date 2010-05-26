
Assembly versus mapping
=======================

Once we have a collection of sequence reads there are two different kinds of analysis. If we do not have any previous genomic information we would have to assemble the reads into a genome or transcriptome. Alternatively, if we had genome already available we could map our reads against that genome. Although both analyses could seem to be similar they are very different. To assemble a genome is computationally much costly than to do a mapping, assembling the human genome was a difficult task, re-sequencing and mapping the reads from a new individual is a more amenable task.

The main computational difference is that the typical software used to assemble requires a time that depends on the total reads length squared while the mapping is just lineal with the reads length. For a review take a look at `Sense from sequence reads <http://www.nature.com/nmeth/journal/v6/n11s/abs/nmeth.1376.html>`_, but the take home message is that the assembly is time and memory consuming while the mapping can be done in standard computer.

Also is it important to notice that the read length is a critical parameter for the assemblies. To assemble a genome with repetitive DNA we need reads longer than the longest repeat, otherwise the assemble would not be able to correctly assembly the repeats. Thus, the shorter the reads the more challenging to do the assembly. If we want to assemble a new genome or transcriptome we would have a much easier time if we are dealing with sanger or 454 reads.


