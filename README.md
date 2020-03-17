# Genome Registration Pipeline

The rapid growth of sequencing data across multiple individuals and populations has driven a paradigm shift 
from linear reference genome assemblies to pangenome representations that incorporate multiple genomic sequences.
Reference genome assemblies play an essential role by providing a coordinate framework for referring to and annotating biological and sequence features. 
However, there are few current standards for assigning coordinates for a graph-based pangenome representations with varying loci.
We propose both a general methodology for genome registration and a graph-based coordinate framework and apply it to an instance of a mouse genetic reference pangenome, the Collaborative Cross Graphical Genome (CCGG). 
Our approach integrates genome assemblies and sequence data from the represented Collaborative Cross population, by
selecting conserved, unique and consistently ordered k-mers as anchors. 
Mapping these anchors back to each linear genome registers and partitions each assembly into disjoint homologous regions consistently.
The anchor sequences establish a sequential ordering of genomic regions. 
Parallel sequence path between anchors reveal the orthogonal sequence structure.
Relative to established anchors, one can assign offsets and designate path names to represent any base position in the graph. 
A general position projection framework between parallel paths is provided based on the anchor coordinate framework.
Anchor sequence based coordinates provide a benefit of isolating genome assembly updates leaving most unaffected annotations unchanged.
Anchor-based genome registration provides a general method for pangenome construction and establishes a graph-based coordinate framework that addresses major issues working with pangenomes
These demoted k-mers help identify assembly errors and potential structural variations in each strain according to their mapping position. 
Joint analysis of k-mer candidates demoted in every linear genome reveals reference-specific genomic features.
