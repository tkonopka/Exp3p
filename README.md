# Exp3p

Exp3p is a suite of java programs for RNA-seq expression analysis.

There exist several programs to evaluate expression from RNA-seq data. Two examples are 
[cufflinks](https://github.com/cole-trapnell-lab/cufflinks) and 
[HTseq](http://www-huber.embl.de/HTSeq/). The Exp3p suite provides similar functionality with some 
additional options that support 3-prime biased RNA-seq.


## Introduction

RNA-seq library preparation protocols typically aim to capture the entire body of mRNA transcripts. 
In such data, the probability of sequencing a read is proportional to the mRNA mass, i.e. it is
higher when there are many transcripts and when a gene is long. An appropriate normalization method
is FPKM - fragments per kilobase of transcript per million reads.

Some library preparation protocols, e.g. [QuantSeq](http://www.nature.com/nmeth/journal/v11/n12/full/nmeth.f.376.html), have different capture properties with the probability of sequencing a read proportional to the number 
of poly-A tails. Expression is expected to be higher when there are many transcripts, but should not 
depend on gene length. An appropriate normalization method therefore includes a division by millions 
of reads, but not by read length. Another consequence of alternative capture is that the distribution 
of reads can be heavily biased toward 3-prime end of transcripts. The Exp3p suite provides options to 
adapt gene and transcript expression calculations to this kind of data.

Features of interest in Exp3p are

* Implementation of FPKM and EPM normalization schemes
* Implementation of 3-prime focused expression estimates
* Support for counting expression on transcript models and gene models defined through refSeq flat files
* Support for counting expression on arbitary regions through bed files
* Modeled low-, mid-, and high- expression estimates for each feature
* Support for multiple samples
* Some support for read filtering
* Fast implementation in java using simple algorithms



## Exp3p programs

The Exp3p suite contains four components. 

### Exp3p eval

The eval tool measures expression on genes, transcripts, or regions of interest by counting 
aligned reads. An example use case is


```
$ java -jar Exp3p.jar eval 
	--bam mybam.bam 
	--sample mysample
	--anno refGene.txt.gz
	--stranded REVERSE 
	--output myoutput.txs.txt.gz
	--outputgene myoutput.gene.txt.gz
```

The program will here scan an input bam file `mybam.bam` and look for expression on genes defined 
in a refGene flat file `refGene.txt.gz`. It will output two files with prefixes `myoutput`, 
one containing expression estimates at the transcript level, the other at the gene level. The output
files will contain columns marked with the sample name, `mysample`. 

*Note:* The functionality of the above command overlaps with [HTseq](http://www-huber.embl.de/HTSeq/) 
and other read-counting programs. Handling of non-unique-mapping reads and output formats are 
slightly different.

Some additional options relevant for QuantSeq data are

```
	--stranded REVERSE
	--tfillnorm 1000
	--filterAR true
```

The `stranded` settings instructs the program to count only reads that align on the strand opposite to
the gene annotation. The `tfillnorm` overrides gene-length normalization and divides expression estimates
by a fixed gene-length of 1000bp instead. The setting 'filterAR' can be used to ignore reads labeled 
with the `AR` code (see below). 



### Exp3p labelArich

The labelArich tool manipulates bam files. It introduces a new boolean tag 'AR' to each read. A nonzero
value of this tag is meant to flag the read as being close to an A-rich or T-rich genomic region, thereby
tagging reads originating by capture via genomic poly-A sequences. An example use case is


```
$ java -jar Exp3p.jar labelArich
	--bam inbam.bam
	--output outbam.bam
	--genome mygenome.fa.gz
	--window 18
	--mincount 12
```

This command will look for reads whose ends are close to genomic regions of length 18 that contain at least 12
Adenine nucleotides (or Thymine nucleotides, depending on strand). The output alignment file will contain tags 'AR' 
associated to each read. This tag can then be used to tune performance in other Exp3p tools, e.g. Exp3p eval
(see above).



### Exp3p callPPA

The callPPA tool is a heuristic peak-finding program. It identifies pileups of reads that may be indications
of transcription end sites. An example use is


```
$ java -jar Exp3p.jar callPPA
	--bam mybam.bam
	--label mysample
	--anno refGene.txt.gz
	--genome mygenome.fa.gz
	--output mysample.ppa.txt.gz
	--depth 5
```

The program will scan an input file and look for peaks of at least five reads. This is not a simple coverage
search as the program will require sharp transitions in coverage, i.e. a steep climb in coverage from one
base to the next.

*Note:* This tool is only relevant for RNA-seq data that is expected to sequence the 3-prime-most ends of mRNA
transcripts, i.e. QuantSeq with T-fill.


### Exp3p richregions

The richregions tool identifies genomic regions with certain nucleotide composition properties. It is a
very simple heuristic tool. An example use case is


```
$ java -jar Exp3p.jar richregions
	--genome mygenome.fa.gz
	--window 20
	--minA 12
	--minG 8
	--output richregions-20-12A-8G.bed
```


This will scan the genome sequence and look for windows of length 20 which contain at least 12 A's and at least 6 G's. The output will consist of a bed table with all the regions of interest.


## Development

The Exp3p suite is written in java using [Netbeans](https://netbeans.org/). The projects uses some third-party libraries:

* [sam-jdk](https://github.com/broadinstitute/picard) Parsing BAM files
* [jopt-simple](http://pholser.github.io/jopt-simple/) Parsing command-line arguments
* [jSeqUtils](https://github.com/tkonopka/jSeqUtils) Toolkit of sequencing- and genomics- related functions.




