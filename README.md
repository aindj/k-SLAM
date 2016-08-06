Please contact david.ainsworth08@imperial.ac.uk with any questions or for updates regarding k-SLAM <br />
A web application (beta) is provided at http://www.sbg.bio.ic.ac.uk/~slam for preliminary analysis of small datasets.
# k-SLAM
Introduction
============
k-SLAM is a program for alignment based metagenomic analysis of large sets of high-throughput sequence data.<br />
k-SLAM uses a k-mer based technique to rapidly find alignments between reads and genomes which are then validated using the Smith-Waterman algorithm. Alignments are chained together into a pesudo-assembly to increase specificity. Taxonomy is inferred using a Lowest Common Ancestor technique. Genes and variants can also be found from the alignments and output in SAM format.<br /> 
k-SLAM is fast and highly parallelisable, with speeds of 5 million reads per minute on 150bp paired end data.

Installation
============
Compilation
------------
k-SLAM requires a modern version of gcc that can compile C++11 and a modern version of the boost libraries.
```
unzip k-SLAM-master.zip
cd SLAM/build
make
cd ..
mv build/SLAM ./
export PATH=$PATH:.
```
Database build
------------
k-SLAM is packaged with an installation script which will download the NCBI taxonomy and bacterial/viral genomes and create k-SLAM's index. <br />
To download taxonomy and bacterial/viral genomes and create k-SLAM database in "database_dir":
```
install_slam.sh database_dir bacteria viruses  (can be bacteria or viruses or both)
```

Usage
============
Metagenomics
------------
k-SLAM takes a FASTQ file as its input (plus an R2 FASTQ file for paired data).<br />
k-SLAM is called from the command line with the following command:
```
SLAM --db=DATABASE_DIR --output-file=OUTFILE R1FILE.fastq (optional R2FILE.fastq) 
```
Note: Paired reads must be split into two FASTQ files with reads in the correct order.<br \>
**Required parameters:** <br />

--------------------------------------
| Command | Arguments | Description | Defaults |
| -------- | ------ | ------ | ----- |
--db | string | k-SLAM database directory which reads will be aligned against | N/A
------------------------------------

**Optional parameters:** <br />

--------------------------------------
| Command | Arguments | Description | Defaults |
| -------- | ------ | ------ | ----- |
--output-file | string | write to this file instead of stdout | stdout  
--min-alignment-score | string | alignment score cutoff | 0  
--score-fraction-threshold | float | For each read, screen alignments with scores less than this fraction of the score of the best alignment | 0.95  
--num-reads | +ve int | Number of reads from R1/R2 File to align | all  
--num-reads-at-once | +ve int | Reduce RAM usage by only analysing "arg" reads at once, this will increase execution time | 10,000,000  
--sam-file | string | write SAM output to this file | none  
--num-alignments | +ve int | Number of alignments (or alignment pairs) to report per read in SAM file | 10 
--sam-xa | N/A | only output primary alignment lines, use XA field for secondary alignments | N/A  
--no-pseudo-assembly | N/A | do not link alignments together | N/A 
--match-score | +ve int | Smith-Waterman match score | 2  
--mismatch-penalty | +ve int | Smith-Waterman mismatch penalty | 3  
--gap-open | +ve int | Smith-Waterman gap opening penalty | 5  
--gap-extend | +ve int | Smith-Waterman gap extend penalty | 2
------------------------------------

Alignment
------------
If metagenomic analysis is not required and only alignment is needed, then the flag "--just-align" can be used. The option "--sam-file" must be used. If a database was built from FASTA files, then the "--just-align" flag must be used.

Custom Database Build
--------------
k-SLAM can parse any genomes in the GenBank flat file format<br />
Create database directory:
```
mkdir custom_db && cd custom_db
```
Build k-SLAM's taxonomy databases using the NCBI taxonomy nodes.dmp and names.dmp
```
SLAM --parse-taxonomy names.dmp nodes.dmp --output-file taxonomy
```
Build k-SLAM's index from any number of Genbank files
```
SLAM --output-file database --parse-genbank file1.gbk file2.gbk file3.gbk ...
```
k-SLAM can also align reads to genomes in FASTA format (only alignment and not metagenomics): <br />
Create database directory
```
mkdir custom_db && cd custom_db
```
Build k-SLAM's index from any number of FASTA files
```
SLAM --output-file database --parse-fasta file1.fa file2.fa file3.fa ...
```
Note: databases produced from FASTA files must be analysed using the "--just-align" flag.

System Requirements
===================
CPU
-------------------
k-SLAM is designed to be run on a parallel platform, ideally with 8 or more cores. It will however run (albeit slower) on a system with fewer cores.
Memory
-------------------
k-SLAM's k-mer list sort algorithm is fast but very memory intensive. For a dataset of 10 million paired reads aligning against the NCBI bacterial genomes, k-SLAM requires around 50GB RAM. For much larger datasets, k-SLAM can split them into smaller subsets (using option --num-reads-at-once see table above) which are analysed sequentially (still producing only one set of output files).

Output
====================
k-SLAM outputs in several formats:
Summary XML
--------------------
Each identified taxon is listed in descending order of the number of reads assigned to it. Because of the use of a Lowest Common Ancestor method, taxons may be at any rank. Each read is assigned to a maximum of one taxon (a genus entry will not contain reads that have been mapped to individual species within that genus but only reads that mapped to several species whose LCA was that genus).<br />
Each taxon has the following tags: abundance (number of reads and percentage of total reads), taxonomyID, lineage (from NCBI taxonomy), name, genes and reads.<br />
Note: Output has to be in XML format therefore any annotations of genes etc will have the characters \textless , \textgreater , & ,'  and " replaced with the relevant entity reference.<br />
For each taxon, the genes found are listed using the gene tag. A maxinum of one gene is inferred for each aligned read, based on its position on the genome. The "count" field describes the number of reads that overlapped with that particular gene. The protein, locus, product, GeneID, reference sequence are listed (using NCBI data) along with the cds range. The same gene may appear in multiple taxons.

Tab Separated Taxonomy
----------------------
Mapped reads are listed with their LCA taxonomy. Unmapped reads are not listed.
Abbreviated Taxonomy
----------------------
Each identified taxon is listed along with the percentage of reads that mapped to it.
SAM
----------------------
Output in the Sequence Alignment/Map format. A Bowtie style output is used (for each read there is a primary line for each reference that it aligned to). A BWA style XA tag can be printed (using the --sam-xa parameter) instead of listing all hits. Unmapped reads are not printed except when they belong to a pair where the other read was mapped. <br />
The following k-SLAM specific tags are used:<br />

1. XS: alignment score assigned by k-SLAM, using pseudo assembly.
2. XO: number of hits for this segment.
3. XT: taxonomy ID of this reference.
4. XG: gene at this position in the reference.
5. XP: protein ID of this gene.
6. XR: product of this gene.
7. XA: BWA style alternate hits in format (chr,pos,CIGAR,NM;)* (this tag is optional).
