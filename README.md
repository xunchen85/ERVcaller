# ERVcaller
1 Introduction
ERVcaller is a tool designed to accurately detect and genotype non-reference unfixed endogenous retroviruses (ERVs) and other transposon elements (TEs) in the human genome using next-generation sequencing (NGS) data. We evaluated the tool using both simulated and benchmark whole-genome sequencing (WGS) datasets. ERVcaller is capable of accurately detecting various TE insertions of any lengths, particularly ERVs. It can be applied to both paired-end and single-end WGS, RNA-Seq or targeted DNA sequencing data. In addition, ERVcaller is capable of detecting the breakpoints at single-nucleotide resolution. It allows for the use of a TE reference library regardless of sequence complexity, such as the entire RepBase database, suggesting that it may be used to detect highly divergent or novel TE insertions. It is easy to install and use from the command line.

Complementary to ERVcaller, other bioinformatics tools designed to detect large deletions, such as Breakdancer, may be used to detect TEs that are present in the human reference genome but not in testing samples.

2 Installation
2.1 Unzip ERVcaller installer
$ tar vxzf ERVcaller_v.1.1.tar.gz

2.2 Installing dependent software
Users need to successfully install the following software separately and make them available in the default path (such as by using the Linux command “export”).
•	BWA-0.7.10
•	Bowtie2
•	Tophat-2.1.1
•	Samtools-1.6 (or later than 1.2)
•	Hydra-0.5.3
•	SE_MEI (Modified version included in the Scripts folder of ERVcaller installer)
•	R-3.3.2 (or higher)

2.3 Databases
•	The human reference genome (hg38 by default; http://hgdownload.soe.ucsc.edu/downloads.html#human) (if you use BAM file as the input, the chromosome IDs of the input BAM file should be identical to those of the human reference genome used, such as “Chr1”, “chr1”, and “1”).
•	One of the following TE reference libraries: 1) a TE reference provided by ERVcaller installer (i.e., the HERV-K reference sequence, the ERV library or the human TE database), or 2) a user-defined TE reference library.

	
2.4 Databases indexing
•	The human reference genome and TE reference library should be indexed separately using BWA “bwa index”.

3 Running ERVcaller 
3.1 Print help page
$ perl user_installed_path/ERVcaller.pl

3.2 TE detection
For detecting and genotyping all ERV insertions: 
$ perl user_installed_path/ERVcaller.pl -i sample_ID -f .bam -H human_genome_hg38.fa -T ERV_library.fa -S -V -G

More examples of ERVcaller running command lines are shown in the help page.

3.3 Parameters
All available parameters are listed below. Four parameters are required, including input_sampleID (-i), file_suffix (-f), Human_reference_genome (-H), and TE_reference_genomes (-T).

Parameter	Full name	Description
-i	input_sampleID	Sample ID (required)

-h	help	Print this help

-t	threads	The number of threads (default: 1)

-f	file_suffix	The suffix of the input data (required; default: .fq.gz)

-d	data_type	Data type, including WGS, and RNA-seq (default: WGS)

-s	sequencing_type	Type of sequencing data, including paired-end, and single-end (default: paired-end)

-H	Human_reference_genome	The FASTA file of the human reference genome (required)

-T	TE_reference_genomes	The TE library (FASTA) used for screening (required)

-l	length_insertsize	Insert size length (bp) (default: 500)

-S	Split	Is the split reads used for detection

-V	Validation	Validation function (input bam file need to be indexed)

-G	Genotyping	Genotyping function (input bam file need to be indexed)

-w	window_size	Window size of selected genomic locations for genotyping (bp) (default: 10,000)
