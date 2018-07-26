#!usr/bin/env perl
#
# Author: 	Xun Chen
# Email: 	Xun.Chen@uvm.edu
# Date:		05/25/18
# Version:	v0.1.0

use Getopt::Long qw(:config no_ignore_case);
use strict;
use warnings;
use Cwd qw();

# ERVcaller v0.1.0:
# 05/25/2018
#
# Updates (v0.1.1):
#   07/08/2018:		Corrected some minor errors;
#   07/08/2018:		Recovered the potential TE insertions in repeat regions as a separate group in the output file;
#   07/08/2018:		Recovered the split reads (< 20 bp) to improve the breakpoint detection;
#   07/14/2018:		Corrected the errors for detecting and genotyping TE insertions using single-end sequencing data;
#   07/14/2018:		Re-organize the output format;
#   07/17/2018:		Corrected bugs about checking input files;
#
# Questions or help:
#
# For questions or support contact Xun Chen (xunchen85@gmail.com or Xun.Chen@uvm.edu);

my $bowtie2_d="";
my $tophat_d="";
my $bwa_d="";
my $samtools_d="";
my $SE_MEI_d="";
my $hydra_d="";

### Except R that you need to load in the enviroment, for other tools, you can either load in the enviroment (export), or specify the path below (remove the "#") e.g.:
# $bwa_d="/ERVcaller/Tools/bwa-master/";
# $samtools_d="/ERVcaller/Tools/samtools-1.6/";
# $tophat_d="/ERVcaller/Tools/tophat-2.1.1.Linux_x86_64/";
# $bowtie2_d="/ERVcaller/Tools/bowtie2-2.2.7/";
# $SE_MEI_d="/ERVcaller/Tools/SE-MEI/";
# $hydra_d="/ERVcaller/Tools/Hydra-Version-0.5.3/";
 
# Parameter variables;
my $help;
my $input_sampleID;
my $file_suffix=".fq.gz";
my $threads=1;
my $data_type="WGS";
my $sequencing_type="paired-end";
my $Human_reference_genome="";
my $TE_reference_genomes="";
my $length_insertsize=500;
my $Split;
my $Validation;
my $Genotyping;
my $window_size=10000;
my $directory=$0;

# Variables
GetOptions(
           "i|input_sampleID=s" => \$input_sampleID,
           "h|help" => \$help,
           "t|threads=i" => \$threads,
           "f|file_suffix=s" => \$file_suffix,
           "d|data_type=s" => \$data_type, 
           "s|sequencing_type=s" => \$sequencing_type,
	   "H|Human_reference_genome=s" => \$Human_reference_genome,
	   "T|TE_reference_genomes=s" =>\$TE_reference_genomes,
           "l|length_insertsize=i" =>\$length_insertsize,
           "S|Split" => \$Split,
	   "V|Validation" => \$Validation,
	   "G|Genotyping" => \$Genotyping,
           "w|window_size=i" =>\$window_size
           );
print "\n";
################################## Step 1
print "                 ########  #######   ##       ##                        ##     ##		\n";
print "                 ########  ##    ##  ##       ##                        ##     ##		\n";
print "                 ##        ##    ##   ##     ##                         ##     ##		\n";
print "                 ##        ##   ##    ##     ##                         ##     ##		\n";
print "                 ########  ######      ##   ##      #####     ####      ##     ##       		\n";
print "                 ########  #####       ##   ##     ##   ##   ##  ##     ##     ##                \n";
print "                 ##        ## ##        ## ##     ##        ##    ##    ##     ##           	\n";
print "                 ##        ##  ##       ## ##     ##        ##    ##    ##     ##           	\n";
print "                 ########  ##   ##       ###       ##   ##   ##   ##    ## ##  ## ##         	\n";
print "                 ########  ##    ##      ###        #####     #### ##    ##     ##           	\n";



print "\n\n# ERVcaller v0.1.1\n";

print "# Please contact Xun Chen Ph.D. for questions and help:\n# Email: xunchen85\@gmail.com or Xun.Chen\@uvm.edu\n\n";

print "\nStep 1: Loading...\n=====================================\n";
if(!defined($input_sampleID) & !defined($help)) {
  prtErr("Error:\n# No samples are provided");
  prtUsa();
  exit;
}

if(defined($help)){
  prtUsa();
  exit;
}

sub prtUsa {
  print "\n=====================================\nperl $0 [arguments]\n";
  print "\n";
  print "# Arguments:\n";
  print "       -i|input_sampleID		Sample ID (required)\n";
  print "       -h|help				Print this help\n";
  print "       -t|threads			The number of threads will be used (default: 1)\n";
  print "       -f|file_suffix			The suffix of the input data, including: zipped FASTQ file (i.e., .fq.gz, and fastq.gz), unzipped FASTQ file (i.e., .fq, and fastq), and BAM file (.bam) (required; default: .fq.gz;)\n";
  print "       -d|data_type			Data type, including: WGS, RNA-seq (default: WGS)\n";
  print "       -s|sequencing_type		Type of sequencing data, including: paired-end, single-end (default: paired-end)\n";
  print "       -H|Human_reference_genome	The FASTA file of the human reference genome (required)\n";
  print "       -T|TE_reference_genomes		The TE library (FASTA) used for screening (required)\n";
  print "       -l|length_insertsize		Insert size length (bp) (default: 500)\n";
  print "       -S|Split			Is the split reads used\n";
  print "       -V|Validation			Validation function (input bam file need to be indexed)\n";
  print "       -G|Genotyping			Genotyping function (input bam file need to be indexed)\n";
  print "       -w|window_size			Window size of selected genomic locations for genotyping (bp) (default: 5000)\n";                 ############# ??????????????
  print "\n\n";
  print "=====================================\nExamples for detecting ERV and other TE insertions:\n";
  print "# WGS data in paired-end fastq format:\n";
  print "       perl ERVcaller.pl -d WGS -i TE_seq.fa_30_100_500 -f .fq.gz -s paired-end -H human_genome_hg38.fa -T ERV_library.fa -t 12 -l 500 -w 10000\n";
  print "# WGS data in single-end fastq format:\n";
  print "       perl ERVcaller.pl -d WGS -i TE_seq.fa_30_100_500_1 -f .fq.gz -s single-end -H human_genome_hg38.fa -T ERV_library.fa -t 12 -l 500 -w 10000\n";
  print "# WGS data in paired-end bam format:\n";
  print "       perl ERVcaller.pl -d WGS -i TE_seq.fa_30_100_500 -f .bam -s paired-end -H human_genome_hg38.fa -T ERV_library.fa -t 12\n\n\n";
  print "# WGS data in paired-end bam format with validation and genotyping functions (indexed input bam file is needed for either validation or genotyping function):\n";
  print "       perl ERVcaller.pl -d WGS -i seq -f .bam -s paired-end -H human_genome_hg38.fa -T ERV_library.fa -S -V -G -t 12\n\n\n";
}

sub prtErr{
  print "\n=====================================\n@_\n\n";
}

# extract the path of the installed ERVcaller software
my @directory=split /\//, $directory;
if($#directory>0){
 $directory=join("/",@directory[0..$#directory-1]);
 $directory=$directory."/";
                 }
else {
  $directory="";
}
 my $chr_name="NA";
 my $human_genome=$Human_reference_genome;
 my $human_genome_tophat=$Human_reference_genome;
 my $order=0;
 my $thread_1=$threads-1;
 my $double_length_insertsize;
 my $cmd="";
 my $header1="Sample_ID Is_Split_mode Is_genotyped Is_validated Human_ref. TE_reference TE_sequence_name Chr. Start End Upstream_breakpoint_on_human Downstream_breakpoint_on_human Upstream_breakpoint_on_TE Downstream_breakpoint_on_TE Information_both_breakpoints Insertion_site Group Average_AS_for_chimeric_and_improper_reads_on_human Average_XS_for_chimeric_and_improper_reads_on_human Maximum_AS_for_chimeric_and_improper_reads_on_human No._supporting_reads Average_AS_supporting_reads\n";
 my $header2="Sample_ID Is_Split_mode Is_genotyped Is_validated Human_ref. TE_reference TE_sequence_name Chr. Start End Upstream_breakpoint_on_human Downstream_breakpoint_on_human Upstream_breakpoint_on_TE Downstream_breakpoint_on_TE Information_both_breakpoints Insertion_site Group Average_AS_for_chimeric_and_improper_reads_on_human Average_XS_for_chimeric_and_improper_reads_on_human Maximum_AS_for_chimeric_and_improper_reads_on_human No._supporting_reads Average_AS_supporting_reads No._chimeric_and_improper_reads V_True_chimeric_and_improper_reads V_False_chimeric_and_improper_reads_PE V_False_chimeric_and_improper_reads_others No._split_reads_(>=20bp) V_True_split_reads_(>=20bp)\n";
 my $header3="Sample_ID Is_Split_mode Is_genotyped Is_validated Human_ref. TE_reference TE_sequence_name Chr. Start End Upstream_breakpoint_on_human Downstream_breakpoint_on_human Upstream_breakpoint_on_TE Downstream_breakpoint_on_TE Information_both_breakpoints Insertion_site Group Average_AS_for_chimeric_and_improper_reads_on_human Average_XS_for_chimeric_and_improper_reads_on_human Maximum_AS_for_chimeric_and_improper_reads_on_human No._supporting_reads Average_AS_supporting_reads No._chimeric_and_improper_reads V_True_chimeric_and_improper_reads V_False_chimeric_and_improper_reads_PE V_False_chimeric_and_improper_reads_others No._split_reads_(>=20bp) V_True_split_reads_(>=20bp) No._split_reads_(<20bp) Geno_No._reads_supporting_non_TE_insertions Geno_Read_depth_of_the_genomic_window_(Mean) Geno_Read_depth_of_the_genomic_window_(SD) Geno_No._random_locations_of_the_genomic_window Geno_Read_depth_of_the_genomic_window_(Quantile=0.1) Genotype\n";
################################# Step 2
# 2.1 Check input file 
 if (-d ${input_sampleID}."_temp") {
  system ("rm -rf ${input_sampleID}_temp");
  system ("mkdir ${input_sampleID}_temp");
 } else {
  system ("mkdir ${input_sampleID}_temp");
 }
print "\nStep 2: Detecting TE insertions...\n=====================================\n";
if(defined($input_sampleID) && -e ${input_sampleID}."_1".${file_suffix} && -e ${input_sampleID}."_2".${file_suffix} && $sequencing_type eq "paired-end" && ($file_suffix =~ "fq" || $file_suffix =~ "fastq")){
  print "~~~~~ paired-end reads in fastq format were loaded\n";
 } elsif (defined($input_sampleID) && -e ${input_sampleID}.${file_suffix} && $sequencing_type eq "single-end" && ($file_suffix =~ "fq" || $file_suffix =~ "fastq")){
  print "~~~~~ single-end read in fastq format was loaded\n";
 } elsif(defined($input_sampleID) && -e ${input_sampleID}.${file_suffix} && $sequencing_type eq "paired-end" && ($file_suffix =~ "bam" || $file_suffix =~ "sam")){
  print "~~~~~ paired-end reads in bam format were loaded\n";
  if ($Genotyping || $Validation){
    if (-e ${input_sampleID}.${file_suffix}.".bai"){
     print "~~~~~ the input bam file was indexed\n";
    } else {
     print "~~~~~ the input bam file was not indexed, please index the bam file using samtools for performing the validation or genotyping function\n";
     exit;
    }
  }
 } elsif (defined($input_sampleID) && -e ${input_sampleID}.${file_suffix} && $sequencing_type eq "single-end" && ($file_suffix =~ "bam" || $file_suffix =~ "sam")){
  print "~~~~~ single-end reads in bam format were loaded\n";
  if ($Genotyping || $Validation){
    if (-e ${input_sampleID}.${file_suffix}.".bai"){
     print "~~~~~ the input bam file was indexed\n";
    } else {
     print "~~~~~ the input bam file was not indexed, please index the bam file using samtools for performing the validation or genotyping function\n";
     exit;
    }
  }
 } else {
  prtErr("# Error: cound not find the input data under the provided sampleID\n");
  prtUsa();
  exit;
 }

# 2.2 Extract supporting reads
 if(${file_suffix} =~ "bam"){
  $order=1;
  convert_bamtofastq(${input_sampleID});
  align_to_hg(${input_sampleID}."_h1",".1fq");
  $order=2;
  convert_bamtofastq(${input_sampleID}."_h1");
  system("mv ${input_sampleID}_h1_sm.bam ${input_sampleID}_sm.bam");
  system("mv ${input_sampleID}_h1_su.bam ${input_sampleID}_su.bam");
  if($Split || $sequencing_type eq "single-end"){
   system("gunzip -c ${input_sampleID}_soft.fastq.gz >${input_sampleID}_1soft.fastq");
   system("gunzip -c ${input_sampleID}_h1_soft.fastq.gz >>${input_sampleID}_1soft.fastq");
  }
  if($sequencing_type eq "paired-end"){
   system("mv ${input_sampleID}_h1_h1_1.1fq ${input_sampleID}_1.1fq");
   system("mv ${input_sampleID}_h1_h1_2.1fq ${input_sampleID}_2.1fq");
  } else {
   system("mv ${input_sampleID}_h1_h1.1fq ${input_sampleID}.1fq");
  }
 } else {
  align_to_hg(${input_sampleID},${file_suffix});
  convert_bamtofastq(${input_sampleID});
  if($Split || $sequencing_type eq "single-end"){
   system("gunzip -c ${input_sampleID}_soft.fastq.gz >${input_sampleID}_1soft.fastq");
  }
  if($sequencing_type eq "paired-end"){
   system("mv ${input_sampleID}_h1_1.1fq ${input_sampleID}_1.1fq");
   system("mv ${input_sampleID}_h1_2.1fq ${input_sampleID}_2.1fq");
  } else {
   system("mv ${input_sampleID}_h1.1fq ${input_sampleID}.1fq");
  }
 }

 if($sequencing_type eq "single-end" || $Split){
  system ("perl ${directory}Scripts/Soft_clipping_filter.pl -length 20 -file ${input_sampleID}_1soft.fastq -o ${input_sampleID}");
  system ("rm ${input_sampleID}_1soft.fastq");
 }

# 2.3 Chimeric reads amd Split reads
print "\nChimeric and split reads...\n=====================================\n";
 if($sequencing_type eq "paired-end"){
  system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $TE_reference_genomes ${input_sampleID}_1.1fq ${input_sampleID}_2.1fq >${input_sampleID}_vsu.sam");
 } else{        
  system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $TE_reference_genomes ${input_sampleID}.1fq >${input_sampleID}_vsu.sam");
 }
 system ("${samtools_d}samtools view -bS -@ $thread_1 ${input_sampleID}_vsu.sam >${input_sampleID}_vsu.bam");
 system ("${samtools_d}samtools sort -@ $thread_1 ${input_sampleID}_vsu.bam -o ${input_sampleID}_vsu.sort.bam");
 system ("${samtools_d}samtools view -@ $thread_1 ${input_sampleID}_vsu.sort.bam >${input_sampleID}_vsu.sort.sam");
 system ("rm ${input_sampleID}_vsu.bam");
 system ("rm ${input_sampleID}_vsu.sam");
 system ("touch ${input_sampleID}_all_breakpoint");

 if($sequencing_type eq "single-end" || ($sequencing_type eq "paired-end" && $Split)){
  system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $TE_reference_genomes ${input_sampleID}_1sf.fastq >${input_sampleID}_vsoft.sam");
  system ("${samtools_d}samtools view -bS -@ $thread_1 ${input_sampleID}_vsoft.sam >${input_sampleID}_vsoft.bam");
  system ("${samtools_d}samtools sort -@ $thread_1 ${input_sampleID}_vsoft.bam -o ${input_sampleID}_vsoft_sort.bam");
  system ("${samtools_d}samtools view -@ $thread_1 ${input_sampleID}_vsoft_sort.bam >${input_sampleID}_vsoft_sort.sam");
  system ("rm ${input_sampleID}_vsoft.bam");
  system ("perl ${directory}Scripts/Soft_clipping_transfer.pl -f ${input_sampleID}_vsoft_sort.sam -o ${input_sampleID}_vsoft_breakpoint");
  system ("cat ${input_sampleID}_vsoft_breakpoint >>${input_sampleID}_all_breakpoint");
  system ("rm ${input_sampleID}_vsoft_breakpoint");
  system ("rm ${input_sampleID}_vsoft_sort.sam");
  system ("rm ${input_sampleID}_vsoft.sam");
 }

 if($sequencing_type eq "paired-end"){ 
  system ("${samtools_d}samtools view ${input_sampleID}_sm.bam >${input_sampleID}_sm.sam"); 
  open SM,"${input_sampleID}_sm.sam";
  open TYPE,">${input_sampleID}.type";
  my $sm_1="";
  my $as="";
  my $xs="";
  while($sm_1=<SM>){
   my @sm_1=split /\s+/, $sm_1;
   $as="NA";$xs="NA";
  unless($sm_1[2] eq "*"){
   for(my $i=11;$i<@sm_1;$i++){
    if($sm_1[$i]=~s/AS:i://){$as=$sm_1[$i];}
    if($sm_1[$i]=~s/XS:i://){$xs=$sm_1[$i];}
                              }
                         }
  if($sm_1[1]%256>=128){print TYPE "$sm_1[0] L $as $xs $sm_1[5]\n";}
  else{print TYPE "$sm_1[0] R $as $xs $sm_1[5]\n";}
 }
 close TYPE;
 close SM;
  system ("perl ${directory}Scripts/Break_point_calling.pl -type ${input_sampleID}.type -position ${input_sampleID}_sm.sam -TE ${input_sampleID}_vsu.sort.sam -o ${input_sampleID}");
  system ("rm ${input_sampleID}_sm.sam");
  system ("cat ${input_sampleID}_breakpoint >>${input_sampleID}_all_breakpoint");
  system ("rm ${input_sampleID}_breakpoint");
 }
 system ("perl ${directory}Scripts/Reads_summary.pl -i ${input_sampleID}_all_breakpoint -o ${input_sampleID}_summary");
 system ("rm ${input_sampleID}_all_breakpoint");
 system ("perl ${directory}Scripts/Order_by_TE_sequence.pl -i ${input_sampleID}_summary -o ${input_sampleID}.integration");
 $double_length_insertsize=$length_insertsize * 2;
 system ("perl ${directory}Scripts/Filtered_fastq.pl $input_sampleID");
 if ($sequencing_type eq "single-end"){
  system ("rm ${input_sampleID}_1.1fuq ${input_sampleID}_2.1fuq");
 }

# 2.4 Improper reads
print "\nImproper reads...\n=====================================\n";
 if($sequencing_type eq "paired-end"){
 if ($file_suffix =~ "bam"){
  system ("${samtools_d}samtools view -F 14 -b -@ $thread_1 ${input_sampleID}$file_suffix >${input_sampleID}_ERV.bam");
 } elsif ($file_suffix =~ "sam") {
  system ("${samtools_d}samtools view -F 14 -b -@ $thread_1 ${input_sampleID}$file_suffix >${input_sampleID}_ERV.bam");
 } else {
  system ("${samtools_d}samtools view -F 14 -b -@ $thread_1 ${input_sampleID}_h.sam >${input_sampleID}_ERV.bam");
 }
 system ("${samtools_d}samtools view ${input_sampleID}_ERV.bam >${input_sampleID}_ERV.sam");
 system ("${samtools_d}samtools sort -n -@ $thread_1 ${input_sampleID}_ERV.bam -o ${input_sampleID}_ERV.sort.bam");
 system ("${samtools_d}samtools fastq -@ $thread_1 -N ${input_sampleID}_ERV.sort.bam -1 ${input_sampleID}_ERV1_1.1fq -2 ${input_sampleID}_ERV1_2.1fq");
 system ("perl ${directory}Scripts/Check_paired_end.pl -s ${input_sampleID}_ERV1 -f .1fq");
 system ("mv ${input_sampleID}_ERV1_1.1fq2 ${input_sampleID}_ERV1_1.1fq");
 system ("mv ${input_sampleID}_ERV1_2.1fq2 ${input_sampleID}_ERV1_2.1fq");
 align_to_hg(${input_sampleID}."_ERV1",".1fq");
 system ("perl ${directory}Scripts/ERV_get_name.pl -f ${input_sampleID}_ERV1_h.sam -o ${input_sampleID}_ERV1_h.sam2 -qc 22");
 system ("perl ${directory}Scripts/ERV_filter_reads.pl ${input_sampleID}_ERV1_h.sam2 >${input_sampleID}_ERV1.bian");
 system ("perl ${directory}Scripts/ERV_get_reads.pl ${input_sampleID}_ERV1.bian ${input_sampleID}_ERV1_1.1fq ${input_sampleID}_ERV1_2.1fq ${input_sampleID}_ERV_1.1fq ${input_sampleID}_ERV_2.1fq");
 system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $TE_reference_genomes ${input_sampleID}_ERV_1.1fq ${input_sampleID}_ERV_2.1fq >${input_sampleID}_ERV_vsu.sam");
 system ("${samtools_d}samtools view -bS -@ $thread_1 ${input_sampleID}_ERV_vsu.sam >${input_sampleID}_ERV_vsu.bam");
 system ("${samtools_d}samtools sort -@ $thread_1 ${input_sampleID}_ERV_vsu.bam -o ${input_sampleID}_ERV_vsu.sort.bam");
 system ("${samtools_d}samtools view -@ $thread_1 ${input_sampleID}_ERV_vsu.sort.bam >${input_sampleID}_ERV_vsu.sort.sam");
 system ("rm ${input_sampleID}_ERV_vsu.bam");
 system ("rm ${input_sampleID}_ERV_vsu.sam");
 system ("perl ${directory}Scripts/ERV_get_name.pl -f ${input_sampleID}_ERV_vsu.sort.sam -b ${input_sampleID}_ERV1.bian -h ${input_sampleID}_ERV.hf -o ${input_sampleID}_ERV_vsu.sort.sam2");
 system ("perl ${directory}Scripts/ERV_organize_reads.pl ${input_sampleID}_ERV1_h.sam2 ${input_sampleID}_ERV_vsu.sort.sam2 ${input_sampleID}_ERV.hf >${input_sampleID}_ERV_breakpoint");
 system ("perl ${directory}Scripts/Reads_summary.pl -i ${input_sampleID}_ERV_breakpoint -o ${input_sampleID}_ERV_summary1");
 system ("cat ${input_sampleID}_summary ${input_sampleID}_ERV_summary1 >${input_sampleID}_ERV_summary");
 } else {
  system ("cat ${input_sampleID}_summary >${input_sampleID}_ERV_summary");
 }
 system ("perl ${directory}Scripts/Order_by_TE_sequence.pl -i ${input_sampleID}_ERV_summary -o ${input_sampleID}_ERV.integration");
 $double_length_insertsize=$length_insertsize * 2;
 system ("perl ${directory}Scripts/Filtered_single_reads.pl -i ${input_sampleID}_ERV.integration -o ${input_sampleID}_ERV.integration2 -r $double_length_insertsize");
 system ("mv ${input_sampleID}_ERV.integration2 ${input_sampleID}_ERV.integration");
 system ("perl ${directory}Scripts/ERV_Filtered_fastq.pl ${input_sampleID}");
 if ($sequencing_type eq "paired-end"){
  system ("cat ${input_sampleID}_ERV_1.2fq ${input_sampleID}_1.1fuq >${input_sampleID}_ERV_1.1fuq");
  system ("cat ${input_sampleID}_ERV_2.2fq ${input_sampleID}_1.1fuq >${input_sampleID}_ERV_2.1fuq");
 }
 system ("cp ${input_sampleID}_1sf.fuq ${input_sampleID}_ERV_1sf.fuq");
 system ("perl ${directory}Scripts/Filtered_integration.pl ${input_sampleID}_ERV >${input_sampleID}_ERV.21");
 system ("perl ${directory}Scripts/Remove_redundancy_split_reads.pl ${input_sampleID}_ERV");
 system ("rm ${input_sampleID}_ERV.21");
 system ("mv ${input_sampleID}_ERV.2 ${input_sampleID}_ERV.3");
 system ("perl ${directory}Scripts/ERV_Assign_reads-filter.pl -i ${input_sampleID}_ERV.3 -o ${input_sampleID}_ERV.f -r $double_length_insertsize -m f");
 system ("perl ${directory}Scripts/Result_filtered.pl ${input_sampleID}_ERV.f >${input_sampleID}_ERV_f");
 system ("perl ${directory}Scripts/ERV_Result_visual.pl ${input_sampleID}_ERV_f ${input_sampleID}_ERV 12 80 >${input_sampleID}_ERV_f2");
 system ("perl ${directory}Scripts/Result_finalize3.pl ${input_sampleID}_ERV_f2 >${input_sampleID}_ERV_f22");
 system ("perl ${directory}Scripts/Results_get.pl ${input_sampleID}_ERV_f22 >${input_sampleID}_ERV.TE_f");
 if ($sequencing_type eq "paired-end") {
  $cmd = q(awk '{if(($15>=50)&&($24>0)&&(($30>2 && $20>=30) || ($30==2 && $20>=50)))print$0}');
 } else {
  $cmd = q(awk '{if(($15>=50)&&($30>2 && $20>=30))print$0}');                            
 }
 system ("$cmd ${input_sampleID}_ERV.TE_f >${input_sampleID}_ERV.TE_f2");
 system ("perl ${directory}Scripts/Extract_specific_loci_final_visualization.pl ${input_sampleID}_ERV_f2 ${input_sampleID}_ERV.TE_f2 >${input_sampleID}_ERV.visualization");
 system ("perl ${directory}Scripts/Fine_mapped.pl ${input_sampleID}_ERV.visualization >${input_sampleID}_ERV.fine_mapped");

# 2.5 Organizing output files
 my @path_current2=split /\//, ${input_sampleID};
 my $path_current2;
 if($#path_current2>0){
  my $path_current2=join("/",@path_current2[0..$#path_current2-1]);
  my $input3=$path_current2[$#path_current2];
                      }
 else {
  $path_current2= Cwd::cwd();
  $path_current2=$path_current2."/";
 }
 out1(${input_sampleID}."_ERV");
# 2.6 Clean files
 if(-e ${input_sampleID}."_ERV_f"){system ("rm ${input_sampleID}_ERV_f");}
 if(-e ${input_sampleID}."_ERV_f2"){system ("rm ${input_sampleID}_ERV_f2");}
 if(-e ${input_sampleID}."_ERV_f22"){system ("rm ${input_sampleID}_ERV_f22");}
 if(-e ${input_sampleID}."_ERV.f"){system ("rm ${input_sampleID}_ERV.f");}
 if(-e ${input_sampleID}."_ERV.integration"){system ("rm ${input_sampleID}_ERV.integration");}
 if(-e ${input_sampleID}."_ERV_summary"){system ("rm ${input_sampleID}_ERV_summary");}
 if(-e ${input_sampleID}."_ERV_summary1"){system ("rm ${input_sampleID}_ERV_summary1");}
 if(-e ${input_sampleID}."_ERV_h.sam"){system ("rm ${input_sampleID}_ERV1_h.sam");}
 if(-e ${input_sampleID}."_ERV_h.sam2"){system ("rm ${input_sampleID}_ERV1_h.sam2");}
 if(-e ${input_sampleID}."_ERV1_1.1fq"){system ("rm ${input_sampleID}_ERV1_1.1fq");}
 if(-e ${input_sampleID}."_ERV1_2.1fq"){system ("rm ${input_sampleID}_ERV1_2.1fq");}
 if(-e ${input_sampleID}."_ERV_fasta"){system ("rm ${input_sampleID}_ERV.fasta");}
 if(-e ${input_sampleID}."_h1_h1_h.sam"){system ("rm ${input_sampleID}_h1_h1_h.sam");}
 if(-e ${input_sampleID}."_h1_h.sam"){system ("rm ${input_sampleID}_h1_h.sam");}
 if(-e ${input_sampleID}.".fasta"){system ("rm ${input_sampleID}.fasta");}
 if(-e ${input_sampleID}."_summary"){system ("rm ${input_sampleID}_summary");}
 if(-e ${input_sampleID}."_vsu.sort.sam"){system ("rm ${input_sampleID}_vsu.sort.sam");}
 if(-e ${input_sampleID}."_soft.fastq.gz"){system ("rm ${input_sampleID}_soft.fastq.gz");}
 if(-e ${input_sampleID}.".integration"){system ("rm ${input_sampleID}.integration");}
 if(-e ${input_sampleID}."_h1_soft.fastq.gz"){system ("rm ${input_sampleID}_h1_soft.fastq.gz");}
 if(-e ${input_sampleID}."_h1_1.1fq"){system ("rm ${input_sampleID}_h1_1.1fq");}
 if(-e ${input_sampleID}."_h1_2.1fq"){system ("rm ${input_sampleID}_h1_2.1fq");}


################################## Step 3. Reciprocal alignemnt to the candidate genomic regions
if ($Validation){
 print "\nStep 3: Validation...\n=====================================\n";
 if (($file_suffix =~ "fq" || $file_suffix =~ "fastq") && !(-e ${input_sampleID}."_h.bam.bai")) {
  print "# Converting SAM to BAM file, and then Sort and index the BAM file......\n";
  system ("${samtools_d}samtools view -b -@ $thread_1 ${input_sampleID}_h.sam >${input_sampleID}_h.bam");
  system ("${samtools_d}samtools sort -@ $thread_1 ${input_sampleID}_h.bam -o ${input_sampleID}_h.sort.bam");
  system ("mv ${input_sampleID}_h.sort.bam ${input_sampleID}_h.bam");
  system ("${samtools_d}samtools index ${input_sampleID}_h.bam");
 }
 open INPUT, "${input_sampleID}_ERV.output";           ##### output file
 open OUTPUT,">${input_sampleID}_ERV.output2";
 print OUTPUT "$header2";
 my $line="";
 my $line2="";
 my $seq="";
 my %genome=();
 my @line=();
 my $title="";
 my %TE_f2=();
 my %countf_c=();
 my %countf_c2=();
 my %countf_sf=();
 my $count_c=();
 my $count_sf=();
 
# 3.1 Read human genome
 open GENOME2, "$Human_reference_genome";
 while($line2=<GENOME2>){ 
  @line=split /\s+/, $line2;
  unless(@line){next;}
  $line=$line[0];
  chomp($line); 
  if($line=~">"){
    $title=$line;
    $title=~s/>//;
    $title=~s/chr//;
    $title=~s/Chr//;
    $genome{$title}="";
                }
  else{
    $genome{$title}.=$line;
      }
                }
 
# 3.2 Read TE_f2 file
 open TE_F2,"${input_sampleID}_ERV.TE_f2";
 while ($line2=<TE_F2>){
  @line=split /\s+/, $line2;
  $line[15]=~s/_ERV//;
  my $name1=$line[15]."%".$line[31]."%".$line[32]."%".$line[33];
  my @array=split /\|/, $line[30];
  @{$TE_f2{$name1}}=@array;
                       }

# 3.3 Read output file
 system ("mkdir ${input_sampleID}_subgenome");
 while($line2=<INPUT>){
  %countf_c=();
  %countf_c2=();
  %countf_sf=();
  $count_c=0;
  $count_sf=0;
  @line=split /\s+/, $line2;
  if ($line[0] =~ "Sample_ID"){next;}
  my $seq_name=$line[0]."%".$line[7]."%".$line[8]."%".$line[9];
  open SUB_REF, ">${input_sampleID}_subgenome/${seq_name}.fa";
  my $length=abs($line[9]-$line[8])+2000;
  $seq=substr $genome{$line[7]},$line[8]-1000,$length;
  print SUB_REF ">$seq_name\n$seq\n";
  system ("${bwa_d}bwa index -a bwtsw ${input_sampleID}_subgenome/${seq_name}.fa");

# 3.4 Extract supporting reads
  my @array = @{$TE_f2{$seq_name}};
  open BIAN2, ">${seq_name}.bian2";
  for(my $i=0;$i<@array;$i+=2){
   print BIAN2 "$array[$i]\n"; 
  }
  if ($sequencing_type eq "paired-end"){
   system ("perl ${directory}Scripts/Extract_fastq.pl -f ${input_sampleID}_ERV_1.1fuq -b ${seq_name}.bian2 -o ${seq_name}_ERV_1.2fuq");
   system ("perl ${directory}Scripts/Extract_fastq.pl -f ${input_sampleID}_ERV_2.1fuq -b ${seq_name}.bian2 -o ${seq_name}_ERV_2.2fuq");
   open FUQ2,"${seq_name}_ERV_1.2fuq";
   while (<FUQ2>){$count_c++;}
   $count_c=($count_c)/4;
  } else {
   $count_c=0;
   system ("touch ${seq_name}_ERV_1.2fuq");
   system ("touch ${seq_name}_ERV_2.2fuq");
  }
  if ($sequencing_type eq "single-end" || $Split) {
   system ("perl ${directory}Scripts/Extract_fastq_sf.pl -f ${input_sampleID}_ERV_1sf.fuq -b ${seq_name}.bian2 -o ${seq_name}_ERV_1sf.2fuq");
   open SF2,"${seq_name}_ERV_1sf.2fuq";
   while (<SF2>){$count_sf++;}
   $count_sf=($count_sf)/4;
  } else {
   $count_sf=0;
   system ("touch ${seq_name}_ERV_1sf.2fuq");
  }
# 3.5 Align to the ref.
  system("${bwa_d}bwa mem -t $threads ${input_sampleID}_subgenome/${seq_name}.fa ${seq_name}_ERV_1.2fuq ${seq_name}_ERV_2.2fuq >${seq_name}_ERV.sam");
  system ("sort ${seq_name}_ERV.sam | uniq > ${seq_name}_ERV.sam2");
  system ("perl ${directory}Scripts/EC_filtered1.pl ${seq_name}_ERV.sam2 >${seq_name}_bp1");
  system ("perl ${directory}Scripts/EC_filtered2.pl ${seq_name}_bp1 >${seq_name}_bp2");
  system ("${bwa_d}bwa mem -t $threads ${input_sampleID}_subgenome/${seq_name}.fa ${seq_name}_ERV_1sf.2fuq >${seq_name}_ERV_sf.sam");
  system ("${samtools_d}samtools view -F 4 ${seq_name}_ERV_sf.sam | uniq > ${seq_name}_ERV_sf.sam2");
  if(-e ${seq_name}."_bp1"){system ("rm ${seq_name}_bp1");}
  if(-e ${seq_name}."_ERV_sf.sam"){system ("rm ${seq_name}_ERV_sf.sam");}
  if(-e ${seq_name}."_ERV_1.2fuq"){system ("rm ${seq_name}_ERV_1.2fuq");}
  if(-e ${seq_name}."_ERV_2.2fuq"){system ("rm ${seq_name}_ERV_2.2fuq");}
  if(-e ${seq_name}."_ERV_1sf.2fuq"){system ("rm ${seq_name}_ERV_1sf.2fuq");}
  if(-e ${seq_name}."_ERV.sam"){system ("rm ${seq_name}_ERV.sam");}

# 3.6 Calculate false chimeric and split reads
  open BP2, "${seq_name}_bp2";
  my @line2=();
  while (my $line2=<BP2>){
   @line2=split /\s+/, $line2;
   if ($line2[0] eq "P") {
    $countf_c{$line2[4]}="P";
   } elsif (!(exists($countf_c{$line2[4]})) && !(exists($countf_c2{$line2[4]})) && $line2[0] ne "P"){
    $countf_c2{$line2[4]}=$line2[0];
   } elsif (!(exists($countf_c{$line2[4]})) && exists($countf_c2{$line2[4]}) && $line2[0] ne "P"){
    if ($countf_c2{$line2[4]} ne $line2[0] && $countf_c2{$line2[4]} ne "P2"){$countf_c2{$line2[4]}="P2";}
   }
  }
  open SF, "${seq_name}_ERV_sf.sam2";
  while ($line2=<SF>){
   @line2=split /\s+/, $line2;
   $countf_sf{$line2[0]}=$line2;
  }
  $line[22]=$count_c;
  $line[26]=$count_sf;
  my @temp=keys %countf_c;
  $line[24]=$#temp+1;                                              ### Total number of false chimeric reads mapped in proper paired;
  @temp=keys %countf_c2;
  my $countf_c2=0;
  my $countf_c3=0;
  for(my $j=0;$j<@temp;$j++){
   if($countf_c2{$temp[$j]} eq "P2"){
    $countf_c2++;
   } else {
    $countf_c3++;
   }
                            }
  $line[25]=$countf_c2;             ### false chimeric with two ends 
  $line[23]=$countf_c3;             ### true chimeric reads;
  @temp=keys %countf_sf;           
  $line[27]=$count_sf-($#temp+1);              ### true split reads
  if($line[22] > 0 && $line[23] eq 0 && $line[26] >0 && $line[27] eq 0){
   next;
  } elsif ($line[22] > 0 && $line[23] eq 0 && $line[26] eq 0) {
   next;
  } elsif ($line[22] > 0 && $line[23] eq 0 && $line[27] > 0) {
   $line[16]=$line[16]."Potential_False";
  } elsif ($line[22] > 0 && $line[23] eq 0 && $line[27] eq 0) {
   $line[16]=$line[16]."Potential_False";
  } elsif ($line[23]<2 && $line[22] eq 2 && $line[26] eq 0) {
   $line[16]=$line[16]."Potential_False";
  } elsif ($line[23]<2 && $line[22] eq 2 && $line[27] > 0) {
   $line[16]=$line[16]."Potential_False";
  } elsif ($line[23]<2 && $line[22] eq 2 && $line[27] eq 0) {
   $line[16]=$line[16]."Potential_False";
  } elsif ($line[23] + $line[27]< 2){
   $line[16]=$line[16]."Potential_False";
  } else {
   $line[16]=$line[16]."Validated";
  }

  print OUTPUT "@line\n";
  if (-e ${seq_name}."_ERV.sam2"){system ("rm ${seq_name}_ERV.sam2");}
  if (-e ${seq_name}."_ERV_sf.sam2"){system ("rm ${seq_name}_ERV_sf.sam2");}
  if (-e ${seq_name}.".bian2"){system ("rm ${seq_name}.bian2");}
  if (-e ${seq_name}."_bp2"){system ("mv ${seq_name}_bp2 ${input_sampleID}_temp/");}
 }
 if (-d ${input_sampleID}."_subgenome"){system ("rm -rf ${input_sampleID}_subgenome");}
}

################################## Step 4. Genotyping
if ($Genotyping){
 print "\nStep 4: Genotyping...\n=====================================\n";
 if (($file_suffix =~ "fq" || $file_suffix =~ "fastq") && !(-e ${input_sampleID}."_h.bam.bai")) {
  print "# Converting SAM to BAM file, and then Sort and index the BAM file......\n";
  system ("${samtools_d}samtools view -b -@ $thread_1 ${input_sampleID}_h.sam >${input_sampleID}_h.bam");
  system ("${samtools_d}samtools sort -@ $thread_1 ${input_sampleID}_h.bam -o ${input_sampleID}_h.sort.bam");
  system ("mv ${input_sampleID}_h.sort.bam ${input_sampleID}_h.bam");
  system ("${samtools_d}samtools index ${input_sampleID}_h.bam");
 }

 if ($Validation) {
  open INPUT, "${input_sampleID}_ERV.output2";
 } else {
  system ("cp ${input_sampleID}_ERV.output ${input_sampleID}_ERV.output2");
  open INPUT, "${input_sampleID}_ERV.output2";
 }           ##### output file
 open OUTPUT,">${input_sampleID}_ERV.output3";
 print OUTPUT "$header3";
 my $line2="";
 my @line=();
 my $Chr;
 my $Position;
 my $R_Position;
 my $Breakpoint;
 my $Number_reads;

# chromosome name
 open GENOME3, "$Human_reference_genome";
 $line2=<GENOME3>;
 if ($line2 =~ ">chr") {
  $chr_name="chr";
 } elsif ($line2 =~ ">Chr"){
  $chr_name="Chr";
 }
 close GENOME3;
 while($line2=<INPUT>){
  @line=split /\s+/, $line2;
  if ($Validation) {
   @line[28..34]=("-") x 7;
  } else{
   @line[22..34]=("-") x 13;
  }
  if ($line[0] =~ "Sample_ID"){next;}
  $Chr=$line[7];
  if ($chr_name ne "NA") {
   $Chr=$chr_name.$Chr;
  }
  $Position=$line[10];
  $R_Position=$line[11];
  $Breakpoint=0;
  $Number_reads=$line[20];
  if ($Position eq "-"){$Position=$R_Position;$Breakpoint=1;}
  if ($R_Position eq "-"){$R_Position=$Position;$Breakpoint=1;}
  if ($Position ne "-" && $R_Position ne "-") {
   $Breakpoint=2;
  }
  my $position1=$Position-int($window_size/2);
  my $position2=$Position+int($window_size/2);
  print "$Chr $position1 $position2 ${input_sampleID} ${file_suffix}\n";
  if ($file_suffix =~ "fq" || $file_suffix =~ "fastq"){
   system ("${samtools_d}samtools view ${input_sampleID}_h.bam ${Chr}:${position1}-${position2} >${input_sampleID}_${Chr}_${Position}_h.sort.sam");
  } else {
   system ("${samtools_d}samtools view ${input_sampleID}${file_suffix} ${Chr}:${position1}-${position2} >${input_sampleID}_${Chr}_${Position}_h.sort.sam");
  }
  system ("sort ${input_sampleID}_${Chr}_${Position}_h.sort.sam | uniq > ${input_sampleID}_${Chr}_${Position}_h.sort.sam2"); 
  system ("perl ${directory}Scripts/EC_filtered1.pl ${input_sampleID}_${Chr}_${Position}_h.sort.sam2 >${input_sampleID}_${Chr}_${Position}_bp1");
  system ("perl ${directory}Scripts/EC_filtered2.pl ${input_sampleID}_${Chr}_${Position}_bp1 >${input_sampleID}_${Chr}_${Position}_bp2");
  system ("Rscript ${directory}Scripts/EC_reads_IS.R ${input_sampleID}_${Chr}_${Position}_bp2");
  system ("perl ${directory}Scripts/EC_filtered3.pl ${input_sampleID}_${Chr}_${Position}_bp2 $Chr $Position >${input_sampleID}_${Chr}_${Position}.bp-F");
  system ("perl ${directory}Scripts/EC_filtered3.pl ${input_sampleID}_${Chr}_${Position}_bp2 $Chr $R_Position >${input_sampleID}_${Chr}_${Position}.bp-R");
  if ($sequencing_type eq "paired-end") {
   system ("perl ${directory}Scripts/EC_genotyping.pl ${input_sampleID}_${Chr}_${Position}.bp-F ${input_sampleID} $Chr $Position $Number_reads $Breakpoint >${input_sampleID}_${Chr}_${Position}.bp2");
   system ("perl ${directory}Scripts/EC_genotyping.pl ${input_sampleID}_${Chr}_${Position}.bp-R ${input_sampleID} $Chr $R_Position $Number_reads $Breakpoint |grep -v ID_Chr_Position >>${input_sampleID}_${Chr}_${Position}.bp2"); 
  } else {
   system ("perl ${directory}Scripts/EC_genotyping_single.pl ${input_sampleID}_${Chr}_${Position}.bp-F ${input_sampleID} $Chr $Position $Number_reads $Breakpoint >${input_sampleID}_${Chr}_${Position}.bp2");
   system ("perl ${directory}Scripts/EC_genotyping_single.pl ${input_sampleID}_${Chr}_${Position}.bp-R ${input_sampleID} $Chr $R_Position $Number_reads $Breakpoint |grep -v ID_Chr_Position >>${input_sampleID}_${Chr}_${Position}.bp2");
  }
  for (my $j=$position1+1000;$j<$position2-1000;$j+=50){
   if(abs($j-$Position)<=$length_insertsize){next;}
   system ("perl ${directory}Scripts/EC_filtered3.pl ${input_sampleID}_${Chr}_${Position}_bp2 $Chr $j >${input_sampleID}_${Chr}_${Position}.bp-B");
   if ($sequencing_type eq "paired-end"){
    system ("perl ${directory}Scripts/EC_genotyping.pl ${input_sampleID}_${Chr}_${Position}.bp-B ${input_sampleID} $Chr $j $Number_reads $Breakpoint |grep -v ID_Chr_Position >>${input_sampleID}_${Chr}_${Position}.bp2-B");
   } else {
    system ("perl ${directory}Scripts/EC_genotyping_single.pl ${input_sampleID}_${Chr}_${Position}.bp-B ${input_sampleID} $Chr $j $Number_reads $Breakpoint |grep -v ID_Chr_Position >>${input_sampleID}_${Chr}_${Position}.bp2-B");
   }
  }
  system ("Rscript ${directory}Scripts/EC_Rscripts_genotyping.R ${input_sampleID}_${Chr}_${Position}.bp2-B ${input_sampleID}_${Chr}_${Position}.bp2");
  if(-e ${input_sampleID}."_".${Chr}."_".${Position}.".bp-F"){system ("rm ${input_sampleID}_${Chr}_${Position}.bp-F");}
  if(-e ${input_sampleID}."_".${Chr}."_".${Position}.".bp-R"){system ("rm ${input_sampleID}_${Chr}_${Position}.bp-R");}
  if(-e ${input_sampleID}."_".${Chr}."_".${Position}.".bp-B"){system ("rm ${input_sampleID}_${Chr}_${Position}.bp-B");}
  if(-e ${input_sampleID}."_".${Chr}."_".${Position}."_bp1"){system ("rm ${input_sampleID}_${Chr}_${Position}_bp1");}
  if(-e ${input_sampleID}."_".${Chr}."_".${Position}."_bp2"){system ("rm ${input_sampleID}_${Chr}_${Position}_bp2");}
  if(-e ${input_sampleID}."_".${Chr}."_".${Position}."_h.sort.sam"){system ("rm ${input_sampleID}_${Chr}_${Position}_h.sort.sam");}
  if(-e ${input_sampleID}."_".${Chr}."_".${Position}."_h.sort.sam2"){system ("rm ${input_sampleID}_${Chr}_${Position}_h.sort.sam2");}
  open BP2, "${input_sampleID}_${Chr}_${Position}.bp2";
  my $line_1=<BP2>;
  $line_1=<BP2>;
  my $line_2=<BP2>;
  my @line_1 = split /\s+/, $line_1;
  my @line_2 = split /\s+/, $line_2;
  if ($line_1[1]<$line_2[1]){
   $line[28]=$line_1[4];
   $line[29]=$line_1[1];            ### No. reads support non ERV
   $line[30]=$line_1[7];            ### Mean
   $line[31]=$line_1[5];            ### SD
   $line[32]=$line_1[6];            ### Count
   $line[33]=$line_1[3];            ### 0.1 quantile;
   if ($line[29]>$line[33]){$line[34]="Heterozygous";}
   else {$line[34]="Homozygous";}         ### genotype
  } else {
   $line[28]=$line_2[4];
   $line[29]=$line_2[1];
   $line[30]=$line_2[7];
   $line[31]=$line_2[5];
   $line[32]=$line_2[6];
   $line[33]=$line_2[3];
   if ($line[29]>$line[33]){$line[34]="Heterozygous";}
   else {$line[34]="Homozygous";}
  }
 print OUTPUT "@line\n";
 if(-e ${input_sampleID}."_".${Chr}."_".${Position}.".bp2"){system ("mv ${input_sampleID}_${Chr}_${Position}.bp2 ${input_sampleID}_temp/");}
 if(-e ${input_sampleID}."_".${Chr}."_".${Position}.".bp2-B"){system ("mv ${input_sampleID}_${Chr}_${Position}.bp2-B ${input_sampleID}_temp/");}
 }
}
 if(-e ${input_sampleID}."_ERV.output3"){system ("mv ${input_sampleID}_ERV.output3 ${input_sampleID}_ERV.output2");}
 if(-e ${input_sampleID}."_ERV_vsu.sort.sam2"){system ("rm ${input_sampleID}_ERV_vsu.sort.sam2");}
 if(-e ${input_sampleID}."_ERV_vsu.sort.sam"){system ("rm ${input_sampleID}_ERV_vsu.sort.sam");}
 if(-e ${input_sampleID}."_ERV.sam"){system ("rm ${input_sampleID}_ERV.sam");}
 if(-e ${input_sampleID}."_ERV_breakpoint"){system ("rm ${input_sampleID}_ERV_breakpoint");}
 if(-e ${input_sampleID}."_ERV_1.1fuq"){system ("rm ${input_sampleID}_ERV_1.1fuq");}
 if(-e ${input_sampleID}."_ERV_2.1fuq"){system ("rm ${input_sampleID}_ERV_2.1fuq");}
 if(-e ${input_sampleID}."_ERV_1sf.fuq"){system ("rm ${input_sampleID}_ERV_1sf.fuq");}
 if(-e ${input_sampleID}."_1.1fuq"){system ("rm ${input_sampleID}_1.1fuq");}
 if(-e ${input_sampleID}."_2.1fuq"){system ("rm ${input_sampleID}_2.1fuq");}
 if(-e ${input_sampleID}."_1sf.fuq"){system ("rm ${input_sampleID}_1sf.fuq");}
 if(-e ${input_sampleID}."_1sf.othu"){system ("rm ${input_sampleID}_1sf.othu");}
 if(-e ${input_sampleID}."_h.sam"){system ("rm ${input_sampleID}_h.sam");}
 if(-e ${input_sampleID}."_ERV.fasta"){system ("rm ${input_sampleID}_ERV.fasta");}
 if(-e ${input_sampleID}."_ERV1_h.sam"){system ("rm ${input_sampleID}_ERV1_h.sam");}
 if(-e ${input_sampleID}."_ERV1_h.sam2"){system ("rm ${input_sampleID}_ERV1_h.sam2");}
 if(-e ${input_sampleID}."_ERV.sort.bam"){system ("rm ${input_sampleID}_ERV.sort.bam");} 

# move temp files to a crated temp folder
 if (-e ${input_sampleID}."_1.1fq"){system ("mv ${input_sampleID}_1.1fq ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_2.1fq"){system ("mv ${input_sampleID}_2.1fq ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_1sf.fastq"){system ("mv ${input_sampleID}_1sf.fastq ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_1sf.others"){system ("mv ${input_sampleID}_1sf.others ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_ERV_1.1fq"){system ("mv ${input_sampleID}_ERV_1.1fq ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_ERV_1.2fq"){system ("mv ${input_sampleID}_ERV_1.2fq ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_ERV1.bian"){system ("mv ${input_sampleID}_ERV1.bian ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_ERV_2.1fq"){system ("mv ${input_sampleID}_ERV_2.1fq ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_ERV_2.2fq"){system ("mv ${input_sampleID}_ERV_2.2fq ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_ERV.3"){system ("mv ${input_sampleID}_ERV.3 ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_ERV.bam"){system ("mv ${input_sampleID}_ERV.bam ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_ERV.error"){system ("mv ${input_sampleID}_ERV.error ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_ERV.fine_mapped"){system ("mv ${input_sampleID}_ERV.fine_mapped ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_ERV.hf"){system ("mv ${input_sampleID}_ERV.hf ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_ERV.soft.bam"){system ("mv ${input_sampleID}_ERV.sort.bam ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_ERV.TE_f"){system ("mv ${input_sampleID}_ERV.TE_f ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_ERV.TE_f2"){system ("mv ${input_sampleID}_ERV.TE_f2 ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_ERV_vsu.sort.bam"){system ("mv ${input_sampleID}_ERV_vsu.sort.bam ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_h1_pe.bam"){system ("mv ${input_sampleID}_h1_pe.bam ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_pe.bam"){system ("mv ${input_sampleID}_pe.bam ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_sm.bam"){system ("mv ${input_sampleID}_sm.bam ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_su.bam"){system ("mv ${input_sampleID}_su.bam ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}.".type"){system ("mv ${input_sampleID}.type ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_vsoft_sort.bam"){system ("mv ${input_sampleID}_vsoft_sort.bam ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_vsu.sort.bam"){system ("mv ${input_sampleID}_vsu.sort.bam ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_h.bam" && $file_suffix ne "_h.bam"){system ("mv ${input_sampleID}_h.bam ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_h.bam.bai" && $file_suffix ne "_h.bam"){system ("mv ${input_sampleID}_h.bam.bai ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_ERV.type"){system ("mv ${input_sampleID}_ERV.type ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}."_u.bam"){system ("mv ${input_sampleID}_u.bam ${input_sampleID}_temp/");}
 if (-e ${input_sampleID}.".1fq"){system ("mv ${input_sampleID}.1fq ${input_sampleID}_temp/");}

############################ subs functions
# Fine_mapped file
sub out1{
 my ($input1)=@_;
 my %output=();
 my $output_header="";
 my @output_line=();

 open FINE, "${input1}.fine_mapped";
 while (<FINE>){
  @output_line=split;
  $output_header=join ("_",(@output_line[0..3],$output_line[5]));
  @{$output{$output_header}} = ("-") x 22;
  $output{$output_header}[0]=$output_line[0];
  $output{$output_header}[0]=~s/_ERV//;
  if($Split){
   $output{$output_header}[1] = "Yes";
  } else{
   $output{$output_header}[1] = "No";
  }
  if($Validation){
   $output{$output_header}[2] = "Yes";
  } else {
   $output{$output_header}[2] = "No";
  }
  if($Genotyping){
   $output{$output_header}[3] = "Yes";
  } else {
   $output{$output_header}[3] = "No";
  }
  ${$output{$output_header}}[4]=$Human_reference_genome;
  ${$output{$output_header}}[5]=$TE_reference_genomes;
  ${$output{$output_header}}[6]=$output_line[5];
  @{$output{$output_header}}[7..9]=@output_line[1..3];
  if(($output_line[8] ne "-" && $output_line[17] eq "-") || ($output_line[8] ne "-" && $output_line[17] ne "-" && ($output_line[22]+$output_line[23]) <= ($output_line[13] + $output_line[14]))){
   if ($output_line[7] eq "-") {
    $output{$output_header}[14]="E(++);";
    $output{$output_header}[10]=$output_line[9];
    $output{$output_header}[13]=$output_line[12];
   } else {
    $output{$output_header}[14]="D(++);";
    $output{$output_header}[10]=$output_line[7];
    $output{$output_header}[13]=$output_line[10];
   }   
  } elsif (($output_line[8] eq "-" && $output_line[17] ne "-") || ($output_line[8] ne "-" && $output_line[17] ne "-" && ($output_line[22]+$output_line[23]) >= ($output_line[13] + $output_line[14]))){
   if ($output_line[16] eq "-") {
    $output{$output_header}[14]="E(+-);";
    $output{$output_header}[10]=$output_line[18];
    $output{$output_header}[12]=$output_line[20];
   } else {
    $output{$output_header}[14]="D(+-);";
    $output{$output_header}[10]=$output_line[16];
    $output{$output_header}[12]=$output_line[19];
   }
  } elsif ($output_line[8] eq "-" && $output_line[17] eq "-") {
    $output{$output_header}[14]="na;";
  }
  
  if(($output_line[26] ne "-" && $output_line[35] eq "-") || ($output_line[26] ne "-" && $output_line[35] ne "-" && ($output_line[40]+$output_line[41]) >= ($output_line[31] + $output_line[32]))){
   if ($output_line[25] eq "-") {
    $output{$output_header}[14]=$output{$output_header}[14]."e(-+)";
    $output{$output_header}[11]=$output_line[26];
    $output{$output_header}[13]=$output_line[30];
   } else {
    $output{$output_header}[14]=$output{$output_header}[14]."d(-+)";
    $output{$output_header}[11]=$output_line[25];
    $output{$output_header}[13]=$output_line[28];
   }
  } elsif (($output_line[26] eq "-" && $output_line[35] ne "-") || ($output_line[26] ne "-" && $output_line[35] ne "-" && ($output_line[40]+$output_line[41]) <= ($output_line[31] + $output_line[32]))){
   if ($output_line[34] eq "-") {
    $output{$output_header}[14]=$output{$output_header}[14]."e(--)";
    $output{$output_header}[11]=$output_line[35];
    $output{$output_header}[12]=$output_line[38];
   } else {
    $output{$output_header}[14]=$output{$output_header}[14]."d(--)";
    $output{$output_header}[11]=$output_line[34];
    $output{$output_header}[12]=$output_line[37];
   }
  } elsif ($output_line[26] eq "-" && $output_line[35] eq "-") {
   $output{$output_header}[14]=$output{$output_header}[14]."na";
  }

  if ($output{$output_header}[10] ne "-" && $output{$output_header}[11] ne "-"){
   if (($output{$output_header}[14] =~ "D" && $output{$output_header}[14] =~ "d") || ($output{$output_header}[14] =~ "E" && $output{$output_header}[14] =~ "e")){
    $output{$output_header}[15]=int(($output{$output_header}[11]+$output{$output_header}[10])/2);
   } elsif ($output{$output_header}[14] =~ "D"){
    $output{$output_header}[15]=$output{$output_header}[10];
   } elsif ($output{$output_header}[14] =~ "d") {
    $output{$output_header}[15]=$output{$output_header}[11];
   }
  } elsif ($output{$output_header}[10] ne "-" && $output{$output_header}[11] eq "-") {
   $output{$output_header}[15]=$output{$output_header}[10];
  } elsif ($output{$output_header}[11] ne "-" && $output{$output_header}[10] eq "-") {
   $output{$output_header}[15]=$output{$output_header}[11];
  }
  $output{$output_header}[14] =~ s/e/E/;
  $output{$output_header}[14] =~ s/d/D/;
 }
 close FINE;

# open SM.bam and ERV1_1.bian files
 my %human_AS=();
 my $file_name=${input1};
 $file_name=~s/_ERV$//;
 if(-e ${input1}.".type"){
   open SM_bam,"${input1}.type";
 } elsif (-e ${file_name}.".type") {
   open SM_bam,"${file_name}.type";
 } else {
   system ("touch ${input1}.type");
   open SM_bam,"${input1}.type";
 }

 while(<SM_bam>){
  my @line=split;
  @{$human_AS{$line[0]}}=@line;
 }
 close SM_bam;

 my @line=();
 if(-e ${input1}."1.bian"){
   open BIAN,"${input1}1.bian"; 
 while(<BIAN>){
  @line=split;
  @{$human_AS{$line[0]}}=@line;
 }
 close BIAN;
                          }
# read TE_f2 file
 my @f2_line=();
 my $f2_header="";
 my %sub_human_AS=();
 open F2,"${input1}.TE_f2";
 while (<F2>){
  %sub_human_AS=();
  @f2_line=split;
  $f2_header=join ("_",($f2_line[15],@f2_line[31..33],$f2_line[35]));
  if(exists($output{$f2_header})){
   $output{$f2_header}[20]=$f2_line[39];
   $output{$f2_header}[21]=$f2_line[41];
# calculate human alignment score;
   my @list_reads=split /\|/, $f2_line[30];
   for(my $i=0;$i<@list_reads;$i+=2){
    if(exists($human_AS{$list_reads[$i]})){
      $sub_human_AS{$list_reads[$i]}=$human_AS{$list_reads[$i]};
                                          }
                                    }
   my @list_reads2=keys %sub_human_AS;
   my @calculate=(0) x 4;
   for(my $i=0;$i<@list_reads2;$i++){
    $calculate[0]++;
    $calculate[1]=$calculate[1]+${$sub_human_AS{$list_reads2[$i]}}[2];
    $calculate[2]=$calculate[2]+${$sub_human_AS{$list_reads2[$i]}}[3];
    if($calculate[3]<(${$sub_human_AS{$list_reads2[$i]}}[2]-${$sub_human_AS{$list_reads2[$i]}}[3])){
     $calculate[3]=${$sub_human_AS{$list_reads2[$i]}}[2]-${$sub_human_AS{$list_reads2[$i]}}[3];
                                                                                                   }
                                    }
   if($calculate[0]>0){
    $output{$f2_header}[17]=$calculate[1]/$calculate[0];
    $output{$f2_header}[18]=$calculate[2]/$calculate[0];
                      }
   $output{$f2_header}[19]=$calculate[3];
   if ($output{$f2_header}[17] eq "-"){
    $output{$f2_header}[17]=0;
   }
   if ($output{$f2_header}[18] eq "-"){
    $output{$f2_header}[18]=0;
   }
   if ($output{$f2_header}[17]-$output{$f2_header}[18]>=30){$output{$f2_header}[16]="Non-repeat;";}
   elsif ($output{$f2_header}[17]-$output{$f2_header}[18]>15 && $output{$f2_header}[17]-$output{$f2_header}[18]<=30){$output{$f2_header}[16]="Likely-repeat_or_false_positive;";}
   elsif ($output{$f2_header}[17]-$output{$f2_header}[18]>5 && $output{$f2_header}[17]-$output{$f2_header}[18]<=15){$output{$f2_header}[16]="Very_likely-repeat_or_false_positive;";}
   elsif ($output{$f2_header}[17]-$output{$f2_header}[18]<=5){$output{$f2_header}[16]="Repeat_or_false_positive;";}
   else{$output{$f2_header}[16]="NA;";}
                               }
            }
 
# output
 my @list=sort{$output{$b}->[19]<=>$output{$a}->[19]} keys %output;
 open OUT, ">${input1}.output";
 print OUT "$header1";
 $" = "\t";
 for(my $i=0;$i<@list;$i++){
  print OUT "@{$output{$list[$i]}}\n";
 }
 close OUT;
} #write to output file end

# Create type file
sub create_type { 
 system ("${samtools_d}samtools view ${input_sampleID}_sm.bam >${input_sampleID}_sm.sam"); 
 open SM,"${input_sampleID}_sm.sam";
 open TYPE,">${input_sampleID}.type";
 my $sm_1="";
 my $as="";
 my $xs="";
 while($sm_1=<SM>){
  my @sm_1=split /\s+/, $sm_1;
# AS
  $as="NA";$xs="NA";
  unless($sm_1[2] eq "*"){
   for(my $i=11;$i<@sm_1;$i++){
    if($sm_1[$i]=~s/AS:i://){$as=$sm_1[$i];}
    if($sm_1[$i]=~s/XS:i://){$xs=$sm_1[$i];}
                              }
                         }
  if($sm_1[1]%256>=128){print TYPE "$sm_1[0] L $as $xs $sm_1[5]\n";}
  else{print TYPE "$sm_1[0] R $as $xs $sm_1[5]\n";}
 }
 close TYPE;
 close SM;
}

# Align to the human reference genome
sub align_to_hg {
 my ($input1,$suffix1)=@_;
 if ($data_type eq "RNA-seq"){
  my $suffix2=$suffix1;
  $suffix2=~s/.gz$//;
  if($sequencing_type eq "paired-end"){
    if($suffix1 =~ "gz"){
     system ("gunzip -c ${input1}_1${suffix1} >${input1}_1${suffix2}");
     system ("gunzip -c ${input1}_2${suffix1} >${input1}_2${suffix2}");
    }
    system ("${tophat_d}tophat2 -p $threads -o ${input1} $human_genome_tophat ${input1}_1${suffix2} ${input1}_2${suffix2}");
  } else {
    if($suffix1 =~ "gz"){
     system ("gunzip -c ${input1}${suffix1} >${input1}${suffix2}");
    }
    system ("${tophat_d}tophat2 -p $threads -o ${input1} $human_genome_tophat ${input1}${suffix2}");
   }
  system ("${samtools_d}samtools merge -f ${input1}_h.bam ${input1}/accepted_hits.bam ${input1}/unmapped.bam");
  system ("${samtools_d}samtools view ${input1}_h.bam >${input1}_h.sam"); 
 } else {
  if($sequencing_type eq "paired-end"){
   system("${bwa_d}bwa mem -t $threads $human_genome ${input1}_1${suffix1} ${input1}_2${suffix1} >${input1}_h.sam");
  } elsif ($sequencing_type eq "single-end") {
   system("${bwa_d}bwa mem -t $threads $human_genome ${input1}${suffix1} >${input1}_h.sam");
  }
 }
}

# Convert bam file to fastq file
sub convert_bamtofastq {
 my ($input1)=@_;
 my $suffix1="";
 if(($file_suffix =~ "bam" || $file_suffix =~ "sam" ) && $order ne 2){
  $suffix1=$file_suffix;
 } else {
  $suffix1="_h.sam";
 }
 if($sequencing_type eq "paired-end"){
  system ("${samtools_d}samtools view -b -f 4 -F 264 -@ $thread_1 ${input1}$suffix1 > ${input1}_su.bam");
  system ("${samtools_d}samtools view -b -f 8 -F 260 -@ $thread_1 ${input1}$suffix1 > ${input1}_sm.bam");
  if ($Split){
   system ("${samtools_d}samtools view -b -f 12 -F 256 -@ $thread_1 ${input1}$suffix1 > ${input1}_pe.bam");
   system ("${samtools_d}samtools view -b -f 2 -@ $thread_1 ${input1}$suffix1 >${input1}_mpe.bam");
   system ("${SE_MEI_d}extractSoftclipped -l 1 ${input1}_mpe.bam >${input1}_soft.fastq.gz");            ##### minumum length of split reads is set as 1 bp;
   if(-e ${input1}."_mpe.bam"){system ("rm ${input1}_mpe.bam");}
   system ("${samtools_d}samtools merge -f ${input1}_h1.bam ${input1}_su.bam ${input1}_sm.bam ${input1}_pe.bam");
  } else{
    system ("${samtools_d}samtools merge -f ${input1}_h1.bam ${input1}_su.bam ${input1}_sm.bam");
  }
  system ("${samtools_d}samtools sort -n -@ $thread_1 ${input1}_h1.bam -o ${input1}_h1$suffix1");
  system ("${samtools_d}samtools fastq -@ $thread_1 -N ${input1}_h1$suffix1 -1 ${input1}_h1_1.1fq -2 ${input1}_h1_2.1fq");
  system ("perl ${directory}Scripts/Check_paired_end.pl -s ${input1}_h1 -f .1fq");
  if(-e ${input1}."_h1_1.1fq2"){system ("mv ${input1}_h1_1.1fq2 ${input1}_h1_1.1fq");}
  if(-e ${input1}."_h1_2.1fq2"){system ("mv ${input1}_h1_2.1fq2 ${input1}_h1_2.1fq");}
  if(-e ${input1}."_h1.bam"){system ("rm ${input1}_h1.bam");}
 } elsif ($sequencing_type eq "single-end"){
  if(-e ${input1}.$suffix1){
   system("${samtools_d}samtools view -b -f 4 -@ $thread_1 ${input1}$suffix1 >${input1}_u.bam");
   system("${samtools_d}samtools view -b -f 4 -@ $thread_1 ${input1}$suffix1 >${input1}_u.bam");       
   system("${samtools_d}samtools view -b -@ $thread_1 ${input1}$suffix1 >${input1}_h.bam");
   system("${SE_MEI_d}extractSoftclipped -l 1 ${input1}$suffix1 >${input1}_soft.fastq.gz");              
   system("${hydra_d}bin/bamToFastq -bam ${input1}_u.bam -fq1 ${input1}_h2.1fq -fq2 ${input1}_h1.1fq");
   if(-e ${input1}."_h2.1fq"){system("rm ${input1}_h2.1fq");}
  }
 }
}

print "                 ########  #######   ##       ##                        ##     ##                \n";
print "                 ########  ##    ##  ##       ##                        ##     ##                \n";
print "                 ##        ##    ##   ##     ##                         ##     ##                \n";
print "                 ##        ##   ##    ##     ##                         ##     ##                \n";
print "                 ########  ######      ##   ##      #####     ####      ##     ##                \n";
print "                 ########  #####       ##   ##     ##   ##   ##  ##     ##     ##                \n";
print "                 ##        ## ##        ## ##     ##        ##    ##    ##     ##                \n";
print "                 ##        ##  ##       ## ##     ##        ##    ##    ##     ##                \n";
print "                 ########  ##   ##       ###       ##   ##   ##   ##    ## ##  ## ##             \n";
print "                 ########  ##    ##      ###        #####     #### ##    ##     ##               \n";
print "\n\n Done!\n";
