#!/usr/bin/env perl
#
# Author: 	Xun Chen
# Email: 	xunchen85@gmail.com
# Date:		02/24/2022
# Version:	v1.4
#
# Updates (v1.4):
#	01/04/2024:		Corrected two bugs of the combining VCF function in the Combine_VCF_files.pl, which led to some missing TE loci in the combined VCF file (i.e., the last or a single TE locus per chromosome)
#       02/24/2022:             Added the Filtering_VCF_v1.4.pl script for filtering out low-quality candidate TE insertions
#	04/23/2019:		Corrected a bug of the genotyping function with the input of a list of BAM files
#	03/12/2019:		Corrected a bug caused by the sample IDs containing the "AS" character
#	02/15/2019:		Re-designed the genotyping process
#	02/10/2019:		Added the scripts to distinguish missing genotypes and none TE insertions genotypes for all samples in the combined VCF file
# 	02/06/2019:		Corrected the output coordinates of TE insertions with TSD
#	02/02/2019:		Further standardized the VCF format for the usage of bcftools
#	02/01/2019:		Add Phred-scale genotype quality and likelihoods
#	01/29/2019:		Adjusted reciprocal-aligned reference genomic region length using the estimated insert size and SD, which significantly reduce the false-positives
#	01/29/2019:		Added a function to estimate insert size and SD
#	01/38/2019:		Corrected several bugs in the main Perl script
#	01/24/2019:		Corrected a bug in the script combing VCF files from multiple samples
#
# Updates (v1.3):
#       11/20/2018:             Added the scripts to merge various samples into a list of known TE loci or TE loci detected from the analyzed samples
# 	11/12/2018:		Updated the Output in VCF_v4.2 format
# 	11/05/2018:		Debugged the support of the BAM files generated using Bowtie2
#
# Updates (v1.2):
#	11/01/2018:		Further optimized the speed of validation steps
#	10/21/2018:		Support the input of multiple bam file
#	10/10/2018:		Optimized the validation steps to increase the specificity
#
# Updates (v1.1):
#	09/02/2018:	    	Optimize the validation steps to significantly increase the speed
#	08/28/2018:         	Updated the parameter of -S to specify the length of split reads used (20 bp by default; >=40 bp is recommend for reads in length of 150 bp)
#	08/10/2018:         	Added component to support BAM files using different chromosome IDs as the reference genome, such as "Chr1", "chr1", "1", and "NC_000001.11"
#	07/17/2018:         	Corrected bugs for checking input files;
#	07/17/2018:         	Corrected the errors for detecting and genotyping TE insertions using single-end sequencing data;
#	07/16/2018:         	Re-formatted the output files
#	07/15/2018:         	Released ERVcaller Version 1.1 and software manual
#
# Release (v1.0):
#	05/27/2018:         	Released ERVcaller Version 1.0  (a testing version) and software manual
#
# Questions or help:
#
# For questions or support contact Xun Chen Ph.D. (xunchen85@gmail.com or Xun.Chen@uvm.edu);

use Getopt::Long qw(:config no_ignore_case);
use strict;
use warnings;
use Cwd qw();
my $tophat_d="";
my $bwa_d="";
my $samtools_d="";
my $SE_MEI_d="";

### Except R that you need to load in the enviroment, for other tools, you can either load in the enviroment (export), or specify the path below (remove the "#") e.g.:
# $bwa_d="/ERVcaller/Tools/bwa-master/";
# $samtools_d="/ERVcaller/Tools/samtools-1.6/";
# $tophat_d="/ERVcaller/Tools/tophat-2.1.1.Linux_x86_64/";
# $SE_MEI_d="/ERVcaller/Tools/SE-MEI/";
 
# Parameter variables;
my $help;
my $input_sampleID;
my $file_suffix;
my $threads;
my $data_type;
my $sequencing_type;
my $Human_reference_genome;
my $TE_reference_genomes;
my $length_insertsize;
my $std_insertsize;
my $L_std_insertsize;
my $read_len;
my $Split;
my $Genotype;
my $directory=$0;
my $number_of_reads;
my $Input_directory;
my $Output_directory;
my $BWA_MEM;
my $multiple_BAM;

# Variables
GetOptions(
           "i|input_sampleID=s" => \$input_sampleID,
           "f|file_suffix=s" => \$file_suffix,
	   "H|Human_reference_genome=s" => \$Human_reference_genome,
	   "T|TE_reference_genomes=s" =>\$TE_reference_genomes,
           "I|Input_directory=s" =>\$Input_directory,
           "O|Output_directory=s" =>\$Output_directory,
	   "n|number_of_reads=i" =>\$number_of_reads,
           "d|data_type=s" => \$data_type,
           "s|sequencing_type=s" => \$sequencing_type,
           "l|length_insertsize=f" =>\$length_insertsize,
           "L|L_std_insertsize=f" =>\$L_std_insertsize,
	   "r|read_len=i" =>\$read_len,
           "t|threads=i" => \$threads,
           "S|Split=i" => \$Split,
	   "m|multiple_BAM" => \$multiple_BAM,
           "B|BWA_MEM" =>\$BWA_MEM,
	   "G|Genotype" => \$Genotype,
           "h|help" => \$help
          );

# default value
if (!defined ($file_suffix)) {
  $file_suffix = ".bam";
}

if (!defined ($read_len)) {
  $read_len = 100;
}

if (!defined ($Split)) {
  $Split = 20;
}

if (defined ($length_insertsize) && defined ($L_std_insertsize)) {
  $std_insertsize = $L_std_insertsize;
} elsif (defined ($length_insertsize) && !(defined ($L_std_insertsize))) {
#  $std_insertsize = $length_insertsize *0.2;
}

if (!defined ($threads)) {
  $threads = 1;
}
if (!defined ($Input_directory)) {
  $Input_directory=Cwd::cwd();
  $Input_directory=$Input_directory."/";
}

if (!defined($Output_directory)) {
  $Output_directory=Cwd::cwd();
  $Output_directory=$Output_directory."/";
}

if (!defined ($number_of_reads)) {
  $number_of_reads = 3;
}

if (!defined($data_type)){
  $data_type = "WGS";
}

if (!defined ($sequencing_type)) {
  $sequencing_type = "paired-end";
}

if ($sequencing_type eq "single-end") {
  $length_insertsize = 500;     ### only for extracting mapped regions;
}

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

print "\n\n# ERVcaller\n";

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

# extract the path of the installed ERVcaller software
my @directory=split /\//, $directory;
if($#directory>0){
  $directory=join("/",@directory[0..$#directory-1]);
  $directory=$directory."/";
} else {
  $directory="";
}
 
 my $TSD_min_len=100;
 my $Alignment_score= 30;
 my $human_genome=$Human_reference_genome;
 my $human_genome_tophat=$Human_reference_genome;
 my %genome=();
 my $order=0;
 my $min_insertsize = 0;
 my $thread_1=$threads-1;
 my $double_length_insertsize;
 my $cmd="";
 my @bp1_tmp3 = ();
 my $header3="Sample_ID Is_Split_mode Is_genotyped Is_validated Human_ref. TE_reference TE_sequence_name Chr. Start End Upstream_breakpoint_on_human Downstream_breakpoint_on_human Upstream_breakpoint_on_TE Downstream_breakpoint_on_TE Information_both_breakpoints Insertion_site Group Average_AS_for_chimeric_and_improper_reads_on_human Average_XS_for_chimeric_and_improper_reads_on_human Maximum_AS_for_chimeric_and_improper_reads_on_human No._supporting_reads Average_AS_supporting_reads No._chimeric_and_improper_reads V_True_chimeric_and_improper_reads V_False_chimeric_and_improper_reads_PE V_False_chimeric_and_improper_reads_others No._split_reads V_True_split_reads No._split_reads_(<20bp) Geno_No._reads_supporting_non_TE_insertions Geno_Read_depth_of_the_genomic_window_(Mean) Geno_Read_depth_of_the_genomic_window_(SD) Geno_No._random_locations_of_the_genomic_window Geno_Read_depth_of_the_genomic_window_(Quantile=0.2) Genotype\n";

########################################
################################# Step 2
 unless (-d $Output_directory) {
   system ("mkdir $Output_directory");
 }
 chdir $Output_directory; 
 ##### 2.1 Check input file 
 unless (-d ${input_sampleID}."_temp") {
  system ("mkdir ${input_sampleID}_temp");
 }
 print "\nStep 2: Detecting TE insertions...\n=====================================\n";
 if(defined($input_sampleID) && -e ${Input_directory}.${input_sampleID}."_1".${file_suffix} && -e ${Input_directory}.${input_sampleID}."_2".${file_suffix} && $sequencing_type eq "paired-end" && ($file_suffix =~ "fq" || $file_suffix =~ "fastq")){
   print "~~~~~ paired-end reads in fastq format were loaded\n";
 } elsif (defined($input_sampleID) && -e ${Input_directory}.${input_sampleID}.${file_suffix} && $sequencing_type eq "single-end" && ($file_suffix =~ "fq" || $file_suffix =~ "fastq")){
   print "~~~~~ single-end read in fastq format was loaded\n";
 } elsif(defined($input_sampleID) && -e ${Input_directory}.${input_sampleID}.${file_suffix} && $sequencing_type eq "paired-end" && ($file_suffix =~ "bam" || $file_suffix =~ "sam")){
   print "~~~~~ paired-end reads in bam format were loaded\n";
   if ($Genotype){
     if (-e ${Input_directory}.${input_sampleID}.${file_suffix}.".bai"){
       print "~~~~~ the input bam file was indexed\n";
     } else {
       print "~~~~~ the input bam file was not indexed, please index the bam file using samtools for performing the validation or genotyping function\n";
       exit;
     }
   }
 } elsif (defined($input_sampleID) && -e ${Input_directory}.${input_sampleID}.${file_suffix} && $sequencing_type eq "single-end" && ($file_suffix =~ "bam" || $file_suffix =~ "sam")){
   print "~~~~~ single-end reads in bam format were loaded\n";
   if ($Genotype){
     if (-e ${Input_directory}.${input_sampleID}.${file_suffix}.".bai"){
       print "~~~~~ the input bam file was indexed\n";
     } else {
       print "~~~~~ the input bam file was not indexed, please index the bam file using samtools for performing the validation or genotyping function\n";
       exit;
     }
   }
 } elsif (defined($input_sampleID) && -e ${Input_directory}.${input_sampleID}.${file_suffix} && $multiple_BAM) {
   print "~~~~~ a list of multiple BAM files were loaded\n";   
 } else {
   prtErr("# Error: cound not find the input data under the provided sampleID\n");
   print "Input: ${Input_directory}${input_sampleID}${file_suffix}\n";
   prtUsa();
   exit;
 }

 ##### 2.2 Extract supporting reads
 if(${file_suffix} =~ "bam" || defined($multiple_BAM)){
  $order=1;
  convert_bamtofastq(${input_sampleID});
  unless ($BWA_MEM) {
    align_to_hg(${input_sampleID}."_h1",".1fq");
    $order=2;
    convert_bamtofastq(${input_sampleID}."_h1");
    if (-e ${input_sampleID}."_h1_sm.bam") { system("mv ${input_sampleID}_h1_sm.bam ${input_sampleID}_sm.bam");}
    if (-e ${input_sampleID}."_h1_su.bam") { system("mv ${input_sampleID}_h1_su.bam ${input_sampleID}_su.bam");}
  }
  if($Split || $sequencing_type eq "single-end"){
    system("gunzip -c ${input_sampleID}_soft.fastq.gz >${input_sampleID}_1sf.fastq");
    unless ($BWA_MEM) {
      system("gunzip -c ${input_sampleID}_h1_soft.fastq.gz >>${input_sampleID}_1sf.fastq");
    }
  }
  if($sequencing_type eq "paired-end"){
    unless ($BWA_MEM) {
      system("mv ${input_sampleID}_h1_h1_1.1fq ${input_sampleID}_1.1fq");
      system("mv ${input_sampleID}_h1_h1_2.1fq ${input_sampleID}_2.1fq");
    } else {
      system("mv ${input_sampleID}_h1_1.1fq ${input_sampleID}_1.1fq");
      system("mv ${input_sampleID}_h1_2.1fq ${input_sampleID}_2.1fq");
    }
  } else {
    unless ($BWA_MEM) {
      system("mv ${input_sampleID}_h1_h1.1fq ${input_sampleID}.1fq");
    } else {
      system("mv ${input_sampleID}_h1.1fq ${input_sampleID}.1fq");
    }
  }
 } else {
   $order=0;
   align_to_hg(${input_sampleID},${file_suffix});
   convert_bamtofastq(${input_sampleID});
   if($Split || $sequencing_type eq "single-end"){
     system("gunzip -c ${input_sampleID}_soft.fastq.gz >${input_sampleID}_1sf.fastq");
   }
   if($sequencing_type eq "paired-end"){
     system("mv ${input_sampleID}_h1_1.1fq ${input_sampleID}_1.1fq");
     system("mv ${input_sampleID}_h1_2.1fq ${input_sampleID}_2.1fq");
   }
 }

 ##### Filter split reads
 open SF1, "${input_sampleID}_1sf.fastq";
 open SF2, ">${input_sampleID}_1sf.fastq2";
 while (my $tmp1=<SF1>) {
   if ($tmp1=~"^\@soft") {
     my @tmp1=split /\|/, $tmp1;
     if ($tmp1[2]%4 >=2 || $sequencing_type eq "single-end") {
       print SF2 "$tmp1";
       $tmp1=<SF1>;
       print SF2 "$tmp1";
       $tmp1=<SF1>;
       print SF2 "$tmp1";
       $tmp1=<SF1>;
       print SF2 "$tmp1";
     } else {
       $tmp1=<SF1>;
       $tmp1=<SF1>;
       $tmp1=<SF1>;
     }
   }
 }
 close SF1;
 close SF2;
 system ("mv ${input_sampleID}_1sf.fastq2 ${input_sampleID}_1sf.fastq");

 ##### 2.3 Chimeric reads amd Split reads
 print "\nChimeric and split reads...\n=====================================\n";
 if($sequencing_type eq "paired-end"){
   system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 30 -h 10000 -a -Y -M $TE_reference_genomes ${input_sampleID}_1.1fq ${input_sampleID}_2.1fq >${input_sampleID}_vsu.sam");
 } else {        
   system ("touch ${input_sampleID}_vsu.sam");
 }
 system ("touch ${input_sampleID}_all_breakpoint");

 if($sequencing_type eq "single-end" || ($sequencing_type eq "paired-end" && $Split)){
   system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $TE_reference_genomes ${input_sampleID}_1sf.fastq >${input_sampleID}_vsoft.sam");
   system ("perl ${directory}Scripts/Soft_clipping_transfer.pl -f ${input_sampleID}_vsoft.sam -o ${input_sampleID}_vsoft_breakpoint");
   system ("cat ${input_sampleID}_vsoft_breakpoint >>${input_sampleID}_all_breakpoint");
   system ("rm ${input_sampleID}_vsoft_breakpoint");
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
     if ($as eq "NA") {
       next;
     } elsif($sm_1[1]%256>=128 && ($as >=$Alignment_score && $as >= 2*$xs)){                        ############ Direction of TE reads; ###############
       print TYPE "$sm_1[0] L $as $xs $sm_1[5]\n";
     } elsif($as >=$Alignment_score && $as >= 2*$xs){
       print TYPE "$sm_1[0] R $as $xs $sm_1[5]\n";
     }
   }
   close TYPE;
   close SM;
   system ("perl ${directory}Scripts/Break_point_calling.pl -type ${input_sampleID}.type -position ${input_sampleID}_sm.sam -TE ${input_sampleID}_vsu.sam -alignment_score $Alignment_score -o ${input_sampleID}");
   system ("cat ${input_sampleID}_breakpoint >>${input_sampleID}_all_breakpoint");
   system ("rm ${input_sampleID}_breakpoint");
 }
 system ("perl ${directory}Scripts/Filtered_fastq.pl $input_sampleID");

 ##### 2.4 Improper reads
 print "\nImproper reads...\n=====================================\n";
 if($sequencing_type eq "paired-end"){
   if ($file_suffix =~ "bam" || $file_suffix =~ "sam" || defined ($multiple_BAM)){
     if ($multiple_BAM) {
       open LIST, "${Input_directory}${input_sampleID}${file_suffix}";
       while (my $list_tmp1=<LIST>){
	 my @list_tmp1=split /\s+/, $list_tmp1;
	 if (-e $list_tmp1[0]) {
	   unless (-e ${input_sampleID}."_m_ERV.bam") {
	     system ("${samtools_d}samtools view -F 14 -b -@ $thread_1 $list_tmp1[0] >${input_sampleID}_m_ERV.bam");
             if (!defined($length_insertsize) || !defined($L_std_insertsize)) {
               ($length_insertsize,$std_insertsize) = estimate_insertsize ($list_tmp1[0],${input_sampleID},0.05);
             }
           } else {
	     system ("${samtools_d}samtools view -F 14 -b -@ $thread_1 $list_tmp1[0] >${input_sampleID}_m1_ERV.bam");
	     system ("${samtools_d}samtools merge -f ${input_sampleID}_m2.bam ${input_sampleID}_m_ERV.bam ${input_sampleID}_m1_ERV.bam");
	     system ("mv ${input_sampleID}_m2.bam ${input_sampleID}_m.bam");
	     system ("rm ${input_sampleID}_m1_ERV.bam");
             if (!defined($length_insertsize) || !defined($L_std_insertsize)) {
               my ($insertsize_estimate_tmp1,$std_estimate_tmp1) = estimate_insertsize ($list_tmp1[0],${input_sampleID},0.05);
               if ($insertsize_estimate_tmp1 >$length_insertsize) {
	         $length_insertsize=$insertsize_estimate_tmp1;
               } 
               if ($std_estimate_tmp1 >$std_insertsize) {
                 $std_insertsize=$std_estimate_tmp1;
               }
             }
	   }
         } else {
	   print "\nCannot read the BAM file in the list of multiple BAM files...\n";
	   last;
	 }
       }
       close LIST;
       system ("mv ${input_sampleID}_m.bam ${input_sampleID}_h_ERV.bam");
     } else {
       system ("${samtools_d}samtools view -F 14 -b -@ $thread_1 ${Input_directory}${input_sampleID}$file_suffix >${input_sampleID}_h_ERV.bam");
       if (!defined($length_insertsize) || !defined($L_std_insertsize)) {
         ($length_insertsize,$std_insertsize) = estimate_insertsize (${Input_directory}.${input_sampleID}.$file_suffix,${input_sampleID},0.05);
       }
     }
     unless ($BWA_MEM) {
       system ("${samtools_d}samtools view -F 14 -b -@ $thread_1 ${input_sampleID}_h1_h.sam >${input_sampleID}_h1_ERV.bam");
       system ("${samtools_d}samtools merge -f ${input_sampleID}_ERV.bam ${input_sampleID}_h_ERV.bam ${input_sampleID}_h1_ERV.bam");
     } else {
       system ("mv ${input_sampleID}_h_ERV.bam ${input_sampleID}_ERV.bam");
     }
   } else {
     system ("${samtools_d}samtools view -F 14 -b -@ $thread_1 ${input_sampleID}_h.sam >${input_sampleID}_ERV.bam");
     if (!defined($length_insertsize) || !defined($L_std_insertsize)) {
       ($length_insertsize,$std_insertsize) = estimate_insertsize (${input_sampleID}."_h.sam",${input_sampleID},0.05);
     }
   }
   if ($sequencing_type eq "paired-end" && !($length_insertsize)) {
     print "# Error: Could not get the insert size\n";
     exit;
   } 
   system ("${samtools_d}samtools fastq -@ $thread_1 -N ${input_sampleID}_ERV.bam -1 ${input_sampleID}_ERV1_1.1fq -2 ${input_sampleID}_ERV1_2.1fq");
   system ("perl ${directory}Scripts/Check_paired_end.pl -s ${input_sampleID}_ERV1 -f .1fq");
   system ("mv ${input_sampleID}_ERV1_1.1fq2 ${input_sampleID}_ERV1_1.1fq");
   system ("mv ${input_sampleID}_ERV1_2.1fq2 ${input_sampleID}_ERV1_2.1fq");
   $order=2;
   align_to_hg(${input_sampleID}."_ERV1",".1fq");
   unless ($BWA_MEM) {
     system ("${samtools_d}samtools view -F 14 -@ $thread_1 ${input_sampleID}_ERV1_h.sam >${input_sampleID}_ERV2_h.sam");          ##### added on 2/6/2019, adjust the false from different aligners
     system ("mv ${input_sampleID}_ERV2_h.sam ${input_sampleID}_ERV1_h.sam");                                                      ##### added on 2/6/2019
   }
   system ("perl ${directory}Scripts/ERV_get_name.pl -f ${input_sampleID}_ERV1_h.sam -o ${input_sampleID}_ERV1_h.sam2 -Alignment_score $Alignment_score");
   my $AS_XS = 20;
   system ("perl ${directory}Scripts/ERV_filter_reads.pl ${input_sampleID}_ERV1_h.sam2 $AS_XS >${input_sampleID}_ERV1.bian");
   system ("perl ${directory}Scripts/ERV_get_reads.pl ${input_sampleID}_ERV1.bian ${input_sampleID}_ERV1_1.1fq ${input_sampleID}_ERV1_2.1fq ${input_sampleID}_ERV_1.1fq ${input_sampleID}_ERV_2.1fq");
   system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 30 -h 10000 -a -Y -M $TE_reference_genomes ${input_sampleID}_ERV_1.1fq ${input_sampleID}_ERV_2.1fq >${input_sampleID}_ERV_vsu.sam");
   system ("perl ${directory}Scripts/ERV_get_name.pl -f ${input_sampleID}_ERV_vsu.sam -b ${input_sampleID}_ERV1.bian -h ${input_sampleID}_ERV.hf -o ${input_sampleID}_ERV_vsu.sam2");
   system ("perl ${directory}Scripts/ERV_organize_reads.pl ${input_sampleID}_ERV1_h.sam2 ${input_sampleID}_ERV_vsu.sam2 ${input_sampleID}_ERV.hf >${input_sampleID}_ERV_breakpoint");
   system ("perl ${directory}Scripts/ERV_Filtered_fastq.pl ${input_sampleID}");
   system ("cat ${input_sampleID}_all_breakpoint >>${input_sampleID}_ERV_breakpoint");
 } else {
   system ("cat ${input_sampleID}_all_breakpoint >${input_sampleID}_ERV_breakpoint");
 }
 system ("perl ${directory}Scripts/Reads_summary.pl -i ${input_sampleID}_ERV_breakpoint -o ${input_sampleID}_ERV_summary"); 
 system ("perl ${directory}Scripts/Order_by_TE_sequence.pl -i ${input_sampleID}_ERV_summary -o ${input_sampleID}_ERV.integration");
 $double_length_insertsize = $length_insertsize * 2;
 system ("perl ${directory}Scripts/Filtered_single_reads.pl -i ${input_sampleID}_ERV.integration -o ${input_sampleID}_ERV.integration2 -r $double_length_insertsize");
 system ("mv ${input_sampleID}_ERV.integration2 ${input_sampleID}_ERV.integration");
 if ($sequencing_type eq "paired-end"){
   system ("cat ${input_sampleID}_ERV_1.2fq ${input_sampleID}_1.1fuq >${input_sampleID}_ERV_1.1fuq");
   system ("cat ${input_sampleID}_ERV_2.2fq ${input_sampleID}_2.1fuq >${input_sampleID}_ERV_2.1fuq");
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

 ##### 2.5 Organizing output files
 my @path_current2=split /\//, ${input_sampleID};
 my $path_current2;
 if($#path_current2>0){
   my $path_current2=join("/",@path_current2[0..$#path_current2-1]);
   my $input3=$path_current2[$#path_current2];
 } else {
   $path_current2= Cwd::cwd();
   $path_current2=$path_current2."/";
 }
 out1(${input_sampleID}."_ERV");
 Clean_tmp_files_detection();

 ### Step 3. Reciprocal alignemnt to the candidate genomic regions
 print "\nStep 3: Validation...\n=====================================\n";
 if (($file_suffix =~ "fq" || $file_suffix =~ "fastq") && !(-e ${input_sampleID}."_h.bam.bai")) {
   print "# Converting SAM to BAM file, and then Sort and index the BAM file......\n";
   system ("${samtools_d}samtools view -b -@ $thread_1 ${input_sampleID}_h.sam >${input_sampleID}_h.bam");
   system ("${samtools_d}samtools sort -@ $thread_1 ${input_sampleID}_h.bam -o ${input_sampleID}_h.sort.bam");
   system ("mv ${input_sampleID}_h.sort.bam ${input_sampleID}_h.bam");
   system ("${samtools_d}samtools index ${input_sampleID}_h.bam");
 }                                                                                       
 open OUTPUT,">${input_sampleID}_ERV.output2";
 print OUTPUT "$header3";
 my $line="";
 my $line2="";
 my $seq="";
 %genome=();
 my @line=();
 my $title="";
 my %TE_f2=();
 my %countf_c=();
 my %countf_c2=();
 my %countf_sf=();
 my $count_c=();
 my $count_sf=();

 ### 3.1 Read human genome 
 open GENOME2, "$Human_reference_genome";
 while($line2=<GENOME2>){ 
   @line=split /\s+/, $line2;
   unless(@line){next;}
   $line=$line[0];
   chomp($line); 
   if($line=~">"){
     $title=$line;
     $title=~s/>//;
     $genome{$title}="";
   } else {
     $genome{$title}.=$line;
   }
 }

 ### 3.2 Read TE_f2 file
 open TE_F2,"${input_sampleID}_ERV.TE_f2";
 while ($line2=<TE_F2>){
   @line=split /\s+/, $line2;
   $line[15]=~s/_ERV//;
   my $name1=$line[15]."%".$line[31]."%".$line[32]."%".$line[33];
   my @array=split /\|/, $line[30];
   @{$TE_f2{$name1}}=@array;
 }
 system ("mkdir ${input_sampleID}_subgenome");
 open SUB_REF, ">${input_sampleID}_subgenome/${input_sampleID}_sub.fa";
 open SUB_1, ">${input_sampleID}_subgenome/${input_sampleID}_ERV_sub_1.2fuq";
 open SUB_2, ">${input_sampleID}_subgenome/${input_sampleID}_ERV_sub_2.2fuq";
 open SUB_SF, ">${input_sampleID}_subgenome/${input_sampleID}_ERV_sub_1sf.2fuq";
 my %bian2=();

 ### 3.3 Load Fuq files
 ### 3.3.1 _1.fuq
 my %fuq_left=();
 if ($sequencing_type eq "single-end") {
   system ("touch ${input_sampleID}_ERV_1.1fuq");
   system ("touch ${input_sampleID}_ERV_2.1fuq");
 }
 open FUQ_LEFT,"${input_sampleID}_ERV_1.1fuq";
 while(my $read_name=<FUQ_LEFT>){
   my $read=$read_name;
   my @temp11=split /\s+/, $read;
   $read_name=$temp11[0];
   $read_name=~s/\/[12]\s+$//;
   $read_name=~s/\/[12]$//;
   $read_name=~s/^@//;
   my $temp01=<FUQ_LEFT>;
   my $temp02=<FUQ_LEFT>;
   my $temp03=<FUQ_LEFT>;
   my $temp04=$read.$temp01.$temp02.$temp03;
   $fuq_left{$read_name}=$temp04;
 }

 ####### 3.3.2 _2.fuq
 my %fuq_right=();
 
 open FUQ_RIGHT,"${input_sampleID}_ERV_2.1fuq";
 while(my $read_name=<FUQ_RIGHT>){
   my $read=$read_name;   
   my @temp11=split /\s+/, $read;
   $read_name=$temp11[0];
   $read_name=~s/\/[12]\s+$//;
   $read_name=~s/\/[12]$//;
   $read_name=~s/^@//;
   my $temp01=<FUQ_RIGHT>;
   my $temp02=<FUQ_RIGHT>;
   my $temp03=<FUQ_RIGHT>;
   my $temp04=$read.$temp01.$temp02.$temp03;
   $fuq_right{$read_name}=$temp04;
 }

 ######## 3.3.3 _sf.fuq
 my %fuq_sf=();
 open FUQ_SF,"${input_sampleID}_ERV_1sf.fuq";
 while(my $read_name=<FUQ_SF>){
   my $read=$read_name;
   my @temp11=split /\s+/, $read;
   $read_name=$temp11[0];
   $read_name=~s/\/[12]\s+$//;
   $read_name=~s/\/[12]$//;
   $read_name=~s/^@//;
   my @read_name=split /\|/, $read_name;
   $read_name=$read_name[0];
   my $temp01=<FUQ_SF>;
   my $temp02=<FUQ_SF>;
   my $temp03=<FUQ_SF>;
   my $temp04=$read.$temp01.$temp02.$temp03;
   $fuq_sf{$read_name}=$temp04;
 }

  ### 3.3 Read output file
  open INPUT, "${input_sampleID}_ERV.output"; 
  while($line2=<INPUT>){
    @line=split /\s+/, $line2;
    if ($line[0] =~ "Sample_ID"){next;}
    my $seq_name=$line[0]."%".$line[7]."%".$line[8]."%".$line[9];
    my $length=abs($line[9]-$line[8])+($length_insertsize+2*$std_insertsize)*4;                  ##### adjusted by the std 2019/01/30, the SD shoud not be too big;
    $seq=substr $genome{$line[7]},$line[8]-($length_insertsize+2*$std_insertsize)*2,$length;     ##### adjusted by the std 2019/01/30, the SD should not be too big;
    print SUB_REF ">$seq_name\n$seq\n";

    ### 3.4 Extract supporting reads
    my @array = @{$TE_f2{$seq_name}};
    for(my $i=0;$i<@array;$i+=2){
      push (@{$bian2{$seq_name}},$array[$i]);
      if ($sequencing_type eq "paired-end") {
        if(exists($fuq_left{$array[$i]})){print SUB_1 "$fuq_left{$array[$i]}";}
        if(exists($fuq_right{$array[$i]})){print SUB_2 "$fuq_right{$array[$i]}";}
      }
      if ($sequencing_type eq "single-end" || $Split) {
        if(exists($fuq_sf{$array[$i]})){print SUB_SF "$fuq_sf{$array[$i]}";}
      }
    } 
  }            
  close INPUT;
  close SUB_1;
  close SUB_2;
  close SUB_SF;

  ### 3.5 Align to the ref. in pool
    system ("${bwa_d}bwa index -a bwtsw ${input_sampleID}_subgenome/${input_sampleID}_sub.fa");
  if ($sequencing_type eq "paired-end") {
    system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 30 -h 10000 -a -Y -M ${input_sampleID}_subgenome/${input_sampleID}_sub.fa ${input_sampleID}_subgenome/${input_sampleID}_ERV_sub_1.2fuq ${input_sampleID}_subgenome/${input_sampleID}_ERV_sub_2.2fuq >${input_sampleID}_subgenome/${input_sampleID}_ERV_sub.sam");
    system ("sort ${input_sampleID}_subgenome/${input_sampleID}_ERV_sub.sam | uniq > ${input_sampleID}_subgenome/${input_sampleID}_ERV_sub.sam2");
    system ("perl ${directory}Scripts/EC_filtered1.pl ${input_sampleID}_subgenome/${input_sampleID}_ERV_sub.sam2 >${input_sampleID}_subgenome/${input_sampleID}_ERV_sub_bp1");
  } else {
    system ("touch ${input_sampleID}_subgenome/${input_sampleID}_ERV_sub_bp1");
  }
  system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 30 -h 10000 -a -Y -M ${input_sampleID}_subgenome/${input_sampleID}_sub.fa ${input_sampleID}_subgenome/${input_sampleID}_ERV_sub_1sf.2fuq >${input_sampleID}_subgenome/${input_sampleID}_ERV_sub_sf.sam");
  system ("${samtools_d}samtools view -F 4 ${input_sampleID}_subgenome/${input_sampleID}_ERV_sub_sf.sam | uniq > ${input_sampleID}_subgenome/${input_sampleID}_ERV_sub_sf.sam2");

  ### bp1 to bp2
  my %bp1_hash=();
  open BP1,"${input_sampleID}_subgenome/${input_sampleID}_ERV_sub_bp1";
  while (my $bp1_tmp1=<BP1>) {
    my @bp1_tmp1=split /\s+/, $bp1_tmp1;
    $bp1_tmp1 =~ s/^\s+//;
    $bp1_tmp1 =~ s/^\t//;
    push (@{$bp1_hash{$bp1_tmp1[3]}},$bp1_tmp1);
  }
  close BP1;
  my @bp1_hash=keys %bp1_hash;

  open BP2, ">${input_sampleID}_subgenome/${input_sampleID}_ERV_sub_bp2";
  for (my $i=0;$i<@bp1_hash;$i++){
    my @bp1_tmp2 = @{$bp1_hash{$bp1_hash[$i]}};
    my @bp1_tmp3 = ();
    @bp1_tmp2 = sort (@bp1_tmp2);
    my @bp1_tmp4 = split /\s+/, $bp1_tmp2[0];
    push (@bp1_tmp3,[@bp1_tmp4]);
    for (my $bp1_i=1;$bp1_i<@bp1_tmp2;$bp1_i++) {
      unless ($bp1_tmp2[$bp1_i] eq $bp1_tmp2[$bp1_i-1]) {
        my @bp1_tmp4 = split /\s+/, $bp1_tmp2[$bp1_i];
        push (@bp1_tmp3,[@bp1_tmp4]);
      }
    }
    @bp1_tmp3 = sort {$a->[1] cmp $b->[1]} @bp1_tmp3;
    reciprocal_validation_filter(@bp1_tmp3); 
  }
  close BP2;

 ### 3.6 split into each candidate region
 my %sub_bp2=();
 open SUB_BP2, "${input_sampleID}_subgenome/${input_sampleID}_ERV_sub_bp2";
 while (my $line_sub=<SUB_BP2>){
   my @line_sub=split /\s+/, $line_sub;
   my $name_sub=$line_sub[6]."|".$line_sub[4];
   $sub_bp2{$name_sub}=$line_sub;
 }
 my %sub_sf=();
 open SUB_SF, "${input_sampleID}_subgenome/${input_sampleID}_ERV_sub_sf.sam2";
 while (my $sub_sf=<SUB_SF>){
   my @sub_sf=split /\s+/, $sub_sf;
   my @sub_sf_header=split /\|/, $sub_sf[0];
   my $name_sub_sf=$sub_sf[2]."|".$sub_sf[0];
   my $name_sub_sf_header=$sub_sf[2]."|".$sub_sf_header[0];
   $sub_sf{$name_sub_sf}=$sub_sf;
   $sub_sf{$name_sub_sf_header}=$sub_sf;
 } 

 ### 3.7 read visualization file;
 my %vis_hash=();
 my $vis_loci="";
 open VISUALIZATION_EACH, "${input_sampleID}_ERV.visualization";
 while (my $vis_line=<VISUALIZATION_EACH>){
   my @vis_line=split /\s+/, $vis_line;
   if (!($vis_line[0]) || $vis_line[0] ne "O1" || $vis_line[0] ne "O2"){
     next;
   }
   elsif ($vis_line[0] eq "O1"){
     $vis_loci=$vis_line[18]."%".$vis_line[19]."%".$vis_line[20];}
   elsif ($vis_line[0] eq "O2"){
     my $vis_loci2=$vis_loci."%".$vis_line[1];
     @{$vis_hash{$vis_loci2}}=@vis_line;
   }
 }
 
 my %sub_bp2_each=();
 my %sub_sf_each=();

 ### 3.8 Calculate false chimeric and split reads
 open INPUT, "${input_sampleID}_ERV.output";
 while($line2=<INPUT>){
   @line=split /\s+/, $line2;
   @line[22..27]= (0) x 6;
   if ($line[0] =~ "Sample_ID"){next;}
   my $seq_name=$line[0]."%".$line[7]."%".$line[8]."%".$line[9];
   %countf_c=();
   %countf_c2=();
   %countf_sf=();
   $count_c=0;
   $count_sf=0;
   #### get bp2 file
   my @sub_bp=@{$bian2{$seq_name}};
   for (my $i=0; $i<@sub_bp;$i++){
     my $name_sub_bp=${seq_name}."|".$sub_bp[$i];
     if (exists ($sub_bp2{$name_sub_bp})){
       push (@{$sub_bp2_each{$seq_name}},$sub_bp2{$name_sub_bp});
                                         }
     if (exists ($fuq_left{$sub_bp[$i]})){$count_c++;}
     if (exists ($sub_sf{$name_sub_bp})){
       push (@{$sub_sf_each{$seq_name}},$sub_sf{$name_sub_bp});
                                        }
     if (exists ($fuq_sf{$sub_bp[$i]})){$count_sf++;}
   }
   my @line2=();
   if (exists ($sub_bp2_each{$seq_name})) {
     for (my $i=0;$i<@{$sub_bp2_each{$seq_name}};$i++){
       @line2=split /\s+/, ${$sub_bp2_each{$seq_name}}[$i];
       if ($line2[0] eq "P") {
         $countf_c{$line2[4]}="P";
       } elsif (!(exists($countf_c{$line2[4]})) && !(exists($countf_c2{$line2[4]})) && $line2[0] ne "P"){
         $countf_c2{$line2[4]}=$line2[0];
       } elsif (!(exists($countf_c{$line2[4]})) && exists($countf_c2{$line2[4]}) && $line2[0] ne "P"){
         if ($countf_c2{$line2[4]} ne $line2[0] && $countf_c2{$line2[4]} ne "P2"){$countf_c2{$line2[4]}="P2";}
       }
     }
   }

   if (exists ($sub_sf_each{$seq_name})) {
     for (my $i=0;$i<@{$sub_sf_each{$seq_name}};$i++){
       @line2=split /\s+/, ${$sub_sf_each{$seq_name}}[$i];
       $countf_sf{$line2[0]}=$line2;
     }
   }

   $line[22]=$count_c;
   $line[26]=$count_sf;
   my @temp=keys %countf_c;
   $line[24]=$#temp+1;
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
   $line[25]=$countf_c2;
   $line[23]=$countf_c3;
   @temp=keys %countf_sf;
   if($count_sf>0){
     $line[27]=$count_sf-($#temp+1);
   } else {
     $line[27]=0;
   }
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
  }
  close OUTPUT;

system ("perl ${directory}Scripts/Filtered_ERVcaller_1.pl ${input_sampleID}_ERV.output2 $read_len $number_of_reads $sequencing_type");

################################## Step 4. Genotyping
if ($Genotype){
  print "\nStep 4: Genotyping...\n=====================================\n";
  if (($file_suffix =~ "fq" || $file_suffix =~ "fastq") && !(-e ${input_sampleID}."_h.bam.bai")) {
    print "# Converting SAM to BAM file, and then Sort and index the BAM file......\n";
    system ("${samtools_d}samtools view -b -@ $thread_1 ${input_sampleID}_h.sam >${input_sampleID}_h.bam");
    system ("${samtools_d}samtools sort -@ $thread_1 ${input_sampleID}_h.bam -o ${input_sampleID}_h.sort.bam");
    system ("mv ${input_sampleID}_h.sort.bam ${input_sampleID}_h.bam");
    system ("${samtools_d}samtools index ${input_sampleID}_h.bam"); 
  }
  open INPUT, "${input_sampleID}_ERV.output2.1";
  open OUTPUT,">${input_sampleID}_ERV.output3";
  print OUTPUT "$header3";
  my $line2="";
  my @line=();
  my $Chr;
  my $Position;
  my $R_Position;
  my $Breakpoint;
  my $Number_reads;
  my %Chr_posi = ();              ##### add 02/24/2019
  my %Each_TE1 = ();
  my %Each_TE2 = ();

  if (-e ${input_sampleID}."_h.sort.sam2") {system ("rm ${input_sampleID}_h.sort.sam2");}
  system ("touch ${input_sampleID}_h.sort.sam2");

  ### read output2.1 file
  while ($line2=<INPUT>){
    @line=split /\s+/, $line2;
    @line[28..34]=("-") x 7;
    if ($line[0] =~ "Sample_ID"){next;}
    $Chr=$line[7];
    $Position=$line[10];
    $R_Position=$line[11];
    $Breakpoint=0;
    $Number_reads=$line[20];
    if ($Position eq "-") {
      $Position=$R_Position;
      $Breakpoint=1;
    }
    if ($R_Position eq "-") {
      $R_Position=$Position;
      $Breakpoint=1;
    }
    if ($Position ne "-" && $R_Position ne "-") {
      $Breakpoint=2;
    }
    if ($sequencing_type eq "single-end") {
      $Position=$line[15];
      $R_Position=$Position;
      $Breakpoint = 1;
    }
    my $position1=$Position-int(($length_insertsize+2*$std_insertsize)*2);
    my $position2=$Position+int(($length_insertsize+2*$std_insertsize)*2);
    my $name_tmp1 = $Chr.":".$Position."-".$R_Position;  ## 02/14/2019
    @{$Chr_posi{$name_tmp1}} = @line;                    ## 02/14/2019
    if (!($Each_TE1{$name_tmp1})) {
      @{$Each_TE1{$name_tmp1}} = (0) x 6;
      ${$Each_TE1{$name_tmp1}}[0] = ${input_sampleID};
      ${$Each_TE1{$name_tmp1}}[1] = $Chr;
      ${$Each_TE1{$name_tmp1}}[2] = $Position;
      ${$Each_TE1{$name_tmp1}}[3] = $Number_reads;
      ${$Each_TE1{$name_tmp1}}[4] = $Breakpoint;
      ${$Each_TE1{$name_tmp1}}[5] = $read_len;
    }
    if ($Position ne $R_Position && !($Each_TE2{$name_tmp1})) {
      @{$Each_TE2{$name_tmp1}} = (0) x 6;
      ${$Each_TE2{$name_tmp1}}[0] = ${input_sampleID};
      ${$Each_TE2{$name_tmp1}}[1] = $Chr;
      ${$Each_TE2{$name_tmp1}}[2] = $R_Position;
      ${$Each_TE2{$name_tmp1}}[3] = $Number_reads;
      ${$Each_TE2{$name_tmp1}}[4] = $Breakpoint;
      ${$Each_TE2{$name_tmp1}}[5] = $read_len;
    }
    my $cmd_1 = q(awk '{print");          ## 02/14/2019
    my $cmd_2 = q(\t"$0}');
    if ($file_suffix =~ "fq" || $file_suffix =~ "fastq"){
      system ("${samtools_d}samtools view ${input_sampleID}_h.bam ${Chr}:${position1}-${position2} | sort | uniq | $cmd_1$name_tmp1$cmd_2 >>${input_sampleID}_h.sort.sam2"); ## 02/14/2019
    } else {
      if ($multiple_BAM) {
        open LIST, "${Input_directory}${input_sampleID}${file_suffix}";
        while (my $list_tmp1=<LIST>){
          my @list_tmp1=split /\s+/, $list_tmp1;
          if (-e $list_tmp1[0]) {
            system ("${samtools_d}samtools view $list_tmp1[0] ${Chr}:${position1}-${position2} | sort | uniq >>${input_sampleID}_h.sort.sam");
          }
        }
	close LIST;
        system ("sort ${input_sampleID}_h.sort.sam | uniq |$cmd_1$name_tmp1$cmd_2 >>${input_sampleID}_h.sort.sam2");             ##### revised 02/14/2019; 04/23/2019
        system ("rm ${input_sampleID}_h.sort.sam");
      } else {
        system ("${samtools_d}samtools view ${Input_directory}${input_sampleID}${file_suffix} ${Chr}:${position1}-${position2} | sort | uniq | $cmd_1$name_tmp1$cmd_2 >>${input_sampleID}_h.sort.sam2");                              ## 02/14/2019
      }
    }
  }

  ### read sam2
  my @sam2 = ();
  my $sam2_id = "";
  my @sam2_tmp3 = ();
  open SAM2, "${input_sampleID}_h.sort.sam2";
  while (my $sam2_tmp1 = <SAM2>) {
    my @sam2_tmp1 = split /\s+/, $sam2_tmp1;
    my @sam2_tmp2 = @sam2_tmp1[1..$#sam2_tmp1];
    if ($sam2_id eq "") {
      $sam2_id = $sam2_tmp1[0];
      push (@sam2,[@sam2_tmp2]);
    } elsif ($sam2_id eq $sam2_tmp1[0] && !(eof(SAM2))) {
      push (@sam2,[@sam2_tmp2]);  
    } elsif ($sam2_id ne $sam2_tmp1[0] || eof (SAM2)) {
      if ($sam2_id eq $sam2_tmp1[0]) {
        push (@sam2,[@sam2_tmp2]);
      }
      @sam2_tmp3 = Genotype_each_TE1(@sam2);
      $min_insertsize = $length_insertsize - 2 * $std_insertsize;
      @sam2_tmp3 = sort {$a->[1] cmp $b->[1] || $b->[0] cmp $a->[0]} @sam2_tmp3;
      my @sam2_tmp4 = Genotype_each_TE2(@sam2_tmp3);
      my @line_1 = ();
      my @line_2 = (); 
      if ($sequencing_type eq "paired-end") {
        if (exists ($Each_TE2{$sam2_id})) {
          @line_1 = Genotype_each_TE3_PE(\@sam2_tmp4,\@{$Each_TE1{$sam2_id}});
          @line_2 = Genotype_each_TE3_PE(\@sam2_tmp4,\@{$Each_TE2{$sam2_id}});
        } else {
          @line_1 = Genotype_each_TE3_PE(\@sam2_tmp4,\@{$Each_TE1{$sam2_id}});
          @line_2 = @line_1;
        }
      } else {
          @line_1 = Genotype_each_TE3_SE(\@sam2_tmp4,\@{$Each_TE1{$sam2_id}});
          @line_2 = @line_1;
      }
      my @sam2_line = @{$Chr_posi{$sam2_id}};
      if ($line_1[1] < $line_2[1]) {
        $sam2_line[28] = $line_1[4];
        $sam2_line[29] = $line_1[1];
        if($line_1[7]){$sam2_line[30] = $line_1[7];}
        if($line_1[5]){$sam2_line[31] = $line_1[5];}
        if($line_1[6]){$sam2_line[32] = $line_1[6];}
        if($line_1[3]){$sam2_line[33] = $line_1[3];}
      } else {
        $sam2_line[28] = $line_2[4];
        $sam2_line[29] = $line_2[1];
        if($line_2[7]){$sam2_line[30] = $line_2[7];}
        if($line_2[5]){$sam2_line[31] = $line_2[5];}
        if($line_2[6]){$sam2_line[32] = $line_2[6];}
        if($line_2[3]){$sam2_line[33] = $line_2[3];}
      } 
      print OUTPUT "@sam2_line\n";
      ### next loci;
      @sam2 = ();
      push (@sam2,[@sam2_tmp2]);
      $sam2_id = $sam2_tmp1[0];
    } 
  }
  close OUTPUT;
  close SAM2;
}

 ##### Convert to VCF format
 print "\nStep 5: Making VCF file...\n=====================================\n";
 if ($Genotype) {
   convert_to_vcf(${input_sampleID}."_ERV.output3");
 } else {
   convert_to_vcf(${input_sampleID}."_ERV.output2.1");
 }
 
 print "\nPrint out:=====================================\n";
 print "# sample_ID: ${input_sampleID}\tInsert size: $length_insertsize\tSD: $std_insertsize\n";
 
 ##### Clean and move tmp files
 Clean_files();
 Move_files();


###############################################################
############################################################### 
############################         Subs functions

##### Fine_mapped file
sub out1{
  my ($input1)=@_;
  my %output=();
  my $output_header="";
  my @output_line=();
  my $RD_left=0;
  my $RD_right=0;
  my $len_left=0;
  my $len_right=0;
  #### fine_mapped file
  open FINE, "${input1}.fine_mapped";
  while (<FINE>) {
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
    $output{$output_header}[2] = "Yes";
    $output{$output_header}[3] = "Yes";
    ${$output{$output_header}}[4]=$Human_reference_genome;
    ${$output{$output_header}}[5]=$TE_reference_genomes;
    ${$output{$output_header}}[6]=$output_line[5];
    @{$output{$output_header}}[7..9]=@output_line[1..3];

    ### upstream breakpoint
    if(($output_line[8] ne "-" && $output_line[17] eq "-") || ($output_line[8] ne "-" && $output_line[17] ne "-" && (($output_line[22]+$output_line[23]) < ($output_line[13] + $output_line[14]) || (($output_line[22]+$output_line[23]) eq ($output_line[13] + $output_line[14]) && $output_line[35] ne "-" && $output_line[26] eq "-")))){
      if ($output_line[7] eq "-") {
        $output{$output_header}[14]="E(++);";
        $output{$output_header}[10]=$output_line[9];
        $output{$output_header}[13]=$output_line[12];
      } else {
        $output{$output_header}[14]="D(++);";
        $output{$output_header}[10]=$output_line[7];
        $output{$output_header}[13]=$output_line[10];                          
      }
      $RD_left=$output_line[13]+$output_line[14];
      $len_left=$output_line[12]-$output_line[11]+1;   
    } elsif (($output_line[8] eq "-" && $output_line[17] ne "-") || ($output_line[8] ne "-" && $output_line[17] ne "-" && (($output_line[22]+$output_line[23]) >= ($output_line[13] + $output_line[14])))){
      if ($output_line[16] eq "-") {
        $output{$output_header}[14]="E(+-);";
        $output{$output_header}[10]=$output_line[18];
        $output{$output_header}[12]=$output_line[20];
      } else {
        $output{$output_header}[14]="D(+-);";
        $output{$output_header}[10]=$output_line[16];
        $output{$output_header}[12]=$output_line[19];
      }
      if ($sequencing_type eq "single-end") {
        $output{$output_header}[10]=$output_line[18];
      }
      $RD_left=$output_line[22]+$output_line[23];
      $len_left=$output_line[21]-$output_line[20]+1;
    } elsif ($output_line[8] eq "-" && $output_line[17] eq "-") {
        $output{$output_header}[14]="na;";
    }

    ### downstream breakpoint
    if(($output_line[26] ne "-" && $output_line[35] eq "-") || ($output_line[26] ne "-" && $output_line[35] ne "-" && (($output_line[40]+$output_line[41]) < ($output_line[31] + $output_line[32]) || (($output_line[40]+$output_line[41]) eq ($output_line[31] + $output_line[32]) && $output_line[17] ne "-" && $output_line[8] eq "-")))){
      if ($output_line[25] eq "-") {
        $output{$output_header}[14]=$output{$output_header}[14]."e(-+)";
        $output{$output_header}[11]=$output_line[26];                             #### right breakpoint
        $output{$output_header}[13]=$output_line[30];
      } else {
        $output{$output_header}[14]=$output{$output_header}[14]."d(-+)";
        $output{$output_header}[11]=$output_line[25];
        $output{$output_header}[13]=$output_line[28];
      }
      if ($sequencing_type eq "single-end") {
        $output{$output_header}[11]=$output_line[26];
      }
      $RD_right=$output_line[31]+$output_line[32];
      $len_right=$output_line[30]-$output_line[29]+1;
    } elsif (($output_line[26] eq "-" && $output_line[35] ne "-") || ($output_line[26] ne "-" && $output_line[35] ne "-" && ($output_line[40]+$output_line[41]) >= ($output_line[31] + $output_line[32]))){
      if ($output_line[34] eq "-") {
        $output{$output_header}[14]=$output{$output_header}[14]."e(--)";
        $output{$output_header}[11]=$output_line[35];
        $output{$output_header}[12]=$output_line[38];
      } else {
        $output{$output_header}[14]=$output{$output_header}[14]."d(--)";
        $output{$output_header}[11]=$output_line[34];
        $output{$output_header}[12]=$output_line[37];
      }
      $RD_right=$output_line[40]+$output_line[41];
      $len_right=$output_line[39]-$output_line[38]+1;
    } elsif ($output_line[26] eq "-" && $output_line[35] eq "-") {
      $output{$output_header}[14]=$output{$output_header}[14]."na";
    }

    ### Filter false TEs
    if ($output{$output_header}[10] ne "-" && $output{$output_header}[11] ne "-") {
      if ($RD_left eq 1 && $RD_right eq 1 && $output{$output_header}[10] - $output{$output_header}[11] >=$TSD_min_len) {
        delete $output{$output_header};
        next;
      } elsif ($RD_left > 1 && $RD_right eq 1 && $output{$output_header}[10] - $output{$output_header}[11] >=$TSD_min_len) {
        $output{$output_header}[11] = "-";
        if ($len_left <= $read_len) {
          delete $output{$output_header};
          next;
        }
      } elsif ($RD_left eq 1 && $RD_right > 1 && $output{$output_header}[10] - $output{$output_header}[11] >=$TSD_min_len) {
        $output{$output_header}[10] = "-";
        if ($len_right <= $read_len) {
          delete $output{$output_header};
          next;
        }
      }  
    } elsif ($output{$output_header}[10] eq "-" && $output{$output_header}[11] ne "-") {
      if ($len_right <= $read_len) {
        delete $output{$output_header};
        next;
      }
    } elsif ($output{$output_header}[10] ne "-" && $output{$output_header}[11] eq "-") {
      if ($len_left <= $read_len) {
        delete $output{$output_header};
        next;
      }
    }

    ### Predicted breakpoint
    if ($output{$output_header}[10] ne "-" && $output{$output_header}[11] ne "-"){
      if (($output{$output_header}[14] =~ "D" && $output{$output_header}[14] =~ "d") || ($output{$output_header}[14] =~ "E" && $output{$output_header}[14] =~ "e")){
        $output{$output_header}[15]=int(($output{$output_header}[11]+$output{$output_header}[10])/2);
      } elsif ($output{$output_header}[14] =~ "D") {
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

  ### read type file
  my %human_AS=();
  my $file_name=${input1};
  $file_name=~s/_ERV$//;
  if (-e ${input1}.".type") {
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

  ### read bian file
  my @line=();
  if(-e ${input1}."1.bian"){
    open BIAN,"${input1}1.bian"; 
    while(<BIAN>){
      @line=split;
      @{$human_AS{$line[0]}}=@line;
    }
    close BIAN;
  }
  ### read TE_f file
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
     for (my $i=0;$i<@list_reads2;$i++){
       $calculate[0]++;
       $calculate[1]=$calculate[1]+${$sub_human_AS{$list_reads2[$i]}}[2];
       $calculate[2]=$calculate[2]+${$sub_human_AS{$list_reads2[$i]}}[3];
       if ($calculate[3]<(${$sub_human_AS{$list_reads2[$i]}}[2]-${$sub_human_AS{$list_reads2[$i]}}[3])){
         $calculate[3]=${$sub_human_AS{$list_reads2[$i]}}[2]-${$sub_human_AS{$list_reads2[$i]}}[3];
                                                                                                       }
                                       }
     if ($calculate[0]>0){
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
     if ($output{$f2_header}[17]-$output{$f2_header}[18]>=30){
       $output{$f2_header}[16]="Non-repeat;";
     } elsif ($output{$f2_header}[17]-$output{$f2_header}[18]>15 && $output{$f2_header}[17]-$output{$f2_header}[18]<=30) {
       $output{$f2_header}[16]="Likely-repeat_or_false_positive;";
     } elsif ($output{$f2_header}[17]-$output{$f2_header}[18]>5 && $output{$f2_header}[17]-$output{$f2_header}[18]<=15) {
       $output{$f2_header}[16]="Very_likely-repeat_or_false_positive;";
     } elsif ($output{$f2_header}[17]-$output{$f2_header}[18]<=5){
       $output{$f2_header}[16]="Repeat_or_false_positive;";
     } else {
       $output{$f2_header}[16]="NA;";
     }
                                 }
             } 
  ### output
  my @list=sort{$output{$b}->[19]<=>$output{$a}->[19]} keys %output;
  open OUT, ">${input1}.output";
  $" = "\t";
  for (my $i=0;$i<@list;$i++){
    print OUT "@{$output{$list[$i]}}\n";
  }
  close OUT;
}


##### Convert from output2 to vcf format
sub convert_to_vcf {
  my ($input1)=@_;
  open OUT3,"$input1";
  open VCF, ">${input_sampleID}.vcf";
  open TAB, ">${input_sampleID}.tab";
  my $vcf_tmp1=0;
  while (my $out3_tmp1=<OUT3>) {
    my @out3_tmp1=split /\s+/, $out3_tmp1;
    if ($out3_tmp1[0] eq "Sample_ID") {
      next;
    } else {
      if ($vcf_tmp1 eq 0) {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	$year=$year+1900;
	print VCF "##fileformat=VCFv4.2\n";
	print VCF "##fileDate=$year$mon$mday\n";
	$directory[$#directory]=~s/.pl$//;
	print VCF "##source=$directory[$#directory]\n";
	print VCF "##reference=file:$out3_tmp1[4]\n";
        if (%genome) {
          my @geno_key = keys %genome;
          for (my $geno_i = 0;$geno_i < @geno_key; $geno_i++) {
            my $geno_len = length ($genome{$geno_key[$geno_i]});
            print VCF "##contig=<ID=$geno_key[$geno_i],length=$geno_len>\n";
          }
        }
	print VCF "##ALT=<ID=INS:MEI:HERVK,Description=\"HERVK insertion\">\n";
	print VCF "##INFO=<ID=TSD,Number=2,Type=String,Description=\"NUCLEOTIDE,LEN, Nucleotides and length of the Target Site Duplication (NULL for unknown)\">\n";
	print VCF "##INFO=<ID=INFOR,Number=6,Type=String,Description=\"NAME,START,END,LEN,DIRECTION,STATUS; NULL for unknown values. Status of detected TE: 0 = Inconsistent direction for the supporting reads; 1 = One breakpoint detected by only chimeric and/or improper reads without split reads; 2 = Only one breakpoint is detected and covered by split reads; 3 = Two breakpoints detected, and both of them are not covered by split reads; 4 = Two breakpoints detected, and one of them are not covered by split reads; 5 = Two breakpoints detected, and both of them are covered by split reads;\">\n";
        print VCF "##INFO=<ID=CR,Number=1,Type=Integer,Description=\"Number of chimeric and improper reads support the TE insertion\">\n";
        print VCF "##INFO=<ID=SR,Number=1,Type=String,Description=\"Number of split reads support TE insertion and the breakpoint\">\n";
        print VCF "##INFO=<ID=GTF,Number=1,Type=String,Description=\"If the detected TE insertions genotyped\">\n";
        print VCF "##INFO=<ID=GR,Number=1,Type=Float,Description=\"The ratio of the number of reads support TE insertions versus the total number of reads at this TE insertion location\">\n";
        print VCF "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        print VCF "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype quality (Phred transformed)\">\n";
        print VCF "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihood\">\n";
        print VCF "##FORMAT=<ID=DPI,Number=1,Type=Integer,Description=\"The number of reads support TE insertions\">\n";
        print VCF "##FORMAT=<ID=DPN,Number=1,Type=Integer,Description=\"The number of reads support non-TE insertions\">\n";
        print VCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$out3_tmp1[0]\n";
      }
      ### extracting information of each TE insertion 
      my $ref_seq=".";
      my $tsdlen="NULL";
      my $tsd="NULL";
      my $pos=$out3_tmp1[15];
      my $svtype=$out3_tmp1[6];
      if ($svtype =~ /ERV/i && !($svtype =~ /LINE1/i ) && !($svtype =~ /L1/i ) && !($svtype =~ /ALU/i ) && !($svtype =~ /SVA/i)) {
        $svtype = "HERV";
      } elsif (!($svtype =~ /ERV/i) && ($svtype =~ /LINE1/i || $svtype =~ /L1/i ) && !($svtype =~ /ALU/i ) && !($svtype =~ /SVA/i)) {
        $svtype = "LINE1";
      } elsif (!($svtype =~ /ERV/i) && !($svtype =~ /LINE1/i) && !($svtype =~ /L1/i ) && $svtype =~ /ALU/i && !($svtype =~ /SVA/i)) {
        $svtype = "ALU";
      } elsif (!($svtype =~ /ERV/i) && !($svtype =~ /LINE1/i) && !($svtype =~ /L1/i ) && !($svtype =~ /ALU/i ) && $svtype =~ /SVA/i) {
        $svtype = "SVA";
      }
      if ($out3_tmp1[10] ne "-" && $out3_tmp1[11] ne "-" && $out3_tmp1[10]>$out3_tmp1[11] && !($out3_tmp1[14] =~ "na") && !($out3_tmp1[14] =~ "E")){
        $tsdlen=$out3_tmp1[10]-$out3_tmp1[11];
        $pos=$out3_tmp1[11];                                             	##### 02/06/2019 adjust the position
        if (%genome) {
          $ref_seq=substr $genome{$out3_tmp1[7]},$pos-1,1;
          $tsd=substr $genome{$out3_tmp1[7]},$pos,$tsdlen;			##### updated on 02/06/2019
        }
      } elsif ($out3_tmp1[10] ne "-" && $out3_tmp1[11] ne "-" && abs($out3_tmp1[10] - $out3_tmp1[11])<=1 && !($out3_tmp1[14] =~ "na") && !($out3_tmp1[14] =~ "E")) {  ##### updated on 02/06/2019;
        $tsdlen=0; 
        if (%genome) {
          $ref_seq=substr $genome{$out3_tmp1[7]},$pos-1,1;
        }
      } else {
        if (%genome) {
          $ref_seq=substr $genome{$out3_tmp1[7]},$pos-1,1;
        }
      }

      my $telen="NULL";
      my $testart="NULL";
      my $teend="NULL";
      my $tedirection="NULL";
      my @tedirection= split /\;/, $out3_tmp1[14];
      my $testatus="NULL";
      if ($tedirection[0] =~ "na" || $tedirection[1] =~ "na") { 
        if ($tedirection[0] =~ /D/i || $tedirection[1] =~ /D/i) {
          $testatus = 2;
        } else {
          $testatus = 1;
        }
        if ($tedirection[0] =~ /\+\-/ || $tedirection[1] =~ /\-\+/) {
           $tedirection="+";
        } else {
           $tedirection="-";
        }
        if ($out3_tmp1[13] ne "-") {
          $teend = $out3_tmp1[13]; 
        } elsif ($out3_tmp1[12] ne "-") {
          $testart = $out3_tmp1[12];
        }
      } else {
        if ($tedirection[0] =~ /D/i && $tedirection[1] =~ /D/i) {
          $testatus = 5;
        } elsif ($tedirection[0] =~ /E/i && $tedirection[1] =~ /E/i) {
          $testatus = 3;
        } else {
          $testatus = 4;
        }
        if ($tedirection[0] =~ /\+\-/ && $tedirection[1] =~ /\-\+/) {
          $tedirection="+";
        } elsif ($tedirection[0] =~ /\+\+/ && $tedirection[1] =~ /\-\-/){
          $tedirection="-";
        } else {
          $tedirection="NULL";
          $testatus = 0;
        }
        ### TE start, end
        if ($out3_tmp1[12] ne "-" && $out3_tmp1[13] ne "-") {
          if ($out3_tmp1[13]-$out3_tmp1[12] > 0) {
            $testart = $out3_tmp1[12];
            $teend = $out3_tmp1[13];
          } else {
            $testart = $out3_tmp1[13];
            $teend = $out3_tmp1[12];
          }
            $telen=abs($out3_tmp1[13]-$out3_tmp1[12])+1;
        } elsif ($out3_tmp1[13] ne "-") {
          $teend = $out3_tmp1[13];
        } elsif ($out3_tmp1[12] ne "-") {
          $testart = $out3_tmp1[12];
        }
      }
      ### genotype
      my $tegenotype="NO";
      my $geno="./.";
      my $genoq=".";
      my $genol=".,.,.";
      my $reads_no_insertion=".";
      my $reads_insertion=$out3_tmp1[22] + $out3_tmp1[26];
      my $genoratio="NULL";
      if ($Genotype) {
        $tegenotype="YES";
        $reads_no_insertion=$out3_tmp1[29];
        $genoratio=$reads_insertion/($reads_insertion + $reads_no_insertion);
        $genoratio=sprintf("%.3f",$genoratio);
      }
      $vcf_tmp1++;
      $pos = $pos-1;                                   ###### adjust the position to start from 0;                      2019/01/20
      print TAB "$out3_tmp1[7]\t$pos\t.\t$ref_seq\t<INS_MEI:$svtype>\t.\t.\tTSD=$tsd,$tsdlen;INFOR=$out3_tmp1[6],$testart,$teend,$telen,$tedirection,$testatus;CR=$out3_tmp1[22];SR=$out3_tmp1[26];GTF=$tegenotype;GR=$genoratio\tGT:GQ:GL:DPN:DPI\t$reads_insertion\t$reads_no_insertion\n";
    }
  }
  close OUT3;
  close TAB;

  ##### 
  if ($Genotype) {
    system ("Rscript ${directory}Scripts/Genotype_likelihood.R ${input_sampleID}.tab ${input_sampleID}.tab2");
    open TAB, "${input_sampleID}.tab2";
    while (my $tmp1_tab = <TAB>) {
      my $geno = "./.";
      my $genoq = ".";
      my $genol = ".,.,.";
      my @tmp1_tab = split /\s+/, $tmp1_tab;
      if ($tmp1_tab[13] eq "NA") {
        $tmp1_tab[13] = 0;
      }
      if ($tmp1_tab[14] eq "NA") {
        $tmp1_tab[14] = 0;
      }  
      if ($tmp1_tab[13] > $tmp1_tab[14]) {
        $geno = "1/1";
        if ($tmp1_tab[16] eq "NA" || $tmp1_tab[16] eq "NaN") {
          $genoq = "0";
        } else {
          $genoq = $tmp1_tab[16];
        }
      } else {
        $geno = "0/1";
        if ($tmp1_tab[17] eq "NA" || $tmp1_tab[17] eq "NaN") {
          $genoq = "0";
        } else {
          $genoq = $tmp1_tab[17];
        }
      }
      if ($tmp1_tab[13] eq "NA" || $tmp1_tab[13] eq "NaN") {
        $tmp1_tab[13] = "0";
      }
      if ($tmp1_tab[14] eq "NA" || $tmp1_tab[14] eq "NaN") {
        $tmp1_tab[14] = "0";
      }
      if ($tmp1_tab[15] eq "NA" || $tmp1_tab[15] eq "NaN") {
        $tmp1_tab[15] = "0";
      }
      $genol = $tmp1_tab[15].",".$tmp1_tab[14].",".$tmp1_tab[13];             ##### No_insertion:Heter_insertion:TE_insertion
      if ($tmp1_tab[3] eq "TRUE") {
        $tmp1_tab[3]="T";
      }
      print VCF "$tmp1_tab[0]\t$tmp1_tab[1]\t$tmp1_tab[2]\t$tmp1_tab[3]\t$tmp1_tab[4]\t$tmp1_tab[5]\t$tmp1_tab[6]\t$tmp1_tab[7]\t$tmp1_tab[8]\t$geno:$genoq:$genol:$tmp1_tab[10]:$tmp1_tab[9]\n";
    }
    close TAB;
  } else {
    open TAB, "${input_sampleID}.tab";
    while (my $tmp1_tab = <TAB>) {
      my @tmp1_tab = split /\s+/, $tmp1_tab;
      if ($tmp1_tab[3] eq "TRUE") {
        $tmp1_tab[3]="T";
      }
      print VCF "$tmp1_tab[0]\t$tmp1_tab[1]\t$tmp1_tab[2]\t$tmp1_tab[3]\t$tmp1_tab[4]\t$tmp1_tab[5]\t$tmp1_tab[6]\t$tmp1_tab[7]\tGT\t1/1\n";
    }
    close TAB;
  }
  close VCF;
}

##### Create type file
sub create_type { 
  system ("${samtools_d}samtools view ${input_sampleID}_sm.bam >${input_sampleID}_sm.sam"); 
  open SM,"${input_sampleID}_sm.sam";
  open TYPE,">${input_sampleID}.type";
  my $sm_1="";
  my $as="";
  my $xs="";
  while ($sm_1=<SM>){
    my @sm_1=split /\s+/, $sm_1;
    $as="NA";$xs="NA";
    unless ($sm_1[2] eq "*"){
      for (my $i=11;$i<@sm_1;$i++){
        if ($sm_1[$i]=~s/AS:i://){$as=$sm_1[$i];}
        if ($sm_1[$i]=~s/XS:i://){$xs=$sm_1[$i];}
      }
    }
    if ($sm_1[1]%256>=128) {
      print TYPE "$sm_1[0] L $as $xs $sm_1[5]\n";
    } else {
      print TYPE "$sm_1[0] R $as $xs $sm_1[5]\n";
    }
  }
  close TYPE;
  close SM;
}

##### Align to the human reference genome
sub align_to_hg {
  my ($input1,$suffix1)=@_;
  if ($data_type eq "RNA-seq" && $order eq 0){
    my $suffix2=$suffix1;
    $suffix2=~s/.gz$//;
    if($sequencing_type eq "paired-end"){
      if($suffix1 =~ "gz"){
        system ("gunzip -c ${Input_directory}${input1}_1${suffix1} >${input1}_1${suffix2}");
        system ("gunzip -c ${Input_directory}${input1}_2${suffix1} >${input1}_2${suffix2}");
        system ("${tophat_d}tophat2 -p $threads -o ${input1} $human_genome_tophat ${input1}_1${suffix2} ${input1}_2${suffix2}");
      } else {
        system ("${tophat_d}tophat2 -p $threads -o ${input1} $human_genome_tophat ${Input_directory}${input1}_1${suffix2} ${Input_directory}.${input1}_2${suffix2}");
      } 
    } else {
      if($suffix1 =~ "gz"){
        system ("gunzip -c ${Input_directory}${input1}${suffix1} >${input1}${suffix2}");
        system ("${tophat_d}tophat2 -p $threads -o ${input1} $human_genome_tophat ${input1}${suffix2}");
      } else {
        system ("${tophat_d}tophat2 -p $threads -o ${input1} $human_genome_tophat ${Input_directory}${input1}${suffix2}");
      }
    }
    system ("${samtools_d}samtools merge -f ${input1}_h.bam ${input1}/accepted_hits.bam ${input1}/unmapped.bam");
    system ("${samtools_d}samtools view ${input1}_h.bam >${input1}_h.sam"); 
  } elsif ($order eq 0) {
    if($sequencing_type eq "paired-end"){
        system("${bwa_d}bwa mem -t $threads $human_genome ${Input_directory}${input1}_1${suffix1} ${Input_directory}${input1}_2${suffix1} >${input1}_h.sam");
    } elsif ($sequencing_type eq "single-end") {
        system("${bwa_d}bwa mem -t $threads $human_genome ${Input_directory}${input1}${suffix1} >${input1}_h.sam");
    }
  } elsif ($order eq 1 || $order eq 2) {
    if($sequencing_type eq "paired-end"){
        system("${bwa_d}bwa mem -t $threads $human_genome ${input1}_1${suffix1} ${input1}_2${suffix1} >${input1}_h.sam");
    } elsif ($sequencing_type eq "single-end") {
        system("${bwa_d}bwa mem -t $threads $human_genome ${input1}${suffix1} >${input1}_h.sam");
    }
  }
}

##### Convert bam file to fastq file
sub convert_bamtofastq {
  my ($input1)=@_;
  my $suffix1="";
  if(($file_suffix =~ "bam" || $file_suffix =~ "sam" ) && $order eq 1){
    $suffix1=$file_suffix;
  } else {
    $suffix1="_h.sam";
  }
  if($sequencing_type eq "paired-end"){
    if ($order eq 1){
      if ($multiple_BAM) {
        print "${Input_directory}${input_sampleID}${file_suffix}\n";
        open LIST, "${Input_directory}${input_sampleID}${file_suffix}";
        while (my $list_tmp1=<LIST>){
          my @list_tmp1=split /\s+/, $list_tmp1;
          unless (-e ${input1}."_su.bam") {
            system ("${samtools_d}samtools view -b -f 4 -F 264 -@ $thread_1 $list_tmp1[0] > ${input1}_su.bam");
            system ("${samtools_d}samtools view -b -f 8 -F 260 -@ $thread_1 $list_tmp1[0] > ${input1}_sm.bam");
            system ("${samtools_d}samtools view -b -f 12 -F 256 -@ $thread_1 $list_tmp1[0] > ${input1}_pe.bam");
          } else { 
            system ("${samtools_d}samtools view -b -f 4 -F 264 -@ $thread_1 $list_tmp1[0] > ${input1}_1su.bam");
            system ("${samtools_d}samtools view -b -f 8 -F 260 -@ $thread_1 $list_tmp1[0] > ${input1}_1sm.bam");
            system ("${samtools_d}samtools view -b -f 12 -F 256 -@ $thread_1 $list_tmp1[0] > ${input1}_1pe.bam");
	    system ("${samtools_d}samtools merge -f ${input1}_2su.bam ${input1}_su.bam ${input1}_1su.bam");
            system ("${samtools_d}samtools merge -f ${input1}_2sm.bam ${input1}_sm.bam ${input1}_1sm.bam");
            system ("${samtools_d}samtools merge -f ${input1}_2pe.bam ${input1}_pe.bam ${input1}_1pe.bam");
            system ("mv ${input1}_2su.bam ${input1}_su.bam");
            system ("mv ${input1}_2sm.bam ${input1}_sm.bam");
            system ("mv ${input1}_2pe.bam ${input1}_pe.bam");
	    system ("rm ${input_sampleID}_1su.bam");
            system ("rm ${input_sampleID}_1sm.bam");
            system ("rm ${input_sampleID}_1pe.bam");
          }
        }
        close LIST;
      } else {
        system ("${samtools_d}samtools view -b -f 4 -F 264 -@ $thread_1 ${Input_directory}${input1}$suffix1 > ${input1}_su.bam");
        system ("${samtools_d}samtools view -b -f 8 -F 260 -@ $thread_1 ${Input_directory}${input1}$suffix1 > ${input1}_sm.bam");
        system ("${samtools_d}samtools view -b -f 12 -F 256 -@ $thread_1 ${Input_directory}${input1}$suffix1 > ${input1}_pe.bam");
      }
    } else {
      system ("${samtools_d}samtools view -b -f 4 -F 264 -@ $thread_1 ${input1}$suffix1 > ${input1}_su.bam");
      system ("${samtools_d}samtools view -b -f 8 -F 260 -@ $thread_1 ${input1}$suffix1 > ${input1}_sm.bam");
      system ("${samtools_d}samtools view -b -f 12 -F 256 -@ $thread_1 ${input1}$suffix1 > ${input1}_pe.bam");
    }
    if ($Split){
      if ($order eq 1){
	if ($multiple_BAM) {
          open LIST, "${Input_directory}${input_sampleID}${file_suffix}";
          while (my $list_tmp1=<LIST>){
            my @list_tmp1=split /\s+/, $list_tmp1;
	    system ("${SE_MEI_d}extractSoftclipped -l $Split $list_tmp1[0] >${input1}_1soft.fastq.gz");
	    system ("gunzip ${input1}_1soft.fastq.gz");
	    system ("cat ${input1}_1soft.fastq >>${input1}_soft.fastq");
	    system ("rm ${input1}_1soft.fastq");
          }
          close LIST;
	  system ("gzip ${input1}_soft.fastq");
        } else {
          system ("${SE_MEI_d}extractSoftclipped -l $Split ${Input_directory}${input1}$suffix1 >${input1}_soft.fastq.gz");
	}
      } else {
        system ("${SE_MEI_d}extractSoftclipped -l $Split ${input1}$suffix1 >${input1}_soft.fastq.gz");
      }
    }
    system ("${samtools_d}samtools merge -f ${input1}_h1.bam ${input1}_su.bam ${input1}_sm.bam ${input1}_pe.bam");
    system ("${samtools_d}samtools sort -n -@ $thread_1 ${input1}_h1.bam -o ${input1}_h1$suffix1");
    system ("${samtools_d}samtools fastq -@ $thread_1 -N ${input1}_h1$suffix1 -1 ${input1}_h1_1.1fq -2 ${input1}_h1_2.1fq");
    system ("perl ${directory}Scripts/Check_paired_end.pl -s ${input1}_h1 -f .1fq");
    if (-e ${input1}."_h1_1.1fq2") {system ("mv ${input1}_h1_1.1fq2 ${input1}_h1_1.1fq");}
    if (-e ${input1}."_h1_2.1fq2") {system ("mv ${input1}_h1_2.1fq2 ${input1}_h1_2.1fq");}
    if (-e ${input1}."_h1.bam") {system ("rm ${input1}_h1.bam");}
  } elsif ($sequencing_type eq "single-end") {
      if ($order eq 1) {
	if ($multiple_BAM) {
          open LIST, "${Input_directory}${input_sampleID}${file_suffix}";
          while (my $list_tmp1=<LIST>){
            my @list_tmp1=split /\s+/, $list_tmp1;
            system ("${SE_MEI_d}extractSoftclipped -l $Split $list_tmp1[0] >${input1}_1soft.fastq.gz");
            system ("gunzip ${input1}_1soft.fastq.gz");
            system ("cat ${input1}_1soft.fastq >>${input1}_soft.fastq");
            system ("rm ${input1}_1soft.fastq");
            unless (-e ${input1}."_u.bam") {
	      system ("${samtools_d}samtools view -b -f 4 -@ $thread_1 $list_tmp1[0] > ${input1}_u.bam");
	    } else {
	      system ("${samtools_d}samtools view -b -f 4 -@ $thread_1 $list_tmp1[0] > ${input1}_1u.bam");
	      system ("${samtools_d}samtools merge -f ${input1}_2u.bam ${input1}_u.bam ${input1}_1u.bam");
	      system ("mv ${input1}_2u.bam ${input1}_u.bam");
	      system ("rm ${input_sampleID}_1u.bam");
	    }
	  }
	  close LIST;
	} else {
          system ("${SE_MEI_d}extractSoftclipped -l $Split ${Input_directory}${input1}$suffix1 >${input1}_soft.fastq.gz");
          system("${samtools_d}samtools view -b -f 4 -@ $thread_1 ${Input_directory}${input1}$suffix1 >${input1}_u.bam");
	}
      } else {
        system ("${SE_MEI_d}extractSoftclipped -l $Split ${input1}$suffix1 >${input1}_soft.fastq.gz");
        system("${samtools_d}samtools view -b -f 4 -@ $thread_1 ${input1}$suffix1 >${input1}_u.bam");
        if ($suffix1 =~ "sam") {
          system("${samtools_d}samtools view -b -@ $thread_1 ${input1}$suffix1 >${input1}_h.bam");
        }
        $suffix1="_h.bam";
      }       
      system ("${samtools_d}samtools fastq -@ $thread_1 -N ${input1}_u.bam -1 ${input1}_h2.1fq -2 ${input1}_h1.2fq >${input1}_h1.1fq");
      if(-e ${input1}."_h2.1fq") {system("rm ${input1}_h2.1fq");}
      if(-e ${input1}."_h2.2fq") {system("rm ${input1}_h2.2fq");}
  }
}

###### calculate insert size                ### newly added 2019/01/29
sub estimate_insertsize {
  my ($input1,$filename,$proportion)=@_;
  my @insert_size = ();
  my $total = 0;
  system ("${samtools_d}samtools view -@ $thread_1 -s $proportion $input1 > ${filename}_$proportion");
  open SUB_SAM, "${filename}_$proportion";
  while (my $tmp1_sam=<SUB_SAM>) {
    my @tmp1_sam = split /\s+/, $tmp1_sam;
    if ($tmp1_sam[1] % 2 >= 1 && $tmp1_sam[1] % 4 >= 2 && $tmp1_sam[1] % 128 >= 64 && $tmp1_sam[1] % 512 < 256 && $tmp1_sam[1] < 256) {
      push (@insert_size, abs($tmp1_sam[8]));
      $total += abs($tmp1_sam[8]);
    }
  }
  close SUB_SAM;
  system ("rm ${filename}_$proportion");
  my $average = $total / ($#insert_size+1);
  $average = sprintf("%.4f",$average);
  my $sqtotal =0;
  foreach(@insert_size) {
    $sqtotal += ($average-$_) ** 2;
  }
  my $std = ($sqtotal / $#insert_size) ** 0.5;
  $std = sprintf("%.4f",$std);
  return ($average,$std);
}

###### reciprocal validation filter 
sub reciprocal_validation_filter {
  my @bp1_tmp4 = @_;
  my $line_bp2="";
  my %read_bp2=();
  my $read_name_bp2="";
  my $name_bp2="";
  my $order_bp2=0;
  for (my $bp_i=0;$bp_i<@bp1_tmp4;$bp_i++) {
    my @line_bp2 = @{$bp1_tmp4[$bp_i]};
    if($line_bp2[2] / 512>=1){next;}
    $line_bp2[5]=$line_bp2[3];
    my $name2_bp2=$line_bp2[1]."_".$line_bp2[3]."_".$line_bp2[8];
    $name_bp2=$line_bp2[1]."_".$line_bp2[3]."_".$line_bp2[4];
    if(!(%read_bp2)){
      @{$read_bp2{$name_bp2}}=@line_bp2;
      $read_name_bp2=$line_bp2[1];
      unshift(@{$read_bp2{$name_bp2}},$line_bp2[11]);
      unshift(@{$read_bp2{$name_bp2}},$line_bp2[3]);
      unshift(@{$read_bp2{$name_bp2}},$line_bp2[0]);
    }
    if($line_bp2[1] eq $read_name_bp2){
      if(exists($read_bp2{$name2_bp2}) && (($line_bp2[0] eq "PF" && ${$read_bp2{$name2_bp2}}[0] eq "PR") || ($line_bp2[0] eq "PR" && ${$read_bp2{$name2_bp2}}[0] eq "PF"))) {
        if($line_bp2[3] eq ${$read_bp2{$name2_bp2}}[6] && $line_bp2[4] eq ${$read_bp2{$name2_bp2}}[11] && $line_bp2[8] eq ${$read_bp2{$name2_bp2}}[7]){
          splice @{$read_bp2{$name2_bp2}},3,0,@line_bp2;
          ${$read_bp2{$name2_bp2}}[0]="P";
          ${$read_bp2{$name2_bp2}}[2]+=$line_bp2[11];
        } else { 
          if(${$read_bp2{$name2_bp2}}[6] ne $line_bp2[3]){
            $name_bp2=$line_bp2[1]."_".$line_bp2[3]."_".$line_bp2[4];
            @{$read_bp2{$name_bp2}}=@line_bp2;
            unshift(@{$read_bp2{$name_bp2}},$line_bp2[11]);
            unshift(@{$read_bp2{$name_bp2}},$line_bp2[3]);
            unshift(@{$read_bp2{$name_bp2}},$line_bp2[0]);
          } else {
            splice @{$read_bp2{$name2_bp2}},3,0,@line_bp2; ${$read_bp2{$name2_bp2}}[0]="P2";
            ${$read_bp2{$name2_bp2}}[2]+=$line_bp2[11];
          }
       }
     } elsif (exists($read_bp2{$name2_bp2}) && (($line_bp2[0] eq "F" && ${$read_bp2{$name2_bp2}}[0] eq "R") || ($line_bp2[0] eq "R" && ${$read_bp2{$name2_bp2}}[0] eq "F")) && $line_bp2[7] eq "="){
       if($line_bp2[3] eq ${$read_bp2{$name2_bp2}}[6] && $line_bp2[4] eq ${$read_bp2{$name2_bp2}}[11] && $line_bp2[8] eq ${$read_bp2{$name2_bp2}}[7]){
         splice @{$read_bp2{$name2_bp2}},3,0,@line_bp2;
         ${$read_bp2{$name2_bp2}}[0]="P2";
         ${$read_bp2{$name2_bp2}}[2]+=$line_bp2[11];
       } else { 
         if(${$read_bp2{$name2_bp2}}[6] ne $line_bp2[3]){
           $name_bp2=$line_bp2[1]."_".$line_bp2[3]."_".$line_bp2[4];
           @{$read_bp2{$name_bp2}}=@line_bp2;
           unshift(@{$read_bp2{$name_bp2}},$line_bp2[11]);
           unshift(@{$read_bp2{$name_bp2}},$line_bp2[3]);
           unshift(@{$read_bp2{$name_bp2}},$line_bp2[0]);
         } else {
           splice @{$read_bp2{$name2_bp2}},3,0,@line_bp2; ${$read_bp2{$name2_bp2}}[0]="P2";
           ${$read_bp2{$name2_bp2}}[2]+=$line_bp2[11];
         }
       }
     } else {
       if($order_bp2 eq 1){next;}
       $name_bp2=$line_bp2[1]."_".$line_bp2[3]."_".$line_bp2[4];
       @{$read_bp2{$name_bp2}}=@line_bp2;
       unshift(@{$read_bp2{$name_bp2}},$line_bp2[11]);
       unshift(@{$read_bp2{$name_bp2}},$line_bp2[3]);
       unshift(@{$read_bp2{$name_bp2}},$line_bp2[0]);
     }
   } else {
     my @order_bp2=sort{$read_bp2{$b}->[2]<=>$read_bp2{$a}->[2]} keys %read_bp2;
     for(my $j=0,my $max=${$read_bp2{$order_bp2[0]}}[2];$j<@order_bp2;$j++){
       ${$read_bp2{$order_bp2[$j]}}[1]="-";
       if(${$read_bp2{$order_bp2[$j]}}[2]<$max) {
         last;
       } else {
         print BP2 "@{$read_bp2{$order_bp2[$j]}}\n";
       }
     }
     %read_bp2=();
     $name_bp2 = $line_bp2[1]."_".$line_bp2[3]."_".$line_bp2[4];
     @{$read_bp2{$name_bp2}}=@line_bp2;
     $read_name_bp2=$line_bp2[1];
     unshift(@{$read_bp2{$name_bp2}},$line_bp2[11]);
     unshift(@{$read_bp2{$name_bp2}},$line_bp2[3]);
     unshift(@{$read_bp2{$name_bp2}},$line_bp2[0]);
   }
 }
 if(%read_bp2) {
   my @order_bp2=sort{$read_bp2{$b}->[2]<=>$read_bp2{$a}->[2]} keys %read_bp2;
   for(my $j=0,my $max=${$read_bp2{$order_bp2[0]}}[2];$j<@order_bp2;$j++){
     ${$read_bp2{$order_bp2[$j]}}[1]="-";
     if (${$read_bp2{$order_bp2[$j]}}[2]<$max) {
       last;
     } else {
       print BP2 "@{$read_bp2{$order_bp2[$j]}}\n";
     }
   }
 }
}

###### Genotyping function 1: filtering
sub Genotype_each_TE1 {
  my (@input) = @_;
  my @return1 = ();
  my $line = "";
  my @line = ();
  my $type1 = "";
  my $k = 0;
  my @len = ();
  my $as = "";
  my $md = "";
  my $posi = "";
  for (my $i = 0; $i < @input; $i++) {
    @line = @{$input[$i]};
    $as = "";
    $md = "";
    $posi = "";
    if($line[0] =~ "#@") {
      next;
    }
    @len = ($line[5]=~/(\d+)/g);
    my $len=0; map {$len+=$_} @len;
    $k=0;
    $type1="";
    if ($line[1]%8<4) { 
      if ($line[1]%4>=2 && $line[1]%256>=128) {
        $type1="PR";
        $k=1;
      } elsif ($line[1]%4>=2 && $line[1]%256<128) {
        $type1="PF";
        $k=1;
      } elsif ($line[1]%256<128) {
        $type1="F";
        $k=1;
      } elsif($line[1]%256>=128) {
        $type1="R";
        $k=1;
      } else {
        $type1="N";
        $k=1;
      }
      for (my $i=11;$i<@line;$i++) {
        if ($line[$i]=~"AS:i:") {
          $as=$line[$i];
          $as=~s/AS\:i\://;
        }
        if ($line[$i]=~"MD:Z:") {
          $md=$line[$i];
          $md=~s/MD\:Z\://; 
        }
        my $hstart = 0;
        my $hend = 0;
        my $hlength = 0;
        my @temp1 = ($line[5]=~/(\d+)/g);
        my @temp2 = ($line[5]=~/([A-Z])/g);
        for (my $j=0;$j<@temp1;$j++) {
          if ($temp2[$j] eq "S") {   
            if ($j==0) {
              $hstart=$temp1[$j]+1;
              $hend=$temp1[$j];
              $hlength=$temp1[$j];
            }
            if ($j==@temp1-1) {
              $hlength+=$temp1[$j];
            }
          } elsif ($temp2[$j] eq "D") {
            next;
          } elsif ($temp2[$j] eq "M" || $temp2[$j] eq "I") {
            if($j == 0) {
              $hstart = 0;
              $hend = $temp1[$j];
              $hlength = $temp1[$j];
            } elsif($j == @temp1-1) {
              $hlength += $temp1[$j];
              $hend += $temp1[$j];
            } else {
              $hend += $temp1[$j];
              $hlength += $temp1[$j];
            }
          }
        }
        $posi=$hstart."_".$hend."_".$hlength."_".$md;
      }
    } else {
      $type1="U";
      next;
    }
    unless ($as) {
      $as=0;
    }
    my @return1_tmp=();
    push (@return1_tmp,$type1);
    push (@return1_tmp,@line[0..8]);
    push (@return1_tmp,$posi);
    push (@return1_tmp,$as);
    push (@return1,[@return1_tmp]);
  }
  return (@return1);
}

##### Genotyping function 2: extract reads corss the breakpoint 
sub Genotype_each_TE2 {
  my (@input2) = @_;
  my $insert_size_geno = $min_insertsize;
  my $line="";
  my @line=();
  my %read=();
  my $read_name="";
  my $name="";
  my $order=0;
  my $order1=0;
  my %human=();
  my @return2=();
  my @tmp1 = (0) x 26;
  push (@return2,[@tmp1]);

  for (my $i = 0; $i < @input2; $i++) {
    @line = @{$input2[$i]};
    if($line[2] / 512>=1){
      next;
    }
    $line[5] = $line[3];
    my $name2=$line[1]."_".$line[3]."_".$line[8];
    $name=$line[1]."_".$line[3]."_".$line[4];
    if(!(%read)) {
      @{$read{$name}} = @line;
      $read_name = $line[1];
      unshift(@{$read{$name}},$line[11]);
      unshift(@{$read{$name}},$line[3]);
      unshift(@{$read{$name}},$line[0]);
    }
    if($line[1] eq $read_name) {
      if(exists($read{$name2}) && (($line[0] eq "PF" && ${$read{$name2}}[0] eq "PR") || ($line[0] eq "PR" && ${$read{$name2}}[0] eq "PF")) && $line[7] eq "="){
        if ($line[3] eq ${$read{$name2}}[6] && $line[4] eq ${$read{$name2}}[11] && $line[8] eq ${$read{$name2}}[7]){        ##### 02/14/2019
          splice @{$read{$name2}},3,0,@line;
          ${$read{$name2}}[0]="P";                                               
          ${$read{$name2}}[2]+=$line[11];
        } else {
          if((${$read{$name2}}[6] ne $line[3])){                                                                            ##### 02/14/2019
            $name=$line[1]."_".$line[3]."_".$line[4];
            @{$read{$name}}=@line;
            unshift(@{$read{$name}},$line[11]);
            unshift(@{$read{$name}},$line[3]);
            unshift(@{$read{$name}},$line[0]); 
          } else {
            splice @{$read{$name2}},3,0,@line;
            ${$read{$name2}}[0]="P2";
            ${$read{$name2}}[2]+=$line[11];
          }
        }
      } else {
        if($order eq 1){next;}
        $name=$line[1]."_".$line[3]."_".$line[4];
        @{$read{$name}}=@line;
        unshift(@{$read{$name}},$line[11]);
        unshift(@{$read{$name}},$line[3]);
        unshift(@{$read{$name}},$line[0]);
      }
    } else {
      my @order=sort{$read{$b}->[2]<=>$read{$a}->[2]} keys %read;
      for(my $j=0,my $max=${$read{$order[0]}}[2];$j<@order;$j++){
        ${$read{$order[$j]}}[1] = $insert_size_geno;
        if(${$read{$order[$j]}}[2]<$max){
          last;
        } else {
          push (@return2,[@{$read{$order[$j]}}]);
        }
      } 
      %read=();
      $name=$line[1]."_".$line[3]."_".$line[4];
      @{$read{$name}}=@line;
      $read_name=$line[1];
      unshift(@{$read{$name}},$line[11]);
      unshift(@{$read{$name}},$line[3]);
      unshift(@{$read{$name}},$line[0]);
    }
  }
  if(%read){
    my @order=sort{$read{$b}->[2]<=>$read{$a}->[2]} keys %read;
    for(my $j=0,my $max=${$read{$order[0]}}[2];$j<@order;$j++) {
      ${$read{$order[$j]}}[1] = $insert_size_geno;
      if(${$read{$order[$j]}}[2]<$max){
        last;
      } else {
        push (@return2,[@{$read{$order[$j]}}]);
      }
    }
  }
  return @return2;
}

##### Genotyping function 3: calculate the reads support non TE
sub Genotype_each_TE3_PE {
  my ($input3,$input4) = @_;
  my @input3 = @{$input3};
  my @input4 = @{$input4};
  my $line="";
  my @line=();
  my %VI_list=();
  my @names=();
  my %split=();
  my %split3=();
  my %split2=();
  ### each PE read
  for (my $i=0;$i < @input3; $i++){
    @line = @{$input3[$i]};
    ### filtering
    if ($line[0] ne "P" && abs($line[7]-$input4[2]) < ($length_insertsize+2*$std_insertsize)*2 && $line[6] eq $input4[1]) {
    } elsif ($line[0] eq "P" && $line[6] eq $input4[1]){
      if ($line[7]>=$line[19] && $input4[2]>=$line[19] && $input4[2]<=($line[7]+$input4[5])) {
   
      } elsif ($line[7]<=$line[19] && $input4[2]<=($line[19]+$input4[5]) && $input4[2]>=$line[7]) {
   
      } else {
        next;
      }
    } else {
      next;
    }
    ### merging PE
    my $name=$input4[0]."%".$input4[1]."%".$input4[2];                      ### corrected 2/15/2019
    unshift (@line,$name);
    if (exists($VI_list{$line[0]."|".$line[5]}) && ${$VI_list{$line[0]."|".$line[5]}}[0] eq $line[0] && ${$VI_list{$line[0]."|".$line[5]}}[7] eq $line[7] && ${$VI_list{$line[0]."|".$line[5]}}[8] eq $line[12] && ${$VI_list{$line[0]."|".$line[5]}}[12] eq $line[8]) {
      ${$VI_list{$line[0]."|".$line[5]}}[1]="P2";
      ${$VI_list{$line[0]."|".$line[5]}}[3]=${$VI_list{$line[0]."|".$line[5]}}[3]+$line[3];
      push(@{$VI_list{$line[0]."|".$line[5]}},@line[4..$#line]);
    } else {
      push(@{$VI_list{$line[0]."|".$line[5]}},@line);
    }
  }
  ### extract PE reads and count
  my %vi=();
  my @reads=keys %VI_list;
  for(my $i=0;$i<@reads;$i++){
    @line=@{$VI_list{$reads[$i]}};
    unless(exists($vi{$line[0]})) {
      @{$vi{$line[0]}}=(0) x 4;
    }
    @names=split /\%/, $line[0];
    if($line[1] ne "P" && $line[1] ne "P2"){
      next;
    }
    ### extract split reads
    my $start1=0;
    my $end=0;
    my $length=0;
    my $start_1=0;
    my $end_1=0;
    my $len_1=0;
    my $start_2=0;
    my $end_2=0;
    my $len_2=0;
    if($line[10] =~ "S" || $line[22] =~ "S"){           ##### split start
      ### first end
      my @temp1=($line[10]=~/(\d+)/g);
      my @temp2=($line[10]=~/([A-Z])/g);
      for(my $i=0;$i<@temp1;$i++) {
        if($temp2[$i] eq "S" || $temp2[$i] eq "H"){   
          if($i==0){
            $start1=$temp1[$i]+1;
            $end=$temp1[$i];
            $length=$temp1[$i];
          }
          if($i==@temp1-1) {
            $length+=$temp1[$i];
          }
        } elsif ($temp2[$i] eq "D") {
          next;
        } elsif ($temp2[$i] eq "M" || $temp2[$i] eq "I"){
          if($i==0) {
            $start1=1;
            $end=$temp1[$i];
            $length=$temp1[$i];
          } elsif($i==@temp1-1) {
            $length+=$temp1[$i];
            $end+=$temp1[$i];
          } else {
            $end+=$temp1[$i];
            $length+=$temp1[$i]
          }
        }
      }
      $start_1=$start1;
      $end_1=$end;
      $len_1=$length;
      ### second end
      @temp1=($line[22]=~/(\d+)/g);
      @temp2=($line[22]=~/([A-Z])/g);
      for(my $i=0;$i<@temp1;$i++){
        if($temp2[$i] eq "S" || $temp2[$i] eq "H"){     
          if($i==0){ 
            $start1=$temp1[$i]+1;
            $end=$temp1[$i];
            $length=$temp1[$i];
          }
          if($i==@temp1-1) {
            $length+=$temp1[$i];
          }
        } elsif ($temp2[$i] eq "D") {
          next;
        } elsif ($temp2[$i] eq "M" || $temp2[$i] eq "I"){
          if($i==0) {
            $start1=1;
            $end=$temp1[$i];
            $length=$temp1[$i];
          } elsif ($i==@temp1-1) {
            $length+=$temp1[$i];
            $end+=$temp1[$i];
          } else {
            $end+=$temp1[$i];
            $length+=$temp1[$i]
          }
        }
      }
      $start_2=$start1;
      $end_2=$end;
      $len_2=$length;
    }                                                        ##### split end
    ### identify split >=20, and split <20 bp (suggestive reads)
    if($names[2] eq $line[8] && $line[10] =~ "S"){
      if($start_1<20){
        $line[1]=$line[1]."_F_".$len_1."_".$start_1."_".$end_1;
        if(exists($split2{$line[0]})) {
          $split2{$line[0]}++;
        } else {
          $split2{$line[0]}=1;
        }
      } else {
        $line[1]=$line[1]."_F_".$len_1."_".$start_1."_".$end_1;
        if (exists($split3{$line[0]})) {
          $split3{$line[0]}++;
        } else {
          $split3{$line[0]}=1;
        }         
      }
      push(@{$split{$line[0]."|".$line[5]}},@line);
    } elsif ($names[2] eq $line[20] && $line[22] =~ "S") {
      if($start_2<20){
        $line[1]=$line[1]."_R_".$len_2."_".$start_2."_".$end_2;
        if(exists($split2{$line[0]})){
          $split2{$line[0]}++;
        } else {
          $split2{$line[0]}=1;
        }                
      } else {
        $line[1]=$line[1]."_R_".$len_2."_".$start_2."_".$end_2;
        if (exists($split3{$line[0]})) {
          $split3{$line[0]}++;
        } else {
          $split3{$line[0]}=1;
        }               
      }
      push(@{$split{$line[0]."|".$line[5]}},@line);
    } elsif (($names[2] - ($line[8]+$end_1))<=5 && $line[10] =~ "S") {
      if ($len_1-$end_1<20) {
        $line[1]=$line[1]."_F2_".$len_1."_".$start_1."_".$end_1;
        if (exists($split2{$line[0]})) {
          $split2{$line[0]}++;
        } else {
          $split2{$line[0]}=1;
        }                
      } else {
        $line[1]=$line[1]."_F2_".$len_1."_".$start_1."_".$end_1;
        if (exists($split3{$line[0]})) {
          $split3{$line[0]}++;
        } else {
          $split3{$line[0]}=1;
        }               
      }
      push(@{$split{$line[0]."|".$line[5]}},@line);
    } elsif (($names[2] - ($line[20]+$end_2))<=5 && $line[22] =~ "S") {
      if($len_2-$end_2<20){
        $line[1]=$line[1]."_R2_".$len_2."_".$start_2."_".$end_2;
        if (exists($split2{$line[0]})) {
          $split2{$line[0]}++;
        } else {
          $split2{$line[0]}=1;
        }                 
      } else {
        $line[1]=$line[1]."_R2_".$len_2."_".$start_2."_".$end_2;
        if (exists($split3{$line[0]})) {
          $split3{$line[0]}++;
        } else {
          $split3{$line[0]}=1;
        }                 
      }
      push(@{$split{$line[0]."|".$line[5]}},@line);
    }
    ### reads support nonTE
    if(exists($vi{$line[0]}) && !(exists($split{$line[0]."|".$line[5]}))){
      my @temp1=($line[10]=~/(\d+)/g);
      my @temp2=($line[10]=~/([A-Z])/g);
      my $hstart="";
      my $hend="";
      my $hlength="";
      for(my $i=0;$i<@temp1;$i++){
        if($temp2[$i] eq "S" || $temp2[$i] eq "H"){                          #### 02/15/2019   
          if($i==0) {
            $hstart=$temp1[$i]+1;
            $hend=$temp1[$i];
            $hlength=$temp1[$i];
          }
          if ($i==@temp1-1) {
            $hlength+=$temp1[$i];
          }
        } elsif ($temp2[$i] eq "D") {
          next;
        } elsif ($temp2[$i] eq "M" || $temp2[$i] eq "I") {
          if ($i==0) {
            $hstart=1;
            $hend=$temp1[$i];
            $hlength=$temp1[$i];
          } elsif ($i==$#temp1) {
            $hlength+=$temp1[$i];
            $hend+=$temp1[$i];
          } else {
            $hend+=$temp1[$i];
            $hlength+=$temp1[$i];
          }
        }
      }
      if (abs(($hend+$line[8]-1)-$input4[2])<=10 && $line[8]>=$line[20]){
        next;
      } elsif (abs(($hstart+$line[8]-1)-$input4[2])<=10 && $line[8]<=$line[20]) {
        next;
      }
      @temp1=($line[22]=~/(\d+)/g);
      @temp2=($line[22]=~/([A-Z])/g);
      $hstart="";
      $hend="";
      $hlength="";

      for (my $i=0;$i<@temp1;$i++) {
        if($temp2[$i] eq "S" || $temp2[$i] eq "H") {                ### 2/20/2019
          if($i==0) {
            $hstart=$temp1[$i]+1;
            $hend=$temp1[$i];
            $hlength=$temp1[$i];
          }
          if ($i==$#temp1) {
            $hlength+=$temp1[$i];
          }
        } elsif ($temp2[$i] eq "D") {
          next;
        } elsif ($temp2[$i] eq "M" || $temp2[$i] eq "I") {
          if ($i==0) {
            $hstart=1;
            $hend=$temp1[$i];
            $hlength=$temp1[$i];
          } elsif($i==@temp1-1) {
            $hlength+=$temp1[$i];
            $hend+=$temp1[$i];
          } else {
            $hend+=$temp1[$i];
            $hlength+=$temp1[$i];
          }
        }
      }
      if (abs(($hstart+$line[20]-1)-$input4[2])<=10 && $line[20]<$line[8]) {
        next;
      } elsif (abs(($hend+$line[20]-1)-$input4[2])<=10 && $line[8]<$line[20]){
        next;
      }
      ${$vi{$line[0]}}[0]++;
    }                                                                    ### 1. here counting number of reads support Non-VI;
  }
  ### number of chimeric and split reads and ratio of the TE;
  my $name=$input4[0]."%".$input4[1]."%".$input4[2];                     # sample_id, chr, position
  unless (exists($split2{$name})){
    $split2{$name}=0;
  }                     # if there is no split reads
  unless (exists($split3{$name})) {
    $split3{$name}=0;
  }                     # if there is no suggestive reads
  unless ($vi{$name}) {
    @{$vi{$name}}=(0) x 3;
    ${$vi{$name}}[2]=1;
  }                     # 0. if there is no reads support non-TE
  ### array structure
  ${$vi{$name}}[1]=$input4[3];                                     # 1. chimeric & split reads 
  ${$vi{$name}}[3]=$split3{$name};                                 # 3. suggestive reads (split reads <20 bp)
  if(${$vi{$name}}[0] eq 0){                                       # if there is no reads support nonTEs
    ${$vi{$name}}[2]=1;
  } elsif(${$vi{$name}}[0] >0 && ${$vi{$name}}[1]>0){
    if($input4[4] eq 2) {
      ${$vi{$name}}[2]=(${$vi{$name}}[1]/2)/((${$vi{$name}}[1]/2)+${$vi{$name}}[0]);
    } elsif ($input4[4] eq 1) {
      ${$vi{$name}}[2]=${$vi{$name}}[1]/(${$vi{$name}}[1]+${$vi{$name}}[0]);
    }
  }
  my @return3 = ();
  push (@return3,$name);
  push (@return3,@{$vi{$name}});
  return @return3;
}

##### Genotyping function 3 for single-end: calculate reads support nonTE
sub Genotype_each_TE3_SE {
  my ($input3,$input4) = @_;
  my @input3 = @{$input3};
  my @input4 = @{$input4};
  my $line="";
  my @line=();
  my %VI_list=();
  my @names=();
  my %split=();
  my %split3=();
  my %split2=();
  my $fully_mapped=0;
  my $split_reads=0;
  my $split_reads_short=0;

  for (my $i = 0; $i < @input3; $i++) {
    @line = @{$input3[$i]};
    ### filtering
    if ($line[0] ne "P" && abs($line[7]-$input4[2]) < ($length_insertsize+2*$std_insertsize)*2 && $line[6] eq $input4[1]) {

    } elsif ($line[0] eq "P" && $line[6] eq $input4[1]){
      if ($line[7]>=$line[19] && $input4[2]>=$line[19] && $input4[2]<=($line[7]+$input4[5])) {
     
      } elsif ($line[7]<=$line[19] && $input4[2]<=($line[19]+$input4[5]) && $input4[2]>=$line[7]) {
    
      } else {
        next;
      }
    } else {
      next;
    }
    my $name=$input4[0]."%".$input4[1]."%".$input4[2];
    unshift (@line,$name);
    push(@{$VI_list{$line[0]."|".$line[5]}},@line);
  }
  my %vi=();
  my @reads=keys %VI_list;
  for(my $i=0;$i<@reads;$i++){
    @line=@{$VI_list{$reads[$i]}};
    @names=split /\%/, $line[0];
    my $start1=0;
    my $end=0;
    my $length=0;
    my @temp1=($line[10]=~/(\d+)/g);
    my @temp2=($line[10]=~/([A-Z])/g);
    for(my $i=0;$i<@temp1;$i++){
      if ($temp2[$i] eq "S" || $temp2[$i] eq "H") {   
        if($i==0){$start1=$temp1[$i]+1;$end=$temp1[$i];$length=$temp1[$i];}
        if($i==@temp1-1){$length+=$temp1[$i];}
      } elsif($temp2[$i] eq "D") {
        next;
      } elsif($temp2[$i] eq "M" || $temp2[$i] eq "I"){
        if($i==0) {
          $start1=1;
          $end=$temp1[$i];
          $length=$temp1[$i];
        } elsif($i==@temp1-1) {
          $length+=$temp1[$i];
          $end+=$temp1[$i];
        } else {
          $end+=$temp1[$i];
          $length+=$temp1[$i];
        }
      }
    }
    if ($line[8] <= $input4[2] && ($line[8] + $length) >= $input4[2] && !($line[10] =~ "S")) {                       ### fully mapped with no split reads
      $fully_mapped++;
    } elsif ($start1 > 0 && abs($line[8] - $input4[2]) <= 5 && $line[10] =~ "S") {
      if ($start1>=20) {
        $split_reads++;
      } else {
        $split_reads_short++;
      }
    } elsif ($end < $length && $line[10] =~ "S" && abs($line[8] + $end - $start1 - $input4[2]) <= 5) {
      if ($length - $end >= 20) {
        $split_reads++;
      } else {
        $split_reads_short++;
      }
    } elsif ((abs($line[8] - $input4[2]) <=5 || abs(($line[8] + $end - $start1) - $input4[2]) <= 5) && $line[10] =~ "S") {
      $fully_mapped++;
    } else {
      next;
    }
  }
  ##### print out
  my $name=$input4[0]."%".$input4[1]."%".$input4[2];
  my @return4 = ();
  push (@return4,$name);
  push (@return4,$fully_mapped);
  push (@return4,$input4[3]);
  my $GR=0;
  if($input4[4] eq 2) {
    $GR = ($input4[3]/2)/(($input4[3]/2) + $fully_mapped);
  } elsif ($input4[4] eq 1) {
    $GR = $input4[3]/($input4[3] + $fully_mapped);
  }
  push (@return4,$GR);
  push (@return4,$split_reads_short);
  return @return4;
}

##### Clean tmp files during detection
sub Clean_tmp_files_detection {
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
  if(-e ${input_sampleID}."_all_breakpoint"){system ("rm ${input_sampleID}_all_breakpoint");}
}

##### Clean files
sub Clean_files {
  if (-e ${input_sampleID}."_ERV_vsu.sam2"){system ("rm ${input_sampleID}_ERV_vsu.sam2");}
  if (-e ${input_sampleID}."_ERV_vsu.sam"){system ("rm ${input_sampleID}_ERV_vsu.sam");}
  if (-e ${input_sampleID}."_ERV.sam"){system ("rm ${input_sampleID}_ERV.sam");}
  if (-e ${input_sampleID}."_ERV_breakpoint"){system ("rm ${input_sampleID}_ERV_breakpoint");}
  if (-e ${input_sampleID}."_1.1fuq"){system ("rm ${input_sampleID}_1.1fuq");}
  if (-e ${input_sampleID}."_2.1fuq"){system ("rm ${input_sampleID}_2.1fuq");}
  if (-e ${input_sampleID}."_1sf.fuq"){system ("rm ${input_sampleID}_1sf.fuq");}
  if (-e ${input_sampleID}."_1sf.othu"){system ("rm ${input_sampleID}_1sf.othu");}
  if (-e ${input_sampleID}."_h.sam"){system ("rm ${input_sampleID}_h.sam");}
  if (-e ${input_sampleID}."_ERV.fasta"){system ("rm ${input_sampleID}_ERV.fasta");}
  if (-e ${input_sampleID}."_ERV1_h.sam"){system ("rm ${input_sampleID}_ERV1_h.sam");}
  if (-e ${input_sampleID}."_ERV1_h.sam2"){system ("rm ${input_sampleID}_ERV1_h.sam2");}
  if (-e ${input_sampleID}."_ERV.sort.bam"){system ("rm ${input_sampleID}_ERV.sort.bam");}
  if (-e ${input_sampleID}."_1.1fq"){system ("rm ${input_sampleID}_1.1fq");}
  if (-e ${input_sampleID}."_2.1fq"){system ("rm ${input_sampleID}_2.1fq");}
  if (-e ${input_sampleID}."_1sf.fastq"){system ("rm ${input_sampleID}_1sf.fastq");}
  if (-e ${input_sampleID}."_1sf.others"){system ("rm ${input_sampleID}_1sf.others");}
  if (-e ${input_sampleID}."_ERV_1.1fq"){system ("rm ${input_sampleID}_ERV_1.1fq");}
  if (-e ${input_sampleID}."_ERV_1.2fq"){system ("rm ${input_sampleID}_ERV_1.2fq");}
  if (-e ${input_sampleID}."_ERV1_h.sam2"){system ("rm ${input_sampleID}_ERV1_h.sam2");}
  if (-e ${input_sampleID}."_ERV_2.1fq"){system ("rm ${input_sampleID}_ERV_2.1fq");}
  if (-e ${input_sampleID}."_ERV_2.2fq"){system ("rm ${input_sampleID}_ERV_2.2fq");}
  if (-e ${input_sampleID}."_ERV.3"){system ("rm ${input_sampleID}_ERV.3");}
  if (-e ${input_sampleID}."_ERV.bam"){system ("rm ${input_sampleID}_ERV.bam");}
  if (-e ${input_sampleID}."_ERV.error"){system ("rm ${input_sampleID}_ERV.error");}
  if (-e ${input_sampleID}."_h1_ERV.bam"){system ("rm ${input_sampleID}_h1_ERV.bam");}
  if (-e ${input_sampleID}."_h1_pe.bam"){system ("rm ${input_sampleID}_h1_pe.bam");}
  if (-e ${input_sampleID}."_h_ERV.bam"){system ("rm ${input_sampleID}_h_ERV.bam");}
  if (-e ${input_sampleID}."_vsu.sam"){system ("rm ${input_sampleID}_vsu.sam");}
  if (-e ${input_sampleID}."_h1_pe.bam"){system ("rm ${input_sampleID}_h1_pe.bam");}
  if (-e ${input_sampleID}."_pe.bam"){system ("rm ${input_sampleID}_pe.bam");}
  if (-e ${input_sampleID}."_sm.bam"){system ("rm ${input_sampleID}_sm.bam");}
  if (-e ${input_sampleID}."_sm.sam"){system ("rm ${input_sampleID}_sm.sam");}
  if (-e ${input_sampleID}."_su.bam"){system ("rm ${input_sampleID}_su.bam");}
  if (-e ${input_sampleID}.".tab"){system ("rm ${input_sampleID}.tab");}
  if (-e ${input_sampleID}.".tab2"){system ("rm ${input_sampleID}.tab2");}
  if (-e ${input_sampleID}."_vsoft.sam"){system ("rm ${input_sampleID}_vsoft.sam");}
  if (-e ${input_sampleID}."_h.soft.sam2"){system ("rm ${input_sampleID}_h.sort.sam2");}
  if (-e ${input_sampleID}."_h1_h1.2fq"){system ("rm ${input_sampleID}_h1_h1.2fq");}
  if (-e ${input_sampleID}.".1fq"){system ("rm ${input_sampleID}.1fq");}
  if (-e ${input_sampleID}."_h1_h.bam"){system ("rm ${input_sampleID}_h1_h.bam");}
  if (-e ${input_sampleID}."_h1_u.bam"){system ("rm ${input_sampleID}_h1_u.bam");}
  if (-e ${input_sampleID}."_u.bam"){system ("rm ${input_sampleID}_u.bam");}
  if (-e ${input_sampleID}."_h.sort.sam2"){system ("rm ${input_sampleID}_h.sort.sam2");}
  if (-d ${input_sampleID}."_subgenome"){system ("rm -rf ${input_sampleID}_subgenome");}
}

##### Move files
sub Move_files {
  if (-e ${input_sampleID}."_ERV_1.1fuq"){system ("mv ${input_sampleID}_ERV_1.1fuq ${input_sampleID}_temp/");}
  if (-e ${input_sampleID}."_ERV_2.1fuq"){system ("mv ${input_sampleID}_ERV_2.1fuq ${input_sampleID}_temp/");}
  if (-e ${input_sampleID}."_ERV_1sf.fuq"){system ("mv ${input_sampleID}_ERV_1sf.fuq ${input_sampleID}_temp/");}
  if (-e ${input_sampleID}."_ERV1.bian"){system ("mv ${input_sampleID}_ERV1.bian ${input_sampleID}_temp/");}
  if (-e ${input_sampleID}."_ERV.fine_mapped"){system ("mv ${input_sampleID}_ERV.fine_mapped ${input_sampleID}_temp/");}
  if (-e ${input_sampleID}."_ERV.hf"){system ("mv ${input_sampleID}_ERV.hf ${input_sampleID}_temp/");}
  if (-e ${input_sampleID}."_ERV.TE_f"){system ("mv ${input_sampleID}_ERV.TE_f ${input_sampleID}_temp/");}
  if (-e ${input_sampleID}."_ERV.TE_f2"){system ("mv ${input_sampleID}_ERV.TE_f2 ${input_sampleID}_temp/");}
  if (-e ${input_sampleID}."_ERV.output"){system ("mv ${input_sampleID}_ERV.output ${input_sampleID}_temp/");}
  if (-e ${input_sampleID}."_ERV.output2"){system ("mv ${input_sampleID}_ERV.output2 ${input_sampleID}_temp/");}
  if (-e ${input_sampleID}."_ERV.output3"){system ("mv ${input_sampleID}_ERV.output3 ${input_sampleID}_temp/");}
  if (-e ${input_sampleID}."_ERV.output2.1"){system ("mv ${input_sampleID}_ERV.output2.1 ${input_sampleID}_temp/");}
  if (-e ${input_sampleID}."_ERV.visualization"){system ("mv ${input_sampleID}_ERV.visualization ${input_sampleID}_temp/");}
  if (-e ${input_sampleID}.".type") {
    system ("gzip ${input_sampleID}.type");
    system ("mv ${input_sampleID}.type.gz ${input_sampleID}_temp/");
  }
  if (-e ${input_sampleID}."_ERV.type") {
    system ("gzip ${input_sampleID}_ERV.type");
    system ("mv ${input_sampleID}_ERV.type.gz ${input_sampleID}_temp/");
  }
}

##### print out Arguments;
sub prtUsa {
  print "\n=====================================\nperl $0 [arguments]\n";
  print "\n";
  print "# Arguments:\n";
  print "       -i|input_sampleID <STR>			Sample ID (required)\n";
  print "       -f|file_suffix <STR>			The suffix of the input data, including: zipped FASTQ file (i.e., .fq.gz, and fastq.gz),\n\t\t\t\t\t\tunzipped FASTQ file (i.e., .fq, and fastq),\n\t\t\t\t\t\tBAM file (.bam), and a bam file list (.list; with \"-multiple_BAM\") (required). Default: .bam\n";
  print "       -H|Human_reference_genome <STR>		The FASTA file of the human reference genome (required)\n";
  print "       -T|TE_reference_genomes <STR>		The TE library (FASTA) used for screening (required)\n";
  print "       -I|Input_directory <STR>			The directory of input data. Default: Not specified (current working directory)\n";
  print "       -O|Output_directory <STR>		The directory for output data. Default: Not specified (current working directory)\n";
  print "       -n|number_of_reads <INT>			The minimum number of reads support a TE insertion. Default: 3\n";
  print "       -d|data_type <STR>			Data type, including: WGS, RNA-seq. Default: WGS\n";
  print "       -s|sequencing_type <STR>			Type of sequencing data, including: paired-end, single-end. Default: paired-end\n";
  print "       -l|length_insertsize <FLOAT>		Insert size length (bp). It will be estimated if it is not specified\n";
  print "       -L|L_std_insertsize <FLOAT>		Standard deviation of insert size length (bp). It will be estimated if it is not specified\n";
  print "       -r|read_len <INT>			Read length (bp), including: 100, 150, and 250 bp. Default: 100\n";
  print "       -t|threads <INT>				The number of threads will be used. Default: 1\n";
  print "       -S|Split <INT>				The minimum length for split reads. A longer length is recommended with longer read length. Default: 20\n";
  print "       -m|multiple_BAM				If multiple BAM files are used as the input (input bam file need to be indexed). Default: not specified\n";
  print "       -B|BWA_MEM				If the bam file is generated using aligner BWA_MEM. Default: Not specified\n";
  print "       -G|Genotype				Genotyping function (input bam file need to be indexed). Default: not specified\n";
  print "       -h|help					Print this help\n";
  print "\n\n";
  print "=====================================\nExamples for detecting ERV and other TE insertions:\n";
  print "# Detecting TE insertions with a BAM file as the input\n";
  print "       perl $0 -i TE_seq -f .bam -H hg38.fa -T TE_consensus.fa -I folder_of_input_data -O folder_for_output_files -t 12 -S 20 -BWA_MEM\n\n";
  print "# Detecting TE insertions with paired-end FASTQ file as the input\n";
  print "       perl $0 -i TE_seq -f .fq.gz -H hg38.fa -T TE_consensus.fa -I folder_of_input_data -O folder_for_output_files -t 12 -S 20\n\n";
  print "# Detecting TE insertions with separated BAM file(s) as the input\n";
  print "       perl $0 -i TE_seq -f .list -H hg38.fa -T TE_consensus.fa -I folder_of_input_data -O folder_for_output_files -t 12 -S 20 -BWA_MEM -m\n\n";
  print "# Detecting and genotyping TE insertions with a BAM file as the input\n";
  print "       perl $0 -i TE_seq -f .bam -H hg38.fa -T TE_consensus.fa -I folder_of_input_data -O folder_for_output_files -t 12 -S 20 -BWA_MEM -G\n\n\n";
}

##### print out error info
sub prtErr{
  print "\n=====================================\n@_\n\n";
}


#######################
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
