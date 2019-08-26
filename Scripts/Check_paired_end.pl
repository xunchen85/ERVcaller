#!/usr/bin/env perl

# This script is to find another end of reads to form PE reads;
# Usage: perl Fastq_PE_format.pl -sample sample
#author: Xun Chen
#email: Xun.Chen@uvm.edu

use strict;
use warnings;
use Getopt::Long;

my $sample="seq";
my $file=".1fq";
GetOptions('sample=s'=>\$sample,
           'file=s'=>\$file
          );

######################### pe files
open PF1, "${sample}_1${file}";
open PR1, "${sample}_2${file}";
open PF2, ">${sample}_1${file}2";
open PR2, ">${sample}_2${file}2";
my %pe2=();
my $title="";
my $line="";

while(<PR1>){
 $title=$_;$title=~s/\/2$//;
 $pe2{$title}=$_;
 $line=<PR1>;$pe2{$title}.=$line;
 $line=<PR1>;$pe2{$title}.=$line;
 $line=<PR1>;$pe2{$title}.=$line;
            }
while(<PF1>){
 $title=$_;
 $title=~s/\/1$//;
 if(exists($pe2{$title})){
  print PR2 "$pe2{$title}";
  print PF2 "$_";
  $line=<PF1>;
  print PF2 "$line";
  $line=<PF1>;
  print PF2 "$line";
  $line=<PF1>;
  print PF2 "$line";
 } else {
  $line=<PF1>;
  $line=<PF1>;
  $line=<PF1>;
 }
            }
