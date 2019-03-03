#!use/bin/perl
use strict;
my $line="";
my @line=();
my $title="";
my %reads=();
my %redundancy=();

################### remove redundant single-end reads;
open REDUNDANCY,"$ARGV[0].redundancy";
while(<REDUNDANCY>){$line=$_;chomp($line);$redundancy{$line}=$line;}

################### update this file;
open TWO,"$ARGV[0].21";
open THREE,">$ARGV[0].2";
open ERROR, ">$ARGV[0].error";
while(<TWO>){
 @line=split;
 if(exists($redundancy{$line[2]})){print ERROR "$_";next;}        # exclude redundant chimeric reads
 if($line[0] =~ "b:S"){                                           # exclude redundant split reads
  $title=$line[4]."|".$line[5]."|".$line[7]."|".$line[9]."|".$line[1]."|".$line[13]."|".$line[14]."|".$line[16];
  if(exists($reads{$title})){print ERROR "$_";}
  else{print THREE "$_";$reads{$title}=0;}
                      }
 else{print THREE "$_";}
            }
close TWO;
close THREE;
close FASTA2;
