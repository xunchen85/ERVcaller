#!/usr/bin/env perl

use strict;

my $line="";
my @line=();
my %redundancy=();
my %type=();
my %integration=();
my $sample=$ARGV[0];
my $directory=$ARGV[1];
my $fasta_title="";
my $directory3=$ARGV[2];

#################
open INTEGRATION,"$ARGV[0]_ERV_breakpoint";            #### ERV.integration file
while(<INTEGRATION>){
 @line=split;
 if(exists($integration{$line[2]})){next;}
 else{$integration{$line[2]}=$line[0];}
                    }
close INTEGRATION;

###############################
open F1,"$ARGV[0]_ERV_1.1fq";
open R1,"$ARGV[0]_ERV_2.1fq";
open F2,">$ARGV[0]_ERV_1.2fq";
open R2,">$ARGV[0]_ERV_2.2fq";

my $f1="";
my $title="";
my $r1="";
#############

open FASTA,">$ARGV[0]_ERV.fasta";
while($f1=<F1>){
 $title=$f1;
 chomp($title);
 $fasta_title=$f1;
 $fasta_title=~s/^@/>/;
 $title=~s/^\@//;
 $title=~s/\/1$//;
 if(exists($integration{$title}) && !(exists($redundancy{$title}))){    # remove redundancy reads and only keep integration reads;
  print F2 "$f1";
  $f1=<F1>;print F2 "$f1";
  print FASTA "$fasta_title$f1";
  $f1=<F1>;print F2 "$f1";
  $f1=<F1>;print F2 "$f1";

  $r1=<R1>;print R2 "$r1";
  $fasta_title=$r1;
  $fasta_title=~s/^@/>/;
  $r1=<R1>;print R2 "$r1";
  print FASTA "$fasta_title$r1";
  $r1=<R1>;print R2 "$r1";
  $r1=<R1>;print R2 "$r1";
                                                                   }
 else{
#  print "d:$f1";
  $f1=<F1>;$f1=<F1>;$f1=<F1>;
  $r1=<R1>;$r1=<R1>;$r1=<R1>;$r1=<R1>;
     }
               }
close F1;
close F2;
close R1;
close R2;
################################################ remove softclip reads;
#####
##### This part is overlap with virome-wide
#####
