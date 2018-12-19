#!usr/bin/perl
use strict;

my %list=();
my @list=();
my @line=();
my $line="";
my $title="";
my $sample=$ARGV[0];
#$sample=~s/_a21//;
my $title2="";
my $ok=0;

##########################
open F1,"$ARGV[1]";
while(<F1>){
 @line=split;
 $line[0]=~s/ad://;
 $title=$line[0]."_".$line[15]."_".$line[35]."_".$line[30]."_".$line[31]."_".$line[32]."_".$line[33];
 $list{$title}=@line;
}

##########################
open F2,"$ARGV[0]";
my @title_O1=();
while(<F2>){
 @line=split;
  if($line[0] eq "O1"){
   $ok=0;
   $title2=$line[2]."_".$line[1]."_".$line[22]."_".$line[17]."_".$line[18]."_".$line[19]."_".$line[20];
   if(exists($list{$title2})){$ok=1;@title_O1=();@title_O1=@line;print "@line\n";}
                      }
   elsif($ok==1){print "$_";}
           }
