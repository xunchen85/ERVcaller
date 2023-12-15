#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
my %file=();
my $temp2="";
my @position="";
my @virus_p="";
my $temp1="";
my $n=0;

my $input="test_summary";
my $output="test.integration";
GetOptions(
           'input=s'=>\$input,
           'output=s'=>\$output
           );

open FILE,"$input";
open C1, "+>${input}_1";
open C2, "+>${input}_2";
open C3, "+>${input}_3";
open C4, "+>${input}_4";
open C5, "+>${input}_5";
open C6, "+>${input}_6";
open C7, "+>${input}_7";
open C8, "+>${input}_8";
open C9, "+>${input}_9";
open C10, "+>${input}_10";
open C11, "+>${input}_11";
open C12, "+>${input}_12";
open C13, "+>${input}_13";
open C14, "+>${input}_14";
open C15, "+>${input}_15";
open C16, "+>${input}_16";
open C17, "+>${input}_17";
open C18, "+>${input}_18";
open C19, "+>${input}_19";
open C20, "+>${input}_20";
open C21, "+>${input}_21";
open C22, "+>${input}_22";
open C23, "+>${input}_23";
open C24, "+>${input}_24";
open C25, "+>${input}_25";
open C26, "+>${input}_26";
open CX, "+>${input}_X";
open CY, "+>${input}_Y";
open CMT, "+>${input}_MT";

#################
while(<FILE>){
  @virus_p=split;
#     $virus_p[4]=~s/Chr//;
#     $virus_p[4]=~s/chr//;
     if($virus_p[4] eq "chr1" || $virus_p[4] eq "Chr1" || $virus_p[4] eq "1" || $virus_p[4] eq "NC_000001.10") {print C1 "@virus_p\n";}
     elsif($virus_p[4] eq "chr2" || $virus_p[4] eq "Chr2" || $virus_p[4] eq "2" || $virus_p[4] eq "NC_000002.11") {print C2 "@virus_p\n";}
     elsif($virus_p[4] eq "chr3" || $virus_p[4] eq "Chr3" || $virus_p[4] eq "3" || $virus_p[4] eq "NC_000003.11") {print C3 "@virus_p\n";}
     elsif($virus_p[4] eq "chr4" || $virus_p[4] eq "Chr4" || $virus_p[4] eq "4" || $virus_p[4] eq "NC_000004.11") {print C4 "@virus_p\n";}
     elsif($virus_p[4] eq "chr5" || $virus_p[4] eq "Chr5" || $virus_p[4] eq "5" || $virus_p[4] eq "NC_000005.9") {print C5 "@virus_p\n";}
     elsif($virus_p[4] eq "chr6" || $virus_p[4] eq "Chr6" || $virus_p[4] eq "6" || $virus_p[4] eq "NC_000006.11") {print C6 "@virus_p\n";}
     elsif($virus_p[4] eq "chr7" || $virus_p[4] eq "Chr7" || $virus_p[4] eq "7" || $virus_p[4] eq "NC_000007.13") {print C7 "@virus_p\n";}
     elsif($virus_p[4] eq "chr8" || $virus_p[4] eq "Chr8" || $virus_p[4] eq "8" || $virus_p[4] eq "NC_000008.10") {print C8 "@virus_p\n";}
     elsif($virus_p[4] eq "chr9" || $virus_p[4] eq "Chr9" || $virus_p[4] eq "9" || $virus_p[4] eq "NC_000009.11") {print C9 "@virus_p\n";}
     elsif($virus_p[4] eq "chr10" || $virus_p[4] eq "Chr10" || $virus_p[4] eq "10" || $virus_p[4] eq "NC_000010.10") {print C10 "@virus_p\n";}
     elsif($virus_p[4] eq "chr11" || $virus_p[4] eq "Chr11" || $virus_p[4] eq "11" || $virus_p[4] eq "NC_000011.9") {print C11 "@virus_p\n";}
     elsif($virus_p[4] eq "chr12" || $virus_p[4] eq "Chr12" || $virus_p[4] eq "12" || $virus_p[4] eq "NC_000012.11") {print C12 "@virus_p\n";}
     elsif($virus_p[4] eq "chr13" || $virus_p[4] eq "Chr13" || $virus_p[4] eq "13" || $virus_p[4] eq "NC_000013.10") {print C13 "@virus_p\n";}
     elsif($virus_p[4] eq "chr14" || $virus_p[4] eq "Chr14" || $virus_p[4] eq "14" || $virus_p[4] eq "NC_000014.8") {print C14 "@virus_p\n";}
     elsif($virus_p[4] eq "chr15" || $virus_p[4] eq "Chr15" || $virus_p[4] eq "15" || $virus_p[4] eq "NC_000015.9") {print C15 "@virus_p\n";}
     elsif($virus_p[4] eq "chr16" || $virus_p[4] eq "Chr16" || $virus_p[4] eq "16" || $virus_p[4] eq "NC_000016.9") {print C16 "@virus_p\n";}
     elsif($virus_p[4] eq "chr17" || $virus_p[4] eq "Chr17" || $virus_p[4] eq "17" || $virus_p[4] eq "NC_000017.10") {print C17 "@virus_p\n";}
     elsif($virus_p[4] eq "chr18" || $virus_p[4] eq "Chr18" || $virus_p[4] eq "18" || $virus_p[4] eq "NC_000018.9") {print C18 "@virus_p\n";}
     elsif($virus_p[4] eq "chr19" || $virus_p[4] eq "Chr19" || $virus_p[4] eq "19" || $virus_p[4] eq "NC_000019.9") {print C19 "@virus_p\n";}
     elsif($virus_p[4] eq "chr20" || $virus_p[4] eq "Chr20" || $virus_p[4] eq "20" || $virus_p[4] eq "NC_000020.10") {print C20 "@virus_p\n";}
     elsif($virus_p[4] eq "chr21" || $virus_p[4] eq "Chr21" || $virus_p[4] eq "21" || $virus_p[4] eq "NC_000021.8") {print C21 "@virus_p\n";}
     elsif($virus_p[4] eq "chr22" || $virus_p[4] eq "Chr22" || $virus_p[4] eq "22" || $virus_p[4] eq "NC_000022.10") {print C22 "@virus_p\n";}
     elsif($virus_p[4] eq "chr23" || $virus_p[4] eq "Chr23" || $virus_p[4] eq "23") {print C23 "@virus_p\n";}
     elsif($virus_p[4] eq "chr24" || $virus_p[4] eq "Chr24" || $virus_p[4] eq "24") {print C24 "@virus_p\n";}
     elsif($virus_p[4] eq "chr25" || $virus_p[4] eq "Chr25" || $virus_p[4] eq "25") {print C25 "@virus_p\n";}
     elsif($virus_p[4] eq "chr26" || $virus_p[4] eq "Chr26" || $virus_p[4] eq "26") {print C26 "@virus_p\n";}
     elsif($virus_p[4] eq "chrX" || $virus_p[4] eq "ChrX" || $virus_p[4] eq "X" || $virus_p[4] eq "NC_000023.10")  {print CX "@virus_p\n";}
     elsif($virus_p[4] eq "chrY" || $virus_p[4] eq "ChrY" || $virus_p[4] eq "Y" || $virus_p[4] eq "NC_000024.9")  {print CY "@virus_p\n";}
     elsif($virus_p[4] eq "chrMT" || $virus_p[4] eq "ChrMT" || $virus_p[4] eq "MT" || $virus_p[4] eq "NC_012920.1") {print CMT "@virus_p\n";}
             }
 open OUT, ">$output";
   seek(C1,0,0);
   my $k=0; 
   my %can=();
   while(<C1>){my @line=split;@{$can{$k}}=@line;$k++;}
   my @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C2,0,0);
   $k=0; 
   %can=();
   while(<C2>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C3,0,0);
   $k=0; 
   %can=();
   while(<C3>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C4,0,0);
   $k=0; 
   %can=();
   while(<C4>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C5,0,0);
   $k=0; 
   %can=();
   while(<C5>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C6,0,0);
   $k=0; 
   %can=();
   while(<C6>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C7,0,0);
   $k=0; 
   %can=();
   while(<C7>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C8,0,0);
   $k=0; 
   %can=();
   while(<C8>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C9,0,0);
   $k=0; 
   %can=();
   while(<C9>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C10,0,0);
   $k=0; 
   %can=();
   while(<C10>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C11,0,0);
   $k=0; 
   %can=();
   while(<C11>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C12,0,0);
   $k=0; 
   %can=();
   while(<C12>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

  seek(C13,0,0);
   $k=0;
   %can=();
   while(<C13>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C14,0,0);
   $k=0;
   %can=();
   while(<C14>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C15,0,0);
   $k=0;
   %can=();
   while(<C15>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C16,0,0);
   $k=0;
   %can=();
   while(<C16>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C17,0,0);
   $k=0;
   %can=();
   while(<C17>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C18,0,0);
   $k=0;
   %can=();
   while(<C18>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C19,0,0);
   $k=0;
   %can=();
   while(<C19>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C20,0,0);
   $k=0;
   %can=();
   while(<C20>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}
 
   seek(C21,0,0);
   $k=0;
   %can=();
   while(<C21>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C22,0,0);
   $k=0;
   %can=();
   while(<C22>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C23,0,0);
   $k=0;
   %can=();
   while(<C23>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(C24,0,0);
   $k=0;
   %can=();
   while(<C24>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}
   
   seek(C25,0,0);
   $k=0;
   %can=();
   while(<C25>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}
   
   seek(C26,0,0);
   $k=0;
   %can=();
   while(<C26>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}
   
   seek(CX,0,0);
   $k=0;
   %can=();
   while(<CX>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(CY,0,0);
   $k=0;
   %can=();
   while(<CY>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   seek(CMT,0,0);
   $k=0;
   %can=();
   while(<CMT>){my @line=split;@{$can{$k}}=@line;$k++;}
   @order=sort{$can{$a}->[5]<=>$can{$b}->[5]} keys %can;
   for (my $i=0;$i<@order;$i++){print OUT "@{$can{$order[$i]}}\n";}

   `rm ${input}_1* ${input}_2* ${input}_3 ${input}_4 ${input}_5 ${input}_6 ${input}_7 ${input}_8 ${input}_9 ${input}_X ${input}_Y ${input}_MT`;
