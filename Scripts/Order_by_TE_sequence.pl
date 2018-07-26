#!usr/bin/perl -w
use strict;
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
open CX, "+>${input}_X";
open CY, "+>${input}_Y";
open CMT, "+>${input}_MT";

#################
while(<FILE>){
  @virus_p=split;
     $virus_p[4]=~s/Chr//;
     $virus_p[4]=~s/chr//;
     if($virus_p[4] eq "1")    {print C1 "@virus_p\n";}
     elsif($virus_p[4] eq "2") {print C2 "@virus_p\n";}
     elsif($virus_p[4] eq "3") {print C3 "@virus_p\n";}
     elsif($virus_p[4] eq "4") {print C4 "@virus_p\n";}
     elsif($virus_p[4] eq "5") {print C5 "@virus_p\n";}
     elsif($virus_p[4] eq "6") {print C6 "@virus_p\n";}
     elsif($virus_p[4] eq "7") {print C7 "@virus_p\n";}
     elsif($virus_p[4] eq "8") {print C8 "@virus_p\n";}
     elsif($virus_p[4] eq "9") {print C9 "@virus_p\n";}
     elsif($virus_p[4] eq "10") {print C10 "@virus_p\n";}
     elsif($virus_p[4] eq "11") {print C11 "@virus_p\n";}
     elsif($virus_p[4] eq "12") {print C12 "@virus_p\n";}
     elsif($virus_p[4] eq "13") {print C13 "@virus_p\n";}
     elsif($virus_p[4] eq "14") {print C14 "@virus_p\n";}
     elsif($virus_p[4] eq "15") {print C15 "@virus_p\n";}
     elsif($virus_p[4] eq "16") {print C16 "@virus_p\n";}
     elsif($virus_p[4] eq "17") {print C17 "@virus_p\n";}
     elsif($virus_p[4] eq "18") {print C18 "@virus_p\n";}
     elsif($virus_p[4] eq "19") {print C19 "@virus_p\n";}
     elsif($virus_p[4] eq "20") {print C20 "@virus_p\n";}
     elsif($virus_p[4] eq "21") {print C21 "@virus_p\n";}
     elsif($virus_p[4] eq "22") {print C22 "@virus_p\n";}
     elsif($virus_p[4] eq "X")  {print CX "@virus_p\n";}
     elsif($virus_p[4] eq "Y")  {print CY "@virus_p\n";}
     elsif($virus_p[4] eq "MT") {print CMT "@virus_p\n";}
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
