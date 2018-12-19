#!usr/bin/perl -w
# Author: Xun Chen
# Email: Xun.Chen@uvm.edu

use strict;
use warnings;
use Getopt::Long;

my $file="test.sam";
my $output="test.sam2";
my $alignment_score=22;
my $bian;
my $hf="";

GetOptions('file=s'=>\$file,
           'output=s'=>\$output,
           'alignment_score=i'=>\$alignment_score,
	   'bian=s'=>\$bian,
           'hf=s'=>\$hf                             ### homolog between human and ERVs
          );

open FILE, "$file";
open OUTPUT, ">$output";
open FILTER, ">$hf";

my @line=();
my @name=();
my $md="";
my $as="";

#####################
my %bian=();
if (defined($bian)) {
 open BIAN, "$bian";
 while(<BIAN>){
  @line=split;
  @{$bian{$line[0]}}=@line;
 }
}
#####################

my $read_direction="NA"; 
while(<FILE>){
 @line=split;
 $as="NA";$md="NA";my $xs="NA";
 unless($line[2] eq "*"){
  for(my $i=11;$i<@line;$i++){
   if($line[$i]=~s/MD:Z://){$md=$line[$i];}
   if($line[$i]=~s/AS:i://){$as=$line[$i];}
   if($line[$i]=~s/XS:i://){$xs=$line[$i];}
                             }
   if($line[1] =~ "[a-zA-Z]" || $line[1]>2048){next;}
   if($line[1]%256>=128){$read_direction="R";}
   else{$read_direction="L";}
   if($as =~ "[a-zA-Z]" || $as<$alignment_score){next;}
   print OUTPUT "$read_direction unknown @line[0..8] $md $as $xs\n";
 #### filter
 if (defined($bian)) {
 if (exists($bian{$line[0]})){
  if ($read_direction ne ${$bian{$line[0]}}[1] && ${$bian{$line[0]}}[2]-$as<=(${$bian{$line[0]}}[2]/2)){     ### updated to half of AS;
   if (${$bian{$line[0]}}[4] eq $line[5]) {
    print FILTER "@{$bian{$line[0]}} False\n";
   } else {
    print FILTER "@{$bian{$line[0]}} Possible\n";
          }
                                                                                     }
                             }
                     }
                       }
           }
