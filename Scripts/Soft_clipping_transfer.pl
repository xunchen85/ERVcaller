#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $file="test_vsoft.sort.sam";
my $output="test_vsoft_breakpoint";

GetOptions('file=s'=>\$file,
           'output=s'=>\$output
          );

open FILE, "$file";
open OUTPUT, ">$output";

my @line=();
my @name=();
my $md="";
my $as="";

#####################
while(<FILE>){
 @line=split;
 @name=split /\|/, $line[0];
 $as="NA";$md="NA";
 unless($line[2] eq "*"){
  for(my $i=11;$i<@line;$i++){
   if($line[$i]=~s/MD:Z://){$md=$line[$i];}
   if($line[$i]=~s/AS:i://){$as=$line[$i];}
                             }
   if($as eq "NA" || $as<20){next;}
    print OUTPUT "S1 unknown @name[1..11] @line[1..8] $md $as\n";
                        }
              }
