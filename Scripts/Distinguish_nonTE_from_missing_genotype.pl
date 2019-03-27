#!/usr/bin/env perl
#
# Author:       Xun Chen
# Email:        Xun.Chen@uvm.edu
# Date:         02/10/19
# 
# Usage: perl Distinguish_nonTE_from_missing_genotype.pl -n nonTE_allsamples -v Combined.VCF -o Combined-final.VCF
#

use Getopt::Long qw(:config no_ignore_case);
use strict;
use warnings;
use Cwd qw();

my $line = "";
my @line = ();
my %nonTE = ();
my @sampleID = ();

my $nonTE;
my $vcf;
my $output;

GetOptions(
  "n|nonTE=s" => \$nonTE,
  "v|vcf=s" =>\$vcf,
  "o|output=s" =>\$output
);


open NON_TE, "$nonTE";          ### Reads among nonTE locations

open VCF, "$vcf";             ### combined VCF file

open OUT, ">$output";


##### load the nonTE locations into hash table
while (my $line = <NON_TE>) {
  @line = split /\s+/, $line;
  $nonTE{$line[0]} = $line[1];
}

##### read VCF
while (my $line = <VCF>) {
  @line = split /\s+/, $line;
  if ($line =~ "^##") {
    print OUT "$line";
  } elsif ($line =~ "#CHROM") {
    @sampleID = @line;
    print OUT join("\t",@line),"\n";
  } else {
    for (my $i = 9; $i < @line; $i++) {
      my $id = $sampleID[$i]."%".$line[0]."%".$line[1];
      if ($line[$i] =~ "0/0") {
        if (exists ($nonTE{$id})) {
          $line[$i] = $nonTE{$id};
        } else {
          $line[$i] = ".";
        }
      } else {
      }
    }
    print OUT join("\t",@line),"\n";
  }
}
