#!usr/bin/perl
# Author:       Xun Chen
# Email:        Xun.Chen@uvm.edu
# Date:         01/07/19
#
use strict;
my $line = "";
my @line = ();
my @info = ();

open VCF,"$ARGV[0]";
open FAL,">$ARGV[0].false";
my $min_depth = $ARGV[1];      		##### threadhold for the number of supporting reads; 3 reads by default;
my $ratio = $ARGV[2];			##### threadhold of maximum ratio for filtering entire loci; 0.2 by default;

while (my $line = <VCF>) {
  @line = split /\s+/, $line;
  if ($line[0] =~ "#") {     		##### header information;
    print "$line";
  } else {                   		##### Each TE insertion;
    @info = split /\;/, $line[7];
    ##### depth
    my $CR = $info[2];
    $CR =~ s/CR=//;
    my $SR = $info[3];
    $SR =~ s/SR=//;
    my $total_depth = 0;
    ##### max ratio
    $info[5] =~ s/GR=//;
    my @ratio = split /\,/, $info[5];
    my $max_ratio = 0;
    for (my $i = 0; $i < @ratio; $i++) {
      if ($ratio[$i] > $max_ratio) {
	$max_ratio = $ratio[$i];
      }
    }
    my $max_depth = 0;
    my $count_MEI = 0;
    my $max_TE_reads = 0;
    my $max_GQ = 0;
    ###### each sample
    for (my $j = 9; $j < @line; $j++){
      my @tmp1 = split /\:/, $line[$j];
      if ($tmp1[1] > $max_GQ) {
        $max_GQ = $tmp1[1];
      }
      if ($tmp1[4] > $max_TE_reads) {
        $max_TE_reads = $tmp1[4];   
      }
      if ($tmp1[0] eq "0/0" || $tmp1[0] eq "./." || $tmp1[0] eq ".") {
        next;
      }

      if ($tmp1[3] + $tmp1[4] >= $max_depth) {
        $max_depth = $tmp1[3] + $tmp1[4];
      }
      $count_MEI++; 
      $total_depth = $total_depth + $tmp1[3] + $tmp1[4];
    }
    my $average_depth =0;
    if ($count_MEI > 0) {
      $average_depth = $total_depth/$count_MEI;
      $average_depth = sprintf("%.2f",$average_depth);
    } else {
      print FAL join("\t",@line),"\n";
      next;
    }
    if ($max_depth < $min_depth) {
      print FAL join("\t",@line),"\n";
      next;
    } else {
      if ($line[6] eq ".") {
        $line[6] = "maxdepth=$max_depth";
      } else {
        $line[6] = $line[6].";maxdepth=$max_depth";
      }
    }
    if ($average_depth < $min_depth) {
      print FAL join("\t",@line),"\n";
      next;
    } else {
      $line[6] = $line[6].",averagedepth=".$average_depth;
    }
    if ($max_ratio < $ratio) {
      print FAL join("\t",@line),"\n";
      next;
    } else {
      $line[6] = $line[6].",maxGR=".$max_ratio;
    }
    if ($max_TE_reads < $ARGV[3]) {
      print FAL join("\t",@line),"\n";
      next;
    }
    if ($max_GQ < $ARGV[4]) {
      print FAL join("\t",@line),"\n";
      next;
    }
    print join("\t",@line),"\n";
  }
}
