#!/usr/bin/env perl

use strict;

my $line="";
my @line=();
my $dir="NA";
my %read=();

open F1,"$ARGV[0]";
my $AS_XS = $ARGV[1];

while(<F1>){
  @line=split;
  my $dir=$line[0];
  if(exists($read{$line[2]})){
    if($dir eq "L") {
      unshift(@{$read{$line[2]}},@line[12..13,7]);
    } else {
      push(@{$read{$line[2]}},@line[12..13,7]);
    }
  } else{
    @{$read{$line[2]}}=@line[12..13,7];
  }
}

my @read_name=keys %read;
for(my $i=0;$i<@read_name;$i++){
  if((${$read{$read_name[$i]}}[0] > ${$read{$read_name[$i]}}[1]+ $AS_XS) && (${$read{$read_name[$i]}}[3] > ${$read{$read_name[$i]}}[4]+$AS_XS)){
    next;
  } elsif((${$read{$read_name[$i]}}[0] - ${$read{$read_name[$i]}}[1] <=$AS_XS) && (${$read{$read_name[$i]}}[3] - ${$read{$read_name[$i]}}[4]>$AS_XS)){
    print "$read_name[$i] L ${$read{$read_name[$i]}}[3] ${$read{$read_name[$i]}}[4] ${$read{$read_name[$i]}}[5]\n";}
  elsif((${$read{$read_name[$i]}}[0] - ${$read{$read_name[$i]}}[1] >$AS_XS) && (${$read{$read_name[$i]}}[3] - ${$read{$read_name[$i]}}[4]<=$AS_XS)){
    print "$read_name[$i] R ${$read{$read_name[$i]}}[0] ${$read{$read_name[$i]}}[1] ${$read{$read_name[$i]}}[2]\n";}
  else{
    next;
  }
}
