#!usr/bin/perl
use strict;

my $line="";
my @line=();
my $dir="NA";
my %read=();

open F1,"$ARGV[0]";

while(<F1>){
 @line=split;
 my $dir=$line[0];
 if(exists($read{$line[2]})){
   if($dir eq "L"){unshift(@{$read{$line[2]}},@line[12..13,7]);}
   else{push(@{$read{$line[2]}},@line[12..13,7]);}
                            }
 else{@{$read{$line[2]}}=@line[12..13,7];}
# print "@{$read{$line[2]}}\n";
           }

my @read_name=keys %read;
for(my $i=0;$i<@read_name;$i++){
# print "@{$read{$read_name[$i]}}\n";
 if((${$read{$read_name[$i]}}[0] > ${$read{$read_name[$i]}}[1]+20) && (${$read{$read_name[$i]}}[3] > ${$read{$read_name[$i]}}[4]+20)){next;}    #### need to updated
 elsif((${$read{$read_name[$i]}}[0] - ${$read{$read_name[$i]}}[1] <=20) && (${$read{$read_name[$i]}}[3] - ${$read{$read_name[$i]}}[4]>20)){  #### updated to use AS - XS >=20
   print "$read_name[$i] L ${$read{$read_name[$i]}}[3] ${$read{$read_name[$i]}}[4] ${$read{$read_name[$i]}}[5]\n";}
 elsif((${$read{$read_name[$i]}}[0] - ${$read{$read_name[$i]}}[1] >20) && (${$read{$read_name[$i]}}[3] - ${$read{$read_name[$i]}}[4]<=20)){
  print "$read_name[$i] R ${$read{$read_name[$i]}}[0] ${$read{$read_name[$i]}}[1] ${$read{$read_name[$i]}}[2]\n";}
 else{next;}
 # print "@{$read{$read_name[$i]}}\n";
 # elsif((${$read{$read_name[$i]}}[12] == ${$read{$read_name[$i]}}[13])&& (${$read{$read_name[$i]}}[26] == ${$read{$read_name[$i]}}[27])){print "PE @{$read{$read_name[$i]}}\n";}
 # print "@{$read{$read_name[$i]}}\n";
                                }