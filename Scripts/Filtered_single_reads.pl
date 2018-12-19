#!usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
$Data::Dumper::Sortkeys=1;

my @line=();
my %reads_id=();
my @all_reads=();
my $start=0;
my $end=0;
my $chr="";

my $input="";
my $output="";
my $region=0;
GetOptions(
 'input=s'=>\$input,
 'output=s'=>\$output,
 'region=s'=>\$region
          );

################
open INPUT,"$input";
open OUTPUT,">$output";
while(<INPUT>){
 @line=split;
 if($line[0] =~ "U"){next;}
 unless($chr){
  $chr=$line[4];
  $start=$line[5];
  push(@all_reads,[@line]);
  $reads_id{$line[2]}=0;
 } elsif(eof(INPUT) || $chr ne $line[4] || ($chr eq $line[4] && abs($line[5]-$start)>($region+500))) {
  if(eof(INPUT)){push(@all_reads,[@line]);$reads_id{$line[2]}=0;}
  my @names=keys %reads_id;
  if($#names>=1){for(my $i=0;$i<@all_reads;$i++){print OUTPUT "@{$all_reads[$i]}\n";}}
  %reads_id=();
  @all_reads=();
  $chr=$line[4];
  $start=$line[5];
  @{$reads_id{$line[2]}}=@line;
 } else {
  push(@all_reads,[@line]);
  $reads_id{$line[2]}=0;
        }
             }


