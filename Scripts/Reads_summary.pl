#!usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
$Data::Dumper::Sortkeys=1;
my $input="test_breakpoint"; 
my $output="test_read_summary";
GetOptions(
  'input=s'=>\$input,
  'output=s'=>\$output
          );

open INPUT, "$input";
open OUTPUT, ">$output";

my $line="";
my @line="";
my @read=();
my @group_number=();
my @TE_sequence_number=();
my @alignment_score=();
my $read_name="";

while($line=<INPUT>){
  chomp($line);
  @line=split /\s+/, $line;
  if ($line[3] % 256 < 128 && $line[0] eq "L") {next;}
  elsif($line[3] % 256 >= 128 && $line[0] eq "R") {next;}
  if(@line eq 23){splice @line,10,1;}
  unless($read_name){$read_name=$line[2];}
  if($read_name eq $line[2]){
    push @read, $line;
    push @group_number, $line[1];
    push @TE_sequence_number, $line[14];
    push @alignment_score, $line[21];
                            }
  else{
    my %count=();$count{$_}++ for @group_number;my @count=%count;my $group_count=@count/2;
    %count=();$count{$_}++ for @TE_sequence_number; @count=%count;my $TE_sequence_count=@count/2;
    %count=();$count{$_}++ for @alignment_score; @count=%count; my $alignment_all=join('|',@count);
    my $len=@read;
    for(my $i=0;$i<$len;$i++){
      print OUTPUT "$read[$i] $group_count $TE_sequence_count $alignment_all\n";
                             }
    @read=();@group_number=();@TE_sequence_number=();@alignment_score=();$read_name=$line[2];
    push @read, $line;
    push @group_number, $line[1];
    push @TE_sequence_number, $line[14];
    push @alignment_score, $line[21];
      }
                   }
