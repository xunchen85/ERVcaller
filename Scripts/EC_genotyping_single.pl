#!usr/bin/perl
use strict;

my $line="";
my @line=();
my %VI_list=();
my @names=();
my %split=();
my %split3=();
my %split2=();

my $fully_mapped=0;
my $split_reads=0;
my $split_reads_short=0;

##########################  No. reads support non-TE insertions
open BP,"$ARGV[0]";            
while(<BP>){
 @line=split;
 my $name=$ARGV[1]."_".$ARGV[2]."_".$ARGV[3];
 unshift (@line,$name);
 push(@{$VI_list{$line[0]."|".$line[5]}},@line);
}

my %vi=();
my @reads=keys %VI_list;
for(my $i=0;$i<@reads;$i++){
 @line=@{$VI_list{$reads[$i]}};
 @names=split /\_/, $line[0];
 my $start1=0;
 my $end=0;
 my $length=0;
  my @temp1=($line[10]=~/(\d+)/g);
  my @temp2=($line[10]=~/([A-Z])/g);
  for(my $i=0;$i<@temp1;$i++){
  if($temp2[$i] eq "S"){   
   if($i==0){$start1=$temp1[$i]+1;$end=$temp1[$i];$length=$temp1[$i];}
   if($i==@temp1-1){$length+=$temp1[$i];}
                       }
  elsif($temp2[$i] eq "D"){next;}
  elsif($temp2[$i] eq "M" || $temp2[$i] eq "I"){
   if($i==0){$start1=1;$end=$temp1[$i];$length=$temp1[$i];}
   elsif($i==@temp1-1){$length+=$temp1[$i];$end+=$temp1[$i];}
   else{$end+=$temp1[$i];$length+=$temp1[$i]}  }
                              }

 if ($line[8] <= $ARGV[3] && ($line[8] + $length) >= $ARGV[3] && !($line[10] =~ "S")){          ### fully mapped with no split reads
  $fully_mapped++;
 } elsif ($start1 > 0 && abs($line[8] - $ARGV[3]) <= 5 && $line[10] =~ "S"){
  if ($start1>=20){
   $split_reads++;
  } else {
   $split_reads_short++;
  }
 } elsif ($end < $length && $line[10] =~ "S" && abs($line[8] + $end - $start1 - $ARGV[3]) <= 5){
  if ($length - $end >= 20){
   $split_reads++;
  } else {
   $split_reads_short++;
  }
 } elsif ((abs($line[8] - $ARGV[3]) <=5 || abs(($line[8] + $end - $start1) - $ARGV[3]) <= 5) && $line[10] =~ "S"){
  $fully_mapped++;
 } else {
  next;
 }
                           }

print "ID_Chr_Position No._reads_support_nonVI No._reads_support_VI Allele_Fraction Suggestive_reads_(split_reads_<20bp)\n";
my $name=$ARGV[1]."_".$ARGV[2]."_".$ARGV[3];
print "$name ";
print "$fully_mapped ";
print "$ARGV[4] ";

############ cellular proportion calculation
my $cellular_proportion=0;
if($ARGV[5] eq 2) {
 $cellular_proportion = ($ARGV[4]/2)/(($ARGV[4]/2) + $fully_mapped);
} elsif ($ARGV[5] eq 1) {
 $cellular_proportion = $ARGV[4]/($ARGV[4] + $fully_mapped);
}
print "$cellular_proportion $split_reads_short\n";
