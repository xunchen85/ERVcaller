#!usr/bin/perl -w

use strict;
my %file=();
my $b1="";
my $read_name="";
my $temp04="";
my $temp05="";
my $file_name="";
my $aa="";
my @temp=();
my $read="";

open F1,"$ARGV[1]"||di("di");
open B1, "$ARGV[3]"||di("di too");
open S1,">$ARGV[5]"||di("di too");

while($read_name=<F1>){
   $read=$read_name;
   my @temp11=split /\s+/, $read;
   $read_name=$temp11[0];
   $read_name=~s/\/[12]\s+$//;
   $read_name=~s/\/[12]$//;
   $read_name=~s/^@//;
   my $temp01=<F1>;
   my $temp02=<F1>;
   my $temp03=<F1>;
   $temp04=$read.$temp01.$temp02.$temp03;
   $file{$read_name}=$temp04;
}
while($temp05=<B1>){
   my $temp06=$temp05;
   $temp05=~s/\s+$//;  
   if(exists($file{$temp05})){print S1 "$file{$temp05}";}

                  }
close F1;
close B1;
close S1;
