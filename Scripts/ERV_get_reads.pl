#!/usr/bin/env perl

use strict;

open F1,"$ARGV[0]";
my %reads=();

my @line=();
while(<F1>){
 @line=split;
 $reads{$line[0]} = $line[0];
}

open F2,"$ARGV[1]";
open F3,"$ARGV[2]";
open F4,">$ARGV[3]";
open F5,">$ARGV[4]";

my $name="";
while(<F2>){

#########################
 my $name=$_;
 $name=~s/@//;
 $name=~s/\/[12]//;
 chomp($name);

############
 if(exists($reads{$name})){
  print F4 "$_";
  my $line=<F2>;
  print F4 "$line";
  my $line=<F2>;
  print F4 "$line";
  my $line=<F2>;
  print F4 "$line";

#############
  my $line=<F3>;
  print F5 "$line";
  my $line=<F3>;
  print F5 "$line";
  my $line=<F3>;
  print F5 "$line";
  my $line=<F3>;
  print F5 "$line";
                         }
  else{
   my $line=<F2>;
   my $line=<F2>;
   my $line=<F2>;

   my $line=<F3>;
   my $line=<F3>;
   my $line=<F3>;
   my $line=<F3>;

   


      }

#########################
}
