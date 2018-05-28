#!usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

################
my $position="test_sm.sam";
my $te="test_vsu.sam";
my $type="test.type";
my $output="test";

GetOptions('position=s'=>\$position,
           'te=s'=>\$te,
           'type=s'=>\$type,
           'output=s'=>\$output);

################
my %position=();
my @position=();
my @te=();
my @type=();
my %type=();
my $md="";
my $as="";

#################
open TYPE, "$type";
while(<TYPE>){
 @type=split;
 if($type[1] eq "L" || $type[1] eq "R"){$type{$type[0]}=$type[1];}
             }

#################    # Read position information;
open POSITION, "$position";
while(<POSITION>){
 @position=split;
 $md="";$as="";
 for(my $i=11;$i<@position;$i++){
  if($position[$i]=~s/MD:Z://){$md=$position[$i];}
  if($position[$i]=~s/AS:i://){$as=$position[$i];}
                                }
 $position{$position[0]}=$position[0]." ".$position[1]." ".$position[2]." ".$position[3]." ".$position[4]." ".$position[5]." ".$position[6]." ".$position[7]." ".$position[8]." ".$md." ".$as;
                 }
##################
open TE, "$te";
open OUT1, ">${output}_breakpoint";
while(<TE>){
 @te=split;
 $md="NA";$as="NA";
 if($te[2] eq "*"){next;}
 for(my $i=11;$i<@te;$i++){
  if($te[$i]=~s/MD:Z://){$md=$te[$i];}
  if($te[$i]=~s/AS:i://){$as=$te[$i];}
                          }
 if($as<22){next;}
 if(exists($type{$te[0]})){
  print OUT1 "$type{$te[0]} unknown $position{$te[0]} @te[1..8] $md $as\n";
                          }
           }
