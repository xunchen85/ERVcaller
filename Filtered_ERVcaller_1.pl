#!usr/bin/perl
use strict;

my $line="";
my @line=();
my $check_No_BPs=0;
my $check_Detected_BP=0;
my $false_or_not="";
my $false_or_not2="";
my $read_len=$ARGV[1];
my $depth = $ARGV[2];
my $sequencing_type = $ARGV[3];
###
open F1,"$ARGV[0]";
open OUT, ">$ARGV[0].1";
while(<F1>){
  @line=split;
  $check_No_BPs=0;
  $check_Detected_BP=0;
  $false_or_not="-";
  $false_or_not2="-";
  ### Check the upstream and downstream;
  if ($line[14] =~ "na") {
    $check_No_BPs=1;
  } else {
    $check_No_BPs=2;
  } 
  ### Check the breakpoints detected
  if ($line[14] =~ "D") {
    if ($line[14] =~ "^D" && $line[14] =~ "\;D") {
      $check_Detected_BP=2;
    } else {
      $check_Detected_BP=1;
    }
  } else {
    $check_Detected_BP=0;
  }

  ### Filtering
  if ($line[23] eq 0 && $sequencing_type eq "paired-end") {
    next;
  } elsif ( $line[23] + $line[27] < $depth ) {
    next;
  } elsif ($line[22] eq $depth && $line[23] eq 1 && $line[27] eq 0 && $sequencing_type eq "paired-end"){
    next;
  } elsif ($line[22] eq $depth && $line[23] eq 1 && $line[27] > 0 && $sequencing_type eq "paired-end") {
    $false_or_not="F_depth";
  } elsif ( $line[24] + $line[25] >= $line[23] && $sequencing_type eq "paired-end") {
    $false_or_not="F_validate";
  } else {
    $false_or_not="T";
  }

  ### parameters
  if ($line[18] eq "NA") {$line[18]=0;}
  if ($read_len <=110 && $read_len >=75 && $line[17] >0) {                     ### read length 75 bp to 110 bp;
    if ($line[21]>=30 && $line[18]/$line[17]<=0.5 && $line[17] >= 50){
      $false_or_not2="T";
    } elsif ($line[21]<30) {
      $false_or_not2="F_AS1";
    } elsif ($line[17] < 50) {
      $false_or_not2="F_AS2";
    } elsif ($line[18]/$line[17]>0.5) {
      $false_or_not2="F_Repeats";
    }
  } elsif ($read_len eq 150) {                                 ### read length 150 bp
    if ($line[21]>=45 && $line[18]/$line[17]<=0.5 && $line[17] >= 75){
      $false_or_not2="T";  
    } elsif ($line[21]<45) {
      $false_or_not2="F_AS1";
    } elsif ($line[17] < 75) {
      $false_or_not2="F_AS2";
    } elsif ($line[18]/$line[17]>0.5) {
      $false_or_not2="F_Repeats";
    }
  } elsif ($read_len eq 250) {                                ### read length 250 bp
    if ($line[21]>=75 && $line[18]/$line[17]<=0.5 && $line[17] >= 125){
      $false_or_not2="T";
    } elsif ($line[21]<75) {
      $false_or_not2="F_AS1";
    } elsif ($line[17] < 125) {
      $false_or_not2="F_AS2";
    } elsif ($line[18]/$line[17]>0.5) {
      $false_or_not2="F_Repeats";
    }    
  }  
  $line[16]=$false_or_not.";".$false_or_not2;
#  print join("\t",@line),"\n";
  print OUT join("\t",@line),"\n";
#  print "@line $false_or_not $false_or_not2 $check_No_BPs $check_Detected_BP\n";
}
