#!usr/bin/perl
use strict;

my $line="";
my @line=();
my $dir="NA";
my %read=();

open F1,"$ARGV[0]";
open F2,"$ARGV[1]";
open HF,"$ARGV[2]";

######## hf file
my %hf=();
while (<HF>) {
 @line=split;
 $hf{$line[0]}=$line[0];
}


########
while(<F1>){
  @line=split;
  my $dir=$line[0];
  if (exists($hf{$line[2]})) {
    next;
  } elsif (exists($read{$line[2]})){
    if($dir eq "L"){
      unshift(@{$read{$line[2]}},@line);
    } else {
      push(@{$read{$line[2]}},@line);
    }
    if((${$read{$line[2]}}[12] > ${$read{$line[2]}}[13]+20) && (${$read{$line[2]}}[26] > ${$read{$line[2]}}[27]+20)){
      $hf{$line[2]}=0;
      delete $read{$line[2]};
    } elsif((${$read{$line[2]}}[12] <= ${$read{$line[2]}}[13]+20)&& (${$read{$line[2]}}[26] > ${$read{$line[2]}}[27]+20)){
      unshift(@{$read{$line[2]}},"PR");
    } elsif((${$read{$line[2]}}[12] > ${$read{$line[2]}}[13]+20) && (${$read{$line[2]}}[26] <= ${$read{$line[2]}}[27]+20)){
      unshift(@{$read{$line[2]}},"PF");
    } else {
      $hf{$line[2]}=0;
      delete $read{$line[2]};
    }
  } else {
    @{$read{$line[2]}}=@line;
  }
}

#my @read_name=keys %read;
#for(my $i=0;$i<@read_name;$i++){
# if((${$read{$read_name[$i]}}[12] > ${$read{$read_name[$i]}}[13]+20) && (${$read{$read_name[$i]}}[26] > ${$read{$read_name[$i]}}[27]+20)){next;}    #### can be updated
# elsif((${$read{$read_name[$i]}}[12] <= ${$read{$read_name[$i]}}[13]+20)&& (${$read{$read_name[$i]}}[26] > ${$read{$read_name[$i]}}[27]+20)){unshift(@{$read{$read_name[$i]}},"PR");}
# elsif((${$read{$read_name[$i]}}[12] > ${$read{$read_name[$i]}}[13]+20) && (${$read{$read_name[$i]}}[26] <= ${$read{$read_name[$i]}}[27]+20)){unshift(@{$read{$read_name[$i]}},"PF");}
# else{next;}
# # print "@{$read{$read_name[$i]}}\n";
# # elsif((${$read{$read_name[$i]}}[12] == ${$read{$read_name[$i]}}[13])&& (${$read{$read_name[$i]}}[26] == ${$read{$read_name[$i]}}[27])){print "PE @{$read{$read_name[$i]}}\n";}
# # print "@{$read{$read_name[$i]}}\n";
#                                }

while(<F2>){
  @line=split;
  if ($line[12]<=30){next;}                       ### can specify the AS of the chimeric reads for ERVs
  if (exists($read{$line[2]}) && !(exists($hf{$line[2]}))){
    if (${$read{$line[2]}}[0] eq "PF"){
      ${$read{$line[2]}}[10]="TE";
      print "$line[0] $line[1] $line[2] @{$read{$line[2]}}[4..13] @line[3..12]\n";
    } elsif (${$read{$line[2]}}[0] eq "PR") {
      ${$read{$line[2]}}[24]="TE";
      print "$line[0] $line[1] $line[2] @{$read{$line[2]}}[18..27] @line[3..12]\n";
    }
  }
}
