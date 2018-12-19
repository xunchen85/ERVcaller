#!usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
no warnings 'deprecated';
$Data::Dumper::Sortkeys=1;

##############
my @line=();
my @integration=();
my %integration_list=();
my @position=();
my $same=0;
my @all_reads=();
my $all_count=0;

###############
my $input="";
my $output="";
my $region=500;
my %redundant_reads=();
my %duplicate_reads=();

############# a: remove reads with all kind of repeats; k: keep all reads; r: remove all reads overlap simple_repeat; f: filter reads based on the length of overlap region;
my $method="f";              # default
GetOptions(
  'input=s'=>\$input,
  'output=s'=>\$output,
  'region=s'=>\$region,
  'method=s'=>\$method
          );

#############
open INPUT, "$input";
my %paired_exclude=();
while(<INPUT>){
  my @line=split;
  if($line[0] =~ "U"){next;}
  my @line2=split /\|/, $line[32];
  if(($line2[2] >=($line[27]-$line[26]+1)) && ($line[31] eq "P1" && $line2[2]>=30) && (($line[0] =~ "R" && $line[25] eq "L") || ($line[0] =~ "L" && $line[25] eq "R"))){
    $paired_exclude{$line[2]}=$line[25];
  }
}
close INPUT;

#############
open INPUT, "$input";
open OUTPUT, ">$output";

#############
while(my $line_input=<INPUT>) {                                  ### read integration file;
  @line=split /\s+/, $line_input;
  $line[0]=~s/^b://;
  $line[0]=~s/^c://;
  if ($line[0] eq "U") {next;}                                   ### exclude paired-end unmapped reads;
  if (exists($paired_exclude{$line[2]})) {next;}

  ### filter out step;
  my @overlap_r=split /\|/, $line[30];
  my @overlap_h=split /\|/, $line[32];
  $line[20]=$line[20]."_".$line[13];
  my $overlap_b=0;
  if($overlap_h[2] eq 0 && $overlap_r[2]) {
    $overlap_b=0;
  } elsif ($overlap_r[2] eq 0) {
    $overlap_b=$overlap_h[2];
  } elsif ($overlap_h[2] eq 0) {
    $overlap_b=$overlap_r[2];
  } elsif ($overlap_h[0] >= $overlap_r[0] && $overlap_h[0] <= $overlap_r[1]) {
      if ($overlap_h[1]>=$overlap_r[1]) {
        $overlap_b=$overlap_h[1]-$overlap_r[0]+1;
      } else {
        $overlap_b=$overlap_r[2];
      }
  } elsif ($overlap_h[1] >= $overlap_r[0] && $overlap_h[1] <= $overlap_r[1]){
      if ($overlap_h[0]>=$overlap_r[0]) {
        $overlap_b=$overlap_r[2];
      } else {
        $overlap_b=$overlap_r[1]-$overlap_h[0];
      }
  } else {
    $overlap_b=$overlap_r[2]+$overlap_h[2];
  }
  if ($line[12]%8>=4) {next;}      ### filter out unmapped reads which cannot mapped to TE database;
  if($method eq "f"){              # filter reads by overlap length;
    if(($line[29] =~ "Simple_repeat" || $line[29] eq "dust") && ($line[27]-$line[26]+1-$overlap_r[2]<20)) {next;}    
    if($line[27]-$line[26]+1<20) {next;}
  }

  if($line[0] =~ "L" && $line[25] =~"R") {
    next;
  } elsif ($line[0] =~ "R" && $line[25] =~"L") {                      ## added on 4/11/16
    next;
  }
  if($line[31] eq "P1" && $overlap_h[2]>0){next;}                     ## added on 5/23/16

##############
  my $read_name_t1=$line[14]+$line[27]-$line[26]+1;
  my $read_name=$line[13]."_".$line[14]."_".$read_name_t1;
  my $read_name2=$line[1]."_".$line[14]."_".$read_name_t1;
  my $read_name3="";
  my $read_name4="";
  if ($line[12]%32>16) {
    $read_name3=$line[13]."_".$line[4]."_".$line[5]."_".$line[7]."_".$read_name_t1;
  } else {
    $read_name3=$line[13]."_".$line[4]."_".$line[5]."_".$line[7]."_".$line[14];
  }
  if ($line[12]%32>16) {
    $read_name4=$line[1]."_".$line[4]."_".$line[5]."_".$line[7]."_".$read_name_t1;
  } else {
    $read_name4=$line[1]."_".$line[4]."_".$line[5]."_".$line[7]."_".$line[14];
  }

##############                           
  if(!(${$integration[0]}[0])){                   ### if this is the first read of each vip;
    $redundant_reads{$read_name}=0;                ### updated on 4/09/16
    $redundant_reads{$read_name2}=0;               ### updated on 4/11/16
    $redundant_reads{$read_name3}=0;
    $redundant_reads{$read_name4}=0;  

    push(@{$integration[0]},($line[1],$line[13],$line[0],@line[2..12],@line[14..32]));
    $integration_list{$line[13]}=$same++;          ### the number of gi number belone to the same vip;
    $position[0]=$line[4];                         ### the chromosome of each vip;
    $position[1]=$line[5];                         ### start physical positions of each vip;
    $position[2]=$line[5];                         ### end physical positions of each vip; 
  } elsif ($line[4] eq $position[0] && abs($line[5]-$position[2])<=$region) {           ### if not the first read and still belong to the same vip;
    ##### remove duplicate reads add on 11/11/18;
    my $tmp1=$line[0]."|".$line[4]."|".$line[5]."|".$line[7];
    if (exists($duplicate_reads{$tmp1})) {
      if ($duplicate_reads{$tmp1} eq $line[2]) {
        next;
      }
    } else {
      $duplicate_reads{$tmp1} = $line[2];
    }
    
    if(exists($redundant_reads{$read_name})){next;}
    elsif(exists($redundant_reads{$read_name2})){next;}
    elsif(exists($redundant_reads{$read_name3})){next;}
    elsif(exists($redundant_reads{$read_name4})){next;}
    else{
      $redundant_reads{$read_name}=0;
      $redundant_reads{$read_name2}=0;
      $redundant_reads{$read_name3}=0;
      $redundant_reads{$read_name4}=0;
        }
    $position[2]=$line[5];                                     ### the end physical positions on chromosome
    push(@all_reads,$line[2]);                                 ### put this read name to all_reads array;

###############                                                ### check if this is the same gi number or not;
    if(exists($integration_list{$line[13]})){                  ### if this is not first read belone to one gi number;     
      my $len=@{$integration[$integration_list{$line[13]}]};   ### the length of this TE seq;
      for(my $i=14;$i<$len;$i+=31){                            ### order them by human position;
        if($line[14]<=${$integration[$integration_list{$line[13]}]}[$i]){
          splice @{$integration[$integration_list{$line[13]}]},$i-12,0,($line[0],@line[2..12],@line[14..32]);
          last;                                                         }
        elsif($line[14]>${$integration[$integration_list{$line[13]}]}[$len-11]){
          push(@{$integration[$integration_list{$line[13]}]},($line[0],@line[2..12],@line[14..32]));
          last;                                                                 }
                                }
    } else{                                         ### if this is the first read belone to the same TE seq; 
      push(@integration,[($line[1],$line[13],$line[0],@line[2..12],@line[14..32])]);
      $integration_list{$line[13]}=$same++;       ### the number of TE seq belone to this TE;
    }
  }
###############                                  ### if the distance larger than regions, it would be different TEs;
 if($line[4] ne $position[0] || abs($line[5]-$position[2])>$region || eof(INPUT)){
#LAST:
   %redundant_reads=();
   my $numb=@integration;                        ### number of total sequence of this TE;
   for(my $i=0;$i<$numb;$i++){                   ### print out each line;
     my @total_reads=();                         
     my $len2=@{$integration[$i]};               ### length of each gi number;
     for(my $j=3;$j<$len2;$j+=31){push @total_reads,${$integration[$i]}[$j];}     ### push all read name belone to each gi;
     my %count=();$count{$_}++ for @total_reads;
     my @count=%count;
     my $read_count=@count/2;                    ### column 3: total number of reads this TE seq;
     my $read_line=join('|',@count);             ### column 4: read name and number of hits;
     my $read_all=@total_reads;                  ### column 2: total number of hits this TE seq;
##############                                   ### total reads each TE;
     %count=();$count{$_}++ for @all_reads;
     @count=%count;
     $all_count=@count/2;                        ### column 1: total number of reads this TE;
     print OUTPUT "$all_count ";
     print OUTPUT "$read_all ";
     print OUTPUT "$read_count ";
     print OUTPUT "$read_line ";
     print OUTPUT "@position ";
     my @tem_inte=@{$integration[$i]};
     for(my $b=0;$b<@tem_inte;$b++){
      if($tem_inte[$b] || $tem_inte[$b] eq 0){print OUTPUT "$tem_inte[$b] ";}
      else{print OUTPUT "q ";}
                                   }   
    print OUTPUT "\n";
                             }
#################                                ### start second TE;
   @integration=();$same=0;%integration_list=(); @all_reads=();$all_count=0;
   if(!(${$integration[0]}[0])){                 ### first read;          # corrected on 1/29/16
     push(@all_reads,$line[2]);
     push(@{$integration[0]},($line[1],$line[13],$line[0],@line[2..12],@line[14..32]));
     $integration_list{$line[13]}=$same++;
     $position[0]=$line[4];$position[1]=$line[5];$position[2]=$line[5];

  ##############
  my $read_name_t1=$line[14]+$line[27]-$line[26]+1;
  my $read_name=$line[13]."_".$line[14]."_".$read_name_t1;
  my $read_name2=$line[1]."_".$line[14]."_".$read_name_t1;
  my $read_name3="";
  my $read_name4="";
  if($line[12]%32>16){$read_name3=$line[13]."_".$line[4]."_".$line[5]."_".$line[7]."_".$read_name_t1;}
  else{$read_name3=$line[13]."_".$line[4]."_".$line[5]."_".$line[7]."_".$line[14];}
  if($line[12]%32>16){$read_name4=$line[1]."_".$line[4]."_".$line[5]."_".$line[7]."_".$read_name_t1;}
  else{$read_name4=$line[1]."_".$line[4]."_".$line[5]."_".$line[7]."_".$line[14];}
   #################

  $redundant_reads{$read_name}=0;
  $redundant_reads{$read_name2}=0;
  $redundant_reads{$read_name3}=0;
  $redundant_reads{$read_name4}=0;
                               }
                                                                                   } 
              }
