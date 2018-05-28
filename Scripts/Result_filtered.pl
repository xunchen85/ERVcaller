#!usr/bin/perl
use strict;

my $line1="";
my @line1=();
my %vip=();
my %virus=();
my $name1="";
my $name2="";
my $name3="";
my $chimeric=0;
my $chimeric_quality=0;
my $split_reads=0;
my $split_reads_quality=0;
my %reads=();

my %all_virus=();
my %all_vip=();
my %all_gi=();

my $unique_virus=0;
my $unique_vip=0;
my $unique_gi=0;

my $average_gi=0;
my $average_unique=0;
my $average_all=0;

open FILE,"$ARGV[0]";
############ Descriptions: 1. average score this one; 2. average score for unique reads; 3. unique reads virus; 4. unique reads vip; 4. unique reads gi; 5. cguneruc number; 6. chimeric quality; 7. split reads number; 8. split reads quality; 9. total reads each vip; 10. total hits each gi; 11. total reads each gi;
print "aver_a\taver_u\tu_virus\tu_vip\tuniq_gi\tchim\tchim_q\tsplit\tsplit_q\tt_vip\tt_hits\tt_gi\n";
#######################
while(<FILE>){                # while1 start;
 my @line=split;
 $name2=$line[4]."|".$line[5]."|".$line[6];

#######################
 if($name1 eq $name2 || !($name1)){              ### the same vip;
  $name1=$name2;
  push(@{$vip{$name1}},[@line]);

######################
  my @temp2=split /\|/, $line[3];
  for(my $i=0;$i<@temp2;$i+=2){
   $name3=$line[8]."|".$temp2[$i];
   $all_gi{$name3}=$temp2[$i+1];
   $name3=$line[7]."|".$temp2[$i];
   if(exists($all_virus{$name3})){
    if($all_virus{$name3}<$temp2[$i+1]){$all_virus{$name3}=$temp2[$i+1];}
                                     }
    else{$all_virus{$name3}=$temp2[$i+1];}

#########
   if(exists($all_vip{$temp2[$i]})){
    if($all_vip{$temp2[$i]}<$temp2[$i+1]){$all_vip{$temp2[$i]}=$temp2[$i+1];}
                                   }
    else{$all_vip{$temp2[$i]}=$temp2[$i+1];}
                              }
   if(eof(FILE)){goto LAST1;}
                                 }
###############
 elsif($name2 ne $name1){                              ### different vip
########## print out results;
LAST1:  my @temp1=@{$vip{$name1}};
  for(my $i=0;$i<@temp1;$i++){
   my @temp3=@{$temp1[$i]};
   for(my $j=9;$j<@temp3;$j+=31){
#    $average_all+=$temp3[$j+19];                     ##### quality
     $average_all+=$temp3[$j+25]-$temp3[$j+24]+1;
    $name3="";
    $name3=$temp3[7]."|".$temp3[$j+1];
    if(exists($all_virus{$name3})){if($all_virus{$name3}<2){$unique_virus++;}}
    if(exists($all_vip{$temp3[$j+1]})){if($all_vip{$temp3[$j+1]}<2){$unique_vip++;}}
    $name3=$temp3[8]."|".$temp3[$j+1];
    if(exists($all_gi{$name3})){
     if($all_gi{$name3}<2){
     $unique_gi++;
#    $average_unique+=$temp3[$j+19];                 ##### quality 
     $average_unique+=$temp3[$j+25]-$temp3[$j+24]+1;
                          }
                               }
    if($temp3[$j]=~"S"){
     $split_reads++;
#    $split_reads_quality+=$temp3[$j+19];            ##### quality
      $split_reads_quality+=$temp3[$j+25]-$temp3[$j+24]+1;
                       }
    elsif($temp3[$j]=~"L" || $temp3[$j]=~"R"){
     $chimeric++;
#     $chimeric_quality+=$temp3[$j+19];               ##### quality
      $chimeric_quality+=$temp3[$j+25]-$temp3[$j+24]+1;
                                             }
                                }
   if($unique_gi){$average_unique=$average_unique/$unique_gi;}
   else{$average_unique=0;}
   $average_all=$average_all/$temp3[1];
   if($chimeric){$chimeric_quality=$chimeric_quality/$chimeric;}
   else{$chimeric_quality=0;}
   if($split_reads){$split_reads_quality=$split_reads_quality/$split_reads;}
   else{$split_reads_quality=0;}
   print "$average_all\t$average_unique\t$unique_virus\t$unique_vip\t$unique_gi\t$chimeric\t$chimeric_quality\t$split_reads\t$split_reads_quality\t$temp3[0]\t$temp3[1]\t$temp3[2]\t@temp3[3..$#temp3]\n";
   $average_unique=0; $average_all=0;
   $unique_vip=0; $unique_virus=0; $unique_gi=0;$chimeric=0;$split_reads=0;
   $chimeric_quality=0;$split_reads_quality=0;
                             }

########## reset the value;
  %vip=();
  $name1=$name2;
  $average_all=0;$average_unique=0;%all_virus=();%all_vip=(); %all_gi=(); $unique_vip=0; $unique_virus=0; $unique_gi=0;
  push(@{$vip{$name1}},[@line]);
  my @temp2=split /\|/, $line[3];

#####
  my @temp2=split /\|/, $line[3];
  for(my $i=0;$i<@temp2;$i+=2){
   $name3=$line[8]."|".$temp2[$i];
   $all_gi{$name3}=$temp2[$i+1];

   $name3=$line[7]."|".$temp2[$i];
   if(exists($all_virus{$name3})){   
    if($all_virus{$name3}<$temp2[$i+1]){$all_virus{$name3}=$temp2[$i+1];}
                                     }
    else{$all_virus{$name3}=$temp2[$i+1];}

#########
   if(exists($all_vip{$temp2[$i]})){
    if($all_vip{$temp2[$i]}<$temp2[$i+1]){$all_vip{$temp2[$i]}=$temp2[$i+1];}
                                   }
    else{$all_vip{$temp2[$i]}=$temp2[$i+1];}
                              }
                        }         # different vip end
             }                # while1 end;
