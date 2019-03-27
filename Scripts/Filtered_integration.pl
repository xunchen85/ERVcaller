#!/usr/bin/env perl

use strict;

my $line="";
my @line=();
my %repeat=();
my %h_unmapped=();
my %h_single=();
my %h_paired=();
my $read_name="";

##########################                   open repeat files;
open REPEAT,"$ARGV[0].repeat2";
$line=<REPEAT>;                   ### remove title;
while(<REPEAT>){
 @line=split;
 @{$repeat{$line[0]}}=@line[1..$#line];
               }
close REPEAT;

##########################                   open realign to human file;
open HUNMAPPED,"$ARGV[0].hunmapped";
while(<HUMAPPED>){
 $read_name="";
 @line=split;
 if($line[0] eq "L"){$read_name=$line[0]."/1";}
 elsif($line[0] eq "R"){$read_name=$line[0]."/2";}
 @{$h_unmapped{$read_name}}=@line;
                 }
close HUNMAPPED;

#####################
my $orderr=0;
my %mapped_temp=();

open HMAPPED,"$ARGV[0].hmapped";
while(<HMAPPED>){
 @line=split;
 my $com2=join(" ",(@line[0..8],$line[17],$line[26]));
 if($line[0] eq "L" || $line[0] eq "R" || $line[0] eq "E"){
  if($line[0] eq "L" && $line[10]%32>=16){my $start_1=$line[5]-$line[4]+1;my $end_1=$line[5]-$line[3]+1;$line[3]=$start_1;$line[4]=$end_1;}
  elsif($line[0] eq "R" && $line[10]%32>=16){my $start_1=$line[8]-$line[7]+1;my $end_1=$line[8]-$line[6]+1;$line[6]=$start_1;$line[7]=$end_1;}
    push(@{$h_single{$line[1]}},[@line]);      
                                                          }
 elsif($line[0] =~ "^P"){
  if($line[10]%32>=16){my $start_1=$line[5]-$line[4]+1;my $end_1=$line[5]-$line[3]+1;$line[3]=$start_1;$line[4]=$end_1;}
  if($line[19]%32>=16){my $start_1=$line[8]-$line[7]+1;my $end_1=$line[8]-$line[6]+1;$line[6]=$start_1;$line[7]=$end_1;}
  if(exists($mapped_temp{$com2})){next;}
  else{
   $mapped_temp{$com2}=$line[1];
   push(@{$h_paired{$line[1]}},[@line[0..8],$line[17],$line[26]]);
      }
                        }
                }
close HMAPPED;
%mapped_temp=();

##############################               end of open human file;
 my @tem_human=();
 my $overlap_h=0;
 my $overlap_total_h=0;
 my $type_h_keep="n";
 my @out_h_keep=(0) x 2;
 my $overlap_total_h_keep=0;

 my $length_h=0;
 my $type_h="N";
 my $vstart=0;my $vend=0;my $vlength=0; my $overlap=0; my $overlap_total=0;
 my $type="";
 my $vlength_align=0;
 my $hlength_align=0;
 my $direction="";
 my $hlength=0;my $hstart=0;my $hend=0;
 my @out_h=(0) x 2;
 my @out_v=(0) x 2;

#############################              # open the file of integration;
open INTEGRATION,"$ARGV[0].integration";
open ERROR,"$ARGV[0].error";
while(<INTEGRATION>){
 @line=split; 
 if($#line<24){print ERROR "@line\n";}
 if($line[0]=~"S1" && $#line eq 25){splice @line,10,1;}  # Remove one colume;
 elsif($#line eq 25 && !($line[0]=~"S1")){splice @line,12,1;}
 if($line[24]=~"^|"){$line[24]="$line[21]"."$line[24]";}
 @tem_human=();
 $overlap_h=0;
 $overlap_total_h=0;
 $length_h=0;
 $type_h="N";
 $vstart=0;$vend=0;$vlength=0;$overlap=0;$overlap_total=0;
 $type="";
 $vlength_align=0;
 $hlength_align=0;
 $direction="";
 $hlength=0;$hstart=0;$hend=0;
 @out_h=(0) x 2;
 @out_v=(0) x 2;

########################                   # get mapping information on human
 my @temp1=($line[7]=~/(\d+)/g);
 my @temp2=($line[7]=~/([A-Z])/g);
 for(my $i=0;$i<@temp1;$i++){
  if($temp2[$i] eq "S"){   
   if($i==0){$hstart=$temp1[$i]+1;$hend=$temp1[$i];$hlength=$temp1[$i];}
   if($i==@temp1-1){$hlength+=$temp1[$i];}
                       }  
  elsif($temp2[$i] eq "D"){next;}
  elsif($temp2[$i] eq "M" || $temp2[$i] eq "I"){
   if($i==0){$hstart=1;$hend=$temp1[$i];$hlength=$temp1[$i];}
   elsif($i==@temp1-1){$hlength+=$temp1[$i];$hend+=$temp1[$i];}
   else{$hend+=$temp1[$i];$hlength+=$temp1[$i];}}
                            }

########################                  # get mapping information on ERV;
 @temp1=($line[16]=~/(\d+)/g);
 @temp2=($line[16]=~/([A-Z])/g);
 for(my $i=0;$i<@temp1;$i++){
  if($temp2[$i] eq "S"){
   if($i==0){$vstart=$temp1[$i]+1;$vend=$temp1[$i];$vlength=$temp1[$i];}
   if($i==@temp1-1){$vlength+=$temp1[$i];}
                       }
  elsif($temp2[$i] eq "D"){next;}
  elsif($temp2[$i] eq "M" || $temp2[$i] eq "I"){
   if($i==0){$vstart=1;$vend=$temp1[$i];$vlength=$temp1[$i];}
   elsif($i==@temp1-1){$vlength+=$temp1[$i];$vend+=$temp1[$i];}
   else{$vend+=$temp1[$i];$vlength+=$temp1[$i];}
                                               }
                            }
  if($line[12]%32>=16){my $start_1=$vlength-$vend+1;my $end_1=$vlength-$vstart+1;$vstart=$start_1;$vend=$end_1;}

#######################                                 # if it is softclip reads;
 if($line[0] eq "S1"){
  if($line[12]%256>=128){$direction="/2";}
  elsif($line[12]%256<128 && $line[12]%256>=64){$direction="/1";}
  else{$direction="/1";}

#######################                                 # to check repeats;
  $read_name=$line[2]."|".$line[3]."|".$line[9];
  if(exists($repeat{$read_name})){
   my @rep_tem1=@{$repeat{$read_name}};
   for(my $i=0;$i<@rep_tem1;$i+=4){
    if($type ne $rep_tem1[$i+3] && $type){$type=$type."|".$rep_tem1[$i+3];}
    elsif($rep_tem1[$i+3]){$type=$rep_tem1[$i+3];}
    else{$type="q";}

   if($vstart < $rep_tem1[$i]){
    if($vend < $rep_tem1[$i]){$overlap=0;}
    elsif($vend >= $rep_tem1[$i] && $vend <= $rep_tem1[$i+1]){$overlap=$vend-$rep_tem1[$i]+1;unless($out_v[0]){$out_v[0]=$rep_tem1[$i];}$out_v[1]=$vend;}
    elsif($vend > $rep_tem1[$i+1]){$overlap=$rep_tem1[$i+1]-$rep_tem1[$i]+1;unless($out_v[0]){$out_v[0]=$rep_tem1[$i];}$out_v[1]=$rep_tem1[$i+1];}
                              }
   elsif($vstart >= $rep_tem1[$i] && $vstart <= $rep_tem1[$i+1]){
    unless($out_v[0]){$out_v[0]=$vstart;}
    if($vend <= $rep_tem1[$i+1]){$overlap=$vend-$vstart+1;$out_v[1]=$vend;}
    elsif($vend > $rep_tem1[$i+1]){$overlap=$rep_tem1[$i+1]-$vstart+1;$out_v[1]=$rep_tem1[$i+1];}
                                                                }
   elsif($vstart > $rep_tem1[$i+1]){$overlap=0;}
   $overlap_total+=$overlap;
                                  }
                                }
  else{$type="n";$overlap=0;$overlap_total=0;}

#####################               # to check realign to human region for softclip reads;
  $type_h="N",$overlap_h=0;$length_h=0;
                    }               # end of if it is softclip reads;

#####################      # if it is single-end mapped reads;
 elsif($line[0] eq "L" || $line[0] eq "R"){             #1 

#####################               # to check repeats;  
  if($line[12]%256>=128){$direction="/2";}
  elsif($line[12]%256<128 && $line[12]%256>=64){$direction="/1";}
  else{$direction="/1";}
  $read_name=$line[2].$direction;
  if(exists($repeat{$read_name})){                      #2 end
   my @rep_tem1=@{$repeat{$read_name}};
   for(my $i=0;$i<@rep_tem1;$i+=4){
    if($type ne $rep_tem1[$i+3] && $type){$type=$type."|".$rep_tem1[$i+3];}
    elsif($rep_tem1[$i+3]){$type=$rep_tem1[$i+3];}
    else{$type="q";}
   if($vstart < $rep_tem1[$i]){
    if($vend < $rep_tem1[$i]){$overlap=0;}
    elsif($vend >= $rep_tem1[$i] && $vend <= $rep_tem1[$i+1]){$overlap=$vend-$rep_tem1[$i]+1;unless($out_v[0]){$out_v[0]=$rep_tem1[$i];}$out_v[1]=$vend;}
    elsif($vend > $rep_tem1[$i+1]){$overlap=$rep_tem1[$i+1]-$rep_tem1[$i]+1;unless($out_v[0]){$out_v[0]=$rep_tem1[$i];}$out_v[1]=$rep_tem1[$i+1];}
                                }
   elsif($vstart >= $rep_tem1[$i] && $vstart <= $rep_tem1[$i+1]){
    unless($out_v[0]){$out_v[0]=$vstart;}
    if($vend <= $rep_tem1[$i+1]){$overlap=$vend-$vstart+1;$out_v[1]=$vend;}
    elsif($vend > $rep_tem1[$i+1]){$overlap=$rep_tem1[$i+1]-$vstart+1;$out_v[1]=$rep_tem1[$i+1];}
                                                                }
   elsif($vstart > $rep_tem1[$i+1]){$overlap=0;}
   $overlap_total+=$overlap;
                                  }
                                 }                     #2 end
  else{$type="N";$overlap=0;$overlap_total=0;}

#####################            # check the overlap with human;
  if(exists($h_paired{$line[2]})){                     # found the pair    #3 start
   my @paired=@{$h_paired{$line[2]}};
   my $row=0;
   my @candidate=();
  $type_h_keep="n";@out_h=(0) x 2;@out_h_keep=(0) x 2;$overlap_total_h_keep=0;$overlap_total_h=0;

   for(my $j=0;$j<@paired;$j++){                       # 4 start
    @tem_human=();
     push(@tem_human,($paired[$j][0],$paired[$j][1],$paired[$j][2],$paired[$j][3],$paired[$j][4],$paired[$j][5],$paired[$j][6],$paired[$j][7],$paired[$j][8],$paired[$j][9],$paired[$j][10],$paired[$j][11]));
   $overlap_total_h=0;
   $overlap_h=0;
   @out_h=(0) x 2;
   my $i=0;
   if($line[12]%256>=128){$i=5;}
   elsif($line[12]%256<128 && $line[12]%256>=64){$i=2;}
   else{$i=2;}
   if($vstart < $tem_human[$i+1]){
    if($vend < $tem_human[$i+1]){$overlap_h=0;}
    elsif($vend >= $tem_human[$i+1] && $vend <= $tem_human[$i+2]){$overlap_h=$vend-$tem_human[$i+1]+1;unless($out_h[0]){$out_h[0]=$tem_human[$i+1];}$out_h[1]=$vend;}
    elsif($vend > $tem_human[$i+2]){$overlap_h=$tem_human[$i+2]-$tem_human[$i+1]+1;unless($out_h[0]){$out_h[0]=$tem_human[$i+1];}$out_h[1]=$tem_human[$i+2];}
                                 }
   elsif($vstart >= $tem_human[$i+1] && $vstart <= $tem_human[$i+2]){
    unless($out_h[0]){$out_h[0]=$vstart;}
    if($vend <= $tem_human[$i+2]){$overlap_h=$vend-$vstart+1;$out_h[1]=$vend;}
    elsif($vend > $tem_human[$i+2]){$overlap_h=$tem_human[$i+2]-$vstart+1;$out_h[1]=$tem_human[$i+2];}
                                                                    }
   elsif($vstart > $tem_human[$i+2]){$overlap_h=0;$overlap_total_h=0;}
   $overlap_total_h+=$overlap_h;
   $type_h=$tem_human[0];
   unless($type_h){$type_h="n";}
   if($type_h_keep eq "P1" && $type_h ne "P1"){next;}
   elsif($type_h_keep ne "P1" && $type_h eq "P1"){$overlap_total_h_keep=$overlap_total_h;@out_h_keep=@out_h;$type_h_keep=$type_h;}
   elsif($overlap_total_h > $overlap_total_h_keep){$overlap_total_h_keep=$overlap_total_h;@out_h_keep=@out_h;$type_h_keep=$type_h;}
                               }                     # 4 end
                                }                    # 3 end
####################                                 # single end reads
  elsif(exists($h_single{$line[2]})){
   $overlap_h=0;$overlap_total_h=0;$type_h="N";
                                    }
  else{$overlap_h=0;$overlap_total_h=0;$type_h="N";}

####################
                                 } # 1 end, it it still can mapped to human;

####################                                # last step, print out;
 my $tem_length=$vend-$vstart+1;
 my $tem_direc="S";
  if($line[12]%256>=128){$tem_direc="R";}
  elsif($line[12]%256<128 && $line[12]%256>=64){$tem_direc="L";}
  else{$tem_direc="q";}
  print "b:@line $tem_direc $vstart $vend $vlength $type $out_v[0]|$out_v[1]|$overlap_total $type_h_keep $out_h_keep[0]|$out_h_keep[1]|$overlap_total_h_keep\n";
 @line=();$overlap=0;$vlength=0;$type="";$overlap_total=0;$overlap_total_h=0;@out_h_keep=();$overlap_total_h_keep=0;
                     }                              # the end of integration file;
