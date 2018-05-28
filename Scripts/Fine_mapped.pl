#!usr/bin/perl
use strict;

my $line="";
my @line=();
my @array=();
my @ff=(0) x 8;
my @fr=(0) x 8;
my @rf=(0) x 8;
my @rr=(0) x 8;
my @O1=();
my $order=0;

##############
open F1,"$ARGV[0]";
while(<F1>){
 @line=split;
 if($line[0] eq "O1"){
  @O1=@line;
  $order=1;
 }
 elsif($line[0] eq "O2" && $order eq 1){

############## fr
  if($line[3] eq "F" && $line[5] eq "R"){
   if($line[2] eq "S1"){
    $rr[7]++;
    $rr[0]=$line[7]+$line[6];
    $rr[3]=$line[9]+$line[8];
   if($rr[1] >0 || $rr[2]>0){
    if($line[7]+$line[6] <= $rr[1]){$rr[1]=$line[7]+$line[6];}
    if($line[7]+$line[6]+($line[10]-$line[6]+1) >= $rr[2]){$rr[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
                            }
   else{$rr[1]=$line[7]+$line[6];$rr[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
   if($rr[4] >0 || $rr[5]>0){
    if($line[9]+$line[8] <= $rr[4]){$rr[4]=$line[9]+$line[8];}
    if($line[9]+$line[8]+($line[12]-$line[8]+1) >= $rr[5]){$rr[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
                            }
   else{$rr[4]=$line[9]+$line[8];$rr[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
                       }
   else{
    $fr[6]++;
    if($fr[1] >0 || $fr[2]>0){
     if($line[7]+$line[6] <= $fr[1]){$fr[1]=$line[7]+$line[6];}
     if($line[7]+$line[6]+($line[10]-$line[6]+1) >= $fr[2]){$fr[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
                             }
    else{$fr[1]=$line[7]+$line[6];$fr[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
    if($fr[4] >0 || $fr[5]>0){
     if($line[9]+$line[8] <= $fr[4]){$fr[4]=$line[9]+$line[8];}
     if($line[9]+$line[8]+($line[12]-$line[8]+1) >= $fr[5]){$fr[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
                             }
   else{$fr[4]=$line[9]+$line[8];$fr[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
       }
                                      }
############## ff
  elsif($line[3] eq "F" && $line[5] eq "F"){
   if($line[2] eq "S1"){
    $rf[7]++;
    $rf[0]=$line[7]+$line[6];
    $rf[3]=$line[9]+$line[8]+($line[12]-$line[8]+1);
   if($rf[1] >0 || $rf[2]>0){
    if($line[7]+$line[6] <= $rf[1]){$rf[1]=$line[7]+$line[6];}
    if($line[7]+$line[6]+($line[10]-$line[6]+1) >= $rf[2]){$rf[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
                            }
   else{$rf[1]=$line[7]+$line[6];$rf[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
   if($rf[4] >0 || $rf[5]>0){
    if($line[9]+$line[8] <= $rf[4]){$rf[4]=$line[9]+$line[8];}
    if($line[9]+$line[8]+($line[12]-$line[8]+1) >= $rf[5]){$rf[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
                            }
   else{$rf[4]=$line[9]+$line[8];$rf[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
                       }
   else{
    $ff[6]++;
   if($ff[1] >0 || $ff[2]>0){
    if($line[7]+$line[6] <= $ff[1]){$ff[1]=$line[7]+$line[6];}
    if($line[7]+$line[6]+($line[10]-$line[6]+1) >= $ff[2]){$ff[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
                            }
   else{$ff[1]=$line[7]+$line[6];$ff[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
   if($ff[4] >0 || $ff[5]>0){
    if($line[9]+$line[8] <= $ff[4]){$ff[4]=$line[9]+$line[8];}
    if($line[9]+$line[8]+($line[12]-$line[8]+1) >= $ff[5]){$ff[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
                            }
   else{$ff[4]=$line[9]+$line[8];$ff[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
       }
                                           }
############## rf
  elsif($line[3] eq "R" && $line[5] eq "F"){
   if($line[2] eq "S1"){
    $fr[7]++;
    $fr[0]=$line[7]+$line[6]+($line[10]-$line[6]+1);
    $fr[3]=$line[9]+$line[8];
   if($fr[1] >0 || $ff[2]>0){
    if($line[7]+$line[6] <= $fr[1]){$fr[1]=$line[7]+$line[6];}
    if($line[7]+$line[6]+($line[10]-$line[6]+1) >= $fr[2]){$fr[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
                            }
   else{$fr[1]=$line[7]+$line[6];$fr[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
   if($fr[4] >0 || $fr[5]>0){
    if($line[9]+$line[8] <= $fr[4]){$fr[4]=$line[9]+$line[8];}
    if($line[9]+$line[8]+($line[12]-$line[8]+1) >= $fr[5]){$fr[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
                            }
   else{$fr[4]=$line[9]+$line[8];$fr[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
                      }
   else{
    $rf[6]++;
   if($rf[1] >0 || $rf[2]>0){
    if($line[7]+$line[6] <= $rf[1]){$rf[1]=$line[7]+$line[6];}
    if($line[7]+$line[6]+($line[10]-$line[6]+1) >= $rf[2]){$rf[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
                            }
   else{$rf[1]=$line[7]+$line[6];$rf[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
   if($rf[4] >0 || $rf[5]>0){
    if($line[9]+$line[8] <= $rf[4]){$rf[4]=$line[9]+$line[8];}
    if($line[9]+$line[8]+($line[12]-$line[8]+1) >= $rf[5]){$rf[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
                            }
   else{$rf[4]=$line[9]+$line[8];$rf[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
       }
                                          }
############## rr
  elsif($line[3] eq "R" && $line[5] eq "R"){
   if($line[2] eq "S1"){
    $ff[7]++;
    $ff[0]=$line[7]+$line[6]+($line[10]-$line[6]+1);
    $ff[3]=$line[9]+$line[8]+($line[12]-$line[8]+1);
    if($ff[1] >0 || $ff[2]>0){
     if($line[7]+$line[6] <= $ff[1]){$ff[1]=$line[7]+$line[6];}
     if($line[7]+$line[6]+($line[10]-$line[6]+1) >= $ff[2]){$ff[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
                             }
    else{$ff[1]=$line[7]+$line[6];$ff[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
    if($ff[4] >0 || $ff[5]>0){
     if($line[9]+$line[8] <= $ff[4]){$ff[4]=$line[9]+$line[8];}
     if($line[9]+$line[8]+($line[12]-$line[8]+1) >= $ff[5]){$ff[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
                             }
   else{$ff[4]=$line[9]+$line[8];$ff[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
                       }
   else{
    $rr[6]++;
   if($rr[1] >0 || $rr[2]>0){
    if($line[7]+$line[6] <= $rr[1]){$rr[1]=$line[7]+$line[6];}
    if($line[7]+$line[6]+($line[10]-$line[6]+1) >= $rr[2]){$rr[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
                            }
   else{$rr[1]=$line[7]+$line[6];$rr[2]=$line[7]+$line[6]+($line[10]-$line[6]+1);}
   if($rr[4] >0 || $rr[5]>0){
    if($line[9]+$line[8] <= $rr[4]){$rr[4]=$line[9]+$line[8];}
    if($line[9]+$line[8]+($line[12]-$line[8]+1) >= $rr[5]){$rr[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
                            }
   else{$rr[4]=$line[9]+$line[8];$rr[5]=$line[9]+$line[8]+($line[12]-$line[8]+1);}
       }
                                         }
#################################
  } elsif ($order eq 1 && $line[0] eq "O3" && $line[1] eq "//End//"){

 my $number_breakpoints=0;
 my $breakpoints_in_region=0;
 my $number_sites=0;
 my $types="";
 my %overlap=();

 if($ff[6]>0 || $ff[7]>0){
  push(@{$overlap{"ff"}},@ff);
  $number_sites++;
  $types=$types."ff|";
  if($ff[7]>0){
   $number_breakpoints++;
   if(($ff[0]>$ff[1] && $ff[0]<$ff[2]) || ($ff[3]>$ff[4] && $ff[3]<$ff[5])){$breakpoints_in_region++;}
              }
 }
 if($fr[6]>0 || $fr[7]>0){
  push(@{$overlap{"fr"}},@fr);
  $number_sites++;
  $types=$types."fr|";
  if($fr[7]>0){
   $number_breakpoints++;
   if(($fr[0]>$fr[1] && $fr[0]<$fr[2]) || ($fr[3]>$fr[4] && $fr[3]<$fr[5])){$breakpoints_in_region++;}
              }
 }
 if($rf[6]>0 || $rf[7]>0){
  push(@{$overlap{"rf"}},@rf);
  $number_sites++;
  $types=$types."rf|";
  if($rf[7]>0){
   $number_breakpoints++;
   if(($rf[0]>$rf[1] && $rf[0]<$rf[2]) || ($rf[3]>$rf[4] && $rf[3]<$rf[5])){$breakpoints_in_region++;}
              }
 }
 if($rr[6]>0 || $rr[7]>0){
  push(@{$overlap{"rr"}},@rr);
  $number_sites++;
  $types=$types."rr|";
  if($rr[7]>0){
   $number_breakpoints++;
   if(($rr[0]>$rr[1] && $rr[0]<$rr[2]) || ($rr[3]>$rr[4] && $rr[3]<$rr[5])){$breakpoints_in_region++;}
              }
 }

 ##################
 for(my $i=0;$i<6;$i++){
  unless($ff[$i]){$ff[$i]="-";}
  unless($fr[$i]){$fr[$i]="-";}
  unless($rf[$i]){$rf[$i]="-";}
  unless($rr[$i]){$rr[$i]="-";}
 }

 my $overlap_h=0;
 my $overlap_v=0;
 my @temp1=sort{$overlap{$a}->[1]<=>$overlap{$b}->[1]} keys %overlap;
 for(my $i=1;$i<@temp1;$i++){
   if(${$overlap{$temp1[$i]}}[1] < ${$overlap{$temp1[$i-1]}}[2]){$overlap_h++;}
 }
 my @temp1=sort{$overlap{$a}->[4]<=>$overlap{$b}->[4]} keys %overlap;
 for(my $i=1;$i<@temp1;$i++){
   unless($temp1[$i]){next;}
   if(${$overlap{$temp1[$i]}}[4] < ${$overlap{$temp1[$i-1]}}[5]){$overlap_v++;}
 }
 print "$O1[1] $O1[18] $O1[19] $O1[20] $O1[21] $O1[22] ff: @ff fr: @fr rf: @rf rr: @rr summary: $types $number_sites $overlap_h $overlap_v $number_breakpoints $breakpoints_in_region\n";
 @ff=(0) x 8;
 @fr=(0) x 8;
 @rf=(0) x 8;
 @rr=(0) x 8;
 $types="";
 $number_sites="";
 $overlap_v="";
 $overlap_h="";
 $number_breakpoints="";
 $breakpoints_in_region="";
}
}

