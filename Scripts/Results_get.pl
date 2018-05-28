#!usr/bin/perl
use strict;

my @VIP=();
my %VIP=();

my $line="";
my @line=();
my $group="";
my $group_keep="";
my $type=0;
############################ open the file;
open F1,"$ARGV[0]";
while(<F1>){
  @line=split;
  if($line[9]+$line[11]==0){$type=0;} #### type of VIPs;
  elsif($line[9]==0 && $line[11]>0){$type=0.5;}
  else{$type=1;}
  $group=$line[15]."_".$line[31]."_".$line[32]."_".$line[33];
  unless($group_keep){                #### first line;
   $group_keep=$group;
   @VIP=();
  if($line[3]+$line[9]>0 && $line[3]+$line[5]+$line[9]+$line[11]>0){
    push(@VIP,[@line,$line[3]+$line[9],$line[3]+$line[5]+$line[9]+$line[11],($line[2]+$line[8])/($line[3]+$line[9]),($line[2]+$line[4]+$line[8]+$line[10])/($line[3]+$line[5]+$line[9]+$line[11]),$type,$group]);
                                                                   }
  elsif($line[3]+$line[9]>0 && $line[3]+$line[5]+$line[9]+$line[11]==0){
push(@VIP,[@line,$line[3]+$line[9],$line[3]+$line[5]+$line[9]+$line[11],($line[2]+$line[8])/($line[3]+$line[9]),0,$type,$group]);
                                                                       }
  elsif($line[3]+$line[9]==0 && $line[3]+$line[5]+$line[9]+$line[11]>0){
   push(@VIP,[@line,$line[3]+$line[9],$line[3]+$line[5]+$line[9]+$line[11],0,($line[2]+$line[4]+$line[8]+$line[10])/($line[3]+$line[5]+$line[9]+$line[11]),$type,$group]);                                                   
                                                                       }
 else{push(@VIP,[@line,$line[3]+$line[9],$line[3]+$line[5]+$line[9]+$line[11],0,0,$type,$group]);}
                     }
 elsif($group eq $group_keep){
   if($line[3]+$line[9]>0 && $line[3]+$line[5]+$line[9]+$line[11]>0){
    push(@VIP,[@line,$line[3]+$line[9],$line[3]+$line[5]+$line[9]+$line[11],($line[2]+$line[8])/($line[3]+$line[9]),($line[2]+$line[4]+$line[8]+$line[10])/($line[3]+$line[5]+$line[9]+$line[11]),$type,$group]);
                                                                    }
  elsif($line[3]+$line[9]>0 && $line[3]+$line[5]+$line[9]+$line[11]==0){
push(@VIP,[@line,$line[3]+$line[9],$line[3]+$line[5]+$line[9]+$line[11],($line[2]+$line[8])/($line[3]+$line[9]),0,$type,$group]);
                                                                       }
  elsif($line[3]+$line[9]==0 && $line[3]+$line[5]+$line[9]+$line[11]>0){
   push(@VIP,[@line,$line[3]+$line[9],$line[3]+$line[5]+$line[9]+$line[11],0,($line[2]+$line[4]+$line[8]+$line[10])/($line[3]+$line[5]+$line[9]+$line[11]),$type,$group]);
                                                                       }
 else{push(@VIP,[@line,$line[3]+$line[9],$line[3]+$line[5]+$line[9]+$line[11],0,0,$type,$group]);}
 
                             }
 if(eof(F1) || $group ne $group_keep){                                #### if it is next VIP region;

################## How to pick candidate and summarize the results;
  my %viruses=();

############# version 1;
   my @VIP_sort = sort { $b->[38] <=> $a->[38] } @VIP;       # according to uniquely mapped reads; 

############################# keep best hit each virus first;
  for(my $i=0;$i<@VIP_sort;$i++){
   my $name=$VIP_sort[$i][30]."_".$VIP_sort[$i][34];
   if(exists($viruses{$name})){
    if($VIP_sort[$i][35] eq ${$viruses{$name}}[35]){
     if($VIP_sort[$i][38]==${$viruses{$name}}[38] && $VIP_sort[$i][40]>=${$viruses{$name}}[40]){     
      if((abs(${$viruses{$name}}[9]- ${$viruses{$name}}[3]) > abs($VIP_sort[$i][9]-$VIP_sort[$i][3])) || ($VIP_sort[$i][9]>0 && $VIP_sort[$i][12]>$VIP_sort[$i][7]) || ($VIP_sort[$i][40]>${$viruses{$name}}[40] && ${$viruses{$name}}[9]+${$viruses{$name}}[11]==0)){@{$viruses{$name}}=@{$VIP_sort[$i]};}
                                                                                               }
                                                   }
    elsif($VIP_sort[$i][38]>=${$viruses{$name}}[38] && $VIP_sort[$i][40]>=${$viruses{$name}}[40] && ((abs(${$viruses{$name}}[9] - ${$viruses{$name}}[3]) > abs($VIP_sort[$i][9]-$VIP_sort[$i][3])) || ( $VIP_sort[$i][9]>0 && $VIP_sort[$i][12]>$VIP_sort[$i][7]) || ($VIP_sort[$i][40]>${$viruses{$name}}[40] && ${$viruses{$name}}[9]+${$viruses{$name}}[11]==0))){@{$viruses{$name}}=@{$VIP_sort[$i]};}
                             }
   else{@{$viruses{$name}}=@{$VIP_sort[$i]};}
                                }

############################# keep best virus each VIP region
  my @names=keys %viruses;
  my @VIP_sort2=();
  my %viruses2=();
  for(my $i=0;$i<@names;$i++){push(@VIP_sort2,[@{$viruses{$names[$i]}}]);}
  my @VIP_sort3=sort { $b->[38] <=> $a->[38] } @VIP_sort2;
  @VIP_sort2=@VIP_sort3;
  for(my $i=0;$i<@VIP_sort2;$i++){
   my $name=$VIP_sort2[$i][43];
   if(exists($viruses2{$name})){
    if($VIP_sort2[$i][35] eq ${$viruses2{$name}}[35]){
     if($VIP_sort2[$i][38]==${$viruses2{$name}}[38] && $VIP_sort2[$i][40]>=${$viruses2{$name}}[40]){
      if((abs(${$viruses2{$name}}[9] - ${$viruses2{$name}}[3]) > abs($VIP_sort2[$i][9]-$VIP_sort2[$i][3])) || ($VIP_sort2[$i][9]>0 && $VIP_sort2[$i][12]>$VIP_sort2[$i][7]) || ($VIP_sort2[$i][40]>${$viruses2{$name}}[40]&& ${$viruses2{$name}}[9]+${$viruses2{$name}}[11]==0)){@{$viruses2{$name}}=@{$VIP_sort2[$i]};}
                                                                                                   }
                                                     }
    elsif($VIP_sort2[$i][38]>=${$viruses2{$name}}[38] && $VIP_sort2[$i][40]>=${$viruses2{$name}}[40] && ((abs(${$viruses2{$name}}[9] - ${$viruses2{$name}}[3]) > abs($VIP_sort2[$i][9]-$VIP_sort2[$i][3])) || ( $VIP_sort2[$i][9]>0 && $VIP_sort2[$i][12]>$VIP_sort2[$i][7])||($VIP_sort2[$i][40]>${$viruses2{$name}}[40]&& ${$viruses2{$name}}[9]+${$viruses2{$name}}[11]==0))){@{$viruses2{$name}}=@{$VIP_sort2[$i]};}
                               }  
   else{@{$viruses2{$name}}=@{$VIP_sort2[$i]};}
                                  }

#############################
 my @names=keys %viruses2;
 my @viruses_all=keys %viruses;
 my $list="";
 for(my $i=0;$i<@viruses_all;$i++){
   $list=$list.$viruses{$viruses_all[$i]}[38]."|".$viruses{$viruses_all[$i]}[39]."|".$viruses{$viruses_all[$i]}[40]."|".$viruses{$viruses_all[$i]}[41]."|".$viruses{$viruses_all[$i]}[42]."|".$viruses{$viruses_all[$i]}[34]."|".$viruses{$viruses_all[$i]}[35]."\$";}
 for(my $i=0;$i<@names;$i++){print "@{$viruses2{$names[$i]}} $list\n";}

################## next VIP region;
  $group_keep=$group;
  @VIP=();
   if($line[3]+$line[9]>0 && $line[3]+$line[5]+$line[9]+$line[11]>0){
    push(@VIP,[@line,$line[3]+$line[9],$line[3]+$line[5]+$line[9]+$line[11],($line[2]+$line[8])/($line[3]+$line[9]),($line[2]+$line[4]+$line[8]+$line[10])/($line[3]+$line[5]+$line[9]+$line[11]),$type,$group]);
                                                                    }
  elsif($line[3]+$line[9]>0 && $line[3]+$line[5]+$line[9]+$line[11]==0){
push(@VIP,[@line,$line[3]+$line[9],$line[3]+$line[5]+$line[9]+$line[11],($line[2]+$line[8])/($line[3]+$line[9]),0,$type,$group]);
                                                                       }
  elsif($line[3]+$line[9]==0 && $line[3]+$line[5]+$line[9]+$line[11]>0){
   push(@VIP,[@line,$line[3]+$line[9],$line[3]+$line[5]+$line[9]+$line[11],0,($line[2]+$line[4]+$line[8]+$line[10])/($line[3]+$line[5]+$line[9]+$line[11]),$type,$group]);
                                                                       }
 else{push(@VIP,[@line,$line[3]+$line[9],$line[3]+$line[5]+$line[9]+$line[11],0,0,$type,$group]);}
     }
         }
