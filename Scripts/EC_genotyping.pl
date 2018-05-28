#!usr/bin/perl
use strict;

my $line="";
my @line=();
my %VI_list=();
my @names=();
my %split=();
my %split3=();
my %split2=();

##########################  Cellular_proportion.results
open BP,"$ARGV[0]";            
#open BP2, ">>$ARGV[0].bp2-2";
while(<BP>){
 @line=split;
 my $name=$ARGV[1]."_".$ARGV[2]."_".$ARGV[3];
 unshift (@line,$name);
 if(exists($VI_list{$line[0]."|".$line[5]}) && ${$VI_list{$line[0]."|".$line[5]}}[0] eq $line[0] && ${$VI_list{$line[0]."|".$line[5]}}[7] eq $line[7] && ${$VI_list{$line[0]."|".$line[5]}}[8] eq $line[12] && ${$VI_list{$line[0]."|".$line[5]}}[12] eq $line[8]){
  ${$VI_list{$line[0]."|".$line[5]}}[1]="P2";
  ${$VI_list{$line[0]."|".$line[5]}}[3]=${$VI_list{$line[0]."|".$line[5]}}[3]+$line[3];
  push(@{$VI_list{$line[0]."|".$line[5]}},@line[4..$#line]);
   }
 else {push(@{$VI_list{$line[0]."|".$line[5]}},@line);}
}

########################## extract PE reads and count;
my %vi=();
my @reads=keys %VI_list;
for(my $i=0;$i<@reads;$i++){
 unless(exists($vi{$line[0]})){@{$vi{$line[0]}}=(0) x 4;}
 @line=@{$VI_list{$reads[$i]}};
 @names=split /\_/, $line[0];
 if($line[1] ne "P" && $line[1] ne "P2"){next;}

########## extract out split reads;
 my $start1=0;
 my $end=0;
 my $length=0;
 my $start_1=0;
 my $end_1=0;
 my $len_1=0;
 my $start_2=0;
 my $end_2=0;
 my $len_2=0;
 if($line[10] =~ "S" || $line[22] =~ "S"){           ##### split start

##### first end
  my @temp1=($line[10]=~/(\d+)/g);
  my @temp2=($line[10]=~/([A-Z])/g);
  for(my $i=0;$i<@temp1;$i++){
  if($temp2[$i] eq "S"){   
   if($i==0){$start1=$temp1[$i]+1;$end=$temp1[$i];$length=$temp1[$i];}
   if($i==@temp1-1){$length+=$temp1[$i];}
                       }
  elsif($temp2[$i] eq "D"){next;}
  elsif($temp2[$i] eq "M" || $temp2[$i] eq "I"){
   if($i==0){$start1=1;$end=$temp1[$i];$length=$temp1[$i];}
   elsif($i==@temp1-1){$length+=$temp1[$i];$end+=$temp1[$i];}
   else{$end+=$temp1[$i];$length+=$temp1[$i]}  }
                              }
  $start_1=$start1;$end_1=$end;$len_1=$length;

###### second end
  my @temp1=($line[22]=~/(\d+)/g);
  my @temp2=($line[22]=~/([A-Z])/g);
  for(my $i=0;$i<@temp1;$i++){
  if($temp2[$i] eq "S"){     
   if($i==0){$start1=$temp1[$i]+1;$end=$temp1[$i];$length=$temp1[$i];}
   if($i==@temp1-1){$length+=$temp1[$i];}
                       }  
  elsif($temp2[$i] eq "D"){next;}
  elsif($temp2[$i] eq "M" || $temp2[$i] eq "I"){
   if($i==0){$start1=1;$end=$temp1[$i];$length=$temp1[$i];}
   elsif($i==@temp1-1){$length+=$temp1[$i];$end+=$temp1[$i];}
   else{$end+=$temp1[$i];$length+=$temp1[$i]}  }
                              }
  $start_2=$start1;$end_2=$end;$len_2=$length;
                                           }         ##### split end

########## identify split >=20, and split <20 bp
 if($names[2] eq $line[8] && $line[10] =~ "S"){
  if($start_1<20){
   $line[1]=$line[1]."_F_".$len_1."_".$start_1."_".$end_1;
   if(exists($split2{$line[0]})){$split2{$line[0]}++;}
   else{$split2{$line[0]}=1;}
                 }
  else{
   $line[1]=$line[1]."_F_".$len_1."_".$start_1."_".$end_1;
   if(exists($split3{$line[0]})){$split3{$line[0]}++;}
   else{$split3{$line[0]}=1;}         
      }
  push(@{$split{$line[0]."|".$line[5]}},@line);
                                              }
 elsif($names[2] eq $line[20] && $line[22] =~ "S"){
  if($start_2<20){
   $line[1]=$line[1]."_R_".$len_2."_".$start_2."_".$end_2;
   if(exists($split2{$line[0]})){$split2{$line[0]}++;}
   else{$split2{$line[0]}=1;}                
                 }
  else{
   $line[1]=$line[1]."_R_".$len_2."_".$start_2."_".$end_2;
   if(exists($split3{$line[0]})){$split3{$line[0]}++;}
   else{$split3{$line[0]}=1;}               
      }
  push(@{$split{$line[0]."|".$line[5]}},@line);
                                                  }

 elsif(($names[2] - ($line[8]+$end_1))<=5 && $line[10] =~ "S"){
  if($len_1-$end_1<20){
   $line[1]=$line[1]."_F2_".$len_1."_".$start_1."_".$end_1;
   if(exists($split2{$line[0]})){$split2{$line[0]}++;}
   else{$split2{$line[0]}=1;}                
                      }
  else{
   $line[1]=$line[1]."_F2_".$len_1."_".$start_1."_".$end_1;
   if(exists($split3{$line[0]})){$split3{$line[0]}++;}
   else{$split3{$line[0]}=1;}               
      }
  push(@{$split{$line[0]."|".$line[5]}},@line);
                                                              }
 elsif(($names[2] - ($line[20]+$end_2))<=5 && $line[22] =~ "S"){
  if($len_2-$end_2<20){
   $line[1]=$line[1]."_R2_".$len_2."_".$start_2."_".$end_2;
   if(exists($split2{$line[0]})){$split2{$line[0]}++;}
   else{$split2{$line[0]}=1;}                 
                 }
  else{
   $line[1]=$line[1]."_R2_".$len_2."_".$start_2."_".$end_2;
   if(exists($split3{$line[0]})){$split3{$line[0]}++;}
   else{$split3{$line[0]}=1;}                 
      }
  push(@{$split{$line[0]."|".$line[5]}},@line);
                                                               }
##########
 if(exists($vi{$line[0]}) && !(exists($split{$line[0]."|".$line[5]}))){
  my @temp1=($line[10]=~/(\d+)/g);
  my @temp2=($line[10]=~/([A-Z])/g);
  my $hstart="";
  my $hend="";
  my $hlength="";

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
   if (abs(($hend+$line[8]-1)-$ARGV[3])<=10 && $line[8]>=$line[20]){next;}
   elsif (abs(($hstart+$line[8]-1)-$ARGV[3])<=10 && $line[8]<=$line[20]){next;}

#####
  my @temp1=($line[22]=~/(\d+)/g);
  my @temp2=($line[22]=~/([A-Z])/g);
  my $hstart="";
  my $hend="";
  my $hlength="";

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
   if (abs(($hstart+$line[20]-1)-$ARGV[3])<=10 && $line[20]<$line[8]) {next;}
   elsif (abs(($hend+$line[20]-1)-$ARGV[3])<=10 && $line[8]<$line[20]){next;}
   ${$vi{$line[0]}}[0]++;
  }                  ##### 1. here counting number of reads support Non-VI;
                           }

######################### number of chimeric and split reads and cellular proportion of the VI;
 
my $name=$ARGV[1]."_".$ARGV[2]."_".$ARGV[3];                     # sample_id, chr, position
 unless(exists($split2{$name})){$split2{$name}=0;}                # if there is no split reads
 unless(exists($split3{$name})){$split3{$name}=0;}                # if there is no suggestive reads
 unless($vi{$name}){@{$vi{$name}}=(0) x 3;${$vi{$name}}[2]=1;}    # 0. if there is no reads support non-VI

############ array structure
 ${$vi{$name}}[1]=$ARGV[4];                                       # 1. chimeric & split reads 
 ${$vi{$name}}[3]=$split3{$name};                                 # 3. suggestive reads (split reads <20 bp)
 if(${$vi{$name}}[0] eq 0){                                       # if there is no reads support nonVIs
  ${$vi{$name}}[2]=1;
 }
############ cellular proportion calculation
 elsif(${$vi{$name}}[0] >0 && ${$vi{$name}}[1]>0){
  if($ARGV[5] eq 2) {
   ${$vi{$name}}[2]=(${$vi{$name}}[1]/2)/((${$vi{$name}}[1]/2)+${$vi{$name}}[0]);
  }
  elsif ($ARGV[5] eq 1) {
   ${$vi{$name}}[2]=${$vi{$name}}[1]/(${$vi{$name}}[1]+${$vi{$name}}[0]);
  }
 }

############ formating the results
print "ID_Chr_Position No._reads_support_nonVI No._reads_support_VI Allele_Fraction Suggestive_reads_(split_reads_<20bp)\n";
print "$name @{$vi{$name}}\n";                  # print out each VI

