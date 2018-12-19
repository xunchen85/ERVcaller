#!/usr/bin/perl
use strict;

my $line="";
my @line=();
my %read=();
my $read_name="";
my $name="";
my $order=0;
my $order1=0;
my %human=();
my %virus=();

#####################
system ("sort -k 2 -k 3 -k 4 -k 5 $ARGV[0] |uniq | sort -k 2 >$ARGV[0]_2");

##################### Average Insert size;

open CANDIDATE,"$ARGV[0]_2";

my @tmp1 = (0) x 26;
print "@tmp1\n";
#####################
while(my $tmp1=<CANDIDATE>){
 $tmp1=~s/\n$//;
 chomp ($tmp1);
 @line=split /\s+/, $tmp1;
 if($line[2] / 512>=1){next;}
 $line[5]=$line[3];
 my $name2=$line[1]."_".$line[3]."_".$line[8];
 $name=$line[1]."_".$line[3]."_".$line[4];
 if(!(%read)){
   @{$read{$name}}=@line;
   $read_name=$line[1];
   unshift(@{$read{$name}},$line[11]);
   unshift(@{$read{$name}},$line[3]);
   unshift(@{$read{$name}},$line[0]);
             }
 if($line[1] eq $read_name){
  if(exists($read{$name2}) && (($line[0] eq "PF" && ${$read{$name2}}[0] eq "PR") || ($line[0] eq "PR" && ${$read{$name2}}[0] eq "PF")) && $line[7] eq "="){
   if((($line[0] eq "PF" && ${$read{$name2}}[0] eq "PR") || ($line[0] eq "PF" && ${$read{$name2}}[0] eq "PR")) && $line[3] eq ${$read{$name2}}[6] && $line[4] eq ${$read{$name2}}[11] && $line[8] eq ${$read{$name2}}[7]){
     splice @{$read{$name2}},3,0,@line;
     ${$read{$name2}}[0]="P";                                               
     ${$read{$name2}}[2]+=$line[11];                                        }
   else{
    if((${$read{$name2}}[6] ne ${$read{$name2}}[18])){
     $name=$line[1]."_".$line[3]."_".$line[4];
     @{$read{$name}}=@line;
     unshift(@{$read{$name}},$line[11]);
     unshift(@{$read{$name}},$line[3]);
     unshift(@{$read{$name}},$line[0]); 
                                                     } 
    else{
      splice @{$read{$name2}},3,0,@line; ${$read{$name2}}[0]="P2";
      ${$read{$name2}}[2]+=$line[11];
        }
       }
                                                                                                                                                          }
  else{
   if($order eq 1){next;}
   $name=$line[1]."_".$line[3]."_".$line[4];
   @{$read{$name}}=@line;
   unshift(@{$read{$name}},$line[11]);
   unshift(@{$read{$name}},$line[3]);
   unshift(@{$read{$name}},$line[0]);
      }
                           }
 else{
  my @order=sort{$read{$b}->[2]<=>$read{$a}->[2]} keys %read;
  for(my $j=0,my $max=${$read{$order[0]}}[2];$j<@order;$j++){
   ${$read{$order[$j]}}[1]="-";
   if(${$read{$order[$j]}}[2]<$max){last;}
   else{print "@{$read{$order[$j]}}\n";}
                                                            } 
  %read=();
  $name=$line[1]."_".$line[3]."_".$line[4];
  @{$read{$name}}=@line;
  $read_name=$line[1];
  unshift(@{$read{$name}},$line[11]);
  unshift(@{$read{$name}},$line[3]);
  unshift(@{$read{$name}},$line[0]);
    }
                  }
 if(%read){
  my @order=sort{$read{$b}->[2]<=>$read{$a}->[2]} keys %read;
  for(my $j=0,my $max=${$read{$order[0]}}[2];$j<@order;$j++){
   ${$read{$order[$j]}}[1]="-";
   if(${$read{$order[$j]}}[2]<$max){last;}
   else{print "@{$read{$order[$j]}}\n";}
                                                            }
          }
  system("rm $ARGV[0]_2");
