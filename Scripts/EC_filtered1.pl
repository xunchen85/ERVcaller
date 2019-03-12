#!usr/bin/perl
use strict;

my $line="";
my @line=();
my $type1="";
my $k=0;
my @len=();
my $as="";
my $md="";
my $posi="";

################
open F1,"$ARGV[0]";

while(<F1>){
 $as="";
 $md="";
 $posi="";
 my @line=split;
 if($line[0]=~"@"){next;}
 @len=($line[5]=~/(\d+)/g);
 my $len=0; map {$len+=$_} @len;
 $k=0;my $type1="";
 if($line[1]%8<4){ 
 if($line[1]%4>=2 && $line[1]%256>=128){$type1="PR";$k=1;}
 elsif($line[1]%4>=2 && $line[1]%256<128){$type1="PF";$k=1;}
 elsif($line[1]%256<128){$type1="F";$k=1;}
 elsif($line[1]%256>=128){$type1="R";$k=1;}
 else{$type1="N";$k=1;}

###############
 for(my $i=11;$i<@line;$i++){
  if($line[$i]=~"AS:i:"){
   $as=$line[$i];
   $as=~s/AS\:i\://;
                     }
  if($line[$i]=~"MD:Z:"){
   $md=$line[$i];
   $md=~s/MD\:Z\://; 
                     }
 my $hstart=0;
 my $hend=0;
 my $hlength=0;

 ### read mapping information
 my @temp1=($line[5]=~/(\d+)/g);
 my @temp2=($line[5]=~/([A-Z])/g);
 for(my $j=0;$j<@temp1;$j++){
  if($temp2[$j] eq "S"){   
   if($j==0){$hstart=$temp1[$j]+1;$hend=$temp1[$j];$hlength=$temp1[$j];}
   if($j==@temp1-1){$hlength+=$temp1[$j];}
                       }  
  elsif($temp2[$j] eq "D"){next;}
  elsif($temp2[$j] eq "M" || $temp2[$j] eq "I"){
   if($j==0){$hstart=0;$hend=$temp1[$j];$hlength=$temp1[$j];}
   elsif($j==@temp1-1){$hlength+=$temp1[$j];$hend+=$temp1[$j];}
   else{$hend+=$temp1[$j];$hlength+=$temp1[$j];}}
                            }
  $posi=$hstart."_".$hend."_".$hlength."_".$md;
                            }    ### each read
                 }
 else{$type1="U";next;}
# if($as<50){next;}
################################
 my $virus_name="";
  unless ($as) {$as=0;}
  print "$type1 @line[0..8] $posi $as\n";
}
