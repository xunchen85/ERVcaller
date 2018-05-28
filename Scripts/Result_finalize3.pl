#!usr/bin/perl
use strict;

my $line="";
my @line=();
my @o1=();
my @o2=();
my @o3_h=();
my @o3_v=();
my @vip_title=();
my @vip_list=();
my @vip_list2=();
my $order=0;
my $vstart=0;
my $hstart=0;
my $vend=0;
my $hend=0;
my @reads=();
my %single=();
my %multiple=();

open INPUT,"$ARGV[0]";
my $seq_len=$ARGV[1];
my $seq_blank1="-" x $seq_len;
my $seq_blank="";
$line=<INPUT>;
while(<INPUT>){
 @line=split;
 if($line[0] eq "O2"){splice @line,1,1;}
#############################
XUN: if($line[1] =~ "End"){                             # organize results and print out the results.
  $order++;
  @vip_list=();
  @vip_list2=();
  my @vip_list_multiple=();
#  print "order:$order\n";
  for(my $j=0;$j<@o2;$j++){                           # single mapped reads;              and also skill first line
   $o2[$j][6]=$o1[$j+1][4];
   $o2[$j][8]=$o1[$j+1][12];
   my $endd=$o2[$j][6]+$o2[$j][9]-$o2[$j][5];
#   print "end:$o2[$j][6] $endd\n";
   unless(@vip_list){push(@vip_list,[$o2[$j][6],$endd,$o1[$j+1][19],$o1[$j+1][1],@{$o2[$j]}[1..12]]);}
   else{
#      print "@{$vip_list[0]}\n";    
#      print "$#vip_list\n";
#      print "@vip_list\n";
##    for(my $k=0;$k<@vip_list;$k++){
#     print "vip_list:\t@{$vip_list[$k]}\n";      
     if($o2[$j][6]>$vip_list[$#vip_list][1]){
      if($o2[$j][6]>$endd){push(@vip_list,[$endd,$o2[$j][6],$o1[$j+1][19],$o1[$j+1][1],@{$o2[$j]}[1..12]]);}
      else{push(@vip_list,[$o2[$j][6],$endd,$o1[$j+1][19],$o1[$j+1][1],@{$o2[$j]}[1..12]]);}
                                            }
     else{
      if($o2[$j][6]<$vip_list[$#vip_list][0]){$vip_list[$#vip_list][0]=$o2[$j][6];}
      if($endd>$vip_list[$#vip_list][1]){$vip_list[$#vip_list][1]=$endd;}
      if($endd<$vip_list[$#vip_list][0]){$vip_list[$#vip_list][1]=$vip_list[$#vip_list][0];$vip_list[$#vip_list][0]=$endd;}
      else{$vip_list[$#vip_list][1]=$endd;}
      push(@{$vip_list[$#vip_list]},($o1[$j+1][19],$o1[$j+1][1],@{$o2[$j]}[1..12]));
         }
#      print "$#vip_list\n";
#                                  }
#      if($endd<$vip_list[$k][0]){$vip_list[$k][1]=$vip_list[$k][0];$vip_list[$k][0]=$endd;}
#      last;
#                                            }
#     else{
#      if($endd<$vip_list[0][0] || ($endd<$vip_list[$k+1][0] && $o2[$j][6]<$vip_list[$k][0])){splice @vip_list,$k,0,[$o2[$j][6],$endd,$o1[$j+1][19],$o1[$j+1][1],@{$o2[$j]}[1..12]];}
#      else{
#       if($vip_list[$k][1]<$endd){$vip_list[$k][1]=$endd;}
#       if($vip_list[$k][0]>$o2[$j][6]){$vip_list[$k][1]=$o2[$j][6];}
#       push(@{$vip_list[$k]},($o1[$j+1][19],$o1[$j+1][1],@{$o2[$j]}[1..12]));
#          }
#         }
##                                  }
       }
                        }
#  print "vip_list:\t@vip_list\n";
#  for(my $dd=0;$dd<@vip_list;$dd++){print "vip_list2:\t@{$vip_list[$dd]}\n";}
  for(my $dd=0;$dd<@vip_list;$dd++){
#   print "vip_list2:\t@{$vip_list[$dd]}\n";
   my @temp11=@{$vip_list[$dd]};
   my @temp22=();
#   print "start:\t@temp11[0..1]\n";
   for(my $dd2=2;$dd2<@temp11;$dd2+=14){       #single
    @temp22=@temp11[$dd2..$dd2+13];
    my $endd2=$temp22[9]+$temp22[12]-$temp22[8];
    if(exists($single{$temp22[1]})){
     unless($vip_list2[$dd]){$vip_list2[$dd][0]=$temp22[9],$vip_list2[$dd][1]=$endd2;}
     else{
      if($endd2<$vip_list2[$dd][1]){next;}
      else{$vip_list2[$dd][1]=$endd2;}
 
      if($temp22[9]>$vip_list2[$dd][0]){next;}
      else{$vip_list2[$dd][0]=$temp22[9];}
         }
                                  }
                                       }
#  for(my $dd3=0;$dd3<@vip_list2;$dd3++){
#   print "a:\t$vip_list2[$dd][0] $vip_list2[$dd][1]\n";
#   print "b:\t$vip_list2[$dd+1][0] $vip_list2[$dd+1][1]\n";
#                                       }
  my $temp_name=""; my $temp_dis=0;my @temp_array=();
  for(my $dd2=2;$dd2<@temp11;$dd2+=14){       #multiple
    @temp22=@temp11[$dd2..$dd2+13];
#    print "temp22:\t@temp22\n";
    my $endd2=$temp22[9]+$temp22[12]-$temp22[8];
    if(exists($multiple{$temp22[1]})){
     unless($vip_list2[$dd]){
      $vip_list2[$dd][0]=$temp22[9];
      $vip_list2[$dd][1]=$endd2;
#      print "vip2:\t@{$vip_list2[$dd]}\n";
                            }
     else{
#      print "vip3:\t@{$vip_list2[$dd]}\n";
#      print "temm:\t@temp_array ddd\n";
      unless(@temp_array){
       if($endd2<$vip_list2[$dd][1]){$temp_array[1]=$vip_list2[$dd][1];}
       else{$temp_array[1]=$endd2;}
 
       if($temp22[9]>$vip_list2[$dd][0]){$temp_array[0]=$vip_list2[$dd][0];}
       else{$temp_array[0]=$temp22[9];}

#       push(@temp_array,($temp22[9],$endd2));
       $temp_dis=$endd2-$vip_list2[$dd][1];
#       print "empty:\t@temp_array\n";
                         }
      elsif($endd2>=$vip_list2[$dd][1] && $temp22[9]>=$vip_list2[$dd][0]){
       if($temp_name eq $temp22[1] && $temp_dis>$endd2-$vip_list2[$dd][1]){
        @temp_array=();push(@temp_array,($vip_list2[$dd][0],$endd2));$temp_dis=$$endd2-$vip_list2[$dd][1];}
       if($temp_name ne $temp22[1] || $dd2+14 >= @temp11){
#        print "temp2:\t@temp_array\n";
        @{$vip_list2[$dd]}=@temp_array;
        @temp_array=();
        push(@temp_array,($vip_list2[$dd][0],$endd2));
        $temp_dis=$endd2-$vip_list2[$dd][1];
                                                        }
                                                                        }
      elsif($endd2<=$vip_list2[$dd][1] && $temp22[9]<=$vip_list2[$dd][0]){
       if($temp_name eq $temp22[1] && $temp_dis>$endd2-$vip_list2[$dd][1]){
        @temp_array=();push(@temp_array,($temp22[9],$vip_list2[$dd][1]));$temp_dis=$$endd2-$vip_list2[$dd][1];
                                                                          }
       elsif($temp_name ne $temp22[1] || $dd2+14 >= @temp11){
#        print "temp3:\t@temp_array\n";
        @{$vip_list2[$dd]}=@temp_array;
        @temp_array=();
        push(@temp_array,($temp22[9],$vip_list2[$dd][1]));
        $temp_dis=$endd2-$vip_list2[$dd][1];
#        print "temp_array2:\t@temp_array\n";
                                                           }           }
#      print "temp_array:\t@temp_array\n";
#      print "vip_list:@{$vip_list2[$dd]}\n************************************\n";
                                                                        
         }
                                      }        # if exist multiple
                                     }         # multiple end
############################
   my @left=(0) x 2;
   my @right=(0) x 2;
   my $rank=-1;my @left_s=(0) x 2;my @left_m=(0) x 2;my @right_s=(0) x 2;my @right_m=(0) x 2;
   my @keep_left=(0) x 4;my @keep_right=(0) x 4; my $keep_rank=0;
    for(my $dd3=0;$dd3<@vip_list2;$dd3++){
        my %vip_list11=();
        my %vip_list21=();
        for(my $id=0;$id<=$dd3;$id++){
         my @temp3=@{$vip_list[$id]};
#         print "temp31\ta:@temp3\t@{$o1[0]}\n";
         for(my $id3=2;$id3<@temp3;$id3+=14){$vip_list11{$temp3[$id3+1]}=$temp3[$id3];}
         if($left[0] eq 0 || $vip_list2[$id][0]<$left[0]){$left[0]=$vip_list2[$id][0];}
         if($left[1] eq 0 || $vip_list2[$id][1]>$left[1]){$left[1]=$vip_list2[$id][1];}
                                     }
        for(my $id2=$dd3+1;$id2<@vip_list2;$id2++){
         my @temp3=@{$vip_list[$id2]};
#         print "temp32\ta:@temp3\t@{$o1[0]}\n";
         for(my $id4=2;$id4<@temp3;$id4+=14){$vip_list21{$temp3[$id4+1]}=$temp3[$id4];}
         if($right[0] eq 0 || $vip_list2[$id2][0]<$right[0]){$right[0]=$vip_list2[$id2][0];}
         if($right[1] eq 0 || $vip_list2[$id2][1]>$right[1]){$right[1]=$vip_list2[$id2][1];}
                                                  }
      if($#vip_list2+1-$dd3>$dd3){$rank=$dd3+1;}
      elsif($#vip_list2-$dd3<=$dd3){$rank=$#vip_list2+1-$dd3;}

####################### calculate the quality score;
        my @temp_left=keys %vip_list11;@left_s=(0) x 2;@left_m=(0) x 2;
        for(my $id5=0;$id5<@temp_left;$id5++){
         if(exists($single{$temp_left[$id5]})){$left_s[0]+=$vip_list11{$temp_left[$id5]};$left_s[1]++;}
         else{$left_m[0]+=$vip_list11{$temp_left[$id5]};$left_m[1]++;}
                                             }
        my @temp_right=keys %vip_list21; @right_s=(0) x 2;@right_m=(0) x 2;
        for(my $id5=0;$id5<@temp_right;$id5++){
         if(exists($single{$temp_right[$id5]})){$right_s[0]+=$vip_list21{$temp_right[$id5]};$right_s[1]++;}
         else{$right_m[0]+=$vip_list21{$temp_right[$id5]};$right_m[1]++;}
                                              }
       print "ad:$order\t$rank\t@left_s\t@left_m\t@left\t@right_s\t@right_m\t@right\t@{$o1[0]}\n";
                                          }
#     print "bid:$rank\t@left_s\t@left_m\t@left\t@right_s\t@right_m\t@right\t@{$o1[0]}\n";
                                   }
  
##############
  $hstart=0;$hend=0;$vstart=0;$vend=0;
  @o3_h=(); @o2=(); @o1=(); @o3_v=();
  %single=();%multiple=();
                          }
############################# first line including vip information;
 if($line[0] eq "O1"){                                 # O1 group;
# push (@o1,[@line[0..23]]);
  push (@o1,[$line[3],@line[1..2],@line[4..24]]);                          ### updated on 4/14/16
  @reads=split /\|/, $line[17];                                            ### updated on 4/14/16
  for(my $i=0;$i<@reads;$i+=2){
   if($reads[$i+1] eq 1){$single{$reads[$i]}=$reads[$i+1];}
   else{$multiple{$reads[$i]}=$reads[$i+1];}
                              }
  for(my $i=25;$i<@line;$i+=31){push(@o1,[@line[$i..$i+30]]);}             ### updated on 4/14/16
                     }
############################# O2 group for each vip;
 elsif($line[0] eq "O2"){
  unless($hstart){$hstart=$line[17];}
  unless($vstart){$vstart=$line[18];}
  push (@o2,[@line]);
                        }
############################# O3 group for each vip;
 elsif($line[0] eq "O3" && $line[1] eq "//Start//"){
  my $seq_type=0;
  my $read_order=0;
  while(<INPUT>){
   @line=split;
   $line[0]=~s/Seq//;
   if($line[1] =~ "End"){
#     print "\naaaaaaaaaaaaaaaaaaaaaaaaaaa\n\n";
     goto XUN;
                        }
   if($line[2] eq $seq_len || !($_=~ "[a-zA-Z0-9]")){next;}
   if($_=~"human"){$seq_type=1;$read_order=0;}
   elsif($_=~"virus"){$seq_type=2;$read_order=0;}
   elsif($seq_type==1){
     unless($o3_h[$line[0]]){if($line[1]>$hstart){$seq_blank="=" x ($line[1]-$hstart);$o3_h[$line[0]]=$seq_blank.$line[2];$hend=$line[3]+1}else{$o3_h[$line[0]]=$line[2];$hend=$line[3]+1;}}
     elsif($line[1]>$hend){my $break="=" x ($line[1]-$hend);$o3_h[$line[0]]=$o3_h[$line[0]].$break.$line[2];}
     elsif($line[1] eq $hend){$o3_h[$line[0]]=$o3_h[$line[0]].$line[2];}
     $read_order++;
     if($read_order eq @o2){$hend=$line[3]+1;$read_order=0;}
                      }
   elsif($seq_type==2){
     unless($o3_v[$line[0]]){if($line[1]>$vstart){$seq_blank="=" x ($line[1]-$vstart);$o3_v[$line[0]]=$seq_blank.$line[2];$vend=$line[3]+1}else{$o3_v[$line[0]]=$line[2];$vend=$line[3]+1;}}
     elsif($line[1]>$vend){my $break="=" x ($line[1]-$vend);$o3_v[$line[0]]=$o3_v[$line[0]].$break.$line[2];}
     elsif($line[1] eq $vend){$o3_v[$line[0]]=$o3_v[$line[0]].$line[2];}
     $read_order++;
     if($read_order eq @o2){$vend=$line[3]+1;$read_order=0;}
                      }
                }
                                                   }
              }
