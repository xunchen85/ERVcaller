#!usr/bin/perl
use strict;
my $line="";
my @line=();
my $seq_name="";
my $seq="";
my %chimeric_reads=();
my %split_reads=();
my @start=();
my @end=();
my @seq=();
my @seq_final=();
my $seq1="";
my $seq2="";
my $region=$ARGV[3];
my $directory="";

########################### open chimeric reads;
open SF,"$ARGV[1]_1.1fuq";
open SR,"$ARGV[1]_2.1fuq";
my $seq_name="";
my $seq="";
while(<SF>){
 $seq_name=$_;
 $seq_name=~s/\/1//;
 $seq_name=~s/^@//;
 chomp($seq_name);
 $seq=<SF>;
 chomp($seq);
 $line=<SF>;
 $line=<SF>;
 $line=<SR>;
 $line=<SR>;
 chomp($line);
 $seq.=" ".$line;
 $line=<SR>;$line=<SR>;
 $chimeric_reads{$seq_name}=$seq;        ### Paired-End
           }

############################# open split reads;
open SOFT,"$ARGV[1]_1sf.fuq";
my $seq_name="";
my $seq="";
while(<SOFT>){
 $line=$_;
 @line=split /\|/;
 if($line=~"soft"){$seq_name=$line[1]."|".$line[2]."|".$line[8];}
 else{$seq_name=$line;$seq_name=~s/^@//;chomp($seq_name);}
 $line=<SOFT>;
 chomp($line);
 $seq=$line;
 $line=<SOFT>;
 $line=<SOFT>;
 $split_reads{$seq_name}=$seq;
             }

########################### open input files;
open INPUT,"$ARGV[0]";
my $col=$ARGV[2];
@start=(0) x 4;
@end=(0) x 4;
@seq=(0) x 4;
@seq_final=(0) x 4;
my @length=(0) x 4;
my $first_line=<INPUT>;
my $chr=$line[13];
my $virus=$line[16];
my $gi_number=$line[17];
my $order=0;
my @vir_position_sort=();
my @hum_position_sort=();
my $seq_name="";
my $seq="";
my @start2=(0) x 2;
my @end2=(0) x 2;
my @final_infor=();
my @alignment_type=(0) x 4;

############################## print out title;
print "Sample\tOrder\tAver_a\tAver_u\tU_virus\tU_vip\tUniq_gi\tChim\tChim_q\tSplit\tSplit_q\tT_vip\tT_hist\tT_gi\tString\tChr.\tH_start\tH_end\tVirus\tGi\tV_start\tV_end\t";
print "Direction\tRead_name\tFlag\tChr.\tPosi.\t0\tMapping\tMate_hchr\tMate_hposi\t0\tAS\tV_flag\tV_posi.\t0\tMapping\tMate_vchr\tMate_vposi.\t0\tAS\tAS2\t0\t0\tAS3\tDirectory2\tStart_r\tEnd_r\tLength_r\tRE\tOverlap_RE\tType_h\tOverlap_h\n";

while(<INPUT>){                        #0      each VIP
 my @line=split;
 unless($line[11]>=2){next;}            # filter
# unless($line[4]>=1){next;}             #
 $order++;

########################               # get the region start and end position on human and virus genome;
 my @vir_position=(); my @hum_position=();
 for(my $j=$col+6;$j<@line;$j+=31){push(@vir_position,$line[$j+12]);push(@hum_position,$line[$j+4]);} 
 @vir_position_sort = sort {$a <=> $b } @vir_position;
 @hum_position_sort = sort {$a <=> $b } @hum_position; 
 for(my $j=$col+6;$j<@line;$j+=31){    #1      each read

########################               # get mapping information on human
 @start=(0) x 4;
 @end=(0) x 4;
 @seq=(0) x 4;
 @start2=(0) x 2;
 @end2=(0) x 2;
 @alignment_type=(0) x 4;

 my @temp1=($line[$j+6]=~/(\d+)/g);
 my @temp2=($line[$j+6]=~/([A-Z])/g);
 for(my $i=0;$i<@temp1;$i++){
  if($temp2[$i] eq "S"){   
   if($i==0){$start[0]=$temp1[$i]+1;$end[0]=$temp1[$i];$length[0]=$temp1[$i];}
   if($i==@temp1-1){$length[0]+=$temp1[$i];}
                       }  
  elsif($temp2[$i] eq "D"){next;}
  elsif($temp2[$i] eq "I"){$end[0]+=$temp1[$i];$length[0]+=$temp1[$i];}
  elsif($temp2[$i] eq "M"){
   if($i==0){$start[0]=1;$end[0]=$temp1[$i];$length[0]=$temp1[$i];}
   elsif($i==@temp1-1){$length[0]+=$temp1[$i];$end[0]+=$temp1[$i];}
   else{$end[0]+=$temp1[$i];$length[0]+=$temp1[$i];}}
                            }
$length[1]=$length[0];$start[1]=$line[$j+4]-$start[0];$end[1]=$start[1]+$length[1];

########################                  # get mapping information on virus;                       # okay
 my @tema1=($line[$j+14]=~/(\d+)/g);
 my @tema2=($line[$j+14]=~/([A-Z])/g);
 for(my $i=0;$i<@tema1;$i++){
  if($tema2[$i] eq "S"){
   if($i==0){$start[2]=$tema1[$i]+1;$end[2]=$tema1[$i];$length[2]=$tema1[$i];}
   if($i==@tema1-1){$length[2]+=$tema1[$i];}
                       }
  elsif($tema2[$i] eq "D"){next;}
  elsif($tema2[$i] eq "I"){$end[2]+=$tema1[$i];$length[2]+=$tema1[$i];}
  elsif($tema2[$i] eq "M"){
   if($i==0){$start[2]=1;$end[2]=$tema1[$i];$length[2]=$tema1[$i];}
   elsif($i==@tema1-1){$length[2]+=$tema1[$i];$end[2]+=$tema1[$i];}
   else{$end[2]+=$tema1[$i];$length[2]+=$tema1[$i];}}
                            }
$length[3]=$length[2];$start[3]=$line[$j+12]-$start[2];$end[3]=$length[3]+$start[3];

#########################                 # extract sequences;
 if($line[$j] =~"S"){                         ### split reads, extract sequences;
  $seq_name=$line[$j+1]."|".$line[$j+2]."|".$line[$j+8];
  my @temp4=split /\s+/, $split_reads{$seq_name};
  my @temp3=();
  @temp3=split //,$temp4[0];
  $seq[0]=join("",@temp3);
  if($line[$j+4] ne $line[$j+8]){
    $seq[0]=("-" x ($start[0]-1)).("@" x ($end[0]-$start[0]+1)).$seq[0];
    if(length($seq[0])<$length[0]){$seq[0]=$seq[0].("+" x ($length[0]-length($seq[0])))}
                                }
  else{
   $seq[0]=$seq[0].("@" x ($end[0]-$start[0]+1)).("-" x ($length[0]-$end[0]));
   if(length($seq[0])<$length[0]){$seq[0]=("+" x ($length[0]-length($seq[0]))).$seq[0];}
      }
  $seq[1]=uc($seq[0]);
  my $len_h1=abs(($hum_position_sort[0]-120)-$start[1]);
  my $len_h2=120+$hum_position_sort[$#hum_position_sort]-$end[1];

############## softclip reads;
  @temp3=();
  @temp3=split //,$temp4[0];
  if($line[$j+11]%32>=16){@temp3=reverse(@temp3);}
  $temp3[$_]=lc($temp3[$_]) for($start[2]-1..$end[2]-1);
  $seq[2]=join("",@temp3);
  $seq[3]=$seq[2];
  my $len_v1=abs(($vir_position_sort[0]-120)-$start[3]);
  my $len_v2=120+$vir_position_sort[$#vir_position_sort]-$end[3];
  $start2[0]=$start[1]-$len_h1;
  $start2[1]=$start[3]-$len_v1;
  $end2[0]=$end[1]+$len_h2;
  $end2[1]=$end[3]+$len_v2;
  if($line[$j+2]%256<128){@alignment_type[0]="S1";}
  elsif($line[$j+2]%256>=128){@alignment_type[0]="S1";}
  if($line[$j+2]%32<16){@alignment_type[1]="F";$seq[1]=("." x $len_h1).$seq[1].">".("." x ($len_h2-1));}          # update on 4/6/2016    # show read direction 
  elsif($line[$j+2]%32>=16){@alignment_type[1]="R";$seq[1]=("." x ($len_h1-1))."<".$seq[1].("." x $len_h2);}      # update on 4/6/2016
  if($line[$j+11]%256<128){@alignment_type[2]="S1";}
  elsif($line[$j+11]%256>=128){@alignment_type[2]="S1";}
  if($line[$j+11]%32<16){@alignment_type[3]="F";$seq[3]=("." x $len_v1).$seq[3].">".("." x ($len_v2-1));}         # update on 4/6/2016
  elsif($line[$j+11]%32>=16){@alignment_type[3]="R";$seq[3]=("." x ($len_v1-1))."<".$seq[3].("." x $len_v2);}     # update on 4/6/2016
                    }
 else{                                        ### chimeric reads, extract sequences
  $seq_name=$line[$j+1];
  my @temp4=split /\s+/, $chimeric_reads{$seq_name};
  my @temp4_2=@temp4;
  my @temp3=();
  if($line[$j+2]%256<128){
   if(length($temp4[0])<$length[0]){$temp4[0]=$temp4[0].("+" x ($length[0]-length($temp4[0])));}
   @temp3=split //,$temp4[0];
   if($line[$j+2]%32>=16){@temp3=reverse(@temp3);}
                         }
  elsif($line[$j+2]%256>=128){
   if(length($temp4[1])<$length[1]){$temp4[1]=$temp4[1].("+" x ($length[1]-length($temp4[1])));}
   @temp3=split //,$temp4[1];
   if($line[$j+2]%32>=16){@temp3=reverse(@temp3);}
                             }
  $temp3[$_]=lc($temp3[$_]) for($start[0]-1..$end[0]-1);
  $seq[0]=join("",@temp3);
  $seq[1]=$seq[0];
  my $len_h1=abs(($hum_position_sort[0]-120)-$start[1]);
  my $len_h2=$hum_position_sort[$#hum_position_sort]+120-$end[1];

########################
  if($line[$j+11]%256<128){
   @temp3=split //,$temp4_2[0];
   if($line[$j+11]%32>=16){@temp3=reverse(@temp3);}
                          }
  elsif($line[$j+11]%256>=128){
   @temp3=split //,$temp4_2[1];
   if($line[$j+11]%32>=16){@temp3=reverse(@temp3);}
                              }
  $temp3[$_]=lc($temp3[$_]) for($start[2]-1..$end[2]-1);
  $seq[2]=join("",@temp3);
  my $len_v1=abs(($vir_position_sort[0]-120)-$start[3]);
  my $len_v2=$vir_position_sort[$#vir_position_sort]+120-$end[3];
  $seq[3]=$seq[2];
  $start2[0]=$start[1]-$len_h1;
  $start2[1]=$start[3]-$len_v1;
  $end2[0]=$end[1]+$len_h2;
  $end2[1]=$end[3]+$len_v2;
  if($line[$j+2]%256<128){@alignment_type[0]="F";}
  elsif($line[$j+2]%256>=128){@alignment_type[0]="S";}
  if($line[$j+2]%32<16){@alignment_type[1]="F";$seq[1]=("." x $len_h1).$seq[1].">".("." x ($len_h2-1));}            # update on 4/6/2016
  elsif($line[$j+2]%32>=16){@alignment_type[1]="R";$seq[1]=("." x ($len_h1-1))."<".$seq[1].("." x $len_h2);}        # update on 4/6/2016
  if($line[$j+11]%256<128){@alignment_type[2]="F";}
  elsif($line[$j+11]%256>=128){@alignment_type[2]="S";}
  if($line[$j+11]%32<16){@alignment_type[3]="F";$seq[3]=("." x $len_v1).$seq[3].">".("." x ($len_v2-1));}           # update on 4/6/2016
  elsif($line[$j+11]%32>=16){@alignment_type[3]="R";$seq[3]=("." x ($len_v1-1))."<".$seq[3].("." x $len_v2);}       # update on 4/6/2016
     }
  push(@final_infor,[@line,$seq_name,@alignment_type,@start,@end,@length,@start2,@end2,@seq]);          # i add read name
                                }      #1     each read

#####################            ##### calculate the length of virus regions;
  my @vir_posi=(); my @non_overlap=(); my $dd=0;
  for(my $o=0;$o<@final_infor;$o++){
   push(@vir_posi,${$final_infor[$o]}[$#{$final_infor[$o]}-16]);
   push(@vir_posi,${$final_infor[$o]}[$#{$final_infor[$o]}-12]);
   unless(@non_overlap){push(@non_overlap, @vir_posi);$dd=1;}
   else{
    for(my $p=0;$p<@non_overlap;$p+=2){
     if($vir_posi[0]>=$non_overlap[$p] && $vir_posi[0]<=$non_overlap[$p+1] && $vir_posi[1]>=$non_overlap[$p+1]){$non_overlap[$p+1]=$vir_posi[1];$dd=1;}
     elsif($vir_posi[1]>=$non_overlap[$p] && $vir_posi[1]<=$non_overlap[$p+1] && $vir_posi[0]<=$non_overlap[$p]){$non_overlap[$p]=$vir_posi[0];$dd=1;}
     elsif($vir_posi[0]>=$non_overlap[$p] && $vir_posi[1]<=$non_overlap[$p+1]){$dd=1;next;}
     elsif($vir_posi[0]<=$non_overlap[$p] && $vir_posi[1]>=$non_overlap[$p+1]){$dd=1;$non_overlap[$p]=$vir_posi[0];$non_overlap[$p+1]=$vir_posi[1];}
                                      }
       }
  if($dd==0){push(@non_overlap, @vir_posi);}
  $dd=0;
  @vir_posi=();
                                   }
  my $total_len=0;
  for(my $o=0;$o<@non_overlap;$o+=2){$total_len+=$non_overlap[$o+1]-$non_overlap[$o];}
  @non_overlap=();
  @vir_posi=();

#####################           ##### Calculate the length of human regions;        #### updated on 4/14/16
  my @vir_posi_h=(); my @non_overlap_h=(); my $dd=0;
  for(my $o=0;$o<@final_infor;$o++){
   push(@vir_posi_h,${$final_infor[$o]}[$#{$final_infor[$o]}-18]+${$final_infor[$o]}[$#{$final_infor[$o]}-19]);
   push(@vir_posi_h,${$final_infor[$o]}[$#{$final_infor[$o]}-18]+${$final_infor[$o]}[$#{$final_infor[$o]}-15]);
   unless(@non_overlap_h){push(@non_overlap_h, @vir_posi_h);$dd=1;}
   else{
    for(my $p=0;$p<@non_overlap_h;$p+=2){
     if($vir_posi_h[0]>=$non_overlap_h[$p] && $vir_posi_h[0]<=$non_overlap_h[$p+1] && $vir_posi_h[1]>=$non_overlap_h[$p+1]){$non_overlap_h[$p+1]=$vir_posi_h[1];$dd=1;}
     elsif($vir_posi_h[1]>=$non_overlap_h[$p] && $vir_posi_h[1]<=$non_overlap_h[$p+1] && $vir_posi_h[0]<=$non_overlap_h[$p]){$non_overlap_h[$p]=$vir_posi_h[0];$dd=1;}
     elsif($vir_posi_h[0]>=$non_overlap_h[$p] && $vir_posi_h[1]<=$non_overlap_h[$p+1]){$dd=1;next;}
     elsif($vir_posi_h[0]<=$non_overlap_h[$p] && $vir_posi_h[1]>=$non_overlap_h[$p+1]){$dd=1;$non_overlap_h[$p]=$vir_posi_h[0];$non_overlap_h[$p+1]=$vir_posi_h[1];}
                                        }
       }
  if($dd==0){push(@non_overlap_h, @vir_posi_h);}
  $dd=0;
  @vir_posi_h=();
                                  }
  my $total_len_h=0;
  for(my $o=0;$o<@non_overlap_h;$o+=2){$total_len_h+=$non_overlap_h[$o+1]-$non_overlap_h[$o];}
  @non_overlap_h=();
  @vir_posi_h=();


#####################           ##### put the sequences into two array;
  my @seq_h=();my @seq_v=();my $len_infor=0; my @seq_temp=();
  my $out_temp1=join("\t",(@{$final_infor[0]}[0..17],$start2[1],$end2[1],@{$final_infor[0]}[18..$#{$final_infor[0]}-24]));
  print "O1\t$ARGV[1]\t$order\t$total_len_h\t$total_len\t$out_temp1\n";
  for(my $m=0;$m<@final_infor;$m++){
  print "O2\t@{$final_infor[$m]}[$#{$final_infor[0]}-24..$#{$final_infor[0]}-4] ${$final_infor[$m]}[$#{$final_infor[0]}-3] ${$final_infor[$m]}[$#{$final_infor[0]}-1]\n";
  $len_infor=@{$final_infor[$m]};
  push(@seq_h,${$final_infor[$m]}[$len_infor-3]);
  push(@seq_v,${$final_infor[$m]}[$len_infor-1]);
                                   }

#####################           ##### organize the results and print out human
   print "\nO3\t\/\/Start\/\/\n\t\t\t\t\t\t********************************** human part\n";
   my $start_temp=${$final_infor[0]}[$#{$final_infor[0]}-7];
   my $end_temp=$start_temp+$region;
   for(my $k=0;$k<length($seq_h[0])+$region;$k+=$region){
    @seq_temp=();
    for(my $n=0;$n<@seq_h;$n++){
     my $seq_temp=substr $seq_h[$n],$k,$region;
     push(@seq_temp,$seq_temp);
                               }
  my $line_temp="." x $region; my $keep=0;
  for(my $d=0;$d<@seq_temp;$d++){if($seq_temp[$d] ne $line_temp && ($seq_temp[$d])){$keep=1;}}
  if($keep==1){
   my $seq_tag=0;
   for(my $d=0;$d<@seq_temp;$d++){
    if(length($seq_temp[$d])<$region){$seq_temp[$d]=$seq_temp[$d].(" " x ($region-length($seq_temp[$d])));}
    print "\t\t\tSeq$seq_tag\t$start_temp\t\t$seq_temp[$d]\t\t$end_temp\n";
    $seq_tag++;
                                 }
   print "\n";
              }
  $start_temp=$end_temp+1;
  $end_temp=$start_temp+$region;
                                                        }
  print "\t\t\t\t\t\t************************************* virus part\n";

#####################                               ##### organize the results and print out virus
   my $start_temp=${$final_infor[0]}[$#{$final_infor[0]}-6];
   my $end_temp=$start_temp+$region;my $tem_blank="";
   for(my $k=0;$k<length($seq_v[0])+$region;$k+=$region){
    @seq_temp=();
    for(my $n=0;$n<@seq_v;$n++){
     my $seq_temp=substr $seq_v[$n],$k,$region;
     push(@seq_temp,$seq_temp);
                               }
    my $line_temp="." x $region; my $keep=0; my $count_interval=0;
    for(my $d=0;$d<@seq_temp;$d++){
     if($seq_temp[$d] ne $line_temp && ($seq_temp[$d])){$keep=1;}
     if($seq_temp[$d] eq $line_temp){$count_interval++;}else{$count_interval=0;}
                                  }
    
    if($keep==1){
     my $seq_tag=0;
     for(my $d=0;$d<@seq_temp;$d++){
      if(length($seq_temp[$d])<$region){$seq_temp[$d]=$seq_temp[$d].(" " x ($region-length($seq_temp[$d])));}          # last line;
      print "\t\t\tSeq$seq_tag\t$start_temp\t\t$seq_temp[$d]\t\t$end_temp\n";
     $seq_tag++;$tem_blank="";
                                   }
     print "\n";
                }
   else{
    my $tem_len=@seq_temp;my $tem_seq="-" x $region;
    unless($tem_blank){
     $tem_blank=$tem_len;
     if($count_interval eq $tem_len && $count_interval){print "\t\t\t\t\t\t$tem_seq\n";}
                      }
       }
   $start_temp=$end_temp+1;
   $end_temp=$start_temp+$region;                        
                                                       }
#####################
  print "O3\t\/\/End\/\/\n\n";
  @final_infor=();
              }                        #0     each VIP
