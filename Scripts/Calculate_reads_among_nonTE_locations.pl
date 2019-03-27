#!/usr/bin/env perl
#
# Author:       Xun Chen
# Email:        Xun.Chen@uvm.edu
# Date:         02/10/19
#
# Usage: perl Calculate_reads_among_nonTE_locations.pl -i Combined.VCF -S SampleID -o Output.nonTE -b bamFile.bam -s paired-end -l Length_insertsize -L Std_insertsize -r Read_length -t Threads
#
#

use Getopt::Long qw(:config no_ignore_case);
use strict;
use warnings;
use Cwd qw();

my $input;
my $output;
my $bamFile;
my $sequencing_type;
my $length_insertsize;
my $L_std_insertsize;
my $std_insertsize;
my $read_len;
my $threads;
my $samtools_d ="";
my $Position;
my $line="";
my $line2="";
my @line=();
my %Chr_posi = ();
my %Each_TE1 = ();
my $Sample_ID;

# Variables
 GetOptions(
   "i|input=s" => \$input,
   "S|Sample_ID=s" =>\$Sample_ID,
   "o|output=s" => \$output,
   "b|bamFile=s" => \$bamFile,
   "s|sequencing_type=s" => \$sequencing_type,
   "l|length_insertsize=f" => \$length_insertsize,
   "L|L_std_insertsize=f" => \$L_std_insertsize,
   "r|read_len=i" => \$read_len,
   "t|threads=i" => \$threads 
 );

if (!defined ($sequencing_type)) {
  $sequencing_type = "paired-end";
}
if (!defined ($read_len)) {
  $read_len = 100;   
}
if (!defined ($threads)) {
  $threads = 1;   
}
if (defined ($L_std_insertsize)) {
  $std_insertsize = $L_std_insertsize;
}

##### estimate insertsize
if (!(defined ($length_insertsize)) || !(defined ($L_std_insertsize))) {
  ($length_insertsize,$std_insertsize) = estimate_insertsize($bamFile,${input},0.05);
}
my $min_insertsize = $length_insertsize - 2 * $std_insertsize;

##### directory
my $directory=$0;
my @directory=split /\//, $directory;
if($#directory>0){
  $directory=join("/",@directory[0..$#directory-1]);
  $directory=$directory."/";
} else {
  $directory="";
}

open INPUT, "${input}";
open OUTPUT,">${output}_tmp";
if (-e ${output}.".sam2") {
  system ("rm ${output}.sam2");
} else {
  system ("touch ${output}.sam2");
}

while($line2=<INPUT>){
  @line=split /\s+/, $line2;
  if ( $line[0]=~ "#" ) {
    next;
  }
  my $name_tmp1 = $line[0].":".$line[1]."-".$line[1];
  @{$Each_TE1{$name_tmp1}} = (0) x 6;
  ${$Each_TE1{$name_tmp1}}[0] = ${input};
  ${$Each_TE1{$name_tmp1}}[1] = $line[0];
  ${$Each_TE1{$name_tmp1}}[2] = $line[1]+1;
  ${$Each_TE1{$name_tmp1}}[3] = 0;
  ${$Each_TE1{$name_tmp1}}[4] = 1;
  ${$Each_TE1{$name_tmp1}}[5] = $read_len;
  @{$Chr_posi{$name_tmp1}} = @line[0..1];
  @{$Chr_posi{$name_tmp1}}[2..8] = ("-") x 7;
  ${$Chr_posi{$name_tmp1}}[2]=$Sample_ID;

  ${$Each_TE1{$name_tmp1}}[1]=$line[0];
  $Position=$line[1]+1;
  my $position1=$Position-int(($length_insertsize+2*$std_insertsize)*2);               ###### needed this script
  my $position2=$Position+int(($length_insertsize+2*$std_insertsize)*2);               ###### need this script
  my $cmd_1 = q(awk '{print");
  my $cmd_2 = q(\t"$0}');
  system ("${samtools_d}samtools view $bamFile ${$Each_TE1{$name_tmp1}}[1]:${position1}-${position2} | sort | uniq | $cmd_1$name_tmp1$cmd_2 >>${output}.sam2");
}

######
 my @sam2 = ();
 my $sam2_id = "";
 my @sam2_tmp3 = ();
 open SAM2, "${output}.sam2";
 while (my $sam2_tmp1 = <SAM2>) {
   my @sam2_tmp1 = split /\s+/, $sam2_tmp1;
   my @sam2_tmp2 = @sam2_tmp1[1..$#sam2_tmp1];
   if ($sam2_id eq "") {
     $sam2_id = $sam2_tmp1[0];
     push (@sam2,[@sam2_tmp2]);
   } elsif ($sam2_id eq $sam2_tmp1[0] && !(eof(SAM2))) {
     push (@sam2,[@sam2_tmp2]);
   } elsif ($sam2_id ne $sam2_tmp1[0] || eof (SAM2)) {
     if ($sam2_id eq $sam2_tmp1[0]) {
       push (@sam2,[@sam2_tmp2]);
     }
     @sam2_tmp3 = Genotype_each_TE1(@sam2);
     @sam2_tmp3 = sort {$a->[1] cmp $b->[1] || $b->[0] cmp $a->[0]} @sam2_tmp3;
     my @sam2_tmp4 = Genotype_each_TE2(@sam2_tmp3);
     my @line_1 = ();
     if ($sequencing_type eq "paired-end") {
       @line_1 = Genotype_each_TE3_PE(\@sam2_tmp4,\@{$Each_TE1{$sam2_id}});
     } else {
       @line_1 = Genotype_each_TE3_SE(\@sam2_tmp4,\@{$Each_TE1{$sam2_id}});
     }
     $line_1[0] =0;
     push (@{$Chr_posi{$sam2_id}},@line_1);
     print OUTPUT "@{$Chr_posi{$sam2_id}}\n";
     @sam2 = ();
     push (@sam2,[@sam2_tmp2]);
     $sam2_id = $sam2_tmp1[0];
   }
 }
 close OUTPUT;
 close SAM2;

 system ("Rscript ${directory}Genotype_likelihood.R ${output}_tmp ${output}_tmp2");
 system ("rm ${output}_tmp");
 system ("rm ${output}.sam2");
 
 open NON_TE, ">${output}";
 open TMP2, "${output}_tmp2";
 while (my $line_TE = <TMP2>) {
  my @line_TE = split /\s+/, $line_TE;
  my $name_TE = $line_TE[2]."%".$line_TE[0]."%".$line_TE[1];
  my $geno = "";
  if ($line_TE[10] eq 0) {
    $geno = "./.:.:.,.,.:.:.";
  } else {
    $geno = "0/0:".$line_TE[18].":".$line_TE[15].",".$line_TE[14].",".$line_TE[13].":".$line_TE[10].":".$line_TE[9];
  }
  print NON_TE "$name_TE\t$geno\n";
}
 close TMP2;
 close NON_TE;
 system ("rm ${output}_tmp2");

######################################## Sub functions
###### calculate insert size
sub estimate_insertsize {
  my ($input1,$filename,$proportion)=@_;
  my @insert_size = ();
  my $total = 0;
  system ("${samtools_d}samtools view -@ $threads -s $proportion $input1 > ${filename}_$proportion");
  open SUB_SAM, "${filename}_$proportion";
  while (my $tmp1_sam=<SUB_SAM>) {
    my @tmp1_sam = split /\s+/, $tmp1_sam;
    if ($tmp1_sam[1] % 2 >= 1 && $tmp1_sam[1] % 4 >= 2 && $tmp1_sam[1] % 128 >= 64 && $tmp1_sam[1] % 512 < 256 && $tmp1_sam[1] < 256) {
      push (@insert_size, abs($tmp1_sam[8]));
      $total += abs($tmp1_sam[8]);
    }
  }
  close SUB_SAM;
  system ("rm ${filename}_$proportion");
  my $average = $total / ($#insert_size+1);
  $average = sprintf("%.4f",$average);
  my $sqtotal =0;
  foreach(@insert_size) {
    $sqtotal += ($average-$_) ** 2;
  }
  my $std = ($sqtotal / $#insert_size) ** 0.5;
  $std = sprintf("%.4f",$std);
  return ($average,$std);
}

###### Genotyping function 1: filtering
sub Genotype_each_TE1 {
  my (@input) = @_;
  my @return1 = ();
  my $line = "";
  my @line = ();
  my $type1 = "";
  my $k = 0;
  my @len = ();
  my $as = "";
  my $md = "";
  my $posi = "";
  for (my $i = 0; $i < @input; $i++) {
    @line = @{$input[$i]};
    $as = "";
    $md = "";
    $posi = "";
    if($line[0] =~ "#@") {
      next;
    }
    @len = ($line[5]=~/(\d+)/g);
    my $len=0; map {$len+=$_} @len;
    $k=0;
    $type1="";
    if ($line[1]%8<4) { 
      if ($line[1]%4>=2 && $line[1]%256>=128) {
        $type1="PR";
        $k=1;
      } elsif ($line[1]%4>=2 && $line[1]%256<128) {
        $type1="PF";
        $k=1;
      } elsif ($line[1]%256<128) {
        $type1="F";
        $k=1;
      } elsif($line[1]%256>=128) {
        $type1="R";
        $k=1;
      } else {
        $type1="N";
        $k=1;
      }
      for (my $i=11;$i<@line;$i++) {
        if ($line[$i]=~"AS:i:") {
          $as=$line[$i];
          $as=~s/AS\:i\://;
        }
        if ($line[$i]=~"MD:Z:") {
          $md=$line[$i];
          $md=~s/MD\:Z\://; 
        }
        my $hstart = 0;
        my $hend = 0;
        my $hlength = 0;
        my @temp1 = ($line[5]=~/(\d+)/g);
        my @temp2 = ($line[5]=~/([A-Z])/g);
        for (my $j=0;$j<@temp1;$j++) {
          if ($temp2[$j] eq "S" || $temp2[$j] eq "H") {   
            if ($j==0) {
              $hstart=$temp1[$j]+1;
              $hend=$temp1[$j];
              $hlength=$temp1[$j];
            }
            if ($j==@temp1-1) {
              $hlength+=$temp1[$j];
            }
          } elsif ($temp2[$j] eq "D") {
            next;
          } elsif ($temp2[$j] eq "M" || $temp2[$j] eq "I") {
            if($j == 0) {
              $hstart = 0;
              $hend = $temp1[$j];
              $hlength = $temp1[$j];
            } elsif($j == @temp1-1) {
              $hlength += $temp1[$j];
              $hend += $temp1[$j];
            } else {
              $hend += $temp1[$j];
              $hlength += $temp1[$j];
            }
          }
        }
        $posi=$hstart."_".$hend."_".$hlength."_".$md;
      }
    } else {
      $type1="U";
      next;
    }
    unless ($as) {
      $as=0;
    }
    my @return1_tmp=();
    push (@return1_tmp,$type1);
    push (@return1_tmp,@line[0..8]);
    push (@return1_tmp,$posi);
    push (@return1_tmp,$as);
    push (@return1,[@return1_tmp]);
  }
  return (@return1);
}

##### Genotyping function 2: extract reads cross the breakpoint 
sub Genotype_each_TE2 {
  my (@input2) = @_;
  my $insert_size_geno = $min_insertsize;
  my $line="";
  my @line=();
  my %read=();
  my $read_name="";
  my $name="";
  my $order=0;
  my $order1=0;
  my %human=();
  my @return2=();
  my @tmp1 = (0) x 26;
  push (@return2,[@tmp1]);

  for (my $i = 0; $i < @input2; $i++) {
    @line = @{$input2[$i]};
    if($line[2] / 512>=1){
      next;
    }
    $line[5] = $line[3];
    my $name2=$line[1]."_".$line[3]."_".$line[8];
    $name=$line[1]."_".$line[3]."_".$line[4];
    if(!(%read)) {
      @{$read{$name}} = @line;
      $read_name = $line[1];
      unshift(@{$read{$name}},$line[11]);
      unshift(@{$read{$name}},$line[3]);
      unshift(@{$read{$name}},$line[0]);
    }
    if($line[1] eq $read_name) {
      if(exists($read{$name2}) && (($line[0] eq "PF" && ${$read{$name2}}[0] eq "PR") || ($line[0] eq "PR" && ${$read{$name2}}[0] eq "PF")) && $line[7] eq "="){
        if ($line[3] eq ${$read{$name2}}[6] && $line[4] eq ${$read{$name2}}[11] && $line[8] eq ${$read{$name2}}[7]){        ##### 02/14/2019
          splice @{$read{$name2}},3,0,@line;
          ${$read{$name2}}[0]="P";                                               
          ${$read{$name2}}[2]+=$line[11];
        } else {
          if((${$read{$name2}}[6] ne $line[3])){                                                                            ##### 02/14/2019
            $name=$line[1]."_".$line[3]."_".$line[4];
            @{$read{$name}}=@line;
            unshift(@{$read{$name}},$line[11]);
            unshift(@{$read{$name}},$line[3]);
            unshift(@{$read{$name}},$line[0]); 
          } else {
            splice @{$read{$name2}},3,0,@line;
            ${$read{$name2}}[0]="P2";
            ${$read{$name2}}[2]+=$line[11];
          }
        }
      } else {
        if($order eq 1){next;}
        $name=$line[1]."_".$line[3]."_".$line[4];
        @{$read{$name}}=@line;
        unshift(@{$read{$name}},$line[11]);
        unshift(@{$read{$name}},$line[3]);
        unshift(@{$read{$name}},$line[0]);
      }
    } else {
      my @order=sort{$read{$b}->[2]<=>$read{$a}->[2]} keys %read;
      for(my $j=0,my $max=${$read{$order[0]}}[2];$j<@order;$j++){
        ${$read{$order[$j]}}[1] = $insert_size_geno;
        if(${$read{$order[$j]}}[2]<$max){
          last;
        } else {
          push (@return2,[@{$read{$order[$j]}}]);
        }
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
    for(my $j=0,my $max=${$read{$order[0]}}[2];$j<@order;$j++) {
      ${$read{$order[$j]}}[1] = $insert_size_geno;
      if(${$read{$order[$j]}}[2]<$max){
        last;
      } else {
        push (@return2,[@{$read{$order[$j]}}]);
      }
    }
  }
  return @return2;
}

##### Genotyping function 3: calculate the reads support non TE
sub Genotype_each_TE3_PE {
  my ($input3,$input4) = @_;
  my @input3 = @{$input3};
  my @input4 = @{$input4};
  my $line="";
  my @line=();
  my %VI_list=();
  my @names=();
  my %split=();
  my %split3=();
  my %split2=();
  ### each PE read
  for (my $i=0;$i < @input3; $i++){
    @line = @{$input3[$i]};
    ### filtering
    if ($line[0] ne "P" && abs($line[7]-$input4[2]) < ($length_insertsize+2*$std_insertsize)*2 && $line[6] eq $input4[1]) {                     ##### can be further optimized
    } elsif ($line[0] eq "P" && $line[6] eq $input4[1]){
      if ($line[7]>=$line[19] && $input4[2]>=$line[19] && $input4[2]<=($line[7]+$input4[5])) {
   
      } elsif ($line[7]<=$line[19] && $input4[2]<=($line[19]+$input4[5]) && $input4[2]>=$line[7]) {
   
      } else {
        next;
      }
    } else {
      next;
    }
    ### merging PE
    my $name=$input4[0]."%".$input4[1]."%".$input4[2];                      ### corrected 2/15/2019
    unshift (@line,$name);
    if (exists($VI_list{$line[0]."|".$line[5]}) && ${$VI_list{$line[0]."|".$line[5]}}[0] eq $line[0] && ${$VI_list{$line[0]."|".$line[5]}}[7] eq $line[7] && ${$VI_list{$line[0]."|".$line[5]}}[8] eq $line[12] && ${$VI_list{$line[0]."|".$line[5]}}[12] eq $line[8]) {
      ${$VI_list{$line[0]."|".$line[5]}}[1]="P2";
      ${$VI_list{$line[0]."|".$line[5]}}[3]=${$VI_list{$line[0]."|".$line[5]}}[3]+$line[3];
      push(@{$VI_list{$line[0]."|".$line[5]}},@line[4..$#line]);
    } else {
      push(@{$VI_list{$line[0]."|".$line[5]}},@line);
    }
  }
  ### extract PE reads and count
  my %vi=();
  my @reads=keys %VI_list;
  for(my $i=0;$i<@reads;$i++){
    @line=@{$VI_list{$reads[$i]}};
    unless(exists($vi{$line[0]})) {
      @{$vi{$line[0]}}=(0) x 4;
    }
    @names=split /\%/, $line[0];
    if($line[1] ne "P" && $line[1] ne "P2"){
      next;
    }
    ### extract split reads
    my $start1=0;
    my $end=0;
    my $length=0;
    my $start_1=0;
    my $end_1=0;
    my $len_1=0;
    my $start_2=0;
    my $end_2=0;
    my $len_2=0;
    if($line[10] =~ "S" || $line[22] =~ "S"){                      ##### split start
      ### first end
      my @temp1=($line[10]=~/(\d+)/g);
      my @temp2=($line[10]=~/([A-Z])/g);
      for(my $i=0;$i<@temp1;$i++) {
        if($temp2[$i] eq "S" || $temp2[$i] eq "H"){   
          if($i==0){
            $start1=$temp1[$i]+1;
            $end=$temp1[$i];
            $length=$temp1[$i];
          }
          if($i==@temp1-1) {
            $length+=$temp1[$i];
          }
        } elsif ($temp2[$i] eq "D") {
          next;
        } elsif ($temp2[$i] eq "M" || $temp2[$i] eq "I"){
          if($i==0) {
            $start1=1;
            $end=$temp1[$i];
            $length=$temp1[$i];
          } elsif($i==@temp1-1) {
            $length+=$temp1[$i];
            $end+=$temp1[$i];
          } else {
            $end+=$temp1[$i];
            $length+=$temp1[$i]
          }
        }
      }
      $start_1=$start1;
      $end_1=$end;
      $len_1=$length;
      ### second end
      @temp1=($line[22]=~/(\d+)/g);
      @temp2=($line[22]=~/([A-Z])/g);
      for(my $i=0;$i<@temp1;$i++){
        if($temp2[$i] eq "S" || $temp2[$i] eq "H"){                ##### 2/20/2019     
          if($i==0){ 
            $start1=$temp1[$i]+1;
            $end=$temp1[$i];
            $length=$temp1[$i];
          }
          if($i==@temp1-1) {
            $length+=$temp1[$i];
          }
        } elsif ($temp2[$i] eq "D") {
          next;
        } elsif ($temp2[$i] eq "M" || $temp2[$i] eq "I"){
          if($i==0) {
            $start1=1;
            $end=$temp1[$i];
            $length=$temp1[$i];
          } elsif ($i==@temp1-1) {
            $length+=$temp1[$i];
            $end+=$temp1[$i];
          } else {
            $end+=$temp1[$i];
            $length+=$temp1[$i]
          }
        }
      }
      $start_2=$start1;
      $end_2=$end;
      $len_2=$length;
    }                                                        ##### split end
    ### identify split >=20, and split <20 bp (suggestive reads)
    if($names[2] eq $line[8] && $line[10] =~ "S"){
      if($start_1<20){
        $line[1]=$line[1]."_F_".$len_1."_".$start_1."_".$end_1;
        if(exists($split2{$line[0]})) {
          $split2{$line[0]}++;
        } else {
          $split2{$line[0]}=1;
        }
      } else {
        $line[1]=$line[1]."_F_".$len_1."_".$start_1."_".$end_1;
        if (exists($split3{$line[0]})) {
          $split3{$line[0]}++;
        } else {
          $split3{$line[0]}=1;
        }         
      }
      push(@{$split{$line[0]."|".$line[5]}},@line);
    } elsif ($names[2] eq $line[20] && $line[22] =~ "S") {
      if($start_2<20){
        $line[1]=$line[1]."_R_".$len_2."_".$start_2."_".$end_2;
        if(exists($split2{$line[0]})){
          $split2{$line[0]}++;
        } else {
          $split2{$line[0]}=1;
        }                
      } else {
        $line[1]=$line[1]."_R_".$len_2."_".$start_2."_".$end_2;
        if (exists($split3{$line[0]})) {
          $split3{$line[0]}++;
        } else {
          $split3{$line[0]}=1;
        }               
      }
      push(@{$split{$line[0]."|".$line[5]}},@line);
    } elsif (($names[2] - ($line[8]+$end_1))<=5 && $line[10] =~ "S") {
      if ($len_1-$end_1<20) {
        $line[1]=$line[1]."_F2_".$len_1."_".$start_1."_".$end_1;
        if (exists($split2{$line[0]})) {
          $split2{$line[0]}++;
        } else {
          $split2{$line[0]}=1;
        }                
      } else {
        $line[1]=$line[1]."_F2_".$len_1."_".$start_1."_".$end_1;
        if (exists($split3{$line[0]})) {
          $split3{$line[0]}++;
        } else {
          $split3{$line[0]}=1;
        }               
      }
      push(@{$split{$line[0]."|".$line[5]}},@line);
    } elsif (($names[2] - ($line[20]+$end_2))<=5 && $line[22] =~ "S") {
      if($len_2-$end_2<20){
        $line[1]=$line[1]."_R2_".$len_2."_".$start_2."_".$end_2;
        if (exists($split2{$line[0]})) {
          $split2{$line[0]}++;
        } else {
          $split2{$line[0]}=1;
        }                 
      } else {
        $line[1]=$line[1]."_R2_".$len_2."_".$start_2."_".$end_2;
        if (exists($split3{$line[0]})) {
          $split3{$line[0]}++;
        } else {
          $split3{$line[0]}=1;
        }                 
      }
      push(@{$split{$line[0]."|".$line[5]}},@line);
    }
    ### reads support nonTE
    if(exists($vi{$line[0]}) && !(exists($split{$line[0]."|".$line[5]}))){
      my @temp1=($line[10]=~/(\d+)/g);
      my @temp2=($line[10]=~/([A-Z])/g);
      my $hstart="";
      my $hend="";
      my $hlength="";
      for(my $i=0;$i<@temp1;$i++){
        if($temp2[$i] eq "S" || $temp2[$i] eq "H"){                          #### 02/15/2019   
          if($i==0) {
            $hstart=$temp1[$i]+1;
            $hend=$temp1[$i];
            $hlength=$temp1[$i];
          }
          if ($i==@temp1-1) {
            $hlength+=$temp1[$i];
          }
        } elsif ($temp2[$i] eq "D") {
          next;
        } elsif ($temp2[$i] eq "M" || $temp2[$i] eq "I") {
          if ($i==0) {
            $hstart=1;
            $hend=$temp1[$i];
            $hlength=$temp1[$i];
          } elsif ($i==$#temp1) {
            $hlength+=$temp1[$i];
            $hend+=$temp1[$i];
          } else {
            $hend+=$temp1[$i];
            $hlength+=$temp1[$i];
          }
        }
      }
      if (abs(($hend+$line[8]-1)-$input4[2])<=10 && $line[8]>=$line[20]){
        next;
      } elsif (abs(($hstart+$line[8]-1)-$input4[2])<=10 && $line[8]<=$line[20]) {
        next;
      }
      @temp1=($line[22]=~/(\d+)/g);
      @temp2=($line[22]=~/([A-Z])/g);
      $hstart="";
      $hend="";
      $hlength="";

      for (my $i=0;$i<@temp1;$i++) {
        if($temp2[$i] eq "S" || $temp2[$i] eq "H") {              ### 2/20/2019
          if($i==0) {
            $hstart=$temp1[$i]+1;
            $hend=$temp1[$i];
            $hlength=$temp1[$i];
          }
          if ($i==$#temp1) {
            $hlength+=$temp1[$i];
          }
        } elsif ($temp2[$i] eq "D") {
          next;
        } elsif ($temp2[$i] eq "M" || $temp2[$i] eq "I") {
          if ($i==0) {
            $hstart=1;
            $hend=$temp1[$i];
            $hlength=$temp1[$i];
          } elsif($i==@temp1-1) {
            $hlength+=$temp1[$i];
            $hend+=$temp1[$i];
          } else {
            $hend+=$temp1[$i];
            $hlength+=$temp1[$i];
          }
        }
      }
      if (abs(($hstart+$line[20]-1)-$input4[2])<=10 && $line[20]<$line[8]) {
        next;
      } elsif (abs(($hend+$line[20]-1)-$input4[2])<=10 && $line[8]<$line[20]){
        next;
      }
      ${$vi{$line[0]}}[0]++;
    }                                                                    ### 1. here counting number of reads support Non-VI;
  }
  ### number of chimeric and split reads and ratio of the TE;
  my $name=$input4[0]."%".$input4[1]."%".$input4[2];                     # sample_id, chr, position
  unless (exists($split2{$name})){
    $split2{$name}=0;
  }                     # if there is no split reads
  unless (exists($split3{$name})) {
    $split3{$name}=0;
  }                     # if there is no suggestive reads
  unless ($vi{$name}) {
    @{$vi{$name}}=(0) x 3;
    ${$vi{$name}}[2]=1;
  }                     # 0. if there is no reads support non-TE
  ### array structure
  ${$vi{$name}}[1]=$input4[3];                                     # 1. chimeric & split reads 
  ${$vi{$name}}[3]=$split3{$name};                                 # 3. suggestive reads (split reads <20 bp)
  if(${$vi{$name}}[0] eq 0){                                       # if there is no reads support nonTEs
    ${$vi{$name}}[2]=1;
  } elsif(${$vi{$name}}[0] >0 && ${$vi{$name}}[1]>0){
    if($input4[4] eq 2) {
      ${$vi{$name}}[2]=(${$vi{$name}}[1]/2)/((${$vi{$name}}[1]/2)+${$vi{$name}}[0]);
    } elsif ($input4[4] eq 1) {
      ${$vi{$name}}[2]=${$vi{$name}}[1]/(${$vi{$name}}[1]+${$vi{$name}}[0]);
    }
  }
  my @return3 = ();
  push (@return3,$name);
  push (@return3,@{$vi{$name}});
  return @return3;
}

##### Genotyping function 3 for single-end: calculate reads support nonTE
sub Genotype_each_TE3_SE {
  my ($input3,$input4) = @_;
  my @input3 = @{$input3};
  my @input4 = @{$input4};
  my $line="";
  my @line=();
  my %VI_list=();
  my @names=();
  my %split=();
  my %split3=();
  my %split2=();
  my $fully_mapped=0;
  my $split_reads=0;
  my $split_reads_short=0;

  for (my $i = 0; $i < @input3; $i++) {
    @line = @{$input3[$i]};
    ### filtering
    if ($line[0] ne "P" && abs($line[7]-$input4[2]) < ($length_insertsize+2*$std_insertsize)*2 && $line[6] eq $input4[1]) {

    } elsif ($line[0] eq "P" && $line[6] eq $input4[1]){
      if ($line[7]>=$line[19] && $input4[2]>=$line[19] && $input4[2]<=($line[7]+$input4[5])) {
     
      } elsif ($line[7]<=$line[19] && $input4[2]<=($line[19]+$input4[5]) && $input4[2]>=$line[7]) {
    
      } else {
        next;
      }
    } else {
      next;
    }
    my $name=$input4[0]."%".$input4[1]."%".$input4[2];
    unshift (@line,$name);
    push(@{$VI_list{$line[0]."|".$line[5]}},@line);
  }
  my %vi=();
  my @reads=keys %VI_list;
  for(my $i=0;$i<@reads;$i++){
    @line=@{$VI_list{$reads[$i]}};
    @names=split /\%/, $line[0];
    my $start1=0;
    my $end=0;
    my $length=0;
    my @temp1=($line[10]=~/(\d+)/g);
    my @temp2=($line[10]=~/([A-Z])/g);
    for(my $i=0;$i<@temp1;$i++){
      if($temp2[$i] eq "S" || $temp2[$i] eq "H") {               ### 2/20/2019   
        if($i==0){$start1=$temp1[$i]+1;$end=$temp1[$i];$length=$temp1[$i];}
        if($i==@temp1-1){$length+=$temp1[$i];}
      } elsif($temp2[$i] eq "D") {
        next;
      } elsif($temp2[$i] eq "M" || $temp2[$i] eq "I"){
        if($i==0) {
          $start1=1;
          $end=$temp1[$i];
          $length=$temp1[$i];
        } elsif($i==@temp1-1) {
          $length+=$temp1[$i];
          $end+=$temp1[$i];
        } else {
          $end+=$temp1[$i];
          $length+=$temp1[$i];
        }
      }
    }
    if ($line[8] <= $input4[2] && ($line[8] + $length) >= $input4[2] && !($line[10] =~ "S")) {                       ### fully mapped with no split reads
      $fully_mapped++;
    } elsif ($start1 > 0 && abs($line[8] - $input4[2]) <= 5 && $line[10] =~ "S") {
      if ($start1>=20) {
        $split_reads++;
      } else {
        $split_reads_short++;
      }
    } elsif ($end < $length && $line[10] =~ "S" && abs($line[8] + $end - $start1 - $input4[2]) <= 5) {
      if ($length - $end >= 20) {
        $split_reads++;
      } else {
        $split_reads_short++;
      }
    } elsif ((abs($line[8] - $input4[2]) <=5 || abs(($line[8] + $end - $start1) - $input4[2]) <= 5) && $line[10] =~ "S") {
      $fully_mapped++;
    } else {
      next;
    }
  }
  ##### print out
  my $name=$input4[0]."%".$input4[1]."%".$input4[2];
  my @return4 = ();
  push (@return4,$name);
  push (@return4,$fully_mapped);
  push (@return4,$input4[3]);
  my $GR=0;
  if($input4[4] eq 2) {
    $GR = ($input4[3]/2)/(($input4[3]/2) + $fully_mapped);
  } elsif ($input4[4] eq 1) {
    $GR = $input4[3]/($input4[3] + $fully_mapped);
  }
  push (@return4,$GR);
  push (@return4,$split_reads_short);
  return @return4;
}

#######################
print "                 ########  #######   ##       ##                        ##     ##                \n";
print "                 ########  ##    ##  ##       ##                        ##     ##                \n";
print "                 ##        ##    ##   ##     ##                         ##     ##                \n";
print "                 ##        ##   ##    ##     ##                         ##     ##                \n";
print "                 ########  ######      ##   ##      #####     ####      ##     ##                \n";
print "                 ########  #####       ##   ##     ##   ##   ##  ##     ##     ##                \n";
print "                 ##        ## ##        ## ##     ##        ##    ##    ##     ##                \n";
print "                 ##        ##  ##       ## ##     ##        ##    ##    ##     ##                \n";
print "                 ########  ##   ##       ###       ##   ##   ##   ##    ## ##  ## ##             \n";
print "                 ########  ##    ##      ###        #####     #### ##    ##     ##               \n";
print "\n\n Done!\n";
