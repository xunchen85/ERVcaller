#!/usr/bin/env perl

#
# Author:       Xun Chen
# Email:        Xun.Chen@uvm.edu
# Date:         11/20/18
#
#
# Updates:	02/03/2019	Use the newer VCF format output by ERVcaller_v1.4.pl;
# Updates:	02/03/2019	Use hash table to speed up the combine process;
#

use strict;
use Getopt::Long;
my $line = "";
my @line = ();
my %sample = ();
my $chr = 0;
my $start = 1;
my $end = 2;
my $region = 500;
my @none = (0) x 10;
my $consensus_TE;
my $consensus_file="";
my $list="sample_list";
my $mode="";
my $output;

GetOptions(
           'consensus_TE=s'=>\$consensus_TE,
           'list=s'=>\$list,
           'output=s'=>\$output
          );

unless (defined ($consensus_TE)) {
  $consensus_file=$list.".consensus";
  create_consensus_loci($list,$consensus_file);
} else {
  $consensus_file=$consensus_TE;
}

### compare to the consensus TE list;
my %TElist = ();
my %TE_position = ();

open CONSENSUS, "$consensus_file";
open OUTPUT, ">$output";

while (my $tmp1 = <CONSENSUS>) {
  chomp($tmp1);
  my @tmp2=split /\t/, $tmp1;
  if ($tmp2[0] =~ "^##" && defined($consensus_TE)) {
    print OUTPUT "$tmp1\n";
  }
  if ($tmp2[0] =~"^#"){
    next;
  }
  @tmp2[$#tmp2+1..8] = (0) x (8-$#tmp2);
  $tmp2[8] = "GT:GQ:GL:DPN:DPI";
  my $name = "";
  if (defined($consensus_TE)) {
    $name = $tmp2[0]."|".$tmp2[1]."|".$tmp2[2];
    $TElist{$name} = join (" ",@tmp2);
#    @{$TElist{$name}} = @tmp2;
  } else {
    my @te_type = split /\;/, $tmp2[4];
    my @te_position = split /\;/, $tmp2[3];
    for (my $m = 0; $m < @te_type; $m++) {
      $name = $tmp2[0]."|".$tmp2[1]."|".$tmp2[2]."|".$te_type[$m]; 
      my @tmp2_2 = @tmp2;
      $tmp2_2[4] = $te_type[$m];
      $TElist{$name} = join (" ",@tmp2_2);
#      @{$TElist{$name}} = @tmp2_2;
      for (my $n = 0; $n < @te_position; $n++){
	my $name2 = $tmp2[0]."|".$te_position[$n]."|".$te_type[$m];
        push(@{$TE_position{$name}},$name2);
      }
    }
  }
}
my @TE_loci = keys %TElist;

### each sample
open LIST, "$list";
my @tmp_list1 = ();
my $header_count = 0;
my %each_sample = ();

while (my $tmp_list1 = <LIST>) {
  print "sample_ID: $tmp_list1\n";
  ##### read each sample
  my @tmp_list1 = split /\s+/, $tmp_list1;
  my $sampleID = $tmp_list1[0];
  open SAMPLE,"${sampleID}.vcf";
  %sample = ();
  %each_sample = ();
  ##### read into a table and print the header
  while (my $tmp1 = <SAMPLE>) {
    @line=split /\s+/, $tmp1;
    if ($header_count eq 0 && $line[0]=~"^##" && !(defined($consensus_TE))) {
      print OUTPUT "$tmp1";
    } elsif ($header_count eq 0 && $line[0]=~"^#CHROM") {
      print OUTPUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    } elsif ($line[0]=~"^#CHROM") { 
      next;
    }
    if ($line[9] =~ "NA" || $line[9] =~ "NaN") {
      $line[9] =~ s/NA/0/g;
      $line[9] =~ s/NaN/0/g;
    }
    push (@{$sample{$line[$chr]}},@line);
    my $name4 = $line[0]."|".$line[1]."|".$line[4];
    @{$each_sample{$name4}} = @line;
  }
  print OUTPUT "\t$sampleID";
  $header_count++;
  my @tmp1=keys %sample;
  my $count1 = 0;
  my $c_start = "";
  my $c_end = "";
  my $c_chr = "";

  ##### compare to the consensus TE loci;
  if (!(defined($consensus_TE))) {
    for (my $k = 0; $k < @TE_loci; $k++) {
      my @tmp1_location = @{$TE_position{$TE_loci[$k]}};
      $count1 = 0;
      my @tmp3 = ();
      for (my $p = 0; $p < @tmp1_location; $p++) {
	if (exists ($each_sample{$tmp1_location[$p]})) {
	  push (@tmp3,@{$each_sample{$tmp1_location[$p]}});
          $count1++;
        }
      }
      if ($count1 eq 0) {
        my $tmp4 = join (" ",@none);
        $TElist{$TE_loci[$k]} = $TElist{$TE_loci[$k]}." ".$sampleID." 0 ".$tmp4; 
#        push (@{$TElist{$TE_loci[$k]}},$sampleID);
#        push (@{$TElist{$TE_loci[$k]}}, 0);
#        push (@{$TElist{$TE_loci[$k]}}, @none);
      } else {
        my $tmp4 = join (" ",@tmp3);
        $TElist{$TE_loci[$k]} = $TElist{$TE_loci[$k]}." ".$sampleID." ".$count1." ".$tmp4;
#        push (@{$TElist{$TE_loci[$k]}}, $sampleID);
#        push (@{$TElist{$TE_loci[$k]}}, $count1);
#        push (@{$TElist{$TE_loci[$k]}}, @tmp3);
      }
    }
  } else {
    for (my $k = 0; $k < @TE_loci; $k++) {
      $count1 = 0;
      my @tmp2 = split /\s+/, $TElist{$TE_loci[$k]};
#      my @tmp2 = @{$TElist{$TE_loci[$k]}};
      $c_chr = $tmp2[0];
      if ($tmp2[1] =~ /^[0-9]+$/ && $tmp2[2] =~ /^[0-9]+$/) {
        $c_start = $tmp2[1];
        $c_end = $tmp2[2];
      } else {
        $c_start = $tmp2[1];
        $c_end = $tmp2[1];
      }
      if (exists($sample{$tmp2[0]})){
        @tmp1 = @{$sample{$tmp2[0]}};
      }
      my @tmp3 = ();
      for(my $j=0;$j<@tmp1;$j+=10){       ####### each event
        if($tmp1[$j] eq $c_chr && $tmp1[$j+1]>=$c_start && $tmp1[$j+1]<=$c_end) {
          $count1++;
          push (@tmp3, @tmp1[$j..$j+9]);
          @tmp1[$j..$j+9] = (0) x 10;
        } elsif($tmp1[$j] eq $c_chr && (abs($tmp1[$j+1]-$c_start)<=$region || abs($tmp1[$j+1]-$c_end)<=$region)) {
          $count1++;
          push (@tmp3, @tmp1[$j..$j+9]);
          @tmp1[$j..$j+9] = (0) x 10;
        }
      }
      if ($count1 eq 0) {
        push (@tmp2, $sampleID);
        push (@tmp2, 0);
        push (@tmp2, @none);
#       @{$TElist{$TE_loci[$k]}} = @tmp2;
        $TElist{$TE_loci[$k]} = join (" ", @tmp2);
      } else {
        push (@tmp2, $sampleID);
        push (@tmp2, $count1);
        push (@tmp2, @tmp3);
        $TElist{$TE_loci[$k]} = join (" ", @tmp2);
#        @{$TElist{$TE_loci[$k]}} = @tmp2;
        @{$sample{$tmp2[0]}} = @tmp1;
      }
    }
  }

}

print OUTPUT "\n";

### Organize the output of genotypes
my %breakpoint_h = ();
my %breakpoint_v = ();
my %TEtype = ();
my %samples_multiple_events = ();
my $count_unique=0;
my $count_none=0;
my @info = ();
my $info_status = 0;


for (my $i = 0; $i < @TE_loci; $i++) {
  %breakpoint_h = ();
  %breakpoint_v = ();
  %TEtype =();
  %samples_multiple_events = ();
  $count_unique=0;
  $count_none=0;
  @info = ();
  $info_status = 0;
#  my @tmp4 = @{$TElist{$TE_loci[$i]}};
  my @tmp4 = split /\s+/, $TElist{$TE_loci[$i]};
  my @tmp4_back = @tmp4;
  for (my $j = 9; $j < @tmp4; $j +=12) {
    if ($tmp4[$j+1] > 1) {                            ### multiple loci
      @{$samples_multiple_events{$tmp4[$j]}} = @tmp4[$j..$j+12+($tmp4[$j+1]-1)*10-1];
      $j = $j + 10 * ($tmp4[$j+1] - 1);
    } elsif ($tmp4[$j+1] eq 0) {
      $count_none++;
    } elsif ($tmp4[$j+1] eq 1) {                      ### unique loci
      $count_unique++;
      if (exists ($breakpoint_h{$tmp4[$j+3]})) {      ### breakpoint_h
        $breakpoint_h{$tmp4[$j+3]}++;
      } else {
        $breakpoint_h{$tmp4[$j+3]}=1;
      }
      if (exists ($TEtype{$tmp4[$j+6]})) {            ### TEtype
        $TEtype{$tmp4[$j+6]}++;
      } else {
        $TEtype{$tmp4[$j+6]}=1;
      }  
    }
  }

  ### extracting infomation
  my @TE_tmp1 = ();
  my @genotype = ();
  my @kept_info = ();
  if (%TEtype) {
    @TE_tmp1 = sort {$TEtype{$a} <=> $TEtype{$b}} keys %TEtype;
    if (@TE_tmp1 >1) {
      $tmp4[4] = "<INS_MEI:Complex>";
    } else {
      $tmp4[4] = $TE_tmp1[0];
    }
  } elsif ($count_unique eq 0 && %samples_multiple_events) {
      $tmp4[4] = "<INS_MEI:Complex>";
  }
  my @breakpoint_h = ();
  if (%breakpoint_h) {
    @breakpoint_h = sort {$breakpoint_h{$b} <=> $breakpoint_h{$a}} keys %breakpoint_h;
    $tmp4[1] = $breakpoint_h[0];
  }
  $tmp4[2] = ".";
  $tmp4[5] = ".";
  $tmp4[6] = ".";
  $tmp4[8] = "GT:GQ:GL:DPN:DPI";

  for (my $j = 9, my $k =0; $j < @tmp4; $j +=12,$k++) {
    @info = ();
    if ($tmp4[$j+1] >1) {
      my @tmp5_info = ();
      my $m = $j;
      for ($m+=2; $m < $j+2+(($tmp4[$j+1])*10); $m+=10) {
        if (@tmp5_info) {
          if ($tmp4[$m+4] eq $tmp4[4]) {
            my $tmp6 = $tmp4[$m+7];
            $tmp6 =~ s/\,/ /g;
            $tmp6 =~ s/\=/ /g;
	    $tmp6 =~ s/\;/ /g;
            @tmp5_info = split /\s+/, $tmp6;
            $genotype[$k] = $tmp4[$m+9];
	    if ($count_unique eq 0) {
	      $tmp4[1] = $tmp4[$m+1];
	      $tmp4[3] = $tmp4[$m+3];
	    }
          } elsif (@TE_tmp1) {
            if ($tmp4[$m+4] eq $TE_tmp1[0]) {
	      my $tmp6 = $tmp4[$m+7];
              $tmp6 =~ s/\,/ /g;
              $tmp6 =~ s/\=/ /g;
	      $tmp6 =~ s/\;/ /g;
              @tmp5_info = split /\s+/, $tmp6;
	      $genotype[$k] = $tmp4[$m+9];
              if ($count_unique eq 0) {
                $tmp4[1] = $tmp4[$m+1];
	        $tmp4[3] = $tmp4[$m+3];
              }
            } 
          }
        } else {
	  my $tmp6 = $tmp4[$m+7];
	  $tmp6 =~ s/\,/ /g;
	  $tmp6 =~ s/\=/ /g;
	  $tmp6 =~ s/\;/ /g;
          @tmp5_info = split /\s+/, $tmp6;
          $genotype[$k] = $tmp4[$m+9];
          if ($count_unique eq 0) {
            $tmp4[1] = $tmp4[$m+1];
	    $tmp4[3] = $tmp4[$m+3];
          }
        }
      }
      $j = $m-12;
      push (@kept_info,[@tmp5_info]); 
    } elsif ($tmp4[$j+1] eq 0) {
      $genotype[$k] = "0/0";
    } elsif ($tmp4[$j+1] eq 1) {
      if ($tmp4[1] eq $tmp4[$j+3]) {
	$tmp4[3] = $tmp4[$j+5];
      }
      $genotype[$k] = $tmp4[$j+11];
      my $tmp6 = $tmp4[$j+9];
      $tmp6 =~ s/\,/ /g;
      $tmp6 =~ s/\=/ /g;
      $tmp6 =~ s/\;/ /g;
      @info = split /\s+/, $tmp6;
      push (@kept_info,[@info]);
    }
  }

  ### organize info 
  my @final_info = ();
  my @kept_info=sort { $b->[9] <=> $a->[9] } @kept_info;
  $final_info[0] = $kept_info[0][0];
  $final_info[1] = $kept_info[0][1];
  $final_info[2] = $kept_info[0][2];
  $final_info[3] = $kept_info[0][3];
  $final_info[4] = $kept_info[0][4];
  $final_info[5] = $kept_info[0][5];
  $final_info[6] = $kept_info[0][6];
  $final_info[7] = $kept_info[0][7];
  $final_info[8] = $kept_info[0][8];
  $final_info[9] = $kept_info[0][9];
  $final_info[10] = $kept_info[0][10];
  $final_info[11] = 0;   ### CR value
  $final_info[12] = $kept_info[0][12];
  $final_info[13] = 0;   ### SR value
  $final_info[14] = $kept_info[0][14];
  $final_info[15] = $kept_info[0][15];
  $final_info[16] = $kept_info[0][16];
  $final_info[17] = 0;   ### GR value
 
  for (my $n =0; $n < @kept_info; $n++) {
    $final_info[11] = $final_info[11] + $kept_info[$n][11];
    $final_info[13] = $final_info[13] + $kept_info[$n][13];
    if ($kept_info[$n][19] ne "NULL" && $final_info[17] ne "NULL") {
      $final_info[17] = $final_info[17].",".$kept_info[$n][17];
    } else {
      $final_info[17] = "NULL";
    }
  }

  my $final_info = $final_info[0]."=".$final_info[1].",".$final_info[2].";".$final_info[3]."=".$final_info[4].",".$final_info[5].",".$final_info[6].",".$final_info[7].",".$final_info[8].",".$final_info[9].";".$final_info[10]."=".$final_info[11].";".$final_info[12]."=".$final_info[13].";".$final_info[14]."=".$final_info[15].";".$final_info[16]."=".$final_info[17];
  $tmp4[7] = $final_info;

  ### print out info and genotype 
  if ($count_unique >0 || %samples_multiple_events) {
    if (defined($consensus_TE)) {
      print OUTPUT join("\t",@tmp4_back[0..8]),"\t";
    } else {
      print OUTPUT join("\t",@tmp4[0..8]),"\t";
    }
    print OUTPUT join("\t",@genotype),"\n";
  }
}
close OUTPUT;

### sub
sub create_consensus_loci {
  my($input1,$input2) = @_;
  my $line = "";
  my @line = ();
  my %vcf_list = ();
  my %location_type = ();
  my @sample_list = ();

  open SAMPLE, "$input1";
  open CONSENSUS, ">$input2";
  while (my $tmp1 = <SAMPLE>) {
    @line=split /\s+/, $tmp1;
    push (@sample_list,$line[0]);
    open F1, "$line[0].vcf";
    while (my $tmp2 = <F1>) {
      my @tmp2 = split /\s+/, $tmp2;
      if ($tmp2[0] =~ "\#") {
        next;
      }
      push (@{$vcf_list{$tmp2[0]}},$tmp2[1]);
      my $tmp1_name = $tmp2[0]."|".$tmp2[1];
      push (@{$location_type{$tmp1_name}},$tmp2[4]);
    }
    close F1;
  }
  close SAMPLE;

  my @chr = keys %vcf_list;
     @chr = sort {$a <=> $b} @chr;
  my @region = ();

  for (my $i = 0; $i < @chr; $i++) {
    my @tmp3 = @{$vcf_list{$chr[$i]}};
    my @tmp3_sort = sort {$a <=> $b} @tmp3;
    @region = ();                                  #### corrected on 02/03/2019
    for (my $j = 0; $j < @tmp3_sort; $j++) {
      if (@region) {
        if (abs($tmp3_sort[$j] - $region[0]) <= $region || abs($tmp3_sort[$j] - $region[1]) <= $region) {
          $region[1] = $tmp3_sort[$j];
          $region[2] = $region[2].";".$tmp3_sort[$j];
          my $tmp1_name = $chr[$i]."|".$tmp3_sort[$j];
          my $tmp1_type = join(";",@{$location_type{$tmp1_name}});
          $region[3] = $region[3].";".$tmp1_type;
	  if ($j == $#tmp3_sort) {                 #### corrected on 01/04/2024
            ##### location
            my @tmp2_location = split /\;/, $region[2];
            @tmp2_location = sort(@tmp2_location);
            my $tmp2_location_sort = "";
            if (@tmp2_location) {
              $tmp2_location_sort = $tmp2_location[0];
            } 
            for (my $k = 1; $k < @tmp2_location; $k++) {
              if ($tmp2_location[$k] eq $tmp2_location[$k-1]) {
                next;
              } else {
                $tmp2_location_sort = $tmp2_location_sort.";".$tmp2_location[$k];
              } 
            } 
          
            ##### type
            my @tmp2_type = split /\;/, $region[3];
            @tmp2_type = sort(@tmp2_type);
            my $tmp2_type_sort = "";
            if (@tmp2_type) {
              $tmp2_type_sort = $tmp2_type[0];
            } 
            for (my $k = 1; $k < @tmp2_type; $k++) {
              if ($tmp2_type[$k] eq $tmp2_type[$k-1]) {
                next;
              } else {
                $tmp2_type_sort = $tmp2_type_sort.";".$tmp2_type[$k];
              } 
            } 
            print CONSENSUS "$chr[$i]\t$region[0]\t$region[1]\t$tmp2_location_sort\t$tmp2_type_sort\n";
	  }
        } else {
          ##### location
          my @tmp2_location = split /\;/, $region[2];
	  @tmp2_location = sort(@tmp2_location);
	  my $tmp2_location_sort = "";
          if (@tmp2_location) {
	    $tmp2_location_sort = $tmp2_location[0];
          }
          for (my $k = 1; $k < @tmp2_location; $k++) {
	    if ($tmp2_location[$k] eq $tmp2_location[$k-1]) {
	      next;
            } else {
	      $tmp2_location_sort = $tmp2_location_sort.";".$tmp2_location[$k];
	    }
          }

          ##### type
          my @tmp2_type = split /\;/, $region[3];
          @tmp2_type = sort(@tmp2_type);
          my $tmp2_type_sort = "";
	  if (@tmp2_type) {
            $tmp2_type_sort = $tmp2_type[0];
          }
          for (my $k = 1; $k < @tmp2_type; $k++) {
	    if ($tmp2_type[$k] eq $tmp2_type[$k-1]) {
	      next;
            } else {
	      $tmp2_type_sort = $tmp2_type_sort.";".$tmp2_type[$k];
            }
          }
          print CONSENSUS "$chr[$i]\t$region[0]\t$region[1]\t$tmp2_location_sort\t$tmp2_type_sort\n";
          @region = ();
	  $region[0] = $tmp3_sort[$j];             #### corrected on 01/24/2019
          $region[1] = $tmp3_sort[$j];
          $region[2] = $tmp3_sort[$j];
          my $tmp1_name = $chr[$i]."|".$tmp3_sort[$j];
          my $tmp1_type = join(";",@{$location_type{$tmp1_name}});
          $region[3] = $tmp1_type;
        }
      } else {
        $region[0] = $tmp3_sort[$j];
        $region[1] = $tmp3_sort[$j];
        $region[2] = $tmp3_sort[$j];
        my $tmp1_name = $chr[$i]."|".$tmp3_sort[$j];
        my $tmp1_type = join(";",@{$location_type{$tmp1_name}});
        $region[3] = $tmp1_type;
	if ($#tmp3_sort == 0) {                   #### corrected on 01/04/2024
          print CONSENSUS "$chr[$i]\t$region[0]\t$region[1]\t$region[2]\t$region[3]\n";
	}
      }
    }
  }
  close SAMPLE;
  close CONSENSUS;
}



