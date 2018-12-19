#!usr/bin/perl
#
# Author:       Xun Chen
# Email:        Xun.Chen@uvm.edu
# Date:         11/20/18
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

GetOptions(
           'consensus_TE=s'=>\$consensus_TE,
           'list=s'=>\$list
          );

unless (defined ($consensus_TE)) {
  $consensus_file="Consensus.TE.loci";
  create_consensus_loci($list,$consensus_file);
} else {
  $consensus_file=$consensus_TE;
}


### compare to the consensus TE list;
my %TElist= ();
open CONSENSUS,"$consensus_file";
while (my $tmp1 = <CONSENSUS>){
  chomp($tmp1);
  my @tmp2=split /\t/, $tmp1;
  if ($tmp2[0] =~ "^##" && defined($consensus_TE)) {
    print "$tmp1\n";
  }
  if ($tmp2[0] =~"^#"){next;}
  @tmp2[$#tmp2+1..8] = (0) x (8-$#tmp2);
  $tmp2[8] = "GT";
  my $name = $tmp2[0]."_".$tmp2[1]."_".$tmp2[2];
  @{$TElist{$name}} = @tmp2;
}
my @TE_loci = keys %TElist;


### each sample
open LIST, "$list";
my @tmp_list1 = ();
my $header_count = 0;

while (my $tmp_list1 = <LIST>) {
  my @tmp_list1 = split /\s+/, $tmp_list1;
  my $sampleID = $tmp_list1[0];
  open SAMPLE,"${sampleID}.vcf";
  %sample = ();
  while (my $tmp1=<SAMPLE>){
    @line=split /\s+/, $tmp1;
    if ($header_count eq 0 && $line[0]=~"^##" && !(defined($consensus_TE))) {
      print "$tmp1";
    } elsif ($header_count eq 0 && $line[0]=~"^#CHROM") {
      print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    } elsif ($line[0]=~"^#CHROM") { 
      next;
    }
    push(@{$sample{$line[$chr]}},@line);
  }
  print "\t$sampleID";
  $header_count++;
  my @tmp1=keys %sample;
  my $count1 = 0;
  my $c_start = "";
  my $c_end = "";
  my $c_chr = "";

  for (my $k = 0; $k < @TE_loci; $k++) {
    $count1 = 0;
    my @tmp2 = @{$TElist{$TE_loci[$k]}};
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
      if($tmp1[$j] eq $c_chr && $tmp1[$j+1]>=$c_start && $tmp1[$j+1]<=$c_end){
      $count1++;
      push (@tmp3, @tmp1[$j..$j+9]);
      @tmp1[$j..$j+9] = (0) x 10;
    } elsif($tmp1[$j] eq $c_chr && (abs($tmp1[$j+1]-$c_start)<=$region || abs($tmp1[$j+1]-$c_end)<=$region)){
      $count1++;
      push (@tmp3, @tmp1[$j..$j+9]);
      @tmp1[$j..$j+9] = (0) x 10;
    }
                               }
    if ($count1 eq 0) {
      push (@tmp2, $sampleID);
      push (@tmp2, 0);
      push (@tmp2, @none);
      @{$TElist{$TE_loci[$k]}} = @tmp2;
    } else {
      push (@tmp2, $sampleID);
      push (@tmp2, $count1);
      push (@tmp2, @tmp3);
      @{$TElist{$TE_loci[$k]}} = @tmp2;
      @{$sample{$tmp2[0]}}=@tmp1;
    }
  }
}

print "\n";

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
#  print "@{$TElist{$TE_loci[$i]}}\n";
  %breakpoint_h = ();
  %breakpoint_v = ();
  %TEtype =();
  %samples_multiple_events = ();
  $count_unique=0;
  $count_none=0;
  @info = ();
  $info_status = 0;
  my @tmp4 = @{$TElist{$TE_loci[$i]}};
  my @tmp4_back = @tmp4;
  for (my $j = 9; $j < @tmp4; $j +=12) {
    if ($tmp4[$j+1] >1) {                             ### multiple loci
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
  $tmp4[8] = "GT";
#  print "\n\naaa: @tmp4\n\n\n";

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
#  print "info: @{$kept_info[0]}\n";
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
  $final_info[11] = 0;   ### DR value
  $final_info[12] = $kept_info[0][12];
  $final_info[13] = 0;   ### SR value
  $final_info[14] = $kept_info[0][14];
  $final_info[15] = $kept_info[0][15];
  $final_info[16] = $kept_info[0][16];
  $final_info[17] = 0;   ### NTEDR
  $final_info[18] = $kept_info[0][18];
  $final_info[19] = 0;   ### GR value
 
  for (my $n =0; $n < @kept_info; $n++) {
    $final_info[11] = $final_info[11] + $kept_info[$n][11];
    $final_info[13] = $final_info[13] + $kept_info[$n][13];
    if ($kept_info[$n][17] ne "NA" && $final_info[17] ne "NA") {
      $final_info[17] = $final_info[17] + $kept_info[$n][17];
    } else {
      $final_info[17] = "NA";
    }
    if ($kept_info[$n][19] ne "NA" && $final_info[17] ne "NA") {
      $final_info[19] = $final_info[19].",".$kept_info[$n][19];
    } else {
      $final_info[19] = "NA";
    }
  }

  my $final_info = $final_info[0]."=".$final_info[1].",".$final_info[2].";".$final_info[3]."=".$final_info[4].",".$final_info[5].",".$final_info[6].",".$final_info[7].",".$final_info[8].",".$final_info[9].";".$final_info[10]."=".$final_info[11].";".$final_info[12]."=".$final_info[13].";".$final_info[14]."=".$final_info[15].";".$final_info[16]."=".$final_info[17].";".$final_info[18]."=".$final_info[19];
  $tmp4[7] = $final_info;

  ### print out info and genotype 
  if ($count_unique >0 || %samples_multiple_events) {
    if (defined($consensus_TE)) {
      print join("\t",@tmp4_back[0..8]),"\t";
    } else {
      print join("\t",@tmp4[0..8]),"\t";
    }
    print join("\t",@genotype),"\n";
  }
}

### sub
sub create_consensus_loci {
  my($input1,$input2) = @_;
  my $line="";
  my @line=();
  my %vcf_list=();
  my @sample_list=();

  open SAMPLE, "$input1";
  open CONSENSUS, ">$input2";
  while (my $tmp1=<SAMPLE>) {
    @line=split /\s+/, $tmp1;
    push (@sample_list,$line[0]);
    open F1, "$line[0].vcf";
    while (my $tmp2=<F1>) {
      my @tmp2 = split /\s+/, $tmp2;
      if ($tmp2[0] =~ "\#") {next;}
      push (@{$vcf_list{$tmp2[0]}},$tmp2[1]);
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
    for (my $j = 0; $j < @tmp3_sort; $j++) {
      if (@region) {
        if (abs ($tmp3_sort[$j] - $region[0]) <= $region || abs ($tmp3_sort[$j] - $region[1]) <= $region) {
          $region[1] = $tmp3_sort[$j];
        } else {
          print CONSENSUS "$chr[$i]\t$region[0]\t$region[1]\n";
          @region = ();
        }
      } else {
        $region[0] = $tmp3_sort[$j];
        $region[1] = $tmp3_sort[$j];
      }
    }
  }
  close SAMPLE;
  close CONSENSUS;
}



