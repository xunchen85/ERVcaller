#!usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

my $length=20;
my $file="test_soft_clipping.fastq"; 
my $output="test";

GetOptions('length=n'=>\$length,
           'file=s'=>\$file,
           'output=s'=>\$output
           );

open FILE, "$file";
open OUT1, ">${output}_1sf.fastq";
open OUT2, ">${output}_1sf.others";

my $line1="";
my $line2="";
my $line3="";
my $line4="";

while($line1=<FILE>){
 my @line=split /\|/, $line1;
 $line[3]=~s/chr//;
 $line1=join("|",@line);
 $line2=<FILE>;
 $line3=<FILE>;
 $line4=<FILE>;
 chomp($line1);
 chomp($line2);
 chomp($line3);
 chomp($line4);
 if(length($line2)>=$length){print OUT1 "$line1\n$line2\n$line3\n$line4\n";}
 else{print OUT2 "$line1 $line2 $line3 $line4\n";}
                    }

