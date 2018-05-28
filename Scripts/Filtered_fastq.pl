#!usr/bin/perl
use strict;

my $line="";
my @line=();
my %type=();
my %integration=();
my $sample=$ARGV[0];
my $fasta_title="";

open TYPE,"$ARGV[0].type";
while(<TYPE>){
 @line=split;
 if(exists($type{$line[0]})){next;}
 else{$type{$line[0]}=$line[1];}
             }
close TYPE;

#################
open INTEGRATION,"$ARGV[0].integration";
while(<INTEGRATION>){
 @line=split;
 if(exists($integration{$line[2]})){next;}
 else{$integration{$line[2]}=$line[0];}
# print "a:$line[2]\t$integration{$line[2]}\n";
                    }
close INTEGRATION;

#################
if (-e $ARGV[0]."_1.1fq") {
 open F1,"$ARGV[0]_1.1fq";
} else {
 open F1,"$ARGV[0].1fq";
}
open R1,"$ARGV[0]_2.1fq";

open F2,">$ARGV[0]_1.1fuq";
open R2,">$ARGV[0]_2.1fuq";

my $f1="";
my $title="";
my $r1="";
#############

open FASTA,">$ARGV[0].fasta";
while($f1=<F1>){
 $title=$f1;
 chomp($title);
 $fasta_title=$f1;
 $fasta_title=~s/^@/>/;
 $title=~s/^\@//;
 $title=~s/\/1$//;
 if(exists($integration{$title})){    # only keep integration reads;
  print F2 "$f1";
  $f1=<F1>;print F2 "$f1";
  print FASTA "$fasta_title$f1";
  $f1=<F1>;print F2 "$f1";
  $f1=<F1>;print F2 "$f1";

  $r1=<R1>;print R2 "$r1";
  $fasta_title=$r1;
  $fasta_title=~s/^@/>/;
  $r1=<R1>;print R2 "$r1";
  print FASTA "$fasta_title$r1";
  $r1=<R1>;print R2 "$r1";
  $r1=<R1>;print R2 "$r1";
                                                                   }
 else{
  $f1=<F1>;$f1=<F1>;$f1=<F1>;
  $r1=<R1>;$r1=<R1>;$r1=<R1>;$r1=<R1>;
     }
               }
close F1;
close F2;
close R1;
close R2;
################################################ remove softclip reads;

$line="";
@line=();
my $title="";
my %fq=();
my %read_n=();
my %other=();
my $new="";
my $read_name="";

#########
open FQ,"$ARGV[0]_1sf.fastq";
open FQ2,">$ARGV[0]_1sf.fuq";
while(<FQ>){
 @line=split /\|/, $_;
 $title=$line[3]."|".$line[4]."|".$line[6]."|".$line[8];
 $read_name=$line[1];
 $fasta_title=">".$line[1]."|".$line[2]."|".$line[8]."\n";
 unless(exists($integration{$read_name})){$line=<FQ>;$line=<FQ>;$line=<FQ>;next;}
 if(exists($fq{$title})){
  $new="@".$line[1]."|".$line[2]."|".$line[8]."\n";
  $line=<FQ>;$new.=$line;
  $fasta_title.=$line;
  my $len1=length($line);
  $line=<FQ>;$new.=$line;
  $line=<FQ>;$new.=$line;
  my @tem=split /\s+/,$fq{$title};
  my $len2=length{$tem[1]};
  if($new && $len1>$len2){
   $read_n{$title}=$read_name;
   $fq{$title}=$new;
   if(exists($integration{$read_name})){print FQ2 "$fq{$title}";print FASTA "$fasta_title";}
                         }
                        }
 else{
  $fq{$title}=$_;
  $fq{$title}="@".$line[1]."|".$line[2]."|".$line[8]."\n";
  $line=<FQ>;$fq{$title}.=$line;
  $fasta_title.=$line;
  $line=<FQ>;$fq{$title}.=$line;
  $line=<FQ>;$fq{$title}.=$line;
  $read_n{$title}=$read_name;
  if(exists($integration{$read_name})){print FQ2 "$fq{$title}";print FASTA "$fasta_title";}
     }
            }
close FQ;
close FQ2;

#########
open OTHER,"$ARGV[0]_1sf.others";
while(<OTHER>){
 my @tmp2=split;
 @line=split /\|/, $tmp2[0];
 my $len1=length($tmp2[1]);
 $title=$line[3]."|".$line[4]."|".$line[6]."|".$line[8];
 if(exists($fq{$title})){
  my @tem=split /\s+/,$fq{$title};
  my $len2=length{$tem[1]};
  if($new && $len1>$len2){$fq{$title}=$new;$other{$title}="y";}
                        }
 else{
  $fq{$title}=$_;
  $other{$title}="y";
     }
               }
close OTHER;

###########
open OTHER2,">$ARGV[0]_1sf.othu";
my @name=keys %other;
for(my $i=0;$i<@name;$i++){print OTHER2 "$fq{$name[$i]}";}
close OTHER2;

