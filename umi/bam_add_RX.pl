#!/usr/bin/perl -w
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;  
use Cwd 'abs_path';
use File::Basename; 
use Data::Dumper; 
use FileHandle; 
my $usage=<<USAGE;
	Usage:
		perl $0 [options]
			-b   --bam	<string>	sort.bam
			-U1  --UMI1	<string>     	umifile1
			-U2  --UMI2	<string>	umifile2
			-o   --out	<string>	outfile
			-h  --help			print help information and exit

	Example:
		perl $0 -b sample.sort.bam -U1 umifile1	-U2 umifile2 -o sample.sort.withUMI.bam

USAGE
#=============global variants=============
my ($bam,$UMI_1,$UMI_2,$out,$help);
my %hash_umi;
#=========================================
GetOptions(
	"bam|b=s"=>\$bam,
	"UMI1|U1=s"=>\$UMI_1,
	"UMI2|U2=s"=>\$UMI_2,
	"out|o=s"=>\$out,
	"help|h:s"=>\$help,

);

#$outdir ||= `pwd`;
if(!$bam || !$UMI_1 || !$UMI_2 || $help ){
        die "$usage";
}


open FQ1, "gunzip -c $UMI_1 |" or die $!;
open FQ2, "gunzip -c $UMI_2 |" or die $!;
#open UMI,">$umi";
while (my $id1 = <FQ1> ){ 
	#my ($seq1, $plus1, $qual1) = (<FQ1>, <FQ1>, <FQ1>);
	#my ($id2, $seq2, $plus2, $qual2) = (<FQ2>, <FQ2>, <FQ2>, <FQ2>);
	my $seq1 = <FQ1> ;
	my $plus1= <FQ1> ;
	my $qual1= <FQ1> ;
	my $id2  = <FQ2> ;
	my $seq2 = <FQ2> ;
        my $plus2= <FQ2> ;
        my $qual2= <FQ2> ;

	chomp ($id1,$seq1,$plus1,$qual1,$id2, $seq2, $plus2, $qual2);
	
	#$id1=~s/(@|\/1)//;
	$id1=~s/@//;
	$id1=~s/\/1//;
	my @id1_tmp=split /\s/,$id1;
	$id1=$id1_tmp[0];
	$hash_umi{$id1}="$seq1";


	#$id2=~s/(@|\/1)//;
        $id2=~s/@//;
        $id2=~s/\/2//;
	my @id2_tmp=split /\s/,$id2;
	$id2=$id2_tmp[0];

	if (exists $hash_umi{$id2}){
                $hash_umi{$id2}.="-$seq2";
         }
	
}


#foreach my $i (keys %hash_umi){
#	print UMI "$i\t$hash_umi{$i}\n";}

open IN, "samtools view --threads 16 -h  $bam | " or die $!;
open OUT,"| samtools view --threads 16 -hbS -o $out" or die $!;
while (<IN>){
      	chomp;
        my $line = $_;	
        my @temp=split /\t/,$_;
	my $readsname=$temp[0];
	if (exists $hash_umi{$readsname}){
		$line.= "\tRX:Z:$hash_umi{$readsname}";
	}
	
	print OUT "$line\n";
}

close FQ1;
close FQ2;
#close UMI;
close IN;
close OUT;
#system ("samtools view -b -S $out -o ZLC210232A.sort.withUMI.bam");
print "####add UMI is OK ###";
