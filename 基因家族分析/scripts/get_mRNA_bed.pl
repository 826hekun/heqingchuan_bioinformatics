#!/usr/bin/env perl
#北京组学生物科技有限公司
#学习perl语言：
#perl入门：https://zzw.xet.tech/s/1MxMKC
#perl高级：https://zzw.xet.tech/s/JKkbr

use Getopt::Long;
my %opts;
use Data::Dumper;
GetOptions( \%opts, "in1=s","in2=s", "out=s", "h" );
if ( !defined( $opts{in1} ) || !defined( $opts{out} ) || defined( $opts{h} ) ) {
	&USAGE;
}

open( IN2, "$opts{in2}" )  || die "open $opts{in2} failed\n";


my %mRNAID=();

while(<IN2>){
	chomp;
	my@tmp=split(/\s+/);
	$mRNAID{$tmp[0]}=1;

}

#print Dumper(\%mRNAID);
close(IN2);




open( IN1, "$opts{in1}" )  || die "open $opts{in1} failed\n";
open( OUT, ">$opts{out}" ) || die "open $opts{out} failed\n";




while (<IN1>) {
	chomp;
	my @a = split /\t/, $_;
	if ($a[2] =~/^mRNA$/i or $a[2] =~/^transcript$/i ) {
		#if ($a[2] eq "mRNA") {
		$a[8] =~ m/ID=([^;]*)/;    #注意这里匹配ID的信息
		$id = $1;
		if (exists $mRNAID{$id}) {
			print OUT "$a[0]\t$a[3]\t$a[4]\t$id\t$a[7]\t$a[6]\n";
		}

	}

}
close OUT;
close IN1;


sub USAGE {
	print "usage: perl $0 -in1  gff  -in2 id.txt -out gene_location.bed ";
	exit;
}
