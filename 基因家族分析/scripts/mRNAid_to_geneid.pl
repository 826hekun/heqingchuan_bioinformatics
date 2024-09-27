#!/usr/bin/env perl
#北京组学生物科技有限公司
#学习perl语言：
#perl入门：https://zzw.xet.tech/s/1MxMKC
#perl高级：https://zzw.xet.tech/s/JKkbr

use strict;
use Cwd qw(abs_path getcwd);
use Getopt::Long;
use Data::Dumper;


die "perl $0 <gff> <outfile>" unless(@ARGV==2);



my$gff=$ARGV[0];
my%gene=();
my%gene_region=();
my%mRNA2Gene=();
open IN,"$gff" or die "$!; can't open file $gff\n";
open OUT ,">$ARGV[1]" or die "$!; can't open file $ARGV[1]\n";
print OUT "#mRNA_ID\tgene_ID\tchr\tstart\tend\tstrand\n";
while(<IN>){
	chomp;
	next if (/^#/);
	my@tmp=split(/\t/);
	if($tmp[2] =~/^gene$/){
		my($id)=($tmp[8]=~/ID=([^;]+)/);
		$gene{$id}=1;
		$gene_region{$id}=[$tmp[0],$tmp[3],$tmp[4],$tmp[6]];
		
		
		#print "gene:$id\n";
		#my$gene_chr->{$id}=$tmp[0];
	}
	if($tmp[2] =~/^mRNA$/i or $tmp[2] =~/^transcript$/i){
		my($id)=($tmp[8]=~/ID=([^;]+)/);
		my($pid)=($tmp[8]=~/Parent=([^;]+)/);
		
		
		if(exists $gene{$pid}){
			print OUT "$id\t$pid\t";
			print OUT "$tmp[0]\t$tmp[3]\t$tmp[4]\t$tmp[6]\n";
		}
		#print "mRNA:$id\n";
	}
}

close(IN);
close(OUT);


	

