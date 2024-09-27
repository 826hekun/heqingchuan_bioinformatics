#!/usr/bin/env perl
#北京组学生物科技有限公司
#学习perl语言：
#perl入门：https://zzw.xet.tech/s/1MxMKC
#perl高级：https://zzw.xet.tech/s/JKkbr

=head1 NAME

get_gene_longest_fa_for_OrthoFinder.pl

=head1 DESCRIPTION

The script aims to get gene longest  cds or pep sequence from gff and fasta file .

=head1 SYNOPSIS

perl get_gene_longest_fa_for_OrthoFinder.pl  --fa cds.fa --gff genome.gff -l 150 -p cds.longest

perl get_gene_longest_fa_for_OrthoFinder.pl  --fa pep.fa --gff genome.gff -l 50 -p pep.longest


=head1 OPTIONS

The options are:

=over 4

=item --help,-h

This help message

=item --fa,-f

set pep.fa or cds.fa file path, gz file supported, required

=item --gff,-g

gff file, gz file support. required

=item --length,-l

seq length cut off , optional

=item  --od,-o

set output directory , optional default current working directory

=item --prefix,-p

set output file name prefix , optional default seq

=back

=head1  CONTACT

huangls
huangls@biomics.com.cn

=head1 COPYRIGHT

Copyright  huangls, www.biomics.com.cn , www.omicsclass.com

=cut


use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Bio::SeqIO;
use Bio::Seq;
use PerlIO::gzip;
use Cwd qw(abs_path getcwd);

my ($fa,$od,$prefix,$gff,$help,$length);
GetOptions(
        "help|h" =>\$help,
        "fa|f=s"=>\$fa,
	"gff|g=s"=>\$gff,
	"length|l=i"=>\$length,
	"prefix|p=s"=>\$prefix,
        "outdir|o=s"=>\$od
                ) or pod2usage(1);

## display documentation
if( $help ){
    pod2usage( {-exitval => 0, -verbose => 2, -output=>\*STDERR} );
}

#if( @ARGV==0 ){
#    pod2usage( {-exitval => 1, -verbose => 2, -message => "Failed to parse command line\n",-output=>\*STDERR} );
#}

unless(  $gff ){
    pod2usage( {-exitval => 1, -verbose => 2, -message => "Must specify  gff3 file (--gff)\n",-output=>\*STDERR} );
}

unless( $fa  ){
    pod2usage( {-exitval => 1, -verbose => 2, -message => "Must specify  fa file (--fa)\n", -output=>\*STDERR} );
}

######################
$od||=getcwd;
$od=abs_path($od);
unless(-d $od){	`mkdir $od`;}
$prefix||="seq";

#####################

if($gff=~/.gz$/){
  open IN,"<:gzip","$gff" or die "$!; can't open file $gff\n";

}else{
  open IN,"$gff" or die "$!; can't open file $gff\n";

}

#读取GFF文件，存储基因ID与转录本ID之间的关系
my%gff;
while(<IN>){
	chomp;
	next if /^#/;
	next if /^\s*$/;
	my@a=split(/\t/);	
	if($a[2] eq "mRNA" || $a[2] eq "transcript"){
	  
	  if($a[8]=~/biotype/i){
	    my($biotype)=($a[8]=~ m/biotype=([^;]*)/i);
	    if($biotype ne "protein_coding"){
	      next;
	    }
	  }
	  
		my($mRNAid)=($a[8]=~ m/ID=([^;]*)/);
		my($geneid)=($a[8]=~ m/Parent=([^;]*)/);
		$gff{$mRNAid}=$geneid;
		}	
}
close IN;


##### 读取蛋白或者cds文件，选取最长的序列作为基因的代表序列
my $in;

if($fa=~/.gz$/){
  open my $FI, "<:gzip", "$fa" or die "$!";
  $in  = Bio::SeqIO->new(-fh => $FI ,
                               -format => 'Fasta');
}else{
  $in  = Bio::SeqIO->new(-file => "$fa" ,
                               -format => 'Fasta');
}

my%fa=();
while ( my $seq = $in->next_seq() ) {
	my($id,$s,$desc,$len)=($seq->id,$seq->seq,$seq->desc,$seq->length);
	my $gene_id="";
	if(exists $gff{$id}){
		$gene_id=$gff{$id};
	}else{
		print STDERR "sequence ID: $id not found in gff file, please check\n";
	next;

	}
	if(exists $fa{$gene_id}){
		if($fa{$gene_id}->length < $len){
			$fa{$gene_id}=$seq;
		}
	}else{
		$fa{$gene_id}=$seq;
	}
}


#输出序列

                               
my $out = Bio::SeqIO->new(-file => ">$od/$prefix.fasta" ,
                               -format => 'Fasta');

for my$id(sort keys %fa){
	my $len=$fa{$id}->length;
	$len=$len-1 if($fa{$id}->seq =~/\*$/);
	
	if ($length){
	  if($len<=$length){
	    next;
	  }
	}
	my$seq_obj=Bio::Seq->new(-seq =>$fa{$id}->seq,
							-desc=>$id."=".$len,
							-id=>$id
							 );
	$out->write_seq($seq_obj);
}
$in->close();
$out->close();
