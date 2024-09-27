#!/usr/bin/env perl
#北京组学生物科技有限公司
#学习perl语言：
#perl入门：https://zzw.xet.tech/s/1MxMKC
#perl高级：https://zzw.xet.tech/s/JKkbr

=head1 NAME

convert_RNAmmer_to_gff3.pl - convert bad gff2 format RNAmmer output to good GFF3

=head1 SYNOPSIS

USAGE: template.pl 
            --input=/path/to/some_file.out 

=head1 OPTIONS

B<--hmmoutfile,-i>
    hmmsearch out file 
B<--fa,-f> 
    fasta file path
B<--evalue,-e> 
    set evalue cut off  ，
B<--domain,-d> 
    This domain’s number to subset seq
B<--outfile,-o> 
    output file path
    
B<--help,-h>
    This help message

=head1  DESCRIPTION

File converter

=head1  INPUT

Input above.

=head1  OUTPUT

GFF3 to STDOUT

=head1  CONTACT

    huangls
    huangls@biomics.com.cn

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Bio::SeqIO;
use Bio::Seq;
my %options = ();
my $results = GetOptions (\%options, 
                          'hmmoutfile|i=s',
                          'fa|f=s',
                          'outfile|o=s',
                          'evalue|e=s',
                           'domain|d=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

#die "perl $0 <hmmoutfile> <fa> <OUT> <E-value>" unless ( @ARGV == 4 );

$in = Bio::SeqIO->new(
	-file   => "$ARGV[1]",
	-format => 'Fasta'
);
$out = Bio::SeqIO->new(
	-file   => ">$ARGV[2]",
	-format => 'Fasta'
);
my %keep = ();
open IN, "$ARGV[0]" or die "$! ; can't open file:  $ARGV[0]\n";

my %domain;

while (<IN>) {
	chomp;
	next if /^#/;

	my @a = split /\s+/;
	next if $a[6] > $ARGV[3];
	
	my @b = ( $a[17], $a[18] );  #结构域在序列中的位置
	my $keys = $a[0];
	$keep{$keys}{$a[9]} = \@b;#结构域在序列中的位置

}
close(IN);
while ( my $seq = $in->next_seq() ) {
	my ( $id, $sequence, $desc ) = ( $seq->id, $seq->seq, $seq->desc );

	if ( exists $keep{$id} ) {
		for my $d(keys %{$keep{$id}}){
			my $subseq = $seq->subseq( $keep{$id}{$d}->[0], $keep{$id}{$d}->[1]); #截取序列
			my $newseqobj = Bio::Seq->new(
				-seq  => $subseq,
				-desc => "domain:$keep{$id}{$d}[0]-$keep{$id}{$d}[1]",
				-id   => "$id.domain$d",
			);

			$out->write_seq($newseqobj);
		}
	}
}
$in->close();
$out->close();
