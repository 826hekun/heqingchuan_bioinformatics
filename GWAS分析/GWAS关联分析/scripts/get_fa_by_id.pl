#!/usr/bin/env perl
#北京组学生物科技有限公司
#学习perl语言：
#perl入门：https://zzw.xet.tech/s/1MxMKC
#perl高级：https://zzw.xet.tech/s/JKkbr


die "perl $0 <idlist> <fa> <OUT>" unless ( @ARGV == 3 );
use Math::BigFloat;
use Bio::SeqIO;
use Bio::Seq;
use warnings;
#读入蛋白序列
$in = Bio::SeqIO->new(
	-file   => "$ARGV[1]",
	-format => 'Fasta'
);

#输出序列：
$out = Bio::SeqIO->new(
	-file   => ">$ARGV[2]",
	-format => 'Fasta'
);

#读取需要提取基因ID
my %keep = ();
open IN, "$ARGV[0]" or die "$!; can't open file $ARGV[0]\n";

while (<IN>) {
	chomp;
	next if /^#/;
	my @a = split /\s+/;
	$keep{$a[0]}=1;
}
close(IN);

#输出想要的基因的序列
while ( my $seq = $in->next_seq() ) {
	my ( $id, $sequence, $desc ) = ( $seq->id, $seq->seq, $seq->desc );

	if ( exists $keep{$id} ) {
		$out->write_seq($seq);
	}
}
$in->close();
$out->close();
