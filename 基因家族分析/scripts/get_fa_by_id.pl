#!/usr/bin/env perl
#北京组学生物科技有限公司
#学习perl语言：
#perl入门：https://zzw.xet.tech/s/1MxMKC
#perl高级：https://zzw.xet.tech/s/JKkbr


die "perl $0 <idlist> <fa> <OUT>" unless ( @ARGV == 3 );
use Math::BigFloat;
use Bio::SeqIO;
use Bio::Seq;

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
	$keep{$a[0]}=0;
}
close(IN);

#输出想要的基因的序列
while ( my $seq = $in->next_seq() ) {
	my ( $id, $sequence, $desc ) = ( $seq->id, $seq->seq, $seq->desc );

	if ( exists $keep{$id} ) {
	  $keep{$id}=1;
		$out->write_seq($seq);
	}
}
$in->close();
$out->close();

#更新IDcheck 2023.1.6 huangls
for my $id (keys %keep){
  if ($keep{$id}==0){
    print "$id : not get sequence from fasta file  please check \n";
  }
}
