#北京组学生物科技有限公司
#学习perl语言：
#perl入门：https://zzw.xet.tech/s/1MxMKC
#perl高级：https://zzw.xet.tech/s/JKkbr

use Getopt::Long;
my %opts;
use Data::Dumper;
GetOptions( \%opts, "in1=s", "in2=s", "out=s", "h" );
if (   !defined( $opts{in1} )
	|| !defined( $opts{in2} )
	|| !defined( $opts{out} )
	|| defined( $opts{h} ) )
{
	&USAGE;
}
open( IN1, "$opts{in1}" )  || die "open $opts{in1} failed\n";
open( IN2, "$opts{in2}" )  || die "open $opts{in2} failed\n";
open( OUT, ">$opts{out}" ) || die "open $opts{out} failed\n";
my %gffs;
while (<IN1>) {
	chomp;
	next if /^#/;
	my @b = split/\s+/, $_;
	$gffs{$b[0]} = 1;
}

#print Dumper(\%gffs);
my %mRNA;
while (<IN2>) {
	chomp;
	next if (/^#/);
	my @a = split /\t/, $_;
	next if $a[2]=~/exon/i;
	if ($a[2] =~/^mRNA$/i or $a[2] =~/^transcript$/i ) {
		($id2) =  ($a[8] =~ m/ID=([^;]*)/);
		($id1) =  ($a[8] =~ m/Parent=([^;]*)/);
		if (exists $gffs{$id1} or exists $gffs{$id2}){
			 $mRNA{$id2}=1;
		}

	}elsif ( $a[2] =~/^CDS$/i or $a[2] =~/utr/i ) {

		($id2) =  ($a[8] =~ m/Parent=([^;]*)/);
	}else{
		next;
	}

	if ( exists $mRNA{$id2} ) {
		print OUT "$_\n";
	}

}
close OUT;
close IN1;
close IN2;

sub USAGE {
	print "usage: perl $0 -in1  mRNA_id.txt -in2  genome.gff3  -out gene_location.txt ";
	exit;
}
