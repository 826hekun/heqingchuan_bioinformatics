#!/usr/bin/env perl
#北京组学生物科技有限公司
#学习perl语言：
#perl入门：https://zzw.xet.tech/s/1MxMKC
#perl高级：https://zzw.xet.tech/s/JKkbr

use Getopt::Long;
my %opts;
use Data::Dumper;
GetOptions (\%opts,"in1=s","in2=s","out=s","h"); 
if (! defined($opts{in1}) ||! defined($opts{in2})||! defined($opts{out}) || defined($opts{h})){
	&USAGE;
}
open (IN1,"$opts{in1}") || die "open $opts{in1} failed\n";
open (IN2,"$opts{in2}") || die "open $opts{in2} failed\n";
open (OUT,">$opts{out}") || die "open $opts{out} failed\n";
my%gffs;
while (<IN1>) {
	next if (/^#/);
	chomp;
	  my@b=split,$_;
	  $keys= $b[0];

      $gffs{$keys} = 0;
 
}

while (<IN2>) {
	 chomp;
          my @a=split /\t/,$_;
		 #if ($a[2] =~/^mRNA$/i or $a[2] =~/^transcript$/i ) {	 
		 if ($a[2] =~/^gene$/i ) {	 
			($m)=($a[8]=~ m/ID=([^;]*)/);             #注意这里匹配的ID信息
		
			if ( exists  $gffs{$m}  ) {
	  	 		print OUT "$m\t$a[0]\t$a[3]\t$a[4]\t$a[6]\n";
				$gffs{$m}=1;
			}
		 }
		 
}
close OUT;
close IN1;
close IN2;

for my $g(keys %gffs){
	if($gffs{$g}==0){
		print "$g failed get location\n";
	}

}


sub USAGE {
       print "usage: perl $0 -in1  gene_id.txt -in2  genome.gff3  -out gene_location.txt ";
	exit;
}
