#!/usr/bin/env perl
#北京组学生物科技有限公司
#学习perl语言：
#perl入门：https://zzw.xet.tech/s/1MxMKC
#perl高级：https://zzw.xet.tech/s/JKkbr

print "perl $0   <id_list>  <data_file> <out_file>\n" and die unless(@ARGV==3);

open IN,"$ARGV[0]" or die "$!; can't open file $ARGV[0]\n";

my%t;
my$head;
while(<IN>){
	chomp;
	my@tmp=split(/\s+/);
	
	
	$t{$tmp[0]}=0;
}

close(IN);

open IN,"$ARGV[1]" or die "$!; can't open file $ARGV[1]\n";

open OUT,">$ARGV[2]" or die "$!; can't open file $ARGV[2]\n";
while(<IN>){
	chomp;
	if (/^#/){
		print OUT "$_\n";
		next ;
	}

	my@tmp=split(/\s+/);

	if(exists $t{$tmp[0]}){
		print OUT "$_\n";
                $t{$tmp[0]}=1;
	}else{
		#print  "$tmp[0]\n";
	}
}
close(IN);

close(OUT);


#更新IDcheck 2023.1.6 huangls
for my $id (keys %t){
  if ($t{$id}==0){
    print "$id : not get data from $ARGV[1] file  please check \n";
  }
  
  
}
