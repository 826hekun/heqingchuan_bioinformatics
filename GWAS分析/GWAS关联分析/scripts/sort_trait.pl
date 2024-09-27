#!/usr/bin/env perl
use strict;
use warnings;
die "perl $0 <tfam> <pheno.table> >pheno.sorted.table\n" unless @ARGV==2;
my $tfam = $ARGV[0];
my $trait = $ARGV[1];

my %trait;

open IN, $trait || die $!;
while(<IN>){
	chomp;
	next if (/^\s*$/);
	my @t = split /\s+/;
	next if ( scalar(@t) != 2);
	if($t[1] =~ /^[0-9.]+$/){
		$trait{$t[0]} = $t[1];
#		print "phe:$_\n";
	}else{
		#warn "skiped ->> trait value not number: $_\n";
	}
}
close IN;

open FM, $tfam || die $!;

while(<FM>){
	chomp;
	next if (/^\s*$/);
	my @t = split /\s+/;
	if(exists $trait{$t[0]}){
		print "$_ $trait{$t[0]}\n";
	}
	else{
		print "$_ NA\n";
	}
}
close FM; 

