#!/usr/bin/env perl
# convert vcf file to hapmap format
# huangls 2020.12.21
#北京组学生物科技有限公司
#学习perl语言：
#perl入门：https://zzw.xet.tech/s/1MxMKC
#perl高级：https://zzw.xet.tech/s/JKkbr

die "perl $0 <vcf> <out> " unless ( @ARGV == 2 );
use strict;
use warnings;
#use PerlIO::gzip; 
my $vcffile = $ARGV[0];
my $outfile = $ARGV[1];

if($vcffile=~/gz$/){
	open INFILE, "gzip -dc $vcffile|" or die "ERROR: could not open $vcffile!\n";
	
}else{
	open INFILE, $vcffile or die "ERROR: could not open $vcffile!\n";
}

if($outfile=~/gz$/){
	open OUTFILE, "| gzip >$outfile" or die "ERROR: could not open $outfile!\n";
	
}else{
	open OUTFILE,">$outfile" or die "ERROR: could not open $outfile!\n";
}

my $line;
my @splitline;
my @headerline;
my @call;
my $count = 0; 
my $i;

while ($line = <INFILE>)
{
    chomp $line;
    if ($line =~ /#CHROM/)
    {
        @headerline = split("\t", $line);
        print OUTFILE "rs\talleles\tchrom\tpos\tstrand\tassembly\tcenter\tprotLSID\tassayLSID\tpanel\tQCcode";
        for ($i = 9; $i < scalar(@headerline); $i = $i + 1)
        {
            print OUTFILE "\t", $headerline[$i]; 
        }
        print OUTFILE "\n";
    }   

    if ($line !~ /^#/)
    {
        
        
        @splitline = split("\t", $line);
        next if(length($splitline[3])!=1 or length($splitline[4])!=1);
        #print OUTFILE "S".$splitline[0] . "_" . $splitline[1], "\t", $splitline[3], "\/", $splitline[4];
        print OUTFILE $splitline[0] . ":" . $splitline[1], "\t", $splitline[3], "\/", $splitline[4];
        
        print OUTFILE"\t", $splitline[0], "\t", $splitline[1], "\t+\tNA\tNA\tNA\tNA\tNA\tNA";
        #print OUTFILE "\n";

        for ($i = 9; $i < scalar(@headerline); $i = $i + 1)
        {
            @call = split(":", $splitline[$i]);
            
            if ($call[0] =~ /0[\/\|]0/)
            {
                print OUTFILE "\t", $splitline[3], $splitline[3];
            }
            elsif ($call[0] =~ /0[\/\|]1/)
            {
                print OUTFILE "\t", $splitline[3], $splitline[4];
            }
            elsif ($call[0] =~ /1[\/\|]1/)
            {
                print OUTFILE "\t", $splitline[4], $splitline[4];
            }
            elsif ($call[0] =~ /\.[\/\|]\./)
            {
                print OUTFILE "\tNN";
            }
            elsif ($call[0] =~ /^\.$/){
            	print OUTFILE "\tNN";
            }else{
                die "ERROR: invalid genotype call at marker position ", $splitline[0], " ", $splitline[1]," ", $splitline[$i] ,"\n";
            }
        }
        print OUTFILE "\n";
    }
}


close(OUTFILE);
