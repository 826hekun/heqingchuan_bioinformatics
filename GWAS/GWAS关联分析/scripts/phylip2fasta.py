#!/usr/bin/env python3
'''
This script was used to  convert from phylip to fasta.

Author:
       haungls
Version:
        1.0;2020-09-16
'''

import sys, os, argparse, os.path
parser = argparse.ArgumentParser(description='This script was used to convert from phylip to fasta')
parser.add_argument('-i', '--phylip', dest='phylip', required=True, help='input  phylip format fasta file')
parser.add_argument('-d','--outDir',dest='outDir',required=False,default=os.getcwd(),help='specify the output file dir,default is current dir')
parser.add_argument('-o','--outFile',dest='outFile',required=False,default='demo.fa',help='specify the output file name,defaul demo.fa')
#parser.add_argument('-W','--width',required=False,default=15,type=int,help='specify the width column index,default is 15')
#parser.add_argument('-H','--height',required=False,default=10,type=int,help='specify the height column index,default is 10')

args = parser.parse_args()



from Bio import SeqIO

records = SeqIO.parse(args.phylip, "phylip")
count = SeqIO.write(records, args.outDir+"/"+args.outFile, "fasta")
print("Converted %i records" % count)
