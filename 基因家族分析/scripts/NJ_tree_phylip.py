#! /usr/bin/env python3
'''
This script is writen to construct NJ tree.
Parts of the codes are writen redundantly and still needed to be improved.

Author:
       huangls@biomics.com.cn
Version:
        1.0;2023-04-11
'''
import os,subprocess,re,sys
import argparse
import matplotlib.pyplot as plt
################################################################################
parser = argparse.ArgumentParser(description="This script is writen to construct NJ tree")
parser.add_argument('-i','--phylip',dest='phylip', required=True,help='Please input complete  aligned  phylip format file path')
parser.add_argument('-t','--seqtype',dest='seqtype', required=True,help='Please set dna or prot ')

parser.add_argument('-o','--outdir',help='Please input complete out_put directory path',default = os.getcwd())
parser.add_argument('-n','--name',help='Please specify the output name')
################################################################################
args = parser.parse_args()
#################################################################################
#Below scripts serve to parse the command line parameters
#################################################################################


if not os.path.exists(args.outdir): os.mkdir(args.outdir)
#name=os.path.abspath(args.outdir)+'/'+args.prefix


if os.path.exists(args.outdir+'/'+"infile"):
  print("infile exist in outdir:"+args.outdir+'/'+"infile")
  sys.exit() 
else:
  os.system("cp %s %s"%(args.phylip,args.outdir+'/'+"infile"))

os.system("cd %s" %args.outdir)
os.system("echo 'R\n100\nI\nY\n9\n'|/share/work/biosoft/Phylip/latest/exe/seqboot")
os.rename('/'+args.outdir.strip()+'/infile','/'+args.outdir.strip()+'/file1')
os.rename('/'+args.outdir.strip()+'/outfile','/'+args.outdir.strip()+'/infile')
#################################################################################
if args.seqtype =="dna":
  os.system("echo 'D\nM\nD\n100\nI\nY\n'|/share/work/biosoft/Phylip/latest/exe/dnadist")
  os.rename('/'+args.outdir.strip()+'/infile','/'+args.outdir.strip()+'/file2')
  os.rename('/'+args.outdir.strip()+'/outfile','/'+args.outdir.strip()+'/infile')
elif(args.seqtype=="prot"):
  os.system("echo 'P\nM\nD\n100\nI\nY\n'|/share/work/biosoft/Phylip/latest/exe/protdist")
  os.rename('/'+args.outdir.strip()+'/infile','/'+args.outdir.strip()+'/file2')
  os.rename('/'+args.outdir.strip()+'/outfile','/'+args.outdir.strip()+'/infile')
else:
  print("unknown seqtype : %s , must be dna or  prot\n"%args.seqtype)
  sys.exit() 

#################################################################################

os.system("echo 'M\n100\n9\nY\n'|/share/work/biosoft/Phylip/latest/exe/neighbor")

os.rename('/'+args.outdir.strip()+'/outtree','/'+args.outdir.strip()+'/intree')
os.rename('/'+args.outdir.strip()+'/outfile','/'+args.outdir.strip()+'/multitree.outfile')
#################################################################################

os.system("echo 'Y\n'|/share/work/biosoft/Phylip/latest/exe/consense")

os.rename('/'+args.outdir.strip()+'/outtree','/'+args.outdir.strip()+'/consense.tre')
os.rename('/'+args.outdir.strip()+'/outfile','/'+args.outdir.strip()+'/consensetree.outfile')
#################################################################################
if args.seqtype =="dna":
  os.rename('/'+args.outdir.strip()+'/infile','/'+args.outdir.strip()+'/file3')
  os.rename('/'+args.outdir.strip()+'/file1','/'+args.outdir.strip()+'/infile')

  os.system("echo 'D\nI\nY\n'|/share/work/biosoft/Phylip/latest/exe/dnadist")

  os.rename('/'+args.outdir.strip()+'/infile','/'+args.outdir.strip()+'/file4')
  os.rename('/'+args.outdir.strip()+'/outfile','/'+args.outdir.strip()+'/infile')
elif(args.seqtype=="prot"):
  os.rename('/'+args.outdir.strip()+'/infile','/'+args.outdir.strip()+'/file3')
  os.rename('/'+args.outdir.strip()+'/file1','/'+args.outdir.strip()+'/infile')
  os.system("echo 'P\nI\nY\n'|/share/work/biosoft/Phylip/latest/exe/protdist")

  os.rename('/'+args.outdir.strip()+'/infile','/'+args.outdir.strip()+'/file4')
  os.rename('/'+args.outdir.strip()+'/outfile','/'+args.outdir.strip()+'/infile')
else:
  print("unknown seqtype : %s , must be dna or  prot\n"%args.seqtype)
  sys.exit() 

#################################################################################

os.system("echo 'Y\n'|/share/work/biosoft/Phylip/latest/exe/neighbor")

#################################################################################
BS=''
DS=''
for i in open('/'+args.outdir.strip()+'/consense.tre'):
     BS+=i.strip()
for i in open('/'+args.outdir.strip()+'/outtree'):
     DS+=i.strip()
BSD={}
DSD={}
#print BS
#print DS
RBS=BS.replace('(','(.')
RBS=RBS.replace(')','.)')
RDS=DS.replace('(','(.')
RDS=RDS.replace(')','.)')
#print RBS
#print RDS
GBS=re.search(RBS,BS)
#print GBS
GDS=re.search(RDS,DS)
#print GDS
for i in range(1,BS.count('(')+1):
    BSD.setdefault(GBS.group(i),set(re.findall('\(([a-zA-Z0-9]+?):',GBS.group(i))+re.findall(',([a-zA-Z0-9]+?):',GBS.group(i))))
for i in range(1,DS.count('(')+1):
    DSD.setdefault(GDS.group(i),set(re.findall('\(([a-zA-Z0-9]+?):',GDS.group(i))+re.findall(',([a-zA-Z0-9]+?):',GDS.group(i))))
#print BSD
#print DSD
kick=lambda x: x.count('(')
ass=DSD.keys()
ass=sorted(ass,key=kick)
#print ass
ass=reversed(ass)
#print DSD.keys()
#print ass
for i in ass:
    for ii in BSD:
        if DSD[i] == BSD[ii]:
#            print 'yes'

            boot=re.findall(re.escape(ii)+':([0-9\.]*)',BS)
#            print ii
#            print boot
            if len(boot)==1 and boot:
#                print 'ok'
                DS=DS.replace(i,i+boot[0])
            elif i==ass[0]:
                DS=DS.replace(i,i+'100')
            break
file=open('/'+args.outdir.strip()+'/'+args.name+'.tre','w')
file.write(DS)
file.close()
 

