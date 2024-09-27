# coding:utf-8
'''
北京组学生物科技有限公司
author huangls
date 2020.04.20
version 1.0
学习python课程：https://zzw.xet.tech/s/2bdd89
'''
from Bio.PDB import *
import sys,os
cwd = os.getcwd()
if(len(sys.argv)==1):
    print("Please input pdb id.")
    print("Use example:")
    print("\tpython %s %s %s"%(__file__,"1FAT","..."))
    sys.exit(0)

pdb_codes=[]

pdb_codes=sys.argv[1:]
print(pdb_codes)
pdbl = PDBList()
#pdbl.retrieve_pdb_file("1FAT")
for i in pdb_codes:
    #print("get %s "%i)
    pdbl.retrieve_pdb_file( i,  pdir=cwd,file_format="pdb")
