#!/usr/bin/env python3
# coding:utf-8
'''
北京组学生物科技有限公司
author huangls
date 2020.04.03
version 1.0
学习python课程：https://zzw.xet.tech/s/2bdd89
'''

import sys, os, argparse, os.path
#from pylab import *
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import math
import re
#sys.path.append('/share/bioCloud/huangls/02.package/python/svgwrite-1.1.6')
#import svgwrite
mpl.rcParams['interactive']=False
from matplotlib.ticker import NullFormatter
from mpl_toolkits.mplot3d import Axes3D


parser = argparse.ArgumentParser(description='This script was used to plot PCA')
#parser.add_argument('-g', '--geneID', dest='geneID', required=False, help='input  geneID list file')
parser.add_argument('-i', '--pca', dest='pca', required=True, help='input  pca  eigenvec file')
parser.add_argument('-v', '--pcav', dest='pcav', required=True, help='input  pca eigenval file')

parser.add_argument('-o','--outDir',dest='outDir',required=False,default=os.getcwd(),help='specify the output file dir,default is current dir')
parser.add_argument('-p','--outFilePrefix',dest='outFilePrefix',required=False,default='pca',help='specify the output file prefix,defaulfstAndpi')
parser.add_argument('-W','--width',required=False,default=15,type=int,help='specify the width column index,default is 15')
parser.add_argument('-H','--height',required=False,default=10,type=int,help='specify the height column index,default is 10')

args = parser.parse_args()






group=[]
sample=[]

sample_color={}
group_color={}

f=open(args.pca,'r')

pc1=0.
pc2=0.
pc3=0.
A = []
B = []
C = []
for i in f:
    i=i.strip()
    if i[0]=='#':
        continue
    #print i
    tmp=re.split("\s+",i)
    s,p1,p2,p3,g=tmp[0],tmp[1],tmp[2],tmp[3],tmp[-1]
    sample.append(s)
    group.append(g)
    sample_color[s]=1
    group_color[g]=1
    
    
    A.append(float(p1))
    B.append(float(p2))
    C.append(float(p3))
    

f.close()

f=open(args.pcav,'r')
enV=[]
for i in f:
    i=i.strip()
    if i[0]=='#':
        continue
    enV.append(float(i))
    
enV=np.array(enV)

pc1_str="PC1(%.2f%%)"%((enV[0]/enV.sum())*100)
pc2_str="PC2(%.2f%%)"%((enV[1]/enV.sum())*100)
pc3_str="PC3(%.2f%%)"%((enV[2]/enV.sum())*100)


# pc1_str="PC1"
# pc2_str="PC2"
# pc3_str="PC3"


# the random data

nullfmt = NullFormatter()         # no labels

# definitions for the axes
left, width = 0.1, 0.4
bottom, height = 0.1, 0.4

gap=0.05
rect_up_left = [left, bottom+height, width, height]
rect_up_right = [left+width+gap, bottom+height, width, height]
rect_down_left = [left, bottom-gap, width, height]
rect_down_right = [left+width+0.05, bottom-0.05, width, height]
# start with a rectangular Figure
fig=plt.figure(1, figsize=(12, 12))








ax_up_left = fig.add_axes(rect_up_left)
ax_up_right = fig.add_axes(rect_up_right,sharey=ax_up_left)
ax_down_left = fig.add_axes(rect_down_left,sharex=ax_up_left)
ax_down_right = fig.add_axes(rect_down_right,projection='3d')
#ax_down_right = fig.add_subplot(2, 2, 1, projection='3d')


ax_down_right.set_xlabel(pc1_str)
ax_down_right.set_ylabel(pc2_str)
ax_down_right.set_zlabel(pc3_str)









# no labels
#ax_up_left.xaxis.set_major_formatter(nullfmt)
#ax_up_right.yaxis.set_major_formatter(nullfmt)



ax_down_left.set_xlabel(pc1_str)
ax_up_left.set_ylabel(pc2_str)
ax_up_left.set_xlabel(pc1_str)
ax_down_left.set_xlabel(pc1_str)
ax_down_left.set_ylabel(pc3_str)
ax_up_right.set_xlabel(pc3_str)

markers = ['s','o','^','P','D','*','p','X','H','>','+','s','o','^','P','D','*','p','X','H','>','+','s','o','^','P','D','*','p','X','H','>','+','s','o','^','P','D','*','p','X','H','>','+','s','o','^','P','D','*','p','X','H','>','+','s','o','^','P','D','*','p','X','H','>','+','s','o','^','P','D','*','p','X','H','>','+','s','o','^','P','D','*','p','X','H','>','+']
cm = plt.get_cmap("Dark2")
n=len(set(group))
mycol=[cm(float(i)/n) for i in range(n)]


if len(set(sample))==len(set(group)):
   
    markers = Line2D.filled_markers
    markers=np.tile(markers,33)
    #cm = plt.get_cmap("Set1")
    
    map={}
    for s,co,mar in zip(sample,mycol,markers):
        map[s]=[co,mar]
    legendPatch=[]
    legendLebel=[]
    for s, x,y,z in zip(sample,A,B,C):
        # the scatter plot:
        c=map[s][0]
        m=map[s][1]
        scatt=ax_up_left.scatter(x, y,c=c, marker=m,label=s)
        legendPatch.append(scatt)
        legendLebel.append(s)
        
        ax_up_left.text(x, y, s,
            horizontalalignment='left',
            verticalalignment='top')
            
        ax_up_right.scatter(z, y,c=c, marker=m,label=s)
        ax_down_left.scatter(x, z,c=c, marker=m,label=s)
        ax_down_right.scatter(x, y, z,c=c, marker=m,label=s)
    
    fig.legend(legendPatch,legendLebel,bbox_to_anchor=(left,bottom+height*2+0.01,width*2,0.049),
                  loc='upper center', shadow=False, fontsize='x-small',framealpha=0.0,ncol=12,bbox_transform=fig.transFigure)
        
else:
    map={}
    for g,co,mar in zip(set(group),mycol,markers):
        map[g]=[co,mar]
        
    sam_to_group={}
    for s,g in zip(sample,group):
        sam_to_group[s]=g
        
    legendPatch=[]
    legendLebel=[]
    has_g={}
    
    for s, x,y,z in zip(sample,A,B,C):
        # the scatter plot:
        c=map[sam_to_group[s]][0]
        m=map[sam_to_group[s]][1]
        scatt=ax_up_left.scatter(x, y,c=c, marker=m,label=sam_to_group[s],alpha=1)
        
        if(sam_to_group[s] in has_g):
            pass
        else:
            legendPatch.append(scatt)
            legendLebel.append(sam_to_group[s])
            has_g[sam_to_group[s]]=1
        ax_up_right.scatter(z, y,c=c, marker=m,label=sam_to_group[s],alpha=1)
        ax_down_left.scatter(x, z,c=c, marker=m,label=sam_to_group[s],alpha=1)
        ax_down_right.scatter(x, y, z,c=c, marker=m,label=sam_to_group[s],alpha=1)
        
    fig.legend(legendPatch,legendLebel,bbox_to_anchor=(left,bottom+height*2+0.01,width*2,0.049),
                  loc='upper center', shadow=False, fontsize='x-small',framealpha=0.0,ncol=12,bbox_transform=fig.transFigure)




plt.savefig(args.outDir.rstrip('/') +'/'+args.outFilePrefix+'.pdf')
plt.savefig(args.outDir.rstrip('/') +'/'+args.outFilePrefix+'.png',dpi=300)
