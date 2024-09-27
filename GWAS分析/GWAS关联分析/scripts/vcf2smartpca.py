#!/usr/bin/env python3
import gzip
import sys, os, argparse, os.path
import re
#!/biosoft/miniconda/bin/python
# coding:utf-8
'''
北京组学生物科技有限公司
author huangls
date 2020.04.03
version 1.0
学习python课程：https://zzw.xet.tech/s/2bdd89
'''


parser = argparse.ArgumentParser(description='This script was used to convert vcf file to smartpca formats')
parser.add_argument("-i",'--vcf', dest='vcf', required=True ,help='set vcf file , gzip supported ,required')
parser.add_argument("-g",'--group', dest='group', required=True ,help='set sample group file,required')
parser.add_argument("-q",'--quality', dest='quality', required=False,type=int ,default=20,help='set genotype quality to keep SNP,[default 20]')
#parser.add_argument("-O",'--omit', dest='omit', required=False,type=int ,default=20,help='set genotype quality to keep SNP,[default 20]')

parser.add_argument("-o",'--outdir',dest='outdir',required=False,default=os.getcwd(),help='specify the output file dir,default is current dir')
parser.add_argument("-p",'--prefix', dest='prefix', required=False, default='all_sample_qc',help='file name prefix,default passing_snps_q20')
args = parser.parse_args()
if not os.path.exists(args.outdir): os.mkdir(args.outdir)
name=os.path.abspath(args.outdir)+'/'+args.prefix


#########################################################################################################################################
#                                                                                                                                       #
#   Written by Jason Ladner                                                                                                             #
#   questions?: jtladner@gmail.com                                                                                                      #
#   Usage:                                                                                                                              #
#   python vcf2smartpca.py   input_vcf_file   prefix   genotype_quality_cutoff  popfile                                                 #
#                                                                                                                                       #                                                                                                                                       #
#   input_vcf_file = the .vcf file that you want to extract genotype information from                                                   #
#   prefix = an arbitrary string of characters that will be added to the output filenames                                              #
#   genotype_quality_cutoff = a floating point number that specifies the minimum quality allowed by genotype information to be used     #
#   popfile = a plain text file with two columns per line, separated by some amount of whitespace                                       #
#             The 1st column is the individual ID (exactly as in the vcf file), the 2nd is the population that the individual belongs to#
#                                                                                                                                       #
# smartpca PCA analysis and Procrustes analysis
# generate smartpca inputs
# vcf2smartpca.py D_magna_RADseq.vcf passing_snps_q20 20 popfile 
# # Conversion script can be obtained from:
# # De Wit P, Pespeni MH, Ladner JT, et al. (2012) The simple fool's guide to population genomics via RNA-Seq: an introduction to high-throughput sequencing data analysis. Molecular Ecology Resources 12, 1058-1067
# 
# # smartpca running commands
# smartpca -p par.passing_snps_q20 > passing_snps_q20_logfile.txt
# 
# # re-run the above steps in order to reduce dataset by MN and IL individuals using
# vcf2smartpca.py D_magna_RADseq.vcf passing_snps_q20 20 popfile INDIVIDUALS_TO_OMIT MN IL
# 
# # and re-run smartpca to generate appropriate PCA eigenvalue and eigenvector files to go into Procrustes analysis




#########################################################################################################################################

'''---------------------------> VCF Format <-----------------------------
    cols[0] = contig name
    cols[1] = SNP position
    cols[2] = ID???? - empty in my vcf files
    cols[3] = Reference base
    cols[4] = Alternative base
    cols[5] = SNP quality (Phred score)
    cols[6] = filter flags
    cols[7] = SNP Info
    cols[8] = genotype format (e.g., GT:AD:DP:GQ:PL)
    cols[9:] = individual genotypes'''

#-------------------------------------------------------------#
#  This function takes a file and creates a dictionary with:  #
#    -one entry for each line                                 #
#    -keys will be first column of each line                  #
#    -values will be the complete line                        #
#-------------------------------------------------------------#
def file_dict(file):
    fin = open(file, 'r')       # open input file
    d={}
    
    for line in fin:
        line=line.rstrip()
        cols=line.split()
        d[cols[0]]=line
    fin.close()
    return d

#---------------------------------------------------------------------------------#
#  This function creates two lists, which will be input to smartpca_input()    :  #
#    -'headers' should be a list, headers[4:] should be individual names          #
#    -'snp_genos' should be a list with one list per indiviual                    #
#        -Each indiv list will contain genotpe data for each SNP                  #
#    -'pops' is a dictionary with indiv names as the keys and their population    #
#      of origin as the values                                                    #
#    - 'contigs' is a list of all the contig names for the snps being output      #
#---------------------------------------------------------------------------------#
def genos4smart(headers, snp_gens, pops, contigs, pos):
    loci_names=[]
    for index, contig in enumerate(contigs):
        loci_names.append(contig + '_' + pos[index])
    lines=[]
    indivcount=0
    for indiv in headers[4:]:
        thisline=[indiv, pops[indiv]]
        for each in snp_gens[indivcount]:
            for hap in each:
                thisline.append(hap)
        indivcount+=1
        lines.append(thisline)
    return loci_names, lines

#---------------------------------------------------------------------------------#
#  This function creates all the necc. input files smartpca:                      #
#     -'loci_names' should be a list that contains the names of all the SNPs      #
#     -'lines' is a list that contains an entry for each indivdual                #
#      -Each entry is also a list, with an entry for the indiv's ID & population  #
#      as well as two entries for each locus saying the # of each allele          #
#      'omit' is an optional parameter                                            #
#---------------------------------------------------------------------------------#
def smartpca_input(loci_names, lines, prefix, omit=[]): 
#Opens an output text files with prefix as specified by user
    GENO = open(prefix + '_Geno', 'w')
    INDIV = open(prefix + '_Indiv', 'w')
    SNP = open(prefix + '_SNP', 'w')
    print ("omitting: " + str(omit))
    pops=[]
    names=[]
    print ('# of loci = %d' % (len(loci_names)))
    
    loci=[]
    for locus in loci_names:
        loci.append([])
    for indiv in lines:
        
        if indiv[0] in omit: continue
        
        INDIV.write('%s\tU\t%s\n' % (indiv[0], indiv[1]))               #Writes an entry in the INDIV file for each individual: 'ID\tU\tPopulation'
        
        names.append(indiv[0])              #Adds each individual's name to the list names
        pops.append(indiv[1])               #Adds each individuals population to the list pops
        
        i=0
        k=0
        while i<len(indiv)-2:                              #the number of columns of actual data
            loci[k].append([indiv[i+2], indiv[i+3]])
            i+=2
            k+=1

#!!!!!!!!!!Can not include data for individuals with only info from one allele, because each snp is only coded with # of ref alleles!!!!!!!!!!!!
    for index, locus in enumerate(loci):
        indcount=0
        SNP.write(loci_names[index] + '\t1\t0.0\t1\n')                             #All SNPs specified to be in the same place on same chrom, I don't think it matters though because I do not utilize this data
        for ind in locus:
            if ind==['.', '.']:                                                     #Identifies individuals without genotype data and excludes them from the GENO file for this locus
                GENO.write('%s\t%s\t%d\n' % (loci_names[index], names[indcount], 9))   #missing data huangls 2020.3.3
                indcount+=1
            else:
                GENO.write('%s\t%s\t%d\n' % (loci_names[index], names[indcount], ind.count('0')))
                indcount+=1

    PAR = open('par.' + prefix, 'w')
    PAR.write('genotypename:    %s_Geno\nsnpname:         %s_SNP\nindivname:       %s_Indiv\nevecoutname:     %s_loci.evec\nevaloutname:     %s_loci.eval\nnumoutevec:      10\nnumoutlieriter:  0\nusenorm:         YES\nsnpweightoutname: %s_snpweights\n' % (prefix, prefix, prefix, prefix, prefix, prefix))
    PAR.close()
    GENO.close()
    SNP.close()
    INDIV.close()

#---------------------------------------------------------------------------------------------------------------
#-------> This function will create a dictionary relating your individuals to their population of origin
#-------> Popfile should be a plain text file with two columns per line, separated by some amount of whitespace
#-------> The first column is the individual ID (exactly as it is in the .vcf file
#-------> The second column is the population that the individual belongs to
#---------------------------------------------------------------------------------------------------------------
def make_pop_dict(popfile):
    pop_dict = {}
    fin = open(popfile, 'r')
    for line in fin:
        cols = line.strip().split()
        pop_dict[cols[0]] = cols[1]
    return pop_dict

#-------------------------------------------------------------#
#  This function takes the genotype description for a SNP     #
#   (e.g., GT:AD:DP:GQ:PL)                                    #
#  and returns the position of the genotype quality (GQ)      #
#  with the first position being 0                            #
#-------------------------------------------------------------#
def get_GQ_pos(genotype_description): 
    if "GQ" in   genotype_description:         
        return genotype_description.split(':').index('GQ')
    else:
        return None

#-----------------------------------------------------------------------------------------------#
#  This function takes a .vcf file and extract genotype info                                    #
#  As inputs it requires a .vcf file and an integer to be used as a genotype quality cutoff     #
#  It returns two variables:                                                                    #
#  'headers' is a list of the headers of interest (i.e., those to be output)                    #
#  'all_gens'is a list of lists with one sublist for each indiv. containing geno info           #
#  'contigs' is a list of all the contig names for the snps being output                        #
#-----------------------------------------------------------------------------------------------#
def extract_geno_info(input_vcf, genotype_quality_cutoff):
    #Opens an infile specified by the user. Should be a .vcf file 
    
    if args.vcf.endswith("gz"):
        IN = gzip.open(input_vcf, mode='rt')
    else:
        IN = open(input_vcf, 'rt')
    quality_cutoff=float(genotype_quality_cutoff)                    #for individual genotypes (not the SNP as a whole)
    
    linecount=0
    all_gens=[]
    locus=[]
    headers=[]
    pos=[]
    contigs=[]
    #Reads through the vcf file line by line
    for line in IN:
        line=line.rstrip()
        cols=line.split('\t')
        gens=[]                           #Will hold individual level genotype data, resets list form previous SNP
        gencount=0                        #Will be used below to keep track of the individuals when appending new SNPs
        
        if cols[0][0] == '#' and cols[0][1] != '#':       #To pull out headers
    
            headers=cols[0:2] + cols[3:5] + cols[9:]      #Creates list of headers of interest
            print (headers[4:])                             #Prints the IDs of the individuals with genotypes in the vcf
    
        if cols[0][0] != '#':                             #Specifies only SNP lines
            linecount+=1                                  #Keeps track of the numbers of SNPs
            
            gens=cols[9:]                  #Creates a list of all the genotype info
            desc=cols[8]
            GQ_pos = get_GQ_pos(desc)               #Pulls out the position of the genotype quality
            
            if linecount==1:                      #will only execute these commands once, when reading the first line of the file
                for gen in gens:
                    all_gens.append([])           #Creates seperate list entry for each individual, will hold genotypes at all SNPs
                    locus.append([])              #Creates seperate list entry for each individual, will hold genoypes at a specific locus 
            
            for gen in gens:                    #Steps through the genotype info for each individual
                parts=gen.split(':')            #Breaks apart the different sections of the genotype info
                genotype=parts[0]               #Most probable genotype for the individual, with the two alleles seperated by '/'
                if genotype==".":
                    genotype="./."
                if "/" in genotype:
                    genotypes=genotype.split('/')
                else:
                    genotypes=genotype.split('|')
                if  genotypes != ['.','.'] and genotypes != ['0','.'] and genotypes != ['.','1']:       #Had to be added because in new version of GATK, .vcf files are written so that if there is no genotype info for an individual, their genotype is './.' as opposed to './.:0:0:0:0', or whatever
                    if GQ_pos==None:
                        gen_quality=0
                    else:
                        try:
                        #print("%s\n"%GQ_pos)
                            gen_quality=float(parts[GQ_pos])
                        except ValueError:
                            print("%s\n"%GQ_pos)
                            print(gen)
                    if gen_quality<quality_cutoff:
                        genotypes=['.','.']
                    
                locus[gencount]=genotypes           #Replace genotype data from last locus with genotype form this locus
                gencount+=1                         #Keeping track of the individual whose data is being examined
        
            #if ['.', '.'] not in locus and locus != []:
            if  locus != []:           # missingdata accept
                contigs.append(cols[0])
                pos.append(cols[1])
                i_count=0
                for each in locus:
                    all_gens[i_count].append(each)                  #Appends genotype data to the all_gens list
                    i_count+=1                                      #keeps track of the number of loci with no missing data

    print ("Number of SNPs in reference contigs = %d" % (linecount))
    return headers, all_gens, contigs, pos

#-------> The Actual Script Execution Starts Here <-------

headers, all_gens, contigs, pos = extract_geno_info(args.vcf, args.quality)
list_headers, list_lines = genos4smart(headers, all_gens, make_pop_dict(args.group), contigs, pos)
smartpca_input(list_headers, list_lines, args.prefix)
#smartpca_input(list_headers, list_lines, args.prefix, sys.argv[5:])
