#安装软件 
conda install -c bioconda Htseq-count fastqc multiqc

conda create -n htseq-count
conda activate htseq-count
conda install -c bioconda htseq-count

#准备一些路径变量，方便后续调用
####################################################

workdir=/u3/hekun/rnaseq  #设置工作路径
refdir=$workdir/ref
datadir=$workdir/data
scriptdir=$workdir/scripts

REF_INDEX=$refdir/Taestivum.dna.toplevel   #hisat2 索引，没有.fa
GFF=$refdir/Taestivum.58.gff3
GTF=$refdir/Taestivum.58.gtf   #这个gtf文件是由gff3文件通过脚本index.sh生成的，不是下载的。
GENE_BED=$refdir/gene.bed
GENE_LENGTH=$refdir/gene_length.txt

######################################################################################
#下载转录组分析docker镜像： 镜像升级v2.0，兼容v1.0,v1.1

#docker pull omicsclass/rnaseq:v2.0
#启动docker容器并交互式进入
#docker desktop方法
#docker run --rm -it -m 8G --cpus 1  -v D:\rnaseq:/work omicsclass/rnaseq:v2.0   # 镜像升级v2.0，兼容v1.0,v1.1
#docker toolbox方法
#docker run --rm -it -m 8G --cpus 1  -v /d/rnaseq:/work omicsclass/rnaseq:v2.0

#linux服务器
#docker run --rm -it -m 32G --cpus 5  -v /share/nas1/huangls/test/rnaseq/:/work omicsclass/rnaseq:v2.0


####################################################
#准备一些路径变量，方便后续调用
####################################################

workdir=/u3/hekun/rnaseq  #设置工作路径
refdir=$workdir/ref
datadir=$workdir/data
scriptdir=$workdir/scripts

###################################################################
# 下载参考基因组数据,并对基因组构建HISAT index
###################################################################
cd $workdir/  ###回到工作目录
cd $refdir  #进入参考基因组ref目录

#hisat2 提供常见物种基因组索引文件：
#下载参考基因组数据到ref文件夹
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/gtf/triticum_aestivum/Triticum_aestivum.IWGSC.58.gtf.gz 
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/gff3/triticum_aestivum/Triticum_aestivum.IWGSC.58.gff3.gz
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/fasta/triticum_aestivum/cdna/Triticum_aestivum.IWGSC.cdna.all.fa
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz

#gff3文件修改，因为ensemble上下载的注释上会有gene：transcript:，后续作分析带着很麻烦，所以提前处理（先提前查看文件看看有没有）
gunzip *.gz
sed 's#gene:##' Taestivum.58.gff3|sed 's#transcript:##'|less -SN
sed 's#gene:##' Taestivum.58.gff3|sed 's#transcript:##' >Taestivum.58.gff3.1 #取个名字
mv Taestivum.58.gff3.1 Taestivum.58.gff3 #重命名，覆盖原文件

#wget -c ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
#wget -c ftp://ftp.ensembl.org/pub/release-99/gff3/homo_sapiens/Homo_sapiens.GRCh38.99.chromosome.22.gff3.gz
#gunzip *gz

#建立索引 常见物种索引下载： http://daehwankimlab.github.io/hisat2/download/
# hisat2 建立索引非常耗内存，如果建立索引失败 被 killed
#这里由由于服务器内存太小，所以其中的--ss和--exon去掉了
#首先提取gtf文件中的剪切位点和外显子位置： 
#extract_splice_sites.py gencode.vM4.annotation.gtf >gencode.vM4.annotation.for.hisat2.ss 
#extract_exons.py gencode.vM4.annotation.gtf >gencode.vM4.annotation.for.hisat2.exon 
#建立索引： 
#hisat2-build -p 30 --ss gencode.vM4.annotation.for.hisat2.ss --exon gencode.vM4.annotation.for.hisat2.exon GRCm38.p3.genome.fa mouseGencodeIndex 
##如果电脑内存<200G，那么可以不用--ss/--exon参数，但是在比对的时候需要加 --known-splicesite-infile参数（但效果依然没有直接加--ss/--exon参数好，如果接下来索引后的比对使用--known-splicesite-infile参数，则如下例子 hisat2 --known-splicesite-infile gencode.vM4.annotation.for.hisat2.ss --dta -t -p 24 -x mouseGencodeIndex -1 samp_1.fq.gz -2 samp_2.fq.gz -S accepted_hits.sam &> alignment_summary.txt）。
#Taestivum.58.gtf没有的话，脚本会根据gff文件自动生成
sh $scriptdir/index.sh Taestivum.dna.toplevel.fa Taestivum.58.gff3 Taestivum.58.gtf



#设置参考基因组相关文件变量，方便后续使用

REF_INDEX=$refdir/Taestivum.dna.toplevel   #hisat2 索引，没有.fa
GFF=$refdir/Taestivum.58.gff3
GTF=$refdir/Taestivum.58.gtf   #这个gtf文件是由gff3文件通过脚本index.sh生成的，不是下载的。
GENE_BED=$refdir/gene.bed
GENE_LENGTH=$refdir/gene_length.txt




#对原始数据进行fastqc质控
#fastqc质控报告查看https://www.omicsclass.com/article/1231


cd $workdir  #回到工作目录
mkdir 1.fastqc

fastqc $datadir/*.gz  -o $workdir/1.fastqc
echo "fastqc $datadir/*.gz  -o $workdir/1.fastqc" #用于查看命令的全称

#对fastqc结果整合与汇总：2023.10.09 更新
multiqc  $workdir/1.fastqc -o $workdir/1.fastqc/multiqc
echo "multiqc  $workdir/1.fastqc -o $workdir/1.fastqc/multiqc" #用于查看命令的全称



#数据质控：对原始序列进行去接头，删除低质量的reads等等

####

cd $workdir #回到工作目录

mkdir 2.data_qc

cd  2.data_qc

for i in /u3/hekun/rnaseq/data/*_1.clean.fq.gz;
do 
samp=`basename ${i} _1.clean.fq.gz`

echo "Processing sample ${samp}"

fastp --thread 4 --qualified_quality_phred 10 \
--unqualified_percent_limit 50 \
--n_base_limit 10 \
-i $datadir/${samp}_1.clean.fq.gz \
-I $datadir/${samp}_2.clean.fq.gz \
-o ${samp}_1.clean.fastp.fq.gz \
-O ${samp}_2.clean.fastp.fq.gz \
--detect_adapter_for_pe \
-h ${samp}.html \
-j ${samp}.json
done
echo "#利用fastp工具去除adapter

#--detect_adapter_for_pe  默认对双端数据则默认不使用自动检测adapter(SE可自动检测)，设置该参数，表示对双端数据也启用自动检测方法；接头序列未知  可设置软件自动识别常见接头

#--detect_adapter_for_pe    # 开启双端的接头检测(应该不同于overlap方法)

#--adapter_fasta    # 指定包含接头序列的fasta文件

                   # 接头序列至少6bp长度，否则将被跳过

                   # 可以指定任何想去除的序列，比如polyA

#对于每一个样本，执行以下操作  文件名如下格式 normal rep3_r2. fastq. gz ，tumor repl_rl. fastq. gz，tumor repl r2. fastq. gz



#使用for关键字开始一个循环，i是一个变量，用于存储每次循环中的当前文件名

#in后面跟着一个通配符表达式，表示要遍历的文件集合

#这里的通配符*表示任意长度的任意字符，_1.clean.fq.gz表示以这个字符串结尾

#所以，这个表达式匹配了data目录下的所有以_1.clean.fq.gz结尾的文件

#例如，ANTF1_FRAS202191140-1r_1.clean.fq.gz，ANTF2_FRAS202191140-1r_1.clean.fq.gz等

#分号表示一条命令的结束

for i in /u3/hekun/rnaseq/data/*_1.clean.fq.gz;

do # do表示循环体的开始，每次循环都会执行循环体中的命令

#使用反引号将basename命令包围，表示执行这个命令并将结果赋值给samp变量

#basename命令用于获取文件的基本名称，即去掉目录和后缀的部分

#basename命令接受两个参数，第一个是文件的完整路径，第二个是要去掉的后缀

#这里，第一个参数是i变量，表示当前循环的文件名，第二个参数是_1.clean.fq.gz，表示要去掉的后缀

#例如，如果i的值是/u3/hekun/rnaseq/data/ANTF1_FRAS202191140-1r_1.clean.fq.gz，那么basename命令的结果就是ANTF1_FRAS202191140-1r

samp=`basename ${i} _1.clean.fq.gz`

#使用echo命令打印出一条信息，表示正在处理哪个样本

#使用双引号将字符串包围，表示字符串中的变量会被替换为实际的值

#使用美元符号和花括号将变量名包围，表示引用变量的值

#例如，如果samp的值是ANTF1_FRAS202191140-1r，那么echo命令会打印出Processing sample ANTF1_FRAS202191140-1r
fastp --thread 1 --qualified_quality_phred 10 \

--unqualified_percent_limit 50 \

--n_base_limit 10 \

-i $datadir/${samp}_1.clean.fq.gz \

-I $datadir/${samp}_2.clean.fq.gz \

-o ${samp}_1.clean.fastp.fq.gz \

-O ${samp}_2.clean.fastp.fq.gz \

#--adapter_fasta $workdir/data/illumina_multiplex.fa \

--detect_adapter_for_pe \

-h ${samp}.html \
-j ${samp}.json

done
"



#质控数据统计汇总：
python $scriptdir/qc_stat.py -d $workdir/2.data_qc/ -o $workdir/2.data_qc/ -p all_sample_qc

#再次fastqc进行质控检查

cd $workdir/2.data_qc/  #回到工作目录
mkdir fastqc
cd fastqc
fastqc $workdir/2.data_qc/*.gz  -o $workdir/2.data_qc/fastqc

echo "fastqc $workdir/2.data_qc/*.gz  -o $workdir/2.data_qc/fastqc" #用于查看命令的全称

#对fastqc结果整合与汇总：2023.10.09 更新
multiqc  $workdir/2.data_qc/fastqc -o $workdir/2.data_qc/fastqc/multiqc
echo "multiqc  $workdir/2.data_qc/fastqc -o $workdir/2.data_qc/fastqc/multiqcc" #用于查看命令的全称

###################################################################
# reads 与基因组进行比对map
##################################################################
cd $workdir  #回到工作目录
mkdir -p $workdir/3.map/hisat2
cd $workdir/3.map/hisat2




# 附加步骤，更改简化clean数据文件名，以后应该在开始分析的时候就简化好文件名，不要中间，将ANTF1_FRAS202191140-1r_1.clean.fastp.fq.gz文件名，改为AN_rep1.fastp.fq.gz
for file in *clean.fastp.fq.gz
do
  # 用sed命令将文件名中的两个下划线及其之间的内容替换为一个下划线,[^r]表示一个非r字符，^在方括号内表示取反的意思。
  newname=$(echo $file | sed 's/_[^r]*r//')

  # 用mv命令将原文件重命名为新文件名
  mv $file $newname
done



# 附加步骤，更改简化clean数据文件名，以后应该在开始分析的时候就简化好文件名，不要中间，将ANTF1_FRAS202191140-1r_1.clean.fastp.fq.gz文件名，改为AN_rep1.fastp.fq.gz
for file in *clean.fastp.fq.gz
do
   newname=$(echo $file | sed 's/clean.//')

  # 用mv命令将原文件重命名为新文件名
  mv $file $newname
done



# 附加步骤，更改简化clean数据文件名，以后应该在开始分析的时候就简化好文件名，不要中间，将ANTF1_FRAS202191140-1r_1.clean.fastp.fq.gz文件名，改为AN_rep1.fastp.fq.gz
for file in *fastp.fq.gz
do
   newname=$(echo $file | sed 's/TF/_rep/')

  # 用mv命令将原文件重命名为新文件名
  mv $file $newname
done
#AN_rep1_1.fastp.fq.gz






cd /u3/hekun/rnaseq/3.map/hisat2/

#hisat2软件，链特异性文库设置
#--rna-strandness RF or FR
#非链特异性文库，不设置此参数

for fn in /u3/hekun/rnaseq/2.data_qc/*_1.fastp.fq.gz;
do 
i=`basename ${fn} _1.fastp.fq.gz`
echo "Processing sample ${i}"
echo "hisat2 -p 20 --rg-id=${i} --rg SM:${i} --rg LB:${i} --rg PL:ILLUMINA \
-x $REF_INDEX --dta \
-1 $workdir/2.data_qc/${i}_1.clean.fastp.fq.gz \
-2 $workdir/2.data_qc/${i}_2.clean.fastp.fq.gz \
-S ${i}.sam 2>${i}.summary"
#--known-splicesite-infile splicesites.tsv只有在构建索引时不添加--ss和--eexon是才用（但若是服务器内存足够，建议使用--ss和--eexon）
hisat2 --known-splicesite-infile splicesites.tsv \
-p 20 --rg-id=${i} --rg SM:${i} --rg LB:${i} --rg PL:ILLUMINA \
-x $REF_INDEX --dta \
-1 $workdir/2.data_qc/${i}_1.fastp.fq.gz \
-2 $workdir/2.data_qc/${i}_2.fastp.fq.gz \
-S ${i}.sam 2>${i}.summary
done


echo 

"hisat2 -p 20 --rg-id=${i} --rg SM:${i} --rg LB:${i} --rg PL:ILLUMINA \
-x $REF_INDEX --dta \
-1 $workdir/2.data_qc/${i}_1.clean.fastp.fq.gz \
-2 $workdir/2.data_qc/${i}_2.clean.fastp.fq.gz \
-S ${i}.sam 2>${i}.summary"
echo "RUN CMD: # 这段代码是用bash shell编写的，它的作用是使用hisat2软件对一些RNA测序数据进行比对和分析
# hisat2是一个用于比对RNA测序数据到参考基因组的软件，它可以处理不同的转录本和剪接变异
# 参考文献：[Hisat2: a fast and sensitive alignment program for mapping next-generation sequencing reads to a population of genomes](https://phoenixnap.com/kb/bash-comment)

# for是一个循环语句，它的作用是对一个列表中的每个元素执行一些命令
# do是一个关键字，它表示循环体的开始
# hisat2是一个命令，它的作用是执行hisat2软件
# -p 8是一个选项，它的作用是指定使用8个线程进行比对
#在sam文件写入id,sample名，建库library id,不知道的话直接写样本id，PL测序平台
# --rg-id=${i}是一个选项，ID，这里使用样本的名称
# --rg SM:${i}是一个选项，样本名，这里使用样本的名称
# --rg LB:${i}是一个选项，文库名，这里使用样本的名称
# --rg PL:ILLUMINA是一个选项，测序平台，这里使用ILLUMINA
# -x $REF_INDEX是一个选项，它的作用是指定参考基因组的索引文件，这里使用一个变量$REF_INDEX，它的值应该在之前定义过
# --dta是一个选项，它的作用是指定输出的SAM文件可以用于StringTie软件进行转录本重建
# --rna-strandness RF是一个选项，它的作用是指定RNA测序的链特异性，这里使用RF，表示第一条链是反向互补的，第二条链是正向的,RF连特异性文库及其类型，**非特异性则删除该参数**
# -1 $workdir/2.data_qc/${i}_1.clean.fq.gz是一个选项，它的作用是指定第一条链的测序数据文件，这里使用一个变量$workdir，它的值应该在之前定义过，表示工作目录的路径，后面跟着样本的名称和文件的后缀
# -2 $workdir/2.data_qc/${i}_2.clean.fq.gz是一个选项，它的作用是指定第二条链的测序数据文件，这里使用一个变量$workdir，它的值应该在之前定义过，表示工作目录的路径，后面跟着样本的名称和文件的后缀
# -S ${i}.sam是一个选项，它的作用是指定输出的SAM文件的名称，这里使用样本的名称和文件的后缀
# 2>${i}.summary是一个重定向符号，它的作用是将标准错误输出重定向到一个文件，这里使用样本的名称和文件的后缀
hisat2 -p 1 --rg-id=${i} --rg SM:${i} --rg LB:${i} --rg PL:ILLUMINA \
-x $REF_INDEX --dta --rna-strandness RF \
-1 $workdir/2.data_qc/${i}_1.clean.fq.gz \
-2 $workdir/2.data_qc/${i}_2.clean.fq.gz \
-S ${i}.sam 2>${i}.summary
--rna-strandness RF
# done是一个关键字，它表示循环体的结束
done"


# sam 转换成bam 格式并排序
for fn in /u3/hekun/rnaseq/2.data_qc/*_1.fastp.fq.gz;
do 
i=`basename ${fn} _1.fastp.fq.gz`
echo "Processing sample ${i}"

echo "RUN CMD: samtools sort --threads 20 -m 10G（给10G 内存） -o ${i}.bam ${i}.sam"
samtools sort  --threads 20 -m 10G -o ${i}.bam ${i}.sam
done








# bam index建立索引
for fn in /u3/hekun/rnaseq/2.data_qc/*_1.fastp.fq.gz;
do 
i=`basename ${fn} _1.fastp.fq.gz`
echo "Processing sample ${i}"

echo "RUN CMD: samtools index ${i}.bam"
samtools index ${i}.bam
done


# 对bam建索引
#samtools index ANrep1.Hisat_aln.sorted.bam ANrep1.Hisat_aln.sorted.bam.bai
#如果报错：bam_sort_core] merging from 42 files and 3 in-memory blocks...
#[E::hts_idx_push] Region 536930333..536930483 cannot be stored in a bai index. Try using a csi index with min_shift = 14, n_lvls >= 6
#samtools index: failed to create index for "/u3/hekun/jinengshu-rnaseq/project/Taestivum-9-state-Trans/Mapping/Hisat2/ANrep3.Hisat_aln.sorted.bam": Numerical result out of range
#这个错误通常是因为BAM文件中的某个区域太大，无法用标准的BAI索引来存储。BAI索引有一个限制，它不能索引超过(2^{29} - 1)的位置，这大约是536870912。当你的BAM文件包含比这更大的位置时，就会出现这个错误。

#解决这个问题的方法是使用CSI索引，它支持更大的区域。你可以通过以下命令来创建CSI索引：
#这里的-c选项是指定创建CSI索引，-m 14是设置最小shift为14，这样就可以支持更大的区域了。使用CSI索引时，不会生成.bam.bai文件。相反，它会生成一个.bam.csi文件，这是一个更先进的索引格式，可以支持更大的基因组区域。这样，你就可以避免因为区域太大而无法存储在BAI索引中的问题
#这里的-m 14意味着索引的最小区间大小将是2^14。如果你需要支持更大的区域，你可以增加这个值。。默认情况下，这个值是2^14。这个参数决定了索引的粒度，也就是说，它决定了索引可以支持的最大区域的大小。
#samtools index -c -m 14 ANrep1.Hisat_aln.sorted.bam
#for i in ./*bam;
#do 
#echo "Processing sample ${i}"
#samtools index -c -m 14  ${i}
#done &





######################################################################
# 只学习到这一步，结果出错了，以后有时间继续学。对map结果进行QC分析
######################################################################
## 采用RSeQC 对比对结果文件进行质控分析
cd $workdir  #回到工作目录
mkdir -p $workdir/3.map/map_QC
cd $workdir/3.map/map_QC

#1.片段inner size，片段选择是否异常（报错了）

for fn in /u3/hekun/rnaseq/2.data_qc/*_1.fastp.fq.gz;
do 
i=`basename ${fn} _1.fastp.fq.gz`
echo "Processing sample ${i}"
echo "RUN CMD: inner_distance.py -i $workdir/3.map/hisat2/${i}.bam  -r $GENE_BED  -o ${i}_inner_size"
inner_distance.py -i $workdir/3.map/hisat2/${i}.bam  -r $GENE_BED  -o ${i}_inner_size
done

#2.基因覆盖情况，RNA是否降解（报错了）
for fn in /u3/hekun/rnaseq/2.data_qc/*_1.fastp.fq.gz;
do 
i=`basename ${fn} _1.fastp.fq.gz`
echo "Processing sample ${i}"

echo "RUN CMD: geneBody_coverage.py -r $GENE_BED -i $workdir/3.map/hisat2/${i}.bam  -o ${i}.genebody"
geneBody_coverage.py -r $GENE_BED -i $workdir/3.map/hisat2/${i}.bam  -o ${i}.genebody
done



###############################################################
#(一直安装不上Htseq-count，所以接下里哀先借用生信技能树的featureCounts定向)
# 基因表达定量及结果展示
## 采用Htseq-count 对已知的基因进行表达定量(一直安装不上Htseq-count，所以接下里哀先借用生信技能树的featureCounts定向)
##################################################################
conda create -n htseq-count
conda activate htseq-count
conda install -c bioconda htseq-count

cd $workdir/
mkdir 4.expression
cd 4.expression

#--stranded yes or no   文库类型设置，是否为链特异性文库
#--mode 包含union，intersection-strict，intersection_nonempty,作者一般推荐union


for fn in /u3/hekun/rnaseq/2.data_qc/*_1.fastp.fq.gz;
do 
i=`basename ${fn} _1.clean.fastp.fq.gz`
echo "Processing sample ${i}"
echo "RUN CMD: htseq-count --format bam --order pos --mode intersection-strict \
--stranded reverse --minaqual 1 --type exon \
--idattr gene_id $workdir/3.map/hisat2/${i}.bam $GTF > ${i}_gene.tsv"

htseq-count --format bam --order pos --mode union \
--stranded no --minaqual 8 --type exon \
--idattr gene_id $workdir/3.map/hisat2/${i}.bam $GTF > ${i}_gene.tsv
done

### 合并不同样品表达定量结果，方便后续基因表达分析，试着写一个循环-l 合并文件的列名
python $scriptdir/merge_gene_count.py -p all_gene_count \
-f  normal_rep1_gene.tsv -l normal_rep1 \
-f  normal_rep2_gene.tsv -l normal_rep2 \
-f  normal_rep3_gene.tsv -l normal_rep3 \
-f  tumor_rep1_gene.tsv -l tumor_rep1 \
-f  tumor_rep2_gene.tsv -l tumor_rep2 \
-f  tumor_rep3_gene.tsv -l tumor_rep3

#下列循环可以批量产生上面的内容
for fn in /u3/hekun/rnaseq/2.data_qc/*_1.clean.fastp.fq.gz;
do 
i=`basename ${fn} _1.clean.fastp.fq.gz`
echo "-f  ${i}_gene.tsv -l ${i} \
" \
>> merge.txt
done



(一直安装不上Htseq-count，所以接下里哀先借用生信技能树的featureCounts定向)
## 定义输入输出文件夹
gtf=/u3/hekun/rnaseq/ref/Taestivum.58.gtf
inputdir=/u3/hekun/rnaseq/3.map/hisat2
# featureCounts对bam文件进行计数
featureCounts -T 6 -p -t exon -g gene_id -a $gtf -o all.id.txt $inputdir/*.bam &
# 对定量结果质控
multiqc all.id.txt.summary &
# 得到表达矩阵
# 处理表头，/home/t_rna/要换成自己的路径
less -S all.id.txt |grep -v '#' |cut -f 1,7- |sed 's#/u3/hekun/rnaseq/3.map/hisat2/##g' |sed 's#.bam##g' >raw_counts.txt










###基因表达定量结果展示
#1.各样本表达量密度 图
#2.各样本表达量box分布图
#4 各样本表达相关性分析热图与聚类图

usr/bin/Rscript $scriptdir/fpkm_and_plot.R -i all_gene_count.tsv  -l $GENE_LENGTH  -o ./


##################################################
## 基因差异表达分析(DESeq2)，并绘制火山图与MA图 ，以及差异基因表达热图
###################################################
cd $workdir/
mkdir 5.deg
cd 5.deg
#创建分组文件normal_vs_tumor.compare.txt，用于查找不同分组之间差异表达基因
#ID	group
#normal_rep1	normal
#normal_rep2	normal
#normal_rep3	normal
#tumor_rep1	tumor
#tumor_rep2	tumor
#tumor_rep3	tumor

Rscript $scriptdir/deseq_analysis.r -i $workdir/4.expression/all_gene_count.tsv -g normal_vs_tumor.compare.txt  -k $workdir/4.expression/all_gene_fpkm.tsv -r normal --fdr 1 --fc 1.1 -p normal_vs_tumor

#差异表达基因热图
Rscript $scriptdir/heatmap.r -i $workdir/4.expression/all_gene_fpkm.tsv -l normal_vs_tumor.DEG.final.tsv -p normal_vs_tumor.deg_gene_heatmap -o ./

############################################
#2023.10.9 更新 差异分析,只需要准备一个meta分组文件即可：-t 设置分组列名，该列名下比较分组 --case tumor --control  normal

#, S1 为免疫侵润高组 DESeq2
Rscript $scriptdir/deseq_analysis1.r  -i $workdir/4.expression/all_gene_count.tsv \
   --fdr 0.01 --fc 2   -m normal_vs_tumor.compare.txt -t group\
   --case tumor --control  normal -p normal_vs_tumor.deseq

#免疫亚型之间差异表达, S1 为免疫侵润高组 edgeR

Rscript $scriptdir/edger_analysis1.r  -i $workdir/4.expression/all_gene_count.tsv \
   --fdr 0.01 --fc 2   -m normal_vs_tumor.compare.txt -t group\
   --case tumor --control  normal -p normal_vs_tumor.edger

##################################################
#差异表达基因功能富集分析
##################################################
cd $workdir/
mkdir 6.enrich
cd 6.enrich

#GO,KEGG 富集分析以下脚本只支持模式物种见表格：https://www.omicsclass.com/article/1244
#非模式物种参考：http://guangchuangyu.github.io/cn/2017/07/clusterprofiler-maize/
Rscript $scriptdir/enrichGO_pip.R --deg.file $workdir/5.deg/normal_vs_tumor.DEG.final.tsv -o GO/ -n normal_vs_tumor --pvalueCutoff 0.5 --ann.db org.Hs.eg.db  --idtype ENSEMBL --totype ENTREZID

Rscript $scriptdir/enrichKEGG_pip.R --deg.file $workdir/5.deg/normal_vs_tumor.DEG.final.tsv -o KEGG -n normal_vs_tumor --pvalueCutoff 1 --ann.db org.Hs.eg.db --organism hsa --idtype ENSEMBL --totype ENTREZID

#GSEA分析只有人类数据
#更多GSEA功能富集数据下载：http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
wget -c https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.1/c2.cp.kegg.v7.1.entrez.gmt
wget -c https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.1/c5.bp.v7.1.entrez.gmt

Rscript $scriptdir/enrichGSEA_pip.R   --all.deg.file $workdir/5.deg/normal_vs_tumor.all.tsv --gmtfile $workdir//6.enrich/c5.bp.v7.1.entrez.gmt -o GSEA -n normal_vs_tumor_BP -p 0.3 --idtype ENSEMBL

Rscript $scriptdir/enrichGSEA_pip.R   --all.deg.file $workdir/5.deg/normal_vs_tumor.all.tsv --gmtfile $workdir//6.enrich/c2.cp.kegg.v7.1.entrez.gmt -o GSEA -n normal_vs_tumor_KEGG -p 0.8 --idtype ENSEMBL




############################################
#非模式物种功能注释与GO KEGG富集分析
###########################################

cd $workdir/
mkdir 7.all_enrich
cd 7.all_enrich

#非模式物种利用eggnog数据库做GO KEGG注释
#eggnog 批量注释：https://www.omicsclass.com/article/1515

#############################################################################
#数据准备  以拟南芥数据为例
########################################################################
#下载拟南芥基因组信息
#wget -c  ftp://ftp.ensemblgenomes.org/pub/plants/release-41/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz
wget -c ftp://ftp.ensemblgenomes.org/pub/plants/release-41/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.41.gff3.gz
wget -c http://ftp.ensemblgenomes.org/pub/plants/release-41/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz


#观察数据，ID保持一致,也就是gff中第9列，ID标签和parent标签与蛋白序列里面的ID是否一致；
#处理GFF 文件里面ID中一些不必要的信息，gene:  transcript: 删除；与蛋白质中的ID保持一致：Arabidopsis_thaliana.TAIR10.pep.all.fa

genome=Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gff=Arabidopsis_thaliana.TAIR10.41.gff3.gz

#保留蛋白编码基因
agat_sp_filter_feature_by_attribute_value.pl --gff  $gff --attribute biotype --value protein_coding -t '!' -o Ath.protein_coding.gff3
#保留最长转录本
agat_sp_keep_longest_isoform.pl --gff Ath.protein_coding.gff3 -o Ath.longest_isoform.gff3
#Ensembl 下载的基因组 会在基因和mRNA ID前面加gene: 和 transcript:这里给去掉
sed -i 's/gene://' Ath.longest_isoform.gff3
sed -i 's/transcript://'  Ath.longest_isoform.gff3
#提取cds序列
gff3_file_to_proteins.pl  --gff3  Ath.longest_isoform.gff3 --fasta $genome  --seqType CDS >Ath.cds.fa
#提取pep序列
gff3_file_to_proteins.pl  --gff3  Ath.longest_isoform.gff3  --fasta  $genome --seqType prot >Ath.pep.fa
#提取cds 转换成OrthoFinder 需要的序列格式 基因ID作为序列ID
perl $scriptdir/get_gene_longest_fa_for_OrthoFinder.pl -l 0 --gff Ath.longest_isoform.gff3 --fa Ath.cds.fa -p Ath.gene.cds
#提取pep 转换成OrthoFinder 需要的序列格式 基因ID作为序列ID
perl $scriptdir/get_gene_longest_fa_for_OrthoFinder.pl -l 0 --gff Ath.longest_isoform.gff3 --fa Ath.pep.fa -p Ath.gene.pep

#######做eggNOG 注释###########################



#使用帮助：https://github.com/eggnogdb/eggnog-mapper/wiki/

#最新数据库下载整理方法：https://www.omicsclass.com/article/1515

#设置数据库路径
eggnogdb=/work/database/eggNOG/latest/  #设置数据库所在路径


#查看 物种分类：emapper.py --list_taxa  ：在线网站：http://eggnog-mapper.embl.de/
emapper.py -i Ath.gene.pep.fasta --itype proteins -o  Ath.eggnog -m diamond \
    --cpu 10 --seed_ortholog_evalue 1e-5 --override --tax_scope Viridiplantae \
    --dmnd_db $eggnogdb/eggnog_proteins.dmnd --data_dir $eggnogdb

# 第二步：注释结果统计

Rscript  $scriptdir/eggNOG/eggnog_stat.r -o Ath.eggnog  \
  --emapper_anno Ath.eggnog.emapper.annotations -f Ath.gene.pep.fasta  \
  --keep.Pathway $scriptdir/eggNOG/pathways_plant.txt


#基因列表

#KEGG富集分析
Rscript $scriptdir/eggNOG/kegg_enrich.r -l gene.list -n Ath.eggnog/kegg.pathway2name.tsv -g Ath.eggnog/kegg.gene2pathway.tsv

#GO富集分析
Rscript $scriptdir/eggNOG/go_enrich.r -l gene.list -d Ath.eggnog/GO_orgdb/

################################
#转录因子分析
################################

cd $workdir/
mkdir 8.TF_ann
cd 8.TF_ann

#下载植物转录因子数据库 ，建议下载对应物种的：http://planttfdb.gao-lab.org/download.php#tf_idseq

wget -c http://planttfdb.gao-lab.org/download/seq/PlantTFDB-all_TF_pep.fas.gz

#解压

gunzip *gz

#建库
makeblastdb -in PlantTFDB-all_TF_pep.fas  -dbtype prot -title PlantTFDB-all_TF_pep.fas

#提取基因列表对应的蛋白序列：
python $scriptdir/get_fa_by_id.py  -i ../7.all_enrich/gene.list \
  -f ../7.all_enrich/Ath.gene.pep.fasta -p query
#转录因子blastp注释
blastp -query query.fa -db PlantTFDB-all_TF_pep.fas  -max_hsps 1  \
  -num_alignments 1 -evalue 1e-5 -num_threads 2   \
  -outfmt "6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle"\
  -out tf.blast.out


#################################
#蛋白互作网络分析
################################


cd $workdir/
mkdir 9.ppi_network
cd 9.ppi_network

python $scriptdir/get_fa_by_id.py  -i ID-info.txt \
  -f ../7.all_enrich/Ath.gene.pep.fasta -p pep






