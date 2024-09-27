### 转录组测序分析

#### 1.工作目录管理

这一部分非常非常非常重要，拥有一个优秀的工作习惯比什么都重要

##### 1.1 目录一览

```bash
## 示例如下：
├── database # 数据库存放目录，包括参考基因组，注释文件，公共数据库等
├── project  # 项目分析目录
    └── Human-16-Asthma-Trans #具体项目
        ├── data # 数据存放目录
        │   ├── cleandata # 过滤后的数据
           	│	├── trim_galore # trim_galore过滤
		   	│	└── fastp	    # fastp过滤
        │   └── rawdata # 原始数据
        ├──  Mapping # 比对目录
        │   ├── Hisat2 # Hisat比对
        │   └── Subjunc # subjunc比对
        └── Expression # 定量
            ├── featureCounts # featureCounts
            └── Salmon # salmon定量
        └── Diff_analysis # 差异分析           
```



##### 1.2 详细命令

```bash
# 进入到个人目录
cd ~

## 1.建立数据库目录：在数据库下建立参考基因组数据库，注意命名习惯：参考基因组版本信息
mkdir -p database/GRCh38.105

## 2.建立项目分析目录
mkdir project
cd project
mkdir Human-16-Asthma-Trans # 注意项目命名习惯：物种-样本数-疾病-分析流程
cd Human-16-Asthma-Trans

# 建立数据存放目录
mkdir -p  data/rawdata  data/cleandata/trim_galore  data/cleandata/fastp
# 建立比对目录
mkdir -p Mapping/Hisat2  Mapping/Subjunc
# 建立定量目录
mkdir -p Expression/featureCounts  Expression/Salmon
# 创建差异分析目录
mkdir Diff_analysis
# 查看整个分析目录准备结构
tree
.
├── data
│   ├── cleandata
│   │   ├── fastp
│   │   └── trim_galore
│   └── rawdata
├── Diff_analysis
├── Expression
│   ├── featureCounts
│   └── Salmon
└── Mapping
    ├── Hisat2
    └── Subjunc
```



##### 课后习题答案

1.统计reads_1.fq文件种共有多少条reads？

```bash
# 答案不只一种，看看你能用集中方法算出来
# NR表示行号
# %表示取余数
zless  SRR1039510_1.fastq.gz | grep "@SRR" -c
zless  SRR1039510_1.fastq.gz | grep '^@SRR' |wc -l
zless -S SRR1039510_1.fastq.gz | paste - - - - |wc -l
zless  SRR1039510_1.fastq.gz |wc -l | awk '{print $0/4}'
zless -S SRR1039510_1.fastq.gz |awk '{ if(NR%4==2) {print} }' |wc -l

# sed 版本 课后习题
```



2.输出reads_1.fq文件中所有的序列ID（即第一行）

```bash
zless  SRR1039510_1.fastq.gz | grep '^@SRR'  |less -S
zless  SRR1039510_1.fastq.gz | paste - - - - |cut -f 1 |less -S
zless -S SRR1039510_1.fastq.gz |awk '{if(NR%4==1){print}}' |less -S
```



3.输出SRR1039510_1.fastq.gz文件中所有的序列（即第二行）

```bash
zless  SRR1039510_1.fastq.gz | paste - - - - |cut -f 2 |less -S
zless -S SRR1039510_1.fastq.gz |awk '{if(NR%4==2){print}}' |less -S 
```



4.统计SRR1039510_1.fastq.gz碱基总数

```bash
# 简单版本
zless -S SRR1039510_1.fastq.gz |paste - - - - |cut -f 2 |tr -d '\n' |wc -m
zless -S SRR1039510_1.fastq.gz |paste - - - - |cut -f 2 |grep -o [ATCGN] |wc -l

# awk的高阶用法：BEGIN END模块
zless -S SRR1039510_1.fastq.gz |awk '{ if(NR%4==2){print} }' | awk 'BEGIN {num=0} {num=num+length($0)}  END{ print "num="num}'
```





#### 2.数据质控

##### 2.1 fastqc

目标：使用fastqc对原始数据进行质量评估

```bash
# 激活conda环境
conda activate rna

# 连接数据到自己的文件夹
cd $HOME/project/Human-16-Asthma-Trans/data/rawdata
ln -s /home/t_rna/data/airway/fastq_raw25000/*gz  ./

# 使用FastQC软件对单个fastq文件进行质量评估，结果输出到qc/文件夹下
nohup fastqc -t 6 -o ./ SRR*.fastq.gz >qc.log &

# 使用MultiQc整合FastQC结果
# 使用绝对路径运行
multiqc=/home/t_rna/miniconda3/envs/rna/bin/multiqc
fastqc=/home/t_rna/miniconda3/envs/rna/bin/fastqc
fq_dir=$HOME/project/Human-16-Asthma-Trans/data/rawdata
outdir=$HOME/project/Human-16-Asthma-Trans/data/rawdata

# 使用绝对路径运行
$fastqc -t 6 -o $outdir ${fq_dir}/SRR*.fastq.gz >${fq_dir}/qc.log

# 报告整合
$multiqc $outdir/*.zip -o $outdir/ >${fq_dir}/multiqc.log
```





#### 3.数据过滤

目标：使用两款软件过滤fq文件，查看过滤前后的区别

##### 3.1 trim_galore过滤

```bash
# 激活小环境
conda activate rna
# 新建文件夹trim_galore
cd $HOME/project/Human-16-Asthma-Trans/data/cleandata/trim_galore

# 先生成一个变量,为样本ID
ls $HOME/project/Human-16-Asthma-Trans/data/rawdata/*_1.fastq.gz | awk -F'/' '{print $NF}' | cut -d'_' -f1 >ID

# 多个样本 vim trim_galore.sh，以下为sh的内容
rawdata=$HOME/project/Human-16-Asthma-Trans/data/rawdata
cleandata=$HOME/project/Human-16-Asthma-Trans/data/cleandata/trim_galore
ID=$HOME/project/Human-16-Asthma-Trans/data/cleandata/trim_galore/ID
cat $ID | while read id
do
  trim_galore -q 20 --length 20 --max_n 3 --stringency 3 --fastqc --paired -o ${cleandata}    ${rawdata}/${id}_1.fastq.gz     ${rawdata}/${id}_2.fastq.gz
done

# 提交任务到后台 可以 用bash或者sh都行 
nohup bash trim_galore.sh >trim_galore.log &

# 使用MultiQc整合FastQC结果
multiqc *.zip


## ==============================================
## 补充技巧：使用掐头去尾获得样本ID
ls $rawdata/*_1.fastq.gz | while read id
do
name=${id##*/}
name=${name%_*}
echo "trim_galore -q 20 --length 20 --max_n 3 --stringency 3 --fastqc --paired -o ${cleandata} ${rawdata}/${name}_1.fastq.gz ${rawdata}/${name}_2.fastq.gz "
done
```



##### 3.2 fastp过滤

```bash
cd $HOME/project/Human-16-Asthma-Trans/data/cleandata/fastp

# 单样本：
fastp -l 20 -q 20 --compression=6 \
  -i ../../rawdata/SRR1039510_1.fastq.gz \
  -I ../../rawdata/SRR1039510_2.fastq.gz  \
  -o ./SRR1039510_clean_1.fq.gz  \
  -O ./SRR1039510_clean_2.fq.gz \
  -R SRR1039510 \
  -h SRR1039510.fastp.html \
  -j ./SRR1039510.fastp.json

# 定义文件夹：vim fastp.sh
cleandata=$HOME/project/Human-16-Asthma-Trans/data/cleandata/fastp/
rawdata=$HOME/project/Human-16-Asthma-Trans/data/rawdata/

## print fastp command
cat ../trim_galore/ID | while read id
do echo "
fastp -l 20 -q 20 --compression=6 \
  -i ${rawdata}/${id}_1.fastq.gz \
  -I ${rawdata}/${id}_2.fastq.gz \
  -o ${cleandata}/${id}_clean_1.fq.gz \
  -O ${cleandata}/${id}_clean_2.fq.gz \
  -R ${cleandata}/${id} \
  -h ${cleandata}/${id}.fastp.html \
  -j ${cleandata}/${id}.fastp.json 
  "
done

## yunxing
cat ../trim_galore/ID | while read id
do
fastp -l 20 -q 20 --compression=6 \
  -i ${rawdata}/${id}_1.fastq.gz \
  -I ${rawdata}/${id}_2.fastq.gz \
  -o ${cleandata}/${id}_clean_1.fq.gz \
  -O ${cleandata}/${id}_clean_2.fq.gz \
  -R ${cleandata}/${id} \
  -h ${cleandata}/${id}.fastp.html \
  -j ${cleandata}/${id}.fastp.json 
done

# 运行fastp脚本
nohup bash fastp.sh >fastp.log &

## 整合
multiqc *json
```



##### 3.3 数据过滤前后比较

```bash
# 进入过滤目录
cd $HOME/project/Human-16-Asthma-Trans/data/cleandata/trim_galore

# 原始数据
zcat $rawdata/SRR1039510_1.fastq.gz | paste - - - - > raw.txt

#  过滤后的数据
zcat SRR1039510_1_val_1.fq.gz |paste - - - - > trim.txt
awk '(length($4)<63){print$1}' trim.txt > Seq.ID

head -n 100 Seq.ID > ID100
grep -w -f ID100 trim.txt | awk '{print$1,$4}' > trim.sm
grep -w -f ID100 raw.txt | awk '{print$1,$4}' > raw.sm
paste raw.sm trim.sm | awk '{print$2,$4}' | tr ' ' '\n' |less -S
```



#### 4.数据比对

目标：使用两个软件对fq数据进行比对，得到比对文件**sam/bam**，并探索比对结果。

##### 4.1 参考基因组

```bash
## 参考基因组准备:注意参考基因组版本信息
# 下载，Ensembl：http://asia.ensembl.org/index.html
# http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/

# 进入到参考基因组目录
mkdir -p $HOME/database/GRCh38.111
cd $HOME/database/GRCh38.111

# 下载基因组序列axel  curl  
nohup axel -n 100 https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz >dna.log &

# 下载转录组序列
nohup axel -n 100  http://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz >rna.log &

# 下载基因组注释文件
nohup axel -n 100 https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz >gtf.log &

nohup axel -n 100 https://ftp.ensembl.org/pub/release-111/gff3/homo_sapiens/Homo_sapiens.GRCh38.111.gff3.gz >gff.log&

# 上述文件下载完整后，再解压；否则文件不完整就解压会报错
# 再次强调，一定要在文件下载完后再进行解压！！！
nohup gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz Homo_sapiens.GRCh38.cdna.all.fa.gz >unzip.log &
```



##### 4.2 课后作业

**1.fastq与fasta文件转换**

```sh
# 答案1
zless -S SRR1039511_1_val_1.fq.gz |awk '{ if(NR%4==1){print">" substr($0,2)} if(NR%4==2){print}  }' | less -S

# 答案2
zless -S SRR1039510_1_val_1.fq.gz |paste - - - - |cut -f 1,2 |tr '@' '>' |tr '\t' '\n' |less -S
```



**2.从gff或者gft文件中获取基因的ID与symbol对应关系，以及biotype类型**

应用：ID与symbol转换本地化，不依赖于第三方工具和软件包，并可以根据biotype类型区分mRNA，lncRNA以及miRNA等信息。

```bash
# 从gff或者gft文件中获取ID与symbol对应关系，以及biotype类型
zless -S Homo_sapiens.GRCh38.104.chr.gtf.gz  |awk -F'\t' '{if($3=="gene"){print$9}}' |awk -F';' '{print$1,$3,$5}' |awk '{print$2"\t"$4"\t"$6}' |sed 's/"//g' |grep 'protein_coding' >protein_coding_id2name.xls
```



##### 4.3 比对:hisat2

````bash
## ----1.构建索引
# 进入参考基因组目录
cd $HOME/database/GRCh38.105
# Hisat2构建索引，构建索引时间比较长，建议提交后台运行，一般会运行20多分钟左右
## 后续索引可直接使用服务器上已经构建好的进行练习
# vim Hisat2Index.sh
mkdir Hisat2Index
fa=Homo_sapiens.GRCh38.dna.primary_assembly.fa
fa_baseName=GRCh38.dna
hisat2-build -p 5 -f ${fa} Hisat2Index/${fa_baseName}

# 运行
nohup bash Hisat2Index.sh >Hisat2Index.sh.log &

## ----2.比对
# 进入比对文件夹
cd $HOME/project/Human-16-Asthma-Trans/Mapping/Hisat2

## 单个样本比对，步骤分解
index=/home/t_rna/database/GRCh38.104/Hisat2Index/GRCh38.dna
inputdir=$HOME/project/Human-16-Asthma-Trans/data/cleandata/trim_galore/
outdir=$HOME/project/Human-16-Asthma-Trans/Mapping/Hisat2

hisat2 -p 5  -x  ${index} \
	   -1 ${inputdir}/SRR1039510_1_val_1.fq.gz \
       -2 ${inputdir}/SRR1039510_2_val_2.fq.gz \
       -S ${outdir}/SRR1039510.Hisat_aln.sam

## ----3.sam转bam
samtools sort -@ 5 -o SRR1039510.Hisat_aln.sorted.bam SRR1039510.Hisat_aln.sam

## ----4.对bam建索引
samtools index SRR1039510.Hisat_aln.sorted.bam SRR1039510.Hisat_aln.sorted.bam.bai


# 多个样本批量进行比对，排序，建索引
# Hisat.sh内容： 注意命令中的-，表示占位符，表示|管道符前面的输出。
## 此处索引直接使用服务器上已经构建好的进行练习
# vim Hisat.sh
index=/home/t_rna/database/GRCh38.104/Hisat2Index/GRCh38.dna
inputdir=$HOME/project/Human-16-Asthma-Trans/data/cleandata/trim_galore/
outdir=$HOME/project/Human-16-Asthma-Trans/Mapping/Hisat2

cat ../../data/cleandata/trim_galore/ID | while read id
do
hisat2 -p 5 -x ${index} -1 ${inputdir}/${id}_1_val_1.fq.gz -2 ${inputdir}/${id}_2_val_2.fq.gz 2>${id}.log  | samtools sort -@ 3 -o ${outdir}/${id}.Hisat_aln.sorted.bam - &&  samtools index ${outdir}/${id}.Hisat_aln.sorted.bam ${outdir}/${id}.Hisat_aln.sorted.bam.bai
done

# 统计比对情况
multiqc -o ./ SRR*log

# 提交后台运行
nohup bash Hisat.sh >Hisat.log &
````



##### 4.4 比对:subjunc

```bash
## ----1.构建索引
# 进入参考基因组目录
cd $HOME/database/GRCh38.105
# subjunc构建索引，构建索引时间比较长大约40分钟左右，建议携程sh脚本提交后台运行
## 后续索引可直接使用服务器上已经构建好的进行练习
# vim SubjuncIndex.sh
mkdir SubjuncIndex
fa=Homo_sapiens.GRCh38.dna.primary_assembly.fa
fa_baseName=GRCh38.dna
subread-buildindex -o SubjuncIndex/${fa_baseName} ${fa}

# 运行
nohup bash SubjuncIndex.sh >SubjuncIndex.sh.log &


## ----比对
# 进入文件夹
cd $HOME/project/Human-16-Asthma-Trans/Mapping/Subjunc

# vim subjunc.sh
index=/home/t_rna/database/GRCh38.104/SubjuncIndex/GRCh38.dna
inputdir=$HOME/project/Human-16-Asthma-Trans/data/cleandata/trim_galore
outdir=$HOME/project/Human-16-Asthma-Trans/Mapping/Subjunc

cat ../../data/cleandata/trim_galore/ID | while read id
do
  subjunc -T 5 -i ${index} -r ${inputdir}/${id}_1_val_1.fq.gz -R ${inputdir}/${id}_2_val_2.fq.gz -o ${outdir}/${id}.Subjunc.bam 1>${outdir}/${id}.Subjunc.log 2>&1 && samtools sort -@ 6 -o ${outdir}/${id}.Subjunc.sorted.bam  ${outdir}/${id}.Subjunc.bam && samtools index ${outdir}/${id}.Subjunc.sorted.bam ${outdir}/${id}.Subjunc.sorted.bam.bai
done

# 运行
nohup bash subjunc.sh >subjunc.log &
```



#### 5.sam/bam应用

##### 5.1 统计比对结果

```bash
# 单个样本
samtools flagstat -@ 3 SRR1039510.Hisat_aln.sorted.bam

# 多个样本，vim flagstat.sh
ls *.sorted.bam | while read id
do
  samtools flagstat -@ 10 ${id} > ${id/bam/flagstat}
done
# 质控
multiqc -o ./  *.flagstat

# 运行
nohup sh flagstat.sh >flagstat.log &
```



##### 5.2 samtools工具使用

```bash
##----view查看bam文件
samtools view SRR1039510.Hisat_aln.sorted.bam
samtools view -H SRR1039510.Hisat_aln.sorted.bam
samtools view -h SRR1039510.Hisat_aln.sorted.bam

##----index对bam文件建索引
samtools index SRR1039510.Hisat_aln.sorted.bam SRR1039510.Hisat_aln.sorted.bam.bai

##----flagstat统计比对结果
samtools flagstat -@ 3 SRR1039510.Hisat_aln.sorted.bam

##----sort排序 sam转bam并排序
samtools sort -@ 3 -o SRR1039510.Hisat_aln.sorted.bam SRR1039510.Hisat_aln.sam

##----depth统计测序深度
# 得到的结果中，一共有3列以指标分隔符分隔的数据，第一列为染色体名称，第二列为位点，第三列为覆盖深度
samtools depth SRR1039510.Hisat_aln.sorted.bam >SRR1039510.Hisat_aln.sorted.bam.depth.txt

##----计算某一个基因的测序深度
# 找到基因的坐标
zless -S Homo_sapiens.GRCh38.95.gff3.gz |awk '{if($3=="gene")print}' |grep 'ID=gene:ENSG00000186092' |awk '{print $1"\t"$4"\t"$5}' >ENSG00000186092.bed
samtools depth  -b  ENSG00000186092.bed SRR1039510.Hisat_aln.sorted.bam >ENSG00000186092.bed.depth

# 如何找到多比对的reads，flag值的理解
# (0x100) 代表着多比对情况，所以直接用samtools view -f 0x100可以提取 multiple比对的 情况
```



#### 6.表达定量

##### 6.1 featureCounts

使用featureCounts对bam文件进行定量。

```bash
cd $HOME/project/Human-16-Asthma-Trans/Expression/featureCounts

## 定义输入输出文件夹
gtf=/home/t_rna/database/GRCh38.104/Homo_sapiens.GRCh38.104.chr.gtf.gz
inputdir=$HOME/project/Human-16-Asthma-Trans/Mapping/Hisat2/

# featureCounts对bam文件进行计数
featureCounts -T 6 -p --countReadPairs -t exon -g gene_id -a $gtf -o all.id.txt $inputdir/*.sorted.bam

# 对定量结果质控
multiqc all.id.txt.summary

# 得到表达矩阵
# 处理表头，/home/t_rna/要换成自己的路径
less -S all.id.txt |grep -v '#' |cut -f 1,7- |sed 's#/home/t_rna/project/Human-16-Asthma-Trans/Mapping/Hisat2//##g' |sed 's#.Hisat_aln.sorted.bam##g' >raw_counts.txt

sed s///
sed s###
sed s%%%


# 列对齐显示
head raw_counts.txt  |column -t
```



##### 6.2 salmon定量

```bash
##----构建索引
## 后续索引可直接使用服务器上已经构建好的进行练习
cd $HOME/database/GRCh38.105
nohup salmon index -t Homo_sapiens.GRCh38.cdna.all.fa -i salmon_index >salmon-index.log &

cd $HOME/project/Human-16-Asthma-Trans/Expression/Salmon

##----运行
# 编写脚本，使用salmon批量对目录下所有fastq文件进行定量
# vim salmon.sh
index=/home/t_rna/database/GRCh38.104/salmon_index
input=$HOME/project/Human-16-Asthma-Trans/data/cleandata/trim_galore
outdir=$HOME/project/Human-16-Asthma-Trans/Expression/Salmon
cat ../../data/cleandata/trim_galore/ID |while read id 
do
  salmon quant -i ${index} -l A -1 ${input}/${id}_1_val_1.fq.gz -2 ${input}/${id}_2_val_2.fq.gz -p 5 -o ${outdir}/${id}.quant
done


##----合并表达矩阵
# 原始count值矩阵
# --quants：ls -d *quant |tr '\n' ',' |sed 's/,$//' |awk '{print "{" $0 "}" }'
# --quants：ls -d *quant |sed ':a;N;s/\n/,/;t a;'|awk '{print "{"$0"}"}'
# --quants：ls -d *quant |xargs  |tr ' ' ',' |awk '{print "{"$0"}"}'
# --names：ls -d *quant |tr '\n' ',' |sed 's/,$//' |awk '{print "{"$0"}"}' |sed 's/.quant//g'

## 完整版：ls -d *quant |tr '\n' ',' |sed 's/,$//' |awk '{print "{"$0"}"}'  | perl -ne 'chomp;$a=$_;$a=~s/\.quant//g;print"salmon quantmerge --quants $_ --names $a --column numreads  --output raw_count.txt \n";'

salmon quantmerge --quants {SRR1039510.quant,SRR1039511.quant,SRR1039512.quant} --names {SRR1039510,SRR1039511,SRR1039512} --column numreads  --output raw_count.txt

# tpm矩阵
salmon quantmerge --quants {SRR1039510.quant,SRR1039511.quant,SRR1039512.quant} --names {SRR1039510,SRR1039511,SRR1039512} --column tpm  --output tpm.txt

# 后台运行脚本
nohup bash salmon.sh >salmon.log &

```



#### 7.差异表达分析

```bash
cd $HOME/project/Human-16-Asthma-Trans/Diff_analysis


#DEG.sh内容
/home/t_rna/miniconda3/envs/R4/bin/Rscript edgeR.R --count filter_count.txt --group group.txt --comp Dex_vs_untreated --fc 1.5 --pvalue 0.05 --od ./Test1
/home/t_rna/miniconda3/envs/R4/bin/Rscript edgeR.R --count filter_count.txt --group group.txt --comp Dex_vs_untreated --fc 2 --pvalue 0.01 --od ./Test2

nohup sh DEG.sh >DEG.sh.log &
```



