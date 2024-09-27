### step1： 使用fastqc对原始数据做质量评估
# 使用绝对路径运行
multiqc=/home/t_rna/miniconda3/envs/rna/bin/multiqc
fastqc=/home/t_rna/miniconda3/envs/rna/bin/fastqc
fq_dir=$HOME/project/Human-16-Asthma-Trans/data/rawdata
outdir=$HOME/project/Human-16-Asthma-Trans/data/rawdata

# 评估
echo "Step1 QC: Start "
echo "$fastqc -t 6 -o $outdir ${fq_dir}/SRR*.fastq.gz >${fq_dir}/qc.log"
$fastqc -t 6 -o $outdir ${fq_dir}/SRR*.fastq.gz >${fq_dir}/qc.log

# 报告整合
echo "$multiqc $outdir/*.zip -o $outdir/ >${fq_dir}/multiqc.log"
$multiqc $outdir/*.zip -o $outdir/ >${fq_dir}/multiqc.log
echo "Step1 QC: Done!"


### step2：  trim_galore过滤
cd $HOME/project/Human-16-Asthma-Trans/data/cleandata/trim_galore

# 先生成一个变量,为样本ID
ls $HOME/project/Human-16-Asthma-Trans/data/rawdata/*_1.fastq.gz | awk -F'/' '{print $NF}' | cut -d'_' -f1 >ID

# 多个样本 vim trim_galore.sh，以下为sh的内容
rawdata=$HOME/project/Human-16-Asthma-Trans/data/rawdata
cleandata=$HOME/project/Human-16-Asthma-Trans/data/cleandata/trim_galore
ID=$HOME/project/Human-16-Asthma-Trans/data/cleandata/trim_galore/ID
# echo 
cat $ID | while read id
do
  echo "trim_galore -q 20 --length 20 --max_n 3 --stringency 3 --fastqc --paired -o ${cleandata}    ${rawdata}/${id}_1.fastq.gz     ${rawdata}/${id}_2.fastq.gz"
done

cat $ID | while read id
do
  trim_galore -q 20 --length 20 --max_n 3 --stringency 3 --fastqc --paired -o ${cleandata}    ${rawdata}/${id}_1.fastq.gz     ${rawdata}/${id}_2.fastq.gz
done



### step3：  hisat比对
index=/home/t_rna/database/GRCh38.104/Hisat2Index/GRCh38.dna
inputdir=$HOME/project/Human-16-Asthma-Trans/data/cleandata/trim_galore/
outdir=$HOME/project/Human-16-Asthma-Trans/Mapping/Hisat2

cat $ID | while read id
do
  echo "hisat2 -p 5 -x ${index} -1 ${inputdir}/${id}_1_val_1.fq.gz -2 ${inputdir}/${id}_2_val_2.fq.gz 2>${id}.log  | samtools sort -@ 3 -o ${outdir}/${id}.Hisat_aln.sorted.bam - &&  samtools index ${outdir}/${id}.Hisat_aln.sorted.bam ${outdir}/${id}.Hisat_aln.sorted.bam.bai"
done

cat $ID | while read id
do
  hisat2 -p 5 -x ${index} -1 ${inputdir}/${id}_1_val_1.fq.gz -2 ${inputdir}/${id}_2_val_2.fq.gz 2>${id}.log  | samtools sort -@ 3 -o ${outdir}/${id}.Hisat_aln.sorted.bam - &&  samtools index ${outdir}/${id}.Hisat_aln.sorted.bam ${outdir}/${id}.Hisat_aln.sorted.bam.bai
done

# 统计比对情况
multiqc -o ./ SRR*log



### step4：  featurecount定量
gtf=/home/t_rna/database/GRCh38.104/Homo_sapiens.GRCh38.104.chr.gtf.gz
inputdir=$HOME/project/Human-16-Asthma-Trans/Mapping/Hisat2/
outdir=$HOME/project/Human-16-Asthma-Trans/Expression/featureCounts

# featureCounts对bam文件进行计数
echo "featureCounts -T 6 -p --countReadPairs -t exon -g gene_id -a $gtf -o $outdir/all.id.txt $inputdir/*.sorted.bam"

featureCounts -T 6 -p --countReadPairs -t exon -g gene_id -a $gtf -o $outdir/all.id.txt $inputdir/*.sorted.bam

# 对定量结果质控
multiqc $outdir/all.id.txt.summary




### step5： 差异表达分析
cd $HOME/project/Human-16-Asthma-Trans/Diff_analysis

/home/t_rna/miniconda3/envs/R4/bin/Rscript edgeR.R --count filter_count.txt --group group.txt --comp Dex_vs_untreated --fc 1.5 --pvalue 0.05 --od ./Test1










