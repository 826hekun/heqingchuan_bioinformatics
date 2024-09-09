

#################################################
#进入docker虚拟机，注意自己安装的docker版本
#################################################

#设置工作目录
Set-Location D:\Docs
Set-Location D:\Docs


#搜索镜像
docker search omicsclass
#下载基因家族分析镜像
#docker pull omicsclass/gene-family:v2.0

#docker desktop windows 中进入 命令Docker Desktop命令进入虚拟机（直接挂载目录到虚拟机）
#-v D:/gene-family:/work 是设置windows中的D:/gene-family映射到虚拟机中的/work 目录，注意中间用“”隔开，并且必须用完整的绝对路径。相当于共享目录设置
#docker run -m 10G --cpus 2 --rm -v D:/gene-family:/work -it omicsclass/gene-family:v2.0
#docker run -m 10G --cpus 2 --rm -v F:/生信修炼手册/生信高手系统修炼一课通/课件及资料/基因家族分析-课程资料/gene_family_demo/gene_family_demo:/work -it omicsclass/gene-family:v2.0
#docker linux 中进入 命令
#docker run -m 10G  --rm -v ~/gene_family/:/work -it omicsclass/gene-family:v2.0
docker run -m 60G  --rm -v /u3/hekun/zuxuedajiangtang_work/gene_family/:/work -it omicsclass/gene-family:v2.0

移除镜像docker rmi -f omicsclass/gene-family:v2.0
Untagged: omicsclass/gene-family:v2.0
Untagged: docker.io/omicsclass/gene-family@sha256:a7d0a7877ebbb1178c5a48e2f421f2af8fd6634f5597f578b67517c2801b6354
########################################################################
#设置环境变量
########################################################################


cd work

mkdir gene_family#分析时只把这个地址改一下
workdir=/work/PF00319_M
scriptsdir=$workdir/scripts

genome=Taestivum.dna.fa.gz
gff=Taestivum.gff3.gz
species=Taestivum
########################################################################
#数据准备
########################################################################

cd $workdir #回到工作路径
mkdir 01.data_prepare
cd 01.data_prepare


#下载小麦基因组信息
wget -c 

#设置基因组文件名字变量
genome=Taestivum.dna.fa.gz
gff=Taestivum.gff3.gz
species=Taestivum
#conda install -c bioconda agat
#保留蛋白编码基因
agat_sp_filter_feature_by_attribute_value.pl --gff  $gff --attribute biotype --value protein_coding -t '!' -o $species.protein_coding.gff3
#保留最长转录本
agat_sp_keep_longest_isoform.pl --gff $species.protein_coding.gff3 -o $species.longest_isoform.gff3
#Ensembl 下载的基因组 会在基因和mRNA ID前面加gene: 和 transcript:这里给去掉
sed -i 's/gene://' $species.longest_isoform.gff3
sed -i 's/transcript://'  $species.longest_isoform.gff3
#提取cds序列
/share/work/biosoft/TransDecoder/latest/util/gff3_file_to_proteins.pl  --gff3  $species.longest_isoform.gff3 --fasta $genome  --seqType CDS >$species.cds.fa
#提取pep序列
/share/work/biosoft/TransDecoder/latest/util/gff3_file_to_proteins.pl  --gff3  $species.longest_isoform.gff3  --fasta  $genome --seqType prot >$species.pep.fa
#提取cds 用基因的ID，作为序列的ID
perl $scriptsdir/get_gene_longest_fa.pl -l 0 --gff $species.longest_isoform.gff3 --fa $species.cds.fa -p $species.gene.cds
#提取pep 用基因的ID，作为序列的ID
perl $scriptsdir/get_gene_longest_fa.pl -l 0 --gff $species.longest_isoform.gff3 --fa $species.pep.fa -p $species.gene.pep




########################################################################
#hmmer搜索结构域
########################################################################


cd $workdir #回到工作路径
mkdir 02.hmmsearch
cd 02.hmmsearch


#下载 hmm模型：https://www.ebi.ac.uk/interpro/entry/pfam/PF03106/curation/

wget -c https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF03106?annotation=hmm -O Target.hmm.gz
gunzip Target.hmm.gz

hmmsearch --domtblout Target_round1_hmmerOut.txt --cut_tc Target.hmm ../01.data_prepare/Taestivum.gene.pep.fasta

##################################################################
#第二次搜索结构域  可选分析

#提取结构域位置，最后的evalue值根据实际情况可调,多个结构域
grep -v '#' Target_round1_hmmerOut.txt|awk 'BEGIN{OFS="\t"}$10==1 && $7 <1.2e-22{print $1,$18,$19 }' >domain1.bed
grep -v '#' Target_round1_hmmerOut.txt|awk 'BEGIN{OFS="\t"}$10==2 && $7 <1.2e-22{print $1,$18,$19 }' >domain2.bed

#截取 结构域序列
seqtk subseq  ../01.data_prepare/Taestivum.gene.pep.fasta domain1.bed > Target_domain1.fa
seqtk subseq  ../01.data_prepare/Taestivum.gene.pep.fasta domain2.bed > Target_domain2.fa

#合并
cat Target_domain1.fa Target_domain2.fa  >Target_domain.fa

###########以下部分为建立物种特异模型再次搜索，可根据自己基因家族情况选做这部分内容#############################
#clusterW或者muscle比对
muscle -align Target_domain.fa -output Target_domain_aligned.fasta

#格式转换，aln
trimal  -in Target_domain_aligned.fasta -out Target_domain_aligned.aln -clustal

#利用比对结果建立物种特异hmm模型
hmmbuild Target_domain_new.hmm Target_domain_aligned.aln

#新建物种特异hmm模型，再次搜索

hmmsearch --domtblout Target_round2_hmmerOut.txt  Target_domain_new.hmm ../01.data_prepare/Taestivum.gene.pep.fasta

#########################################

#hmmerOut文件，筛选EValue  <0.001 hmm搜索结果，可以用excel手动筛选

#如果只想用hmmer搜索一次，可将下面的文件：Target_round2_hmmerOut.txt 替换成第一次hmmer搜索生成的文件：Target_round1_hmmerOut.txt
#grep -v "^#" Target_round2_hmmerOut.txt|awk '$7<0.001 {print}' >Target_hmmerOut_cutEValue.txt

#跳过二次搜索
grep -v "^#" Target_round2_hmmerOut.txt|awk '$7<0.1 {print}' >Target_hmmerOut_cutEValue.txt

#获得初步基因家族基因列表
cut -f 1 -d " " Target_hmmerOut_cutEValue.txt |sort|uniq >Target_raw_IDlist.txt
#这一步我们应该也可以再blast获取基因ID并于 hmm搜索的基因id文件Target_raw_IDlist.txt合并，然后在提取fasta序列并用于数据库检验，接下来就跟着hmm的分析走就可以
blast
#利用脚本得到对应基因的蛋白序列：
#perl $scriptsdir/get_fa_by_id.pl Target_raw_IDlist.txt ../01.data_prepare/Taestivum.gene.pep.fasta Target_pep_need_to_confirm.fa

seqkit grep -f Target_raw_IDlist.txt ../01.data_prepare/Taestivum.gene.pep.fasta -o Target_pep_need_to_confirm.fa

##################################################
#结构域在线数据库中确认
#################################################

echo "需要手动到各个数据库确认结构域，删除一些假阳性基因，得到可靠的基因列表：Target_IDlist_final.txt"




#将上面Target_pep_need_to_confirm.fa文件中的蛋白序列，再手动验证一下，把不需要的ID删除，最终确认的ID存成新文件：Target_IDlist_final.txt
#手动确认结构域，CDD，SMART，PFAM
CDD:https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi
SMART:https://smart.embl.de/smart/batch.pl
PFAM:https://www.ebi.ac.uk/interpro/     #下载tsv格式文件

#三大数据库网站，筛选之后去除一些不确定的基因ID，最终得到可靠的基因家族基因列表,存储在文件：Target_IDlist_final.txt ; 
#这里将Target_raw_IDlist.txt 中的在数据库中不符合的基因删除，得到Target_IDlist_final.txt
Target_IDlist_final.txt


#利用脚本得到对应基因的蛋白序列：
seqkit grep -f Target_IDlist_final.txt ../01.data_prepare/Taestivum.gene.pep.fasta -o Target_pep_final.fa
#确定分子量大小：http://web.expasy.org/protparam/
#引用可以写bioperl
perl $scriptsdir/stat_protein_fa.pl Target_pep_final.fa Target_pep_final.MW.txt

##################################################################
#筛选整理数据，用于后续分析；脚本提取hmm结果文件，重新筛选一下hmm的结果,提取结构域序列，蛋白全长，cds全长，用于后续分析
#################################################################
#grep读取-f 文件的第一列ID，然后到输入文件中去筛选对应ID的数据输出到第三个文件中
grep -f Target_IDlist_final.txt Target_hmmerOut_cutEValue.txt >Target_hmmerOut_final.txt

#截取得到序列上的保守结构域序列，注意多个结构域分开提取
#结构域位置
grep -v '#' Target_hmmerOut_final.txt|awk 'BEGIN{OFS="\t"}$10==1 {print $1,$18,$19 }' >domain1_final.bed
grep -v '#' Target_hmmerOut_final.txt|awk 'BEGIN{OFS="\t"}$10==2 {print $1,$18,$19 }' >domain2_final.bed
#截取序列
seqtk subseq  ../01.data_prepare/Taestivum.gene.pep.fasta domain1_final.bed > Target_domain1_final.fa
seqtk subseq  ../01.data_prepare/Taestivum.gene.pep.fasta domain2_final.bed > Target_domain2_final.fa


#序列ID换回基因ID，这里如果一个基因存在两个相同结构域，则会有相同的基因ID重复，如果要合并，可以在基因ID末尾加-1和-2做区分，做MADS一个基因只有一个M结构域，则Target_domain1_geneID_final.fa直接就是Target_domain_final.fa，所以不用加-1和-2

#sed  's/:/-1 /' Target_domain1_final.fa >Target_domain1_geneID_final.fa
#sed  's/:/-2 /' Target_domain2_final.fa >Target_domain2_geneID_final.fa
sed  's/:/ /' Target_domain1_final.fa >Target_domain1_geneID_final.fa
sed  's/:/ /' Target_domain2_final.fa >Target_domain2_geneID_final.fa
#合并，注意序列ID是否有重复,这里如果一个基因存在两个相同结构域，则会有相同的基因ID重复，如果要合并，可以在基因ID末尾加-1和-2做区分，做MADS一个基因只有一个M结构域，则Target_domain1_geneID_final.fa直接就是Target_domain_final.fa
#cat Target_domain1_geneID_final.fa Target_domain2_geneID_final.fa  >Target_domain_final.fa
cat Target_domain1_geneID_final.fa >Target_domain_final.fa

#得到对应基因的蛋白序列全长，脚本会读取第一个文件的第一列的ID信息，根据ID提取相应的序列：
seqkit grep -f Target_IDlist_final.txt  ../01.data_prepare/Taestivum.gene.pep.fasta -o  Target_pep_final.fa
sed -i 's/*$//' Target_pep_final.fa #终止密码子*删除


#得到对应基因的cds序列，脚本会读取第一个文件的第一列的ID信息，根据ID提取相应的序列：
seqkit grep -f Target_IDlist_final.txt  ../01.data_prepare/Taestivum.gene.cds.fasta -o  Target_cds_final.fa












########################################################################
#Blast 搜索鉴定基因家族
########################################################################

#blastp比对寻找基因家族成员，Target部分

#NCBI上搜索Target蛋白序列protein数据库：搜索条件：(mads[title] NOT putative[title] AND plants[filter]) AND "Arabidopsis thaliana"[porgn]。搜索条件：(mads[title] NOT putative[title] AND plants[filter]) AND "Oryza sativa"[porgn]
# 然后点击send to，下载fasta序列
#参考文献：https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4955-8

cd $workdir #回到工作路径
#回到工作路径  my_gene_family
mkdir 16.blast_search #创建文件夹不要再docker环境下，因为在root下创建，自己没办法编辑读写
cd 16.blast_search

#blast比对首先建库  #mads蛋白质序列库
makeblastdb -in Target_NCBI_pep.fasta -dbtype prot -title Target_NCBI_pep.fasta   

#blastp比对,将整个小麦蛋白组序列分别于自己所建的MADS结构域库进行比对
nohup blastp -query  ../01.data_prepare/Taestivum.gene.pep.fasta -db Target_NCBI_pep.fasta -out ncbi_Target_blast.out -outfmt 6 -evalue 1e-10  -num_threads 10 &

###这一步我们应该先获取基因ID并于 hmm搜索的基因id文件Target_raw_IDlist.txt合并，然后在提取fasta序列并用于数据库检验，接下来就跟着hmm的分析走就可以


#获得比对上的候选基因，存储在Target.fa文件中，get_fa_by_id.pl脚本可以读取ncbi_Target_blast.out第一列并在Taestivum.gene.pep.fasta中搜索，并得到fasta序列
#perl $scriptsdir/get_fa_by_id.pl ncbi_Target_blast.out ../01.data_prepare/Taestivum.gene.pep.fasta Target.fa

#然后将 Target.fa 提交到 NCBI CDD  SMART 上确认，把不存在结构域的基因ID可以删除


























########################################################################
#进化树分析
########################################################################

cd $workdir #回到工作路径
mkdir 03.gene_tree_analysis
cd 03.gene_tree_analysis
cp ../02.hmmsearch/Target_domain_final.fa .
cp ../02.hmmsearch/Target_pep_final.fa .
cp ../02.hmmsearch/Target_cds_final.fa .
cp ../02.hmmsearch/Target_hmmerOut_final.txt .

##############################################
##构建进化树所用的序列有两种，可以用结构域序列或者蛋白质全长序列
##构建进化树有两种方法最大似然法和NJ
#1 结构域序列多序列比对使用muscle
muscle -align Target_domain_final.fa -output Target_domain_final_aligned.fasta


#最大似然法构建进化树,-m MFP自动查找最优的氨基酸替换模型，-T 10十个线程,-B,即Bootstrap，树的可靠性检验次数1000
iqtree2 -s Target_domain_final_aligned.fasta -m MFP -B 1000  --bnni -T 10 --prefix  iqtree2.pep.domain
#绘图
Rscript $scriptsdir/ggtree.r -t iqtree2.pep.domain.treefile  -p iqtree2.pep.domain

################################################
#2 蛋白质全长多序列比对
nohup muscle -align Target_pep_final.fa -output Target_pep_final_aligned.fasta &

#剪切掉不保守的部分（手动的话看视频课）
trimal -in Target_pep_final_aligned.fasta -out Target_pep_final_aligned_trimed.fasta -automated1

#（手动的话看视频课）手动设置剪切：剪切掉不保守的部分 http://trimal.cgenomics.org/_media/manual.b.pdf 命令行帮助：http://trimal.cgenomics.org/getting_started_with_trimal_v1.2
#trimal  -in Target_pep_final_aligned.fasta -out Target_pep_final_aligned_trimed.fasta  -htmlout Target_pep_final_aligned_trimed_result.html -gt 0.8 -st 0.001 -cons 60

#最大似然法构建进化树
nohup iqtree2 -s Target_pep_final_aligned_trimed.fasta  -m MFP -B 1000  --bnni -T 10 --prefix  iqtree2.pep.fullLength &

Rscript $scriptsdir/ggtree.r -t iqtree2.pep.fullLength.treefile  -p iqtree2.pep.fullLength

#NJ树，采用 mega(看视频课)



########################################################################
#利用meme软件做motif分析，其实motif可以简单地理解为保守的结构域，只不过有些motif研究的多，功能重要，给他起了个名字，比如mads
########################################################################


cd $workdir #回到工作路径
mkdir 04.meme_motif_analysis
cd 04.meme_motif_analysis
#搜索结构域：
#-nmotifs 15  搜索motif的总个数
#-minw 6   motif的最短长度
#-maxw 100   motif的最大长度（如果有文章报告过motif特征，也可以调整）
#-p 10线程数10

cp ../02.hmmsearch/Target_pep_final.fa .

meme Target_pep_final.fa -protein -oc ./ -p 10 \
  -nostatus -time 18000 -maxsize 6000000 -mod anr \
  -nmotifs 15 -minw 6 -maxw 100 &
 
#有个html网页文件，打开看详情
#然后也可以是用TBtools展示，看pdf
########################################################################
#基因结构分析structure
########################################################################


cd $workdir #回到工作路径
#回到工作路径 my_gene_family

mkdir 05.gene_structure_analysis
cd 05.gene_structure_analysis
cp ../02.hmmsearch/Target_IDlist_final.txt .

#获得基因的在染色体上的外显子，CDS，UTR位置信息，用于绘制基因结构图
#注意脚本读取的是第一个（-in1）文件第一列信息，里面是转录本ID；然后把转录本对应的位置cds utr等信息筛选出来
perl $scriptsdir/get_gene_exon_from_gff.pl -in1 Target_IDlist_final.txt -in2 ../01.data_prepare/Taestivum.longest_isoform.gff3 -out gene_exon_info.gff
#然后使用GSDS绘制，。看视频课以及教程：https://www.omicsclass.com/article/63
#然后也可以进化树+motif+基因结构三个一起显示: https://www.omicsclass.com/article/1269 


######################################
这一步没用R代码，请看视频和教程
蛋白序列多序列比对图展示genedoc绘图参考:  https://www.omicsclass.com/article/502



########################################################################
#基因定位到染色体
########################################################################

cd $workdir #回到工作路径
mkdir 06.map_to_chr
cd 06.map_to_chr
cp ../02.hmmsearch/Target_IDlist_final.txt .    #Target基因家族ID列表文件


#获得基因的在染色体上的位置信息，用于绘制染色体位置图,注意提取的是基因位置还是mRNA位置,以下代码是提取的 mRNA位置
perl $scriptsdir/get_gene_weizhi.pl -in1 Target_IDlist_final.txt -in2 ../01.data_prepare/Taestivum.longest_isoform.gff3 -out gene_location.txt

#获得基因组染色体长度：
seqkit fx2tab --length --name  ../01.data_prepare/Taestivum.dna.toplevel.fa.gz|awk '{print $1,$NF}' >genome.len

然后使用Mapchart绘图，参考：https://www.omicsclass.com/article/397，看视频




########################################################################
#基因上游顺势作用原件分析#
########################################################################

#回到工作路径 my_gene_family
cd $workdir #回到工作路径
mkdir 07.gene_promoter
cd 07.gene_promoter
cp ../02.hmmsearch/Target_IDlist_final.txt .

#得到基因在染色体上的位置，此脚本会把基因组所有的序列读入内存，如果基因组较大，可能因为内存不足使脚本运行不成功，可以分染色体分开分析：
perl $scriptsdir/get_gene_weizhi.pl -in1 Target_IDlist_final.txt -in2 ../01.data_prepare/Taestivum.longest_isoform.gff3 -out gene_location.txt
#根据位置信息提取，promoter序列 1500
perl $scriptsdir/get_promoter.pl ../01.data_prepare/Taestivum.dna.toplevel.fa.gz gene_location.txt promoter.fa 2000



########################################################################
#共线性分析，基因加倍与串联重复分析MCScanX
########################################################################

#回到工作路径  my_gene_family
cd $workdir #回到工作路径
mkdir 08.MCScanX
cd 08.MCScanX

cp ../02.hmmsearch/Target_IDlist_final.txt .
#获取基因对应的蛋白序列信息，注意：基因的第一个转录本或者最长的转录本为代表序列；

ln -s ../01.data_prepare/Taestivum.gene.pep.fasta pep.fa

#blast建库
makeblastdb -in pep.fa  -dbtype prot -title pep.fa

#blast+  比对，所有基因对所有基因，-max_target_seqs 5 一个基因最多比对目标数，e只也可以调节，这一步比较久
nohup blastp -query  pep.fa -db pep.fa -out Taestivum.blast -max_target_seqs 10 -outfmt 6 -evalue 1e-10  -num_threads 10 &

#从gff中提取MCScanX需要的基因位置信息，生成文件并不是典型的gff文件，而是符合MCScanX的输入要求的gff
awk -F '\t|;' '$3=="gene" {print $1"\t"$9"\t"$4"\t"$5}' ../01.data_prepare/Taestivum.longest_isoform.gff3| sed 's/ID=//' > Taestivum.gff


#运行MCSCANX  查找共线性，Taestivum是指以Taestivum为前缀搜索blast和gff文件并运行
MCScanX -e 1e-10 Taestivum

#查看染色体800是指图片的分辨率1,2,3,4,5值得染色体号，可以由cut -f 1 Taestivum.gff|sort|uniq|tr '\n' "," 得到
cut -f 1 Taestivum.gff|sort|uniq|tr '\n' "," #查看物种染色体

echo "800
1,2,3,4,5" >Taestivum.circle.ctl

#生成基因家族基因列表，一行
cat Target_IDlist_final.txt|tr "\n" "\t" |sed -r 's/^/Target\t/g;s/\t$//g;' >Target_family.txt

#绘图
java family_circle_plotter -g Taestivum.gff -s Taestivum.collinearity -c  Taestivum.circle.ctl -f Target_family.txt  -o Target.circle.PNG

#找到我们分析的基因家族的共线性分析关系（大片段复制基因对）：
perl /share/work/biosoft/MCScanX/downstream_analyses/detect_collinearity_within_gene_families.pl -i Target_family.txt -d Taestivum.collinearity -o Target.collinear.pairs

#从MCSCAN分析的结果文件：Taestivum.tandem 提取串联重复基因对；串联重复基因对介绍：http://www.omicsclass.com/article/399

perl $scriptsdir/get_tandem_gene.pl -id Target_IDlist_final.txt -tandem Taestivum.tandem -od ./ -name Target



########################################################################
#基因加倍分析绘图，circos
########################################################################


cd $workdir #回到工作路径
mkdir 09.circos
cd 09.circos


cp ../08.MCScanX/Taestivum.collinearity .   
cp ../08.MCScanX/Target.collinear.pairs .


## 获得染色体长度
seqkit fx2tab --length --name   ../01.data_prepare/Taestivum.dna.toplevel.fa.gz|awk '{print $1,$NF}' >genome.len

## 生成染色体文件 7列，  karyotype.txt 文件行的顺序决定了染色体顺序，可以手动调整
cat genome.len|sort -k1,1| awk '{print "chr\t-\t"$1"\t"$1"\t0\t"$2"\tchr"NR}'  > karyotype.txt

#准备circos绘图数据文件，脚本从gff里面获得位置信息并整理出circos需要的link文件
perl $scriptsdir/colline_v3.pl -gff ../08.MCScanX/Taestivum.gff -list Target.collinear.pairs -colline Taestivum.collinearity  -name Target

#绘图，主要准备config.txt配置文件，以及染色体长度文件等等。
#chromosomes_units=1000000单位刻度，越小换上刻度越密集，拟南芥基因组小可以为10的六次方
echo "chromosomes_units=10000000    #刻度单位Mb
karyotype=karyotype.txt  #染色体信息配置文件
show_tick_labels=yes
show_ticks=yes
spacing=10u
<ticks>  #设置染色体刻度
    color=black
    format=%d
    multiplier=1e-6
    radius=1r
    thickness=2p
    <tick>
        size=10p
        spacing=5u
    </tick>
    <tick>
        color=black
        format=%d
        label_offset=10p
        label_size=25p
        show_label=yes
        size=15p
        spacing=25u
        thickness=4p
    </tick>
</ticks>

<ideogram> #染色体绘制设置
    fill=yes   #是否填充颜色
    label_font=default
    label_parallel=yes
    label_radius=dims(image,radius)-60p
    label_size=45
    radius=0.6r  #设置半径，以免基因名称过长超出显示范围
    show_label=yes
    <spacing>
        default=0.005r
    </spacing>
    stroke_color=dgrey
    stroke_thickness=2p
    thickness=0.03r
</ideogram>

<links>  #设置连线
    bezier_radius=0r
    bezier_radius_purity=0.75
    color=black
    crest=0.5
    <link>  #基因家族共线性
        bezier_radius=0r
        bezier_radius_purity=0.75
        color=set2-8-qual-1
        crest=0.5
        file=./Target.link.txt    #输入文件
        radius=0.98r
        <rules>
            <rule>
                color=red
                condition=var(intrachr)
            </rule>
            <rule>
                color=red
                condition=var(interchr)
            </rule>
        </rules>
        thickness=8
        z=20   #层数
    </link>
    <link>  #全基因组共线性
        bezier_radius=0r
        bezier_radius_purity=0.75
        color=230,230,230,0.2
        crest=0.5
        ribbon=yes
        file=./genome.blocklink.txt   #输入文件
        radius=0.98r
        thickness=1
        z=15
    </link>
    radius=0.40r
    thickness=1
</links>
<plots>
    <plot>
       color=black
        label_snuggle=yes  #如多个文本文字距离过近，避免重叠  参考：https://www.omicsclass.com/article/678
        file=./Target.text.txt  #输入文件
        label_font=condensed
        label_size=24p
        link_color=red
        #link_dims=10p,10p,50p,20p,10p
        link_dims=2p,2p,30p,2p,50p #设置方法 https://www.omicsclass.com/article/678
        link_thickness=2p
        r0=1r
        r1=1r+500p
        rpadding=0p
        padding=0p
        show_links=yes
        type=text
    </plot>
    type=histogram
</plots>


<colors>
<<include etc/colors.conf>>
<<include etc/brewer.conf>>
#<<include etc/colors_fonts_patterns.conf>>
#<<include colors.ucsc.conf>>
#<<include colors.hsv.conf>>
</colors>

<fonts>
<<include etc/fonts.conf>>
</fonts>

<image>
<<include etc/image.conf>>
</image>
<<include etc/housekeeping.conf>>
">config.txt 


circos -conf config.txt -outputdir ./ -outputfile Target



########################################################################
#复制基因kaks分析
########################################################################

cd $workdir #回到工作路径
mkdir 10.gene_duplication_kaks_batch
cd 10.gene_duplication_kaks_batch
cp ../02.hmmsearch/Target_IDlist_final.txt .
cp ../01.data_prepare/Taestivum.gene.cds.fasta .
cp ../01.data_prepare/Taestivum.gene.pep.fasta .


#串联重复与大片段重复基因对都可以计算
cp  ../08.MCScanX/Target.collinear.pairs .
cp ../08.MCScanX/Target.tandem .

#处理成两列ID,一行为一对复制基因，并删除第一行
cat Target.collinear.pairs | sed 's#\t#\n#g' |sed 's#:#\t#g'|sed  '1d' >collinearity.pairs

#合并串联重复基因对和大片段复制基因对
cat Target.tandem  collinearity.pairs >homolog.pairs

#ParaAT批量比对 参考文献：https://doi.org/10.1016/j.bbrc.2012.02.101，-p线程
perl $scriptsdir/ParaAT/ParaAT.pl -h  homolog.pairs  -a Taestivum.gene.pep.fasta  -n Taestivum.gene.cds.fasta -p 3  -o  align_out  -m muscle -f axt 

#合并axt文件，批量计算kaks，-m YN指定模型，YN比较常用。也可以不指定，则会求很多模型的平均值，会慢一些
cat align_out/*.axt >all_merge_align.axt
KaKs_Calculator -m YN -i ./all_merge_align.axt  -o  all_kaks_result.txt   #-m YN 指定一种算法



########################################################################
#复制基因查找blast分析方法 基因组草图，没有到染色体水平推荐
########################################################################

#回到工作路径  my_gene_family
cd $workdir #回到工作路径
mkdir 11.gene_duplication_blast
cd 11.gene_duplication_blast
cp ../02.hmmsearch/Target_IDlist_final.txt .
cp ../02.hmmsearch/Target_cds_final.fa .

#blast建库
makeblastdb -in Target_cds_final.fa  -dbtype nucl -title  Target_cds_final.fa

#获取基因cds长度

samtools faidx Target_cds_final.fa
#blast+  比对，所有基因对所有基因
blastn -query  Target_cds_final.fa -db Target_cds_final.fa -out blast.out  -outfmt 6 -evalue 1e-10  -num_threads 6


#根据比对筛选结果,比对长度和比对率0.7 70
perl $scriptsdir/KAKS_SHAIXUAN.pl  -in1  Target_cds_final.fa.fai  -in2  blast.out -out  dup.txt
 
#注意这里没有筛选位置信息，因此需要自己再看一下基因在染色体上的位置确定是否为串联重复基因


########################################################################
#不同物种之间的共线性分析
########################################################################

# 回到工作路径  my_gene_family
cd $workdir #回到工作路径
mkdir 12.mcscan_between_genome
cd 12.mcscan_between_genome


#准备白菜基因组
wget -c ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/fasta/brassica_rapa/dna/Brassica_rapa.Brapa_1.0.dna.toplevel.fa.gz
wget -c ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/gff3/brassica_rapa/Brassica_rapa.Brapa_1.0.56.gff3.gz


genome=Brassica_rapa.Brapa_1.0.dna.toplevel.fa.gz
gff=Brassica_rapa.Brapa_1.0.56.gff3.gz
species=Brassica_rapa

#保留蛋白编码基因
agat_sp_filter_feature_by_attribute_value.pl --gff  $gff --attribute biotype --value protein_coding -t '!' -o $species.protein_coding.gff3
#保留最长转录本
agat_sp_keep_longest_isoform.pl --gff $species.protein_coding.gff3 -o $species.longest_isoform.gff3
#Ensembl 下载的基因组 会在基因和mRNA ID前面加gene: 和 transcript:这里给去掉
sed -i 's/gene://' $species.longest_isoform.gff3
sed -i 's/transcript://'  $species.longest_isoform.gff3
#提取cds序列
/share/work/biosoft/TransDecoder/latest/util/gff3_file_to_proteins.pl  --gff3  $species.longest_isoform.gff3 --fasta $genome  --seqType CDS >$species.cds.fa
#提取pep序列
/share/work/biosoft/TransDecoder/latest/util/gff3_file_to_proteins.pl  --gff3  $species.longest_isoform.gff3  --fasta  $genome --seqType prot >$species.pep.fa
#提取cds 转换成OrthoFinder 需要的序列格式，用基因的ID，作为序列的ID
perl $scriptsdir/get_gene_longest_fa.pl -l 0 --gff $species.longest_isoform.gff3 --fa $species.cds.fa -p $species.gene.cds
#提取pep 转换成OrthoFinder 需要的序列格式，用基因的ID，作为序列的ID
perl $scriptsdir/get_gene_longest_fa.pl -l 0 --gff $species.longest_isoform.gff3 --fa $species.pep.fa -p $species.gene.pep



#基因位置信息 bed文件准备
python3 -m jcvi.formats.gff bed --type=gene Brassica_rapa.longest_isoform.gff3 -o Brassica_rapa.bed
python3 -m jcvi.formats.gff bed --type=gene ../01.data_prepare/Taestivum.longest_isoform.gff3 -o Taestivum.bed

#基因的cds序列文件准备
ln -s Brassica_rapa.gene.cds.fasta Brassica_rapa.cds
ln -s ../01.data_prepare/Taestivum.gene.cds.fasta Taestivum.cds

#共线性分析
python3 -m jcvi.compara.catalog ortholog Taestivum Brassica_rapa --cscore=0.7 --no_strip_names

#对共线性区域进行过滤
python3 -m jcvi.compara.synteny screen --minsize=0 --minspan=30 --simple Taestivum.Brassica_rapa.anchors   Taestivum.Brassica_rapa.anchors.new

#############################################################
#绘制共线性图片：
#准备两个配置文件： mcscan_seqid  mcscan_layout

#查看染色体名称，然后将查看的结果复制到下列的配置文件中
cat Brassica_rapa.bed|cut -f1|sort|uniq|tr '\n' ','
cat Taestivum.bed|cut -f1|sort|uniq|tr '\n' ','

#写出配置文件
echo "1,2,3,4,5
A01,A02,A03,A04,A05,A06,A07,A08,A09,A10" >mcscan_seqid

echo "# y, xstart, xend, rotation, color, label, va,  bed
 .7,     .2,    .9,       30,     red, ATH, top, Taestivum.bed
 .3,     .2,    .9,       0,      blue, BRA, bottom, Brassica_rapa.bed
# edges
e, 0, 1, Taestivum.Brassica_rapa.anchors.simple" >mcscan_layout

#绘图

python3 -m jcvi.graphics.karyotype  -o Taestivum.Brassica_rapa.karyotype.pdf --format=pdf  --figsize=10x5   mcscan_seqid mcscan_layout


#基因家族基因对颜色添加，可以用excel vlookup




########################################################################
#结合转录组分析基因家族成员表达量绘制热图
########################################################################

cd $workdir #回到工作路径
#回到工作路径  my_gene_family

mkdir 13.rna_seq_heatmap
cd 13.rna_seq_heatmap

cp ../02.hmmsearch/Target_IDlist_final.txt .
然后再Target_IDlist_final.txt的第一行添加上geneID，这样后续就可以直接得到表头
不知道能成不，试一试sed -i '1i geneID\t' Target_IDlist_final.txt

#提取Target家族基因表达量：grep读取-f 文件的第一列ID，然后到输入文件中去筛选对应ID的数据输出到第三个文件中
#数据地址：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121407
grep -f  Target_IDlist_final.txt GSE121407_fpkm.txt >gene_family_fpkm.txt
#添加表头
#sed -i '1i geneID\tWT_1\tWT_2\tWT_4\tH1_1\tH1_2\tH1_4' gene_family_fpkm.txt

#热图绘制详解：https://www.omicsclass.com/article/1162
Rscript $scriptsdir/heatmap.r -i gene_family_fpkm.txt -s -n heatmap
#也可以使用在线工具对gene_family_fpkm.txt进行标准化，报错的话记得删除表达量都为零的基因


########################################################################
#蛋白互作网络分析
########################################################################

cd $workdir #回到工作路径
#回到工作路径  my_gene_family

mkdir 14.PPI_network
cd 14.PPI_network
cp ../02.hmmsearch/Target_pep_final.fa .

#白菜的蛋白序列
seqkit grep -f Target_bra_IDlist.txt ../12.mcscan_between_genome/Brassica_rapa.gene.pep.fasta -o Target_bra_pep.fa



########################################################################
#蛋白3D结构预测
########################################################################


cd $workdir #回到工作路径
#回到工作路径  my_gene_family
mkdir 15.3D_structure
cd 15.3D_structure

#下载PDB数据库中蛋白序列：https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
wget -c https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
#解压
gunzip *gz
#改名
mv pdb_seqres.txt pdb_seqres.fa
#建库
makeblastdb -in pdb_seqres.fa -dbtype prot -title pdb_seqres.fa

#PSI-blast 比对，
psiblast -query ../02.hmmsearch/Target_pep_final.fa  -db pdb_seqres.fa -out psi.txt -num_iterations 3 -out_ascii_pssm psi.pssm  -outfmt 6  -max_target_seqs 1

python3 $scriptsdir/get_pdb_file.py 6ir8

