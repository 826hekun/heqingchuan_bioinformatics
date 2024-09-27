######################################################################################
#下载基因组GWAS分析docker镜像
#docker pull omicsclass/pop-evol-gwas:v2.0  #本shell适用pop-evol-gwas:v2.0的镜像，以此为准 (视频有时更新不及时)  
#启动docker容器并交互式进入
#docker desktop方法
#docker run --rm -it -m 4G --name gwas   -v D:\gwas:/work omicsclass/pop-evol-gwas:v2.0
#docker linux方法
#docker run --rm -it -m 4G --name gwas  -v ~/gwas:/work omicsclass/pop-evol-gwas:v2.0

####################################################
# 创建目录，以及准备一些路径变量，方便后续调用
####################################################
#docker run --rm -it -v /share/work/huangls/piplines/omicsclass/pop-evol-gwas/scripts/:/work/scripts -v /share/nas1/huangls/test/gwas_rice:/work omicsclass/pop-evol-gwas:v2.0

workdir=/work/my_gwas  #设置工作路径
refdir=$workdir/ref
datadir=$workdir/data
scriptdir=$workdir/scripts
tmpdir=$workdir/tmp
export TMPDIR=$tmpdir      # 设置临时目录 防止容器中默认临时目录写满报错
export PATH=$scriptdir:$PATH    #组学大讲堂提供的脚本 添加到环境变量中方便使用

#一些关键文件所在的路径变量
GFF=$refdir/all.gff3
REF=$refdir/all.con.fa
FAI=$refdir/all.con.fa.fai
GROUP=$datadir/pop_group.txt

#水稻GWAS数据来源：https://www.nature.com/articles/ng.3596

#分组：根据自己的分组名称进行拆分

#cat $GROUP |grep GROUP1 |cut -f 1 >$datadir/pop1.txt
#cat $GROUP |grep GROUP2 |cut -f 1 >$datadir/pop2.txt

####################################################################
#数据过滤
#####################################################################
cd $workdir  #回到工作目录
mkdir 00.filter
cd 00.filter

#对数据进行过滤
#
#过滤1：vcfutils.pl  过滤掉indel附近的snp
# -w INT    SNP within INT bp around a gap to be filtered [3]
# -W INT    window size for filtering adjacent gaps [10]
#
vcfutils.pl varFilter -w 5 -W 10 "gzip -dc $datadir/rice.raw.vcf.gz|" |gzip - >all.varFilter.vcf.gz
#
#过滤2：vcftools
#--max-missing Exclude sites on the basis of the proportion of missing data 
#(defined to be between 0 and 1, where 0 allows sites that are completely missing 
#and 1 indicates no missing data allowed).

vcftools --gzvcf all.varFilter.vcf.gz --recode --recode-INFO-all --stdout     --maf 0.05  --max-missing 0.9  --minDP 2  --maxDP 1000      \
    --minQ 30 --minGQ 0 --min-alleles 2  --max-alleles 2 --remove-indels |gzip - > clean.vcf.gz  

#转换成hapmap格式，方便后续tassel和Gapit读取
perl $scriptdir/vcf2hmp.pl clean.vcf.gz hapmap.txt.gz

#设置变量方便后续命令调用
vcf=$workdir/00.filter/clean.vcf.gz
hmp=$workdir/00.filter/hapmap.txt.gz
##########################################################
#进化树分析
###########################################################
cd $workdir  #回到工作目录
mkdir 01.phylo_tree
cd 01.phylo_tree
#文件格式转换

run_pipeline.pl -Xmx30G  -SortGenotypeFilePlugin -inputFile $vcf  \
    -outputFile clean.sorted.vcf.gz -fileType VCF

run_pipeline.pl  -Xmx5G -importGuess  clean.sorted.vcf.gz  \
    -ExportPlugin -saveAs supergene.phy -format Phylip_Inter

#最大似然法构建进化树
#方法1：fasttree 构建进化树
fasttree -nt -gtr  supergene.phy   >  fasttree.nwk

#进化树绘制
ggtree.r -t fasttree.nwk -f $GROUP -g Group -n fasttree
#进化树手动美化：https://www.omicsclass.com/article/671

########################################################
#PCA分析
#########################################################
cd $workdir  #回到工作目录
mkdir 02.PCA
cd 02.PCA
## plink分析PCA
plink --vcf  $vcf --pca 10 --out  plink_pca   \
    --allow-extra-chr --set-missing-var-ids @:#    --vcf-half-call missing
#绘图
pca_plink_plot.r -i plink_pca.eigenvec -f $GROUP -g Group --name plink_pca


#########################################################
#群体结构STRUCTURE分析
######################################################
cd $workdir  #回到工作目录
mkdir 03.STRUCTURE
cd 03.STRUCTURE
## filter LD 
#50 10 0.2   50个SNP 窗口  step 10个  r2 0.2  ; 50 5 0.4
plink --vcf  $vcf  --indep-pairwise 50 10 0.2 --out ld   \
    --allow-extra-chr --set-missing-var-ids @:# 
plink --vcf  $vcf  --make-bed --extract ld.prune.in  \
    --out LDfiltered --recode vcf-iid  --keep-allele-order  --allow-extra-chr --set-missing-var-ids @:#  


#转换成admixture要求的ped格式
plink --noweb --vcf LDfiltered.vcf  --recode12 --out admixture \
     --allow-extra-chr  --keep-allele-order

#admixture 群体结构分析
for k in {2..10};do
    echo "admixture -j8 -C 0.01 --cv admixture.ped $k >admixture.log$k.out"
    admixture -j8 -C 0.01 --cv admixture.ped $k >admixture.log$k.out
done
#绘图展示 
structure_plot.r  -d ./ -s admixture.nosex 
#确定最佳K，CV值最小时对应的K值为最佳K
grep "CV error" *out
#
##########################################################
#连锁不平衡分析 LDdecay
########################################################
#
cd $workdir  #回到工作目录
mkdir 04.LDdecay
cd 04.LDdecay

PopLDdecay -InVCF  ../00.filter/all.varFilter.vcf.gz \
        -SubPop  popid.txt -MaxDist 500 -OutStat ld.stat

#绘图输出
Plot_OnePop.pl -inFile ld.stat.gz -output ld

#Plot_OnePop.pl -inFile ld.bin.gz -output ld.bin   #平滑处理

###################################################################
#性状数据分析
###################################################################
cd $workdir  #回到工作目录
mkdir 05.trait
cd 05.trait

#注意，下方两个R脚本输入格式要求为tassel格式：https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load

Rscript $scriptdir/trait_plot.r  -i $datadir/traits_gapit.txt 

#如果有多年或者多点的数据可以计算BLUP值再做GWAS分析：https://www.omicsclass.com/article/1380
Rscript $scriptdir/blup_year_or_site.r  -i heading_days.rm_outers.txt

#查看 性状数据分布，去掉一些极端值样本，或者标记为NA

###################################################################
#GWAS分析  数据来源：doi:10.1038/ng.3596
###################################################################
#基因型填充
cd $workdir  #回到工作目录
mkdir 06.impute
cd 06.impute


#v2.0 镜像更新，2023.11.6 使用以下命令

beagle  -Xmx4g -Djava.io.tmpdir=$TMPDIR gt=$vcf   out=all.impute impute=true window=10 nthreads=2


###################################################################
#利用tassel进行GWAS分析
###################################################################

cd $workdir  #回到工作目录
mkdir 07.gwas_tassel
cd 07.gwas_tassel

#对变异结果文件排序
#注意：tassel 做关联分析，可能会报排序错误，因此需要预先排序；
#tassel 排序之后染色体编号会改变，需要注意
run_pipeline.pl -Xmx30G  -SortGenotypeFilePlugin -inputFile  $hmp \
    -outputFile hapmap.ordered -fileType Hapmap


#计算PCA和kinship, 也可以使用前面的群体遗传进化分析的结果
run_pipeline.pl -Xmx30G -fork1 -h hapmap.ordered.hmp.txt \
    -PrincipalComponentsPlugin  -ncomponents 3  -covariance true -endPlugin \
    -export pca -runfork1
    
run_pipeline.pl -Xmx30G -h hapmap.ordered.hmp.txt \
    -KinshipPlugin -method Centered_IBS -endPlugin -export kinship.txt -exportType SqrMatrix 

#数据路径变量设置
tassel_trait=$workdir/05.trait/heading_days_BLUP_value.tsv  #多年数据，使用BLUP值
kinship=$workdir/07.gwas_tassel/kinship.txt
structure=$workdir/07.gwas_tassel/pca1.txt

#GLM模型
run_pipeline.pl -Xmx30G -fork1 -h hapmap.ordered.hmp.txt \
    -fork2 -importGuess $tassel_trait   \
   -combine3 -input1 -input2  -intersect -FixedEffectLMPlugin \
   -endPlugin -export hd_GLM

#GLM模型，加入 群体结构，Q矩阵或者PCA  
#也可以用admixture计算的Q矩阵，但是需要去掉最后一列  -excludeLastTrait  可以自动去除

run_pipeline.pl -Xmx30G -fork1 -h hapmap.ordered.hmp.txt \
    -fork2 -t $tassel_trait -fork3 -q $structure -excludeLastTrait \
   -combine5 -input1 -input2 -input3 -intersect -FixedEffectLMPlugin \
   -endPlugin -export hd_GLM_structure

#可视化
cut -f 2,6 hd_GLM1.txt|sed 's/:/\t/' >glm_pvalue.txt
sed -i '1s/Marker/chr\tpos/' glm_pvalue.txt
Rscript $scriptdir/gwas_manhattan_plot.r -i glm_pvalue.txt -F $FAI -T heading_days -n heading_days_glm_manhattan -c 1.3e-6
Rscript $scriptdir/qq_plot.r  -i glm_pvalue.txt -n heading_days_glm_qq

cut -f 2,6 hd_GLM_structure1.txt|sed 's/:/\t/'>glm_structure_pvalue.txt
sed -i '1s/Marker/chr\tpos/' glm_structure_pvalue.txt
Rscript $scriptdir/gwas_manhattan_plot.r -i glm_structure_pvalue.txt -F $FAI -T heading_days -n heading_days_glm_structure_manhattan -c 1.3e-6
Rscript $scriptdir/qq_plot.r  -i glm_structure_pvalue.txt -n heading_days_glm_structure_qq

#MLM模型  加入随机效应亲缘关系  混合线性模型
#单个性状数据关联 
run_pipeline.pl -Xmx30g -fork1 -h hapmap.ordered.hmp.txt \
    -fork2 -r $tassel_trait -fork3 -q $structure -excludeLastTrait \
    -fork4 -k $kinship -combine5 -input1 -input2 -input3 \
    -intersect -combine6 -input5 -input4 -mlm -mlmVarCompEst P3D \
    -mlmCompressionLevel None  -export hd_MLM -runfork1 -runfork2 -runfork3 -runfork4

#结果可视化
cut -f 2,7 hd_MLM2.txt|sed 's/:/\t/'|grep -v "NaN" >mlm_pvalue.txt
sed -i '1s/Marker/chr\tpos/' mlm_pvalue.txt
Rscript $scriptdir/gwas_manhattan_plot.r -i mlm_pvalue.txt -F $FAI -T heading_days -n heading_days_mlm_manhattan -c 1.3e-6
Rscript $scriptdir/qq_plot.r  -i mlm_pvalue.txt -n heading_days_mlm_qq


#只用亲缘关系矩阵的MLM
#run_pipeline.pl -Xmx30g -fork1 -h hapmap.ordered.hmp.txt \
#    -fork2 -r $tassel_trait   \
#    -fork3 -k $kinship -combine4 -input1 -input2  \
#    -intersect -combine5 -input4 -input3 -mlm -mlmVarCompEst P3D \
#    -mlmCompressionLevel None  -export hd_MLM_onlyK -runfork1 -runfork2 -runfork3 
#
#结果可视化
#cut -f 2,7 hd_MLM_onlyK2.txt|sed 's/:/\t/'|grep -v "NaN" >mlm_onlyK_pvalue.txt
#sed -i '1s/Marker/chr\tpos/' mlm_onlyK_pvalue.txt
#Rscript $scriptdir/gwas_manhattan_plot.r -i mlm_onlyK_pvalue.txt -F $FAI -T heading_days -n heading_days_mlm_onlyK_manhattan -c 1.3e-6
#Rscript $scriptdir/qq_plot.r  -i mlm_onlyK_pvalue.txt -n heading_days_mlm_onlyK_qq


#CMLM 压缩混合线性模型
run_pipeline.pl -Xmx30g -fork1 -h hapmap.ordered.hmp.txt \
    -fork2 -r $tassel_trait -fork3 -q $structure -excludeLastTrait \
    -fork4 -k $kinship -combine5 -input1 -input2 -input3 \
    -intersect -combine6 -input5 -input4 -mlm -mlmVarCompEst P3D \
    -mlmCompressionLevel Optimum  -export hd_CMLM -runfork1 -runfork2 -runfork3 -runfork4

#结果可视化
cut -f 2,7 hd_CMLM2.txt|sed 's/:/\t/'|grep -v "NaN" >cmlm_pvalue.txt
sed -i '1s/Marker/chr\tpos/' cmlm_pvalue.txt
Rscript $scriptdir/gwas_manhattan_plot.r -i cmlm_pvalue.txt -F $FAI -T heading_days -n heading_days_cmlm_manhattan -c 1.3e-6
Rscript $scriptdir/qq_plot.r  -i cmlm_pvalue.txt -n heading_days_cmlm_qq


###################################################################
#利用gapit进行GWAS分析
###################################################################
cd $workdir  #回到工作目录
mkdir 08.gwas_gapit
cd 08.gwas_gapit


#关联分析 , GAPIT 会自动根据相应模型计算 PCA与亲缘关系矩阵
for m in GLM MLM MLMLM CMLM ECMLM SUPER FarmCPU Blink FaST EMMA EMMAx;do
    echo "Rscript $scriptdir/gapit_gwas.r -i $hmp -t $workdir/05.trait/heading_days_BLUP_value.tsv  -m $m  -n hd.$m.GWAS.Results"
    Rscript $scriptdir/gapit_gwas.r -i $hmp -t $workdir/05.trait/heading_days_BLUP_value.tsv  -m $m  -n hd.$m.GWAS.Results
done


#可视化绘图
for m in GLM MLM MLMLM CMLM ECMLM SUPER FarmCPU Blink FaST EMMA EMMAx;do

    cut  -f 1,4 hd.$m.GWAS.Results.txt |sed 's/:/\t/'>pvalue_$m.txt
    sed -i '1s/SNP/Chr\tpos/' pvalue_$m.txt
    Rscript $scriptdir/gwas_manhattan_plot.r -i pvalue_$m.txt -F $FAI -T heading_days_$m -n heading_days_${m}_manhattan -c 1.3e-6
    Rscript $scriptdir/qq_plot.r  -i pvalue_${m}.txt -n heading_days_${m}_qq
done


## GAPIT： GLM模型不添加群体结构
Rscript $scriptdir/gapit_gwas.r -i $hmp -t $workdir/05.trait/heading_days_BLUP_value.tsv  -m GLM  --PCA.total 0 -n hd.GLM0.GWAS.Results
cut  -f 1,4 hd.GLM0.GWAS.Results.txt  |sed 's/:/\t/'>pvalue_GLM0.txt
sed -i '1s/SNP/Chr\tpos/' pvalue_GLM0.txt
Rscript $scriptdir/gwas_manhattan_plot.r -i pvalue_GLM0.txt -F $FAI -T heading_days_GLM0 -n heading_days_GLM0_manhattan -c 1.3e-6
Rscript $scriptdir/qq_plot.r  -i pvalue_GLM0.txt -n heading_days_GLM0_qq



###################################################################
#利用emmax进行GWAS分析  官方文档： https://genome.sph.umich.edu/wiki/EMMAX
###################################################################
cd $workdir  #回到工作目录
mkdir 09.gwas_emmax
cd 09.gwas_emmax

##  vcf  转换格式
#--maf 0.05 --geno 0.1  

plink --vcf $vcf \
     --recode12  --output-missing-genotype 0 \
    --transpose --out snp_var   --set-missing-var-ids @:#  --allow-extra-chr

##性状数据和SNP基因型样本保持一致  并转换成emmax识别格式
perl $scriptdir/sort_trait.pl snp_var.tfam $workdir/05.trait/heading_days_BLUP_value.tsv |awk '{print $1,$2,$NF}' >emmax_heading_days_BLUP_value.tsv

#kinship 亲缘关系矩阵计算 ibs-kinship
emmax-kin-intel64 -v  -s -d 10  -o ./pop.kinship.IBS  snp_var

#kinship 亲缘关系矩阵计算 BN-kinship
emmax-kin-intel64 -v  -d 10  -o ./pop.kinship.BN  snp_var

#pca 群体结构  前4轴
awk  '{print $1" "$2" 1 "$3" "$4" "$5" "$6}' ../02.PCA/plink_pca.eigenvec > pca.cov

#关联分析 emmax模型, pca作为协变量
emmax-intel64 -v -d 10 -t snp_var  \
    -p emmax_heading_days_BLUP_value.tsv \
    -k pop.kinship.IBS -c pca.cov  -o emmax.out 



#整理数据，可视化
#第一列ID，第二列回归系数(beta),第三列回归系数的标准差，第四列为P值。
awk -F"[:\t]" 'BEGIN{print "CHR\tPOS\tP"}{print $1"\t"$2"\t"$NF}' emmax.out.ps > emmax.pvalue
Rscript $scriptdir/gwas_manhattan_plot.r -i emmax.pvalue -F $FAI -T heading_days -n heading_days_emmax_manhattan -c 1.3e-6
Rscript $scriptdir/qq_plot.r  -i emmax.pvalue -n heading_days_emmax_qq


###########################emmax +Q################################## 
#PCA也可也用 Q矩阵代替
#注意Q矩阵要去掉最后一列，Q矩阵作为协变量
paste ../03.STRUCTURE/LDfiltered.nosex  ../03.STRUCTURE/admixture.5.Q |awk '{print $1" "$2" 1 "$3" "$4" "$5" "$6}'  > pop.Qmatrix
emmax-intel64 -v -d 10 -t snp_var  -p emmax_heading_days_BLUP_value.tsv \
     -k pop.kinship.IBS -c pop.Qmatrix   -o emmax.out.Q 
awk -F"[:\t]" 'BEGIN{print "CHR\tPOS\tP"}{print $1"\t"$2"\t"$NF}' emmax.out.Q.ps  > emmaxQ.pvalue
Rscript $scriptdir/gwas_manhattan_plot.r -i emmaxQ.pvalue -F $FAI -T heading_days -n heading_days_emmaxQ_manhattan -c 1.3e-6
Rscript $scriptdir/qq_plot.r  -i emmaxQ.pvalue -n heading_days_emmaxQ_qq




###################################################################
#利用GEMMA进行GWAS分析:  

# 参考：https://github.com/frankvogt/vcf2gwas/blob/main/MANUAL.md
#http://www.xzlab.org/software/GEMMAmanual.pdf
###################################################################
cd $workdir  #回到工作目录
mkdir 10.gwas_gemma
cd 10.gwas_gemma

#基因型数据准备：格式转换

plink --vcf $vcf   --make-bed  --out snp_var   \
    --set-missing-var-ids @:#  --allow-extra-chr

#表型数据准备：
#因为GEMMA是通过识别plink格式的fam文件中的性状来进行关联分析的，所以我们需要预处理一下，
#把性状数据添加进去。添加到fam文件的6，7，8...等等列。
#性状排序并合并到*fam文件中

perl $scriptdir/sort_trait.pl snp_var.fam $workdir/05.trait/heading_days_BLUP_value.tsv  |cut -d " " -f 1-5,7 >snp_var.fam1
mv -f snp_var.fam1 snp_var.fam

#获得kinship矩阵，使用混合线性模型进行分析。
gemma -bfile snp_var -gk 1 -o kin  -miss 1  -outdir ./

################################################
#GEMMA关联分析
###############################################
#lmm 混合线性模型
#-n 1  fam文件中第几个性状做关联分析
# -miss  1    屏蔽缺失过滤
# -outdir  指定输出目录
# -o 指定输出文件名前缀

####lm 线性模型  会自动计算PCA，不用提供PCA也行
gemma -bfile snp_var   -lm 4 -o hd_lm -n 1 -miss 1 -outdir ./

#也可提供PCA和上一步结果一样
gemma -bfile snp_var  -d  ../02.PCA/plink_pca.eigenval \
    -u ../02.PCA/plink_pca.eigenvec  \
    -lm 4 -miss 1  -n 1  -o hd_lm_pca -outdir ./

#或者Q矩阵，注意这里我们计算Q的vcf和snp_var是相同的vcf因此样本顺序是一致的，可也直接使用；
cat ../03.STRUCTURE/admixture.4.Q|cut -d " " -f 1-3 >4.Q
gemma -bfile snp_var  -c 4.Q -n 1 -lm 4 -miss 1    -o hd_lm_Q  -outdir ./


####lmm混合线性模型 K
gemma -bfile snp_var  -k kin.cXX.txt -lmm 4 -miss 1 -n 1    -o hd_lmm_k -outdir ./

###Q+K矩阵
gemma -bfile snp_var -c 4.Q  -k kin.cXX.txt -lmm 4 -miss 1 -n 1    -o hd_lmm_QK -outdir ./

#其中：
#chr：SNP所在染色体号
#rs： SNP名称
#ps： SNP物理位置
#n_miss： SNP缺失个体数
#allele1： 次等位基因
#allele0： 主等位基因
#af：SNP频率
#beta： SNP效应值
#se： beta估计标准误
#l_remle： 计算该SNP效应时对应的lamda的remle估计值。
#p_wald ：wald检验P值
#其中，我们最关心的三个结果是chr, ps, p_wald，我们可以借助这三个结果画曼哈顿图和QQ图。l_remle比较难理解，需要懂模型才知道它的含义，但对分析来说，不是很重要。

################
cat hd_lm.assoc.txt |awk -F "[:\t]" 'BEGIN{OFS="\t"}{if(NR>1){print $2,$3,$12}}'>pvalue_lm.txt
sed -i '1ichr\tpos\tp' pvalue_lm.txt
Rscript $scriptdir/gwas_manhattan_plot.r -i pvalue_lm.txt -F $FAI -T heading_days -n heading_days_lm_manhattan -c 1.3e-6
Rscript $scriptdir/qq_plot.r  -i pvalue_lm.txt -n heading_days_lm_qq


cat hd_lmm_k.assoc.txt |awk -F "[\t:]" 'BEGIN{OFS="\t"}{if(NR>1){print $2,$3,$14}}'>pvalue_lmm_K.txt
sed -i '1ichr\tpos\tp' pvalue_lmm_K.txt
Rscript $scriptdir/gwas_manhattan_plot.r -i pvalue_lmm_K.txt -F $FAI -T heading_days -n heading_days_lmm_K_manhattan -c 1.3e-6
Rscript $scriptdir/qq_plot.r  -i pvalue_lmm_K.txt -n heading_days_lmm_qq

cat hd_lmm_QK.assoc.txt |awk -F "[\t:]" 'BEGIN{OFS="\t"}{if(NR>1){print $2,$3,$14}}'>pvalue_lmm_QK.txt
sed -i '1ichr\tpos\tp' pvalue_lmm_QK.txt
Rscript $scriptdir/gwas_manhattan_plot.r -i pvalue_lmm_QK.txt -F $FAI -T heading_days -n heading_days_lmm_QK_manhattan -c 1.3e-6
Rscript $scriptdir/qq_plot.r  -i pvalue_lmm_QK.txt -n heading_days_lmm_QK_qq



#############################################
##获得关联区域里面的基因
################################################
cd  $workdir
mkdir 11.get_gene
cd 11.get_gene

#1.根据P值获取关联区域

ln -s ../07.gwas_tassel/glm_pvalue.txt  pvalue.txt

#根据阈值筛选关联的SNP，并生成bed文件，排序
awk 'BEGIN{OFS="\t"}NR>1 && $3<1e-7{print $1,$2,$2,$3}' pvalue.txt |sort -k1,1 -k2,2n >pvalue.bed

#生成显著性P值附近10k区域
bedtools flank  -i pvalue.bed  -b 10000 -g $FAI  |sort -k1,1 -k2,2n>region.bed

#合并P值 防止，SNP 所在位置也合并 2023.04.10 更新
cat pvalue.bed  region.bed |sort -k1,1 -k2,2n >all.region.bed

#合并区域
bedtools merge -i all.region.bed >region_merged.bed

#2.获取基因的位置信息
#gff中获取所有基因位置,注意染色体ID保持一致
awk '$3=="gene"' $GFF > gene.gff

#3.筛选关联区域内基因
bedtools  intersect -wo -nonamecheck -F 0.1  -a  region_merged.bed -b gene.gff | awk -F "\t"  '{print $4"\t"$7"\t"$8"\t"$10"\t"$12}' | sort -u  > gene.txt

#4.筛选关联区域SNP变异信息
bedtools  intersect -wo -nonamecheck -F 0.1  -a  region_merged.bed -b $vcf >var.txt

##############################################################
#感兴趣的区域LD热图绘制
LDBlockShow  -InVCF $vcf -OutPut LD -Region Chr1:36000000-37000000 \
    -OutPng -SeleVar 2   

#关联区域添加上 Pvalue值  ,更多绘图参数参照：https://github.com/bgi-shenzhen/ldblockshow

awk '$1=="Chr1" && $2>36000000 && $2<37000000{print}' pvalue.txt >ld.pvalue.txt
LDBlockShow  -Cutline  7 -InVCF $vcf -OutPut LD.p \
    -Region Chr1:36000000-37000000 -OutPng -SeleVar 2 -InGWAS ld.pvalue.txt




