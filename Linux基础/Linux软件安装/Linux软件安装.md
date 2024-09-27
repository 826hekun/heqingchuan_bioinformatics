## conda 的安装和使用

### conda的安装

下载安装包 -- bash 安装 -- 接受协议 -- 选择默认安装路径（回车） -- 重新激活环境 -- 调用帮助文档

```sh
## 下载，其实不需要，只是为了让大家了解一下，服务器已经下载好了，直接cp或者软链接即可
# wget -c https://mirrors.bfsu.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh

## 软链接即可
cd  ~
ln  -s  /home/t_linux/Miniconda3-latest-Linux-x86_64.sh  ./
## 安装，安装过程只需要输入 yes 或者按 Enter
bash  Miniconda3-latest-Linux-x86_64.sh 

## 重新激活环境
source  ~/.bashrc 

## 查看 conda 的帮助文档
conda --help

```

### 配置镜像

我们使用 conda 安装软件时，conda 会去 channel 中搜索软件，如果使用的服务器是在国内，channel 就选择国内的，推荐清华，如果清华镜像出问题，再选择其他。

```shell
## 配置镜像

# 下面四行配置北京外国语大学的conda的channel地址（首选）
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/main/ 
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge/ 
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda/ 
conda config --set show_channel_urls yes 

# 下面这四行配置清华大学的conda的channel地址（首选北外，如果体验不好再换成清华）
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes

# 如果需要官方频道，可以添加下面这两行配置官网的channel地址（不推荐）
conda config --add channels conda-forge 
conda config --add channels bioconda

# 删除defaults频道
sed -i '/defaults/d' ~/.condarc

## 配置镜像成功
# 查看配置结果
cat ~/.condarc
```



### 创建小环境

```sh
# 创建名为rna的软件环境来安装转录组学分析的生物信息学软件
conda create -y -n  rna  python=3.7
# 创建小环境成功，并成功安装python3版本
# 每建立一个小环境，安装一个python=3的软件作为依赖

# 查看当前conda环境
conda info -e
conda env list

# 每次运行前，激活创建的小环境rna
conda activate rna

# 退出小环境
conda deactivate
```

### 在小环境中安装生信软件

注：软件都要安装在小环境中，不要安装在 base 

```sh
# 激活环境
conda activate rna
# 安装 fastqc 软件
conda  install  fastqc

# 调出帮助文档
fastqc --help

# 可以指定软件版本
conda install -y samtools=1.14 

# aspera 
conda install -y -c hcc aspera-cli
ascp --help

# 可以一次安装多个软件
conda install -y python=3.7 libstdcxx-ng=9.1.0 trim-galore  hisat2  subread  multiqc  samtools=1.14  salmon=1.4.0 fastp fastqc
# mamba install -y python=3.7 libstdcxx-ng=9.1.0 trim-galore  hisat2  subread  multiqc  samtools=1.14  salmon=1.4.0 fastp fastqc

## 不是通过软件名来调用帮助文档，而是软件的命令
# sra-tools
prefetch --help
fastq-dump --help
which prefetch

#  trim-galore
trim_galore --help

# hisat2
hisat2 --help

# subread
featureCounts --help

# multiqc
multiqc --help

# samtools
samtools --help

# salmon
salmon --help

# fastp
fastp --help

# 
```

### R4环境

```bash
# 创建R4环境
conda create -y -n R4 python=3.8

# 激活R4环境
conda activate R4

# (可选步骤：在R4里安装mamba)
# conda install mamba

# 安装R语言本体
conda install -y r-base=4.1.2
## 或者使用mamba安装： mamba install -y r-base=4.1.2

# 安装R语言软件包
conda install -y r-getopt r-tidyverse r-ggplot2=3.3.5 bioconductor-limma bioconductor-edger bioconductor-deseq2 bioconductor-clusterprofiler bioconductor-org.hs.eg.db=3.13.0
## 或者使用mamba安装：mamba install -y r-getopt r-tidyverse r-ggplot2=3.3.5 bioconductor-limma bioconductor-edger bioconductor-deseq2 bioconductor-clusterprofiler bioconductor-org.hs.eg.db=3.13.0
```

Q: 如何验证R语言的包安装情况？

A: 进入R语言环境中用`library()` 

激活R：

```bash
# 输入R进入R语言的交互
R
```

在R语言里验证安装包：

```R
library(getopt)
library(tidyverse)
library(ggplot2)
library(limma)
library(edgeR)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
```

Q: 如何知道这些包到底叫啥？哪儿该大写哪儿该小写？

A: [Bioconductor - Home](https://bioconductor.org/) 在Bioconductor的官网搜索即可。

### 直接导入配置文件以安装软件

下载钉钉群里的yaml文件

``` bash
conda env create -n rna -f rna.yaml
# 如果有mamba的话可以用mamba安装
# mamba env create -n rna -f rna.yaml
conda env create -n R4 -f R4.yaml
```





### 其他用法

```r
卸更新软件：conda  update  软件名
载软件：conda  remove  软件名
删除环境：conda  remove  -n   环境名
克隆环境：conda  create –n 新环境名 –clone  旧环境名
查找软件：conda  search  软件名
```

查找软件常用的链接：

- https://anaconda.org/search

- [https://bioconda.github.io](https://bioconda.github.io/)[/](https://bioconda.github.io/)

## 环境变量 PATH

### 常见的环境变量

`$HOME` 记录了用户的家目录所在的路径

`PS1` 命令行配色

```shell
$ echo  $HOME
/trainee2/vip28

$ echo  $PS1
\[\033]2;\h:\u \w\007\033[33;1m\]\u \033[35;1m\t\033[0m \[\033[36;1m\]\w\[\033[0m\]\n\[\e[32;1m\]$ \[\e[0m\]

$ echo  $PATH
/trainee2/vip28/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin

```

### 修改命令行配色

感兴趣的自行搜索

```shell
echo  'export PS1="\[\033]2;\h:\u \w\007\033[33;1m\]\u \033[35;1m\t\033[0m \[\033[36;1m\]\w\[\033[0m\]\n\[\e[32;1m\]$ \[\e[0m\]" ' >> ~/.bashrc
source  ~/.bashrc
```

`~/.bashrc`：系统配置文件，包含专用于你的 bash shell 的bash信息、设置，每次登录或打开新的 shell 时，该文件会被自动读取和执行。



### PATH

`$PATH`：输入命令时Linux会去查找PATH里面记录的路径，如果命令存在某一个路径中，就可以成功调用。

`<PATH1>:<PATH2>:<PATH3>:------:<PATHN>`

打个比方，PATH 是一个工具箱，有很多层（对应很多个路径），每一层放着各式各样的工具（对应各种命令）。

```shell
$ echo $PATH
/trainee2/vip28/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin

# 可以把 : 替换成换行符 \n 
$ echo $PATH | tr ':'  '\n'
/trainee2/vip28/miniconda3/condabin
/usr/local/sbin
/usr/local/bin
/usr/sbin
/usr/bin
/sbin
/bin
/usr/games
/usr/local/games
/snap/bin

# 比如 ls 命令存在
$ ls  
$ which ls 
/bin/ls
```



### 如何管理 PATH

如何管理 `$PATH`：理解环境变量 `$PATH`  是非常重要的，对后续的环境和软件管理都非常重要。

推荐方法：在自己家目录下创建一个 `~/bin/` 文件夹并将其添加到环境变量，后续安装软件，就将软件的可执行文件拷贝或软链接（绝对路径）到这个 bin 文件夹：

```shell
mkdir  ~/bin 
echo  'export "PATH=~/bin:$PATH" ' >> ~/.bashrc 
source  ~/.bashrc
```



## 其他软件安装

```shell
mkdir ~/biosoft 
cd ~/biosoft
# wget -c https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download -O hisat2-2.2.1-Linux_x86_64.zip
ln  -s  /teach/software/hisat2-2.2.1-Linux_x86_64.zip  ./
unzip hisat2-2.2.1-Linux_x86_64.zip
cd hisat2-2.2.1/
./hisat2 --help
# echo 'export PATH="${HOME}/biosoft/hisat2-2.2.1/:$PATH" ' >> ~/.bashrc 
ln  -s  ~/biosoft/hisat2-2.2.1/hisat2*   ~/bin/
```



































