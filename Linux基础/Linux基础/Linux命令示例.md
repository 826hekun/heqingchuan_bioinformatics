### Mac电脑登录前需要先设置

对于 **Mac pro** 来说，经常出现一种情况是 ：登录没问题，但是超过 **5 分钟不操作**就会出现问题，没法输入，不得不重启 终端 或者 iTerm 。解决方法是：

【1】在mac，打开终端，不要登录服务器

【2】然后在本地运行下面命令

```sh
cat  >~/.ssh/config
Host *
    ServerAliveInterval 120
    TCPKeepAlive no
^C
```



### 登录服务器

打开邮箱找到曾老师给大家发的邮件，里面有用户名、密码和ip地址，登录方式为：`ssh   用户名@ip地址`，如：

```sh
ssh  vip28@94.191.82.93
```

回车，然后输入密码

### 修改命令行配色

复制粘贴下面两行代码：

```sh
echo  'export PS1="\[\033]2;\h:\u \w\007\033[33;1m\]\u \033[35;1m\t\033[0m \[\033[36;1m\]\w\[\033[0m\]\n\[\e[32;1m\]$ \[\e[0m\]"' >> ~/.bashrc
source  ~/.bashrc
```

### 查看帮助文档

man 命令，help 命令，或者某个命令的  --help  参数

```sh
man  ls		## 用 man 命令查看 ls 命令的帮助文档
help  ls	## 用 help 命令查看 ls 命令的帮助文档	
ls  --help	## 用 --help 参数查看 ls 命令的帮助文档
```

### ls 命令

列出目录文件情况：

```sh
ls				## 列出当前目录的文件
ls  ./			## 同上，‘.’号代表当前目录
ls  ./*txt		## 列出当前目录下以 txt 结尾的文件
ls  ../ 		## 列出上层目录的文件
ls  -a			## 列出当前目录下的所有文件，包括隐藏文件
ls  -l			## 列出当前目录下文件的详细信息
ll				## ls  -la 的简写
ls  -lh 		## 加上 -h 参数，以 K、M、G 的形式显示文件大小
ls  -lh  /		## 列出根目录下文件的详细信息
```

### cd 命令

切换工作目录

```sh
cd  ..       ## 切换到上层目录，相对路径
cd  /        ## 切换到根目录
cd  /teach/  ## 切换到根目录下的teach，绝对路径
cd  -        ## 返回上一次的工作目录
cd  ~        ## 回到用户家目录
cd           ## 同上，回到用户家目录
```

### mkdir

```sh
# 创建目录
mkdir dir0
ls
mkdir dir0/sub1/sub2
ls
ls dir0
mkdir -p dir0/sub1/sub2
ls dir0
ls dir0/sub1/
mkdir -p  test{1..3}/test{1..3}
tree
```

### touch

```sh
ls
touch  file.txt  new.txt
ls
touch  file{1..5}
ls
```

### rm

```
rm  -i  file.txt
ls  file*
rm  file*
rm  -r  test1
```

### mv

```sh
mv  file1   Data/file2
```

### cp

```sh
cp   readme.txt   Data/
mkdir  dir0
cp  -r  dir0  Data/

```

### ln

```sh
ln -s /teach/software/Miniconda3-latest-Linux-x86_64.sh  ./

```

### tar

```sh
## 解压
tar  -zxvf  Data.tar.gz
## 压缩
tar  -zcvf  Data.tar.gz    Data  ...

```

### cat

```sh
cat  readme.txt
cat  -n  readme.txt
## 写入文件
cat >file
Welcome to Biotrainee() !
^C			## 这里是按Crtl  C
## 查看
cat file
Welcome to Biotrainee() !

```

### head、tail

```sh
head  -n  20  Data/example.fq
## 查看 .bashrc 的最后 10 行
tail  ~/.bashrc
## 查看第20行
head  -n  20  Data/example.fq | tail -1

```

### less

按 q 退出

```sh
less  Data/example.fq
less -S Data/example.fq
less -N Data/example.fq
zless -N Data/reads.1.fq.gz

```

### wc

```sh
cat -n readme.txt
cat readme.txt | wc 
wc -l readme.txt
```

### cut

```sh
less -S Data/example.gtf | cut -f 1,3-5
less -S Data/example.gtf | cut -d 'h' -f 1
```

### sort

```sh
less -S Data/example.gtf | sort -k 4 | less -S
less -S Data/example.gtf | sort -n -k 4 | less -S

```

### uniq

```sh
less -S Data/example.gtf | cut -f 3 | sort | uniq -c

```

### paste

```sh
less -S Data/example.fq | paste - - - | less -S
paste file1 file2
```

### tr

```sh
cat readme.txt | tr 'e' 'E'
cat readme.txt | tr '\n' '\t'
cat readme.txt | tr -d 'e' 

```

### grep

```sh
grep Biotrainee -r ./
less -S Data/example.fq | grep 'gene'
less -S Data/example.fq | grep -w 'gene'
less -S Data/example.fq | grep -v -w 'gene'

```

### 正则表达式

```sh
cat readme.txt  | grep '^T'
cat readme.txt  | grep ')$'
cat readme.txt  | grep 'f.ee'
cat readme.txt  | grep 'f\?ee'
cat readme.txt  | grep 're\+'
cat readme.txt  | grep [bB]

```



### sed

```sh
cat readme.txt | sed '1i Welcome to Biotrainee() '
cat readme.txt | sed '1a Welcome to Biotrainee() '
cat readme.txt | sed '1c Welcome to Biotrainee() 
cat readme.txt | sed 's/is/IS/g'
cat readme.txt | sed '/^$/d'
cat readme.txt | sed  'y/abc/ABC/'
```

### awk

```sh
less -S  Data/example.gtf | awk '{print $9}' | less -S
less -S  Data/example.gtf | awk '{print $9,$10}' | less -S
less -S  Data/example.gtf | awk -F '\t' '{print $9}' | less -S
less -S  Data/example.gtf | awk '{if($3=="gene") print $0}' | less -S
less -S  Data/example.gtf | awk '{if($3=="gene") {print $0} else{print $3 " is not gene "}}' | less -S
less -S Data/example.gtf | awk '/gene/{print $0}' | less -S
less -S Data/example.gtf | awk 'BEGIN{print "find UTR feature"} /UTR/{print $0} END{print "end"}'
less -S Data/example.gtf | awk 'BEGIN{FS="\t"} {print $9}' | less -S
less -S Data/example.gtf | awk 'BEGIN{FS="\t";OFS="\t"} {gsub("gene","Gene",$3);print $0}' | less -S

```

