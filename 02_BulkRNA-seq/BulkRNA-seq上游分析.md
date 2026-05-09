# BulkRNA-seq上游分析

## 0.准备阶段

0.1下载数据

https://ngdc.cncb.ac.cn/

从网站上输入项目找到数据（下载的数据为双端测序的数据）

项目号：PRJCA034153

数据名称：CRR1435510

0.2参考基因组

[Ensembl genome browser 115](http://asia.ensembl.org/index.html?)

下载基因组和注释基因组 （罗非鱼参考序列）

0.3BulkRNA分析环境

```java
BulkRNA_seq环境配置
conda create -n BulkRNA_seq
conda activate BulkRNA_seq
conda install -c bioconda fastp

mamba install python=3.8
#multiqc不能在python<3.8下运行
pip install --upgrade pip
pip install multiqc -i https://pypi.tuna.tsinghua.edu.cn/simple some-package

#安装featurecounts
mamba install -c bioconda subread
#安装tpmcalculator
mamba install -c bioconda tpmcalculator
```

## 1.质控

```java
生成list.txt
cd /home/xuanwenxin/xwx_data/re_jgg/1.rawdata
ls *.gz|sed 's/_R[12]_001\.fastq\.gz//'| sort -u > /home/xuanwenxin/xwx_data/re_jgg/list.txt 
```

```java
fastp

cd /home/xuanwenxin/xwx_data/re_jgg/1.rawdata
less /home/xuanwenxin/xwx_data/re_jgg/list.txt | while read id;
do
echo ${id}
fastp -i ${id}_R1_001.fastq.gz \
-o /home/xuanwenxin/xwx_data/re_jgg/2.cleandata/${id}_1_clean.fq.gz \
-I ${id}_R2_001.fastq.gz \
-O /home/xuanwenxin/xwx_data/re_jgg/2.cleandata/${id}_2_clean.fq.gz \
--thread=16
done
```

```java
fastqc

cd /home/xuanwenxin/xwx_data/re_jgg/2.cleandata
mkdir ./fastq

less /home/xuanwenxin/xwx_data/re_jgg/list.txt|while read id;
do
echo ${id}
fastqc ${id}_1_clean.fq.gz -o ./fastq
fastqc ${id}_2_clean.fq.gz -o ./fastq
done

mkdir ./multiqc
multiqc ./fastq -o ./multiqc -n test
```

## 2.比对

```java
cd /home/xuanwenxin/xwx_data/re_jgg/Tilapia_gene

hisat2-build -p 10 Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.fa Tilapia
#建立索引这步一定要去到索引文件夹中

mkdir /home/xuanwenxin/xwx_data/re_jgg/3.bam
less /home/xuanwenxin/xwx_data/re_jgg/list.txt | while read id;
do
echo ${id}
hisat2 -p 50 -x /home/xuanwenxin/xwx_data/re_jgg/Tilapia_gene/Tilapia \
-1 /home/xuanwenxin/xwx_data/re_jgg/2.cleandata/${id}_1_clean.fq.gz \
-2 /home/xuanwenxin/xwx_data/re_jgg/2.cleandata/${id}_2_clean.fq.gz \
-S /home/xuanwenxin/xwx_data/re_jgg/3.bam/${id}.sam > ${id}.hisat2.log 2>&1
done
```

## 3.SAM转BAM

```java
cd /home/xuanwenxin/xwx_data/re_jgg
mkdir ./4.bam/sort
bam=/home/xuanwenxin/xwx_data/re_jgg/4.bam
cd /home/xuanwenxin/xwx_data/re_jgg/3.bam
less /home/xuanwenxin/xwx_data/re_jgg/list.txt | while read id;
do
echo ${id}
samtools view -b -S ${id}.sam > ${bam}/${id}.bam
samtools sort ${bam}/${id}.bam -o ${bam}/sort/${id}_sorted.bam
samtools index ${bam}/sort/${id}_sorted.bam
#rm ${id}.bam
#rm ${id}.sam
done
```

## 4.reads定量

```java
mkdir -p ./5.count/
cd /home/xuanwenxin/xwx_data/re_jgg

gtf=/home/xuanwenxin/xwx_data/re_jgg/Tilapia_gene/Oreochromis_niloticus.O_niloticus_UMD_NMBU.115.chr.gtf

less /home/xuanwenxin/xwx_data/re_jgg/list.txt | while read id; 
do     
echo ${id}
     featureCounts \
         -a ${gtf}  \
        -T 16 -p --countReadPairs -t exon -g gene_id \
         -o /home/xuanwenxin/xwx_data/re_jgg/5.count/${id}_counts_id.txt \
         /home/xuanwenxin/xwx_data/re_jgg/4.bam/sort/${id}_sorted.bam
 done

cd ./5.count
less /home/xuanwenxin/xwx_data/re_jgg/list.txt | while read id; 
do
     echo ${id}
     cat ${id}_counts_id.txt | cut -f 1,7 > ${id}_featureCounts.txt
done

paste *featureCounts.txt > countstable
awk '{for(i=3;i<=NF;i+=2) $i=""; print}' countstable > RNA_Seq_countstable.txt

```

5.tpm计算

tpm calculator介绍https://github.com/NLM-DIR/TPMCalculator/

```java
cd /home/xuanwenxin/xwx_data/re_jgg
mkdir ./6.tpm
echo "=== 使用TPMCalculator计算TPM ==="
while read id; do
    echo "计算TPM: ${id}"
    TPMCalculator \
        -g Tilapia_gene/Oreochromis_niloticus.O_niloticus_UMD_NMBU.115.chr.gtf \
        -p \
        -b 4.bam/sort/${id}_sorted.bam \
        -o 6.tpm/${id}
done < list.txt
echo "TPM计算完成！"

#生成TPM最终结果矩阵列表
TPM_results.txt
```