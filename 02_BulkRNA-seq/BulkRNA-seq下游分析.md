1.load RNAdata


```R
#读取txt文件（rawdata）
```


```R
library(dplyr)
```


```R
BulkRNA_df <- read.table("/storage/public/xuanwenxin/xwx_data/re_jgg/5.count/RNA-Seq_Practice_countstable")
rownames(BulkRNA_df) <- BulkRNA_df$V1
BulkRNA_df <- BulkRNA_df[,2:19]
```


```R
colnames(BulkRNA_df) <- scan("/storage/public/xuanwenxin/xwx_data/re_jgg/list.txt", what = "character", sep = "\n")
BulkRNA_df <- BulkRNA_df[-1,]
BulkRNA_df[] <- lapply(BulkRNA_df, as.numeric)
head(BulkRNA_df[,1:18],6)
```


<table class="dataframe">
<caption>A data.frame: 6 × 18</caption>
<thead>
	<tr><th></th><th scope=col>Cold_HK1</th><th scope=col>Cold_HK2</th><th scope=col>Cold_HK3</th><th scope=col>Cold_SP1</th><th scope=col>Cold_SP2</th><th scope=col>Cold_SP3</th><th scope=col>Cold_TH1</th><th scope=col>Cold_TH2</th><th scope=col>Cold_TH3</th><th scope=col>Ctrl_HK1</th><th scope=col>Ctrl_HK2</th><th scope=col>Ctrl_HK3</th><th scope=col>Ctrl_SP1</th><th scope=col>Ctrl_SP2</th><th scope=col>Ctrl_SP3</th><th scope=col>Ctrl_TH1</th><th scope=col>Ctrl_TH2</th><th scope=col>Ctrl_TH3</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ENSONIG00000018423</th><td>37</td><td>42</td><td>36</td><td>304</td><td>364</td><td>281</td><td>21</td><td>22</td><td>40</td><td>16</td><td> 97</td><td>19</td><td>271</td><td>59</td><td>403</td><td>34</td><td> 38</td><td>20</td></tr>
	<tr><th scope=row>ENSONIG00000032772</th><td> 0</td><td> 3</td><td> 1</td><td>  5</td><td>  4</td><td>  3</td><td> 4</td><td> 3</td><td> 5</td><td> 0</td><td>  4</td><td> 0</td><td>  5</td><td> 5</td><td> 15</td><td> 3</td><td>  1</td><td> 1</td></tr>
	<tr><th scope=row>ENSONIG00000007536</th><td> 0</td><td> 1</td><td> 0</td><td>  1</td><td>  1</td><td>  0</td><td> 0</td><td> 0</td><td> 0</td><td> 0</td><td>  0</td><td> 0</td><td>  0</td><td> 0</td><td>  4</td><td> 0</td><td>  0</td><td> 0</td></tr>
	<tr><th scope=row>ENSONIG00000029674</th><td> 0</td><td> 0</td><td> 1</td><td>  0</td><td>  1</td><td>  0</td><td> 1</td><td> 0</td><td> 0</td><td> 0</td><td>  0</td><td> 0</td><td>  0</td><td> 1</td><td>  1</td><td> 4</td><td>  2</td><td> 0</td></tr>
	<tr><th scope=row>ENSONIG00000032019</th><td> 1</td><td> 0</td><td> 1</td><td>  0</td><td>  0</td><td>  2</td><td> 0</td><td> 1</td><td> 0</td><td> 0</td><td>  2</td><td> 0</td><td>  0</td><td> 0</td><td>  4</td><td> 4</td><td>  5</td><td> 0</td></tr>
	<tr><th scope=row>ENSONIG00000033498</th><td>59</td><td>60</td><td>49</td><td> 51</td><td> 75</td><td> 39</td><td>48</td><td>49</td><td>92</td><td>13</td><td>100</td><td>14</td><td> 85</td><td>69</td><td>385</td><td>88</td><td>122</td><td>94</td></tr>
</tbody>
</table>




```R
featureCounts_meta <- BulkRNA_df[,1:18]
BulkRNA_df <- BulkRNA_df[rowSums(BulkRNA_df)> 0,] 
prefix <- "Tilapia"
```


```R
RNA_coldata <- data.frame(id = colnames(BulkRNA_df),
                         group = gsub("[0-9]", "", colnames(BulkRNA_df)),
                         condition = stringr::str_split(colnames(BulkRNA_df),pattern = "_",simplify = TRUE)[,1] )

rownames(RNA_coldata) <- RNA_coldata$id
head(RNA_coldata,2)
```


<table class="dataframe">
<caption>A data.frame: 2 × 3</caption>
<thead>
	<tr><th></th><th scope=col>id</th><th scope=col>group</th><th scope=col>condition</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>Cold_HK1</th><td>Cold_HK1</td><td>Cold_HK</td><td>Cold</td></tr>
	<tr><th scope=row>Cold_HK2</th><td>Cold_HK2</td><td>Cold_HK</td><td>Cold</td></tr>
</tbody>
</table>




```R

```


```R

```


```R
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = BulkRNA_df,
                              colData = RNA_coldata,
                              design = ~ condition)
```

    Loading required package: S4Vectors
    
    Loading required package: stats4
    
    Loading required package: BiocGenerics
    
    
    Attaching package: ‘BiocGenerics’
    
    
    The following objects are masked from ‘package:dplyr’:
    
        combine, intersect, setdiff, union
    
    
    The following objects are masked from ‘package:stats’:
    
        IQR, mad, sd, var, xtabs
    
    
    The following objects are masked from ‘package:base’:
    
        anyDuplicated, aperm, append, as.data.frame, basename, cbind,
        colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
        get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
        match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
        Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
        table, tapply, union, unique, unsplit, which.max, which.min
    
    
    
    Attaching package: ‘S4Vectors’
    
    
    The following object is masked from ‘package:clusterProfiler’:
    
        rename
    
    
    The following objects are masked from ‘package:dplyr’:
    
        first, rename
    
    
    The following object is masked from ‘package:utils’:
    
        findMatches
    
    
    The following objects are masked from ‘package:base’:
    
        expand.grid, I, unname
    
    
    Loading required package: IRanges
    
    
    Attaching package: ‘IRanges’
    
    
    The following object is masked from ‘package:clusterProfiler’:
    
        slice
    
    
    The following objects are masked from ‘package:dplyr’:
    
        collapse, desc, slice
    
    
    Loading required package: GenomicRanges
    
    Loading required package: GenomeInfoDb
    
    Loading required package: SummarizedExperiment
    
    Loading required package: MatrixGenerics
    
    Loading required package: matrixStats
    
    
    Attaching package: ‘matrixStats’
    
    
    The following object is masked from ‘package:dplyr’:
    
        count
    
    
    
    Attaching package: ‘MatrixGenerics’
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
        colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
        colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
        colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
        colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
        colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
        colWeightedMeans, colWeightedMedians, colWeightedSds,
        colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
        rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
        rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
        rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
        rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
        rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
        rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
        rowWeightedSds, rowWeightedVars
    
    
    Loading required package: Biobase
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    
    Attaching package: ‘Biobase’
    
    
    The following object is masked from ‘package:MatrixGenerics’:
    
        rowMedians
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        anyMissing, rowMedians
    
    
    converting counts to integer mode
    
    Warning message in DESeqDataSet(se, design = design, ignoreRank):
    “some variables in design formula are characters, converting to factors”



```R
tissue <- c("HK","SP","TH")
```


```R
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
```


```R
condition_level <- c("Ctrl","Cold")
```


```R

```


```R
DESeq_find_diffgenes <- function(sample_id,Col_data,condition_level){
    print(paste("run",sample_id,sep = ":"))
    sample.id <- grep(sample_id,rownames(Col_data),value=T)
    # For sample_id
    dds.sample <- dds[,sample.id]
    colData(dds.sample) <- droplevels(colData(dds.sample))
    dds.sample$condition <- factor(dds.sample$condition, levels = condition_level )
    design(dds.sample) <- formula(~ condition)
    
    dds.sample <- DESeq(dds.sample)
    #res.sample <- results(dds.sample)
    
    # select the genes which showing differentially expressed in any one of HT time points
    temp <- c("condition_Cold_vs_Ctrl")
    lapply(1:length(temp),function(x){
        res1 <- results(dds.sample, name=temp[x])
        res1.ordered <- as.data.frame(res1[order(res1$pvalue),])
        res1.ordered$gene <- rownames(res1.ordered)
        #res1.ordered$condition <- paste0("HT",x)
        res1.ordered$tissue <- sample_id
        return(res1.ordered)
    }) -> res.sample
    res.sample <- do.call(rbind,res.sample)

    res.list <- list(dds.sample,res.sample)
    return(res.list)
}
```


```R
tissue_deg_list <- lapply(tissue,DESeq_find_diffgenes,Col_data = RNA_coldata,condition_level)
names(tissue_deg_list) <- tissue
```

    [1] "run:HK"


    estimating size factors
    
    estimating dispersions
    
    gene-wise dispersion estimates
    
    mean-dispersion relationship
    
    final dispersion estimates
    
    fitting model and testing
    


    [1] "run:SP"


    estimating size factors
    
    estimating dispersions
    
    gene-wise dispersion estimates
    
    mean-dispersion relationship
    
    final dispersion estimates
    
    fitting model and testing
    


    [1] "run:TH"


    estimating size factors
    
    estimating dispersions
    
    gene-wise dispersion estimates
    
    mean-dispersion relationship
    
    final dispersion estimates
    
    fitting model and testing
    



```R
DEG.RNA <- do.call(rbind,lapply(tissue_deg_list,'[[',2))
DEG.RNA_2 <- na.omit(DEG.RNA) %>% filter(padj < 0.05) %>% group_by(tissue) %>% arrange(desc(log2FoldChange),.by_group = TRUE)
head(DEG.RNA_2)
```


<table class="dataframe">
<caption>A grouped_df: 6 × 8</caption>
<thead>
	<tr><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th><th scope=col>gene</th><th scope=col>tissue</th></tr>
	<tr><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>15.349371</td><td>6.732761</td><td>1.221477</td><td>5.511984</td><td>3.548114e-08</td><td>1.113293e-06</td><td>ENSONIG00000014479</td><td>HK</td></tr>
	<tr><td> 9.273605</td><td>5.997890</td><td>1.269832</td><td>4.723374</td><td>2.319638e-06</td><td>4.484633e-05</td><td>ENSONIG00000042238</td><td>HK</td></tr>
	<tr><td> 8.578649</td><td>5.883817</td><td>1.290609</td><td>4.558946</td><td>5.141091e-06</td><td>8.950031e-05</td><td>ENSONIG00000019199</td><td>HK</td></tr>
	<tr><td>11.521892</td><td>5.544759</td><td>1.148675</td><td>4.827091</td><td>1.385420e-06</td><td>2.839257e-05</td><td>ENSONIG00000011052</td><td>HK</td></tr>
	<tr><td> 4.971530</td><td>5.080199</td><td>1.393481</td><td>3.645688</td><td>2.666774e-04</td><td>2.580710e-03</td><td>ENSONIG00000008590</td><td>HK</td></tr>
	<tr><td> 4.958700</td><td>5.077875</td><td>1.386471</td><td>3.662446</td><td>2.498187e-04</td><td>2.442672e-03</td><td>ENSONIG00000035414</td><td>HK</td></tr>
</tbody>
</table>




```R
DEG.RNA_2.list <- split(DEG.RNA_2$gene,DEG.RNA_2$tissue)
```


```R
geneID_name <- readRDS("/storage/public/xuanwenxin/xwx_data/re_jgg/eggnog_geneID2name.rds")
geneID_name <- geneID_name %>% select(gene = Gene.stable.ID,everything())
geneID_name$eggnog_Preferred_name[is.na(geneID_name$eggnog_Preferred_name)] <- "-"
head(geneID_name,2)
```


<table class="dataframe">
<caption>A data.frame: 2 × 4</caption>
<thead>
	<tr><th></th><th scope=col>gene</th><th scope=col>Gene.name</th><th scope=col>eggnog_Preferred_name</th><th scope=col>eggnog_Description</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>ENSONIG00000000001</td><td>ppfibp1a</td><td>PPFIBP1</td><td>PTPRF interacting protein, binding protein 1 (liprin beta 1)</td></tr>
	<tr><th scope=row>2</th><td>ENSONIG00000000002</td><td>syt10   </td><td>SYT10  </td><td>Synaptotagmin X                                             </td></tr>
</tbody>
</table>




```R
DEG.RNA_3 <- DEG.RNA_2 %>% group_by(gene) %>% left_join(geneID_name,by = "gene") %>% distinct()
```

    Warning message in left_join(., geneID_name, by = "gene"):
    “[1m[22mDetected an unexpected many-to-many relationship between `x` and `y`.
    [36mℹ[39m Row 86 of `x` matches multiple rows in `y`.
    [36mℹ[39m Row 4545 of `y` matches multiple rows in `x`.
    [36mℹ[39m If a many-to-many relationship is expected, set `relationship = "many-to-many"` to silence this warning.”



```R

```

Venn Plot


```R
library(ggvenn)
```

    Loading required package: ggplot2
    



```R
options(repr.plot.width = 6,repr.plot.height = 6)
RNA_venn <- ggvenn(
  data = DEG.RNA_2.list[c("HK","TH","SP")],         # 数据列表
  columns = c("HK","TH","SP"),           # 对选中的列名绘图，最多选择4个，NULL为默认全选
  show_elements = F,        # 当为TRUE时，显示具体的交集情况，而不是交集个数
  label_sep = "\n",         # 当show_elements = T时生效，分隔符 \n 表示的是回车的意思
  show_percentage = T,      # 显示每一组的百分比
  digits = 1,               # 百分比的小数点位数
  fill_color = c("#A5C496", "#C7988C", "#8891DB"), # 填充颜色
  fill_alpha = .6,         # 填充透明度
  stroke_color = "black",   # 边缘颜色
  stroke_alpha = 0.9,       # 边缘透明度
  stroke_size = 0.2,        # 边缘粗细
  stroke_linetype = "solid", # 边缘线条 # 实线：solid  虚线：twodash longdash 点：dotdash dotted dashed  无：blank
  set_name_color = "black", # 组名颜色
  set_name_size = 6,        # 组名大小
  text_color = "black",     # 交集个数颜色
  text_size = 5             # 交集个数文字大小
)
```


```R
RNA_venn
```


    
![png](output_24_0.png)
    


PCA


```R
#stabilize varaince within libraries
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
vsd <- vst(dds)
vsd_filtered <- assay(vsd)
```

    gene-wise dispersion estimates
    
    mean-dispersion relationship
    
    final dispersion estimates
    



```R
options(repr.plot.width = 6,repr.plot.height = 5)
RNA_PCA_plot <- plotPCA(vsd[,grep(colnames(vsd),pattern = "SP",invert = T)], intgroup='group', ntop=1000) + 
    ggsci::scale_color_npg() + 
    theme_classic(base_size = 15) + 
    theme(axis.line = element_blank(),
          panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
```

    using ntop=1000 top features by variance
    
    Warning message:
    “[1m[22m`aes_string()` was deprecated in ggplot2 3.0.0.
    [36mℹ[39m Please use tidy evaluation idioms with `aes()`.
    [36mℹ[39m See also `vignette("ggplot2-in-packages")` for more information.
    [36mℹ[39m The deprecated feature was likely used in the [34mDESeq2[39m package.
      Please report the issue to the authors.”
    Warning message:
    “[1m[22mThe `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    [36mℹ[39m Please use the `linewidth` argument instead.”



```R
temp <- RNA_PCA_plot$data
temp$treat <- stringr::str_split(temp$group,pattern = "_",simplify = T)[,1]
temp$tissue <- stringr::str_split(temp$group,pattern = "_",simplify = T)[,2]
temp$treat <- factor(temp$treat,levels = c("Ctrl","Cold"))
```


```R
options(repr.plot.width = 6,repr.plot.height = 4)
ggplot(temp,aes(x = PC1 , y = PC2 ,color = treat,shape = tissue)) + 
    geom_point(size = 4) + scale_color_manual(values = c("#D6AFB9","#7E9BB7")) +
    labs(x = "PC1: 92% variance",y = "PC2: 4% variance" , color = NULL) +
    theme_classic(base_size = 15) + 
    theme(axis.line = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
```


    
![png](output_29_0.png)
    



```R
options(repr.plot.width = 6,repr.plot.height = 5.5)
plotPCA(vsd[,grep(colnames(vsd),pattern = "HK")], intgroup='group', ntop=1000) + ggsci::scale_color_npg() + theme_bw(base_size = 15)
```

    using ntop=1000 top features by variance
    



    
![png](output_30_1.png)
    



```R
pca_function <- function(data){
    mydata<-as.data.frame(t(data))
    mydata=vegan::decostand(mydata,method = "hellinger")
    otu_pca<- prcomp(mydata,scal=F)
    pc12 <- as.data.frame(otu_pca$x[,1:2])*100
    pc12$samples<-rownames(pc12)
    return(pc12)
                             }
```


```R
library(ggrepel)
```


```R
pc12 <- pca_function(data = vsd_filtered)
pc12$treat <- gsub("[1-9]","",pc12$samples)
pc12$tissue <- stringr::str_split(pc12$samples,pattern = "_",simplify = T)[,2]
```


```R
options(repr.plot.width = 7.5,repr.plot.height = 5.5)
ggplot(data = pc12,aes(x = PC1,y = PC2,col = treat)) + 
    geom_point(size = 4) + 
    geom_text_repel(data = pc12,aes(x = PC1,y = PC2,label = tissue),show.legend = F,direction = "both",max.overlaps = 15) +
    ggsci::scale_color_npg() +
    ylim(c(-29,25)) +
    theme_bw(base_size = 17) + labs(x = "PC1 : [61%]", y = "PC2 : [16%]") 
```


    
![png](output_34_0.png)
    


### Heatmap


```R
library(ComplexHeatmap)
```

    Loading required package: grid
    
    ========================================
    ComplexHeatmap version 2.22.0
    Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
    Github page: https://github.com/jokergoo/ComplexHeatmap
    Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
    
    If you use it in published research, please cite either one:
    - Gu, Z. Complex Heatmap Visualization. iMeta 2022.
    - Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
        genomic data. Bioinformatics 2016.
    
    
    The new InteractiveComplexHeatmap package can directly export static 
    complex heatmaps into an interactive Shiny app with zero effort. Have a try!
    
    This message can be suppressed by:
      suppressPackageStartupMessages(library(ComplexHeatmap))
    ========================================
    
    



```R
DEGs <- DEG.RNA_2 %>% filter(tissue != "SP") %>% .$gene %>% unique
```


```R
BulkRNA_df2 <- BulkRNA_df[,grep(names(BulkRNA_df),pattern = "SP",invert = T)] %>% dplyr::select(7:9,1:3,10:12,4:6)
```


```R
scale_row <- function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}
```


```R
BulkRNA_df.scale <- BulkRNA_df2[DEGs,] %>% na.omit
BulkRNA_df.scale[,1:6] <- scale_row(BulkRNA_df.scale[,1:6])
BulkRNA_df.scale[,7:12] <- scale_row(BulkRNA_df.scale[,7:12])
```


```R
BulkRNA_df.scale2 <- BulkRNA_df.scale[complete.cases(BulkRNA_df.scale), ]
```


```R
mat <- vst(dds, blind=FALSE) %>% assay
```


```R
mat2 <- mat[,grep(colnames(mat),pattern = "SP",invert = T)]
```


```R
mat2 <- mat2[,c(7:9,1:3,10:12,4:6)]
```


```R
mat2.scale <- mat2[DEGs,]
mat2.scale[,1:6] <- scale_row(mat2.scale[,1:6])
mat2.scale[,7:12] <- scale_row(mat2.scale[,7:12])
```


```R
mat2.scale2 <- mat2.scale[complete.cases(mat2.scale), ]
```


```R
vst.heatmap <- Heatmap(mat2.scale2,
        col = colorRampPalette(c("navy","white","firebrick3"))(100),
        clustering_distance_rows = function(m) dist(m),
        clustering_method_rows = "ward.D2",
        row_km = 3,
        use_raster = TRUE,
        border = TRUE,
        column_split = factor( rep(c("HK","TH"), each = 6)),
        heatmap_legend_param = list(
             title='Z score',
             title_gp = gpar(fontsize = 14),
             labels_gp = gpar(fontsize = 14),
             border='black'),
        show_row_names = FALSE)
vst.heatmap2 <- draw(vst.heatmap)
```

    'magick' package is suggested to install to give better rasterization.
    
    Set `ht_opt$message = FALSE` to turn off this message.
    



    
![png](output_47_1.png)
    



```R
cluster_gene <- do.call(rbind, lapply(names( row_order(vst.heatmap2) ) , function(name) {
      data.frame(cluster = name, gene = rownames(mat2.scale2)[ row_order(vst.heatmap2)[[name]] ], stringsAsFactors = FALSE)
}))
```


```R
#table(cluster_gene$gene) 获取row_cluster相关信息
row_cluster <- c( rep("cluster3",times = 932),
  rep("cluster2",times = 3317),
  rep("cluster1",times = 3611))
row_cluster <- factor(row_cluster,levels = c("cluster1","cluster2","cluster3"))
```


```R
cluster_gene
```


<table class="dataframe">
<caption>A data.frame: 7879 × 2</caption>
<thead>
	<tr><th scope=col>cluster</th><th scope=col>gene</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>1</td><td>ENSONIG00000005810</td></tr>
	<tr><td>1</td><td>ENSONIG00000043105</td></tr>
	<tr><td>1</td><td>ENSONIG00000039702</td></tr>
	<tr><td>1</td><td>ENSONIG00000032392</td></tr>
	<tr><td>1</td><td>ENSONIG00000010264</td></tr>
	<tr><td>1</td><td>ENSONIG00000015915</td></tr>
	<tr><td>1</td><td>ENSONIG00000039801</td></tr>
	<tr><td>1</td><td>ENSONIG00000020173</td></tr>
	<tr><td>1</td><td>ENSONIG00000028785</td></tr>
	<tr><td>1</td><td>ENSONIG00000035343</td></tr>
	<tr><td>1</td><td>ENSONIG00000000490</td></tr>
	<tr><td>1</td><td>ENSONIG00000030984</td></tr>
	<tr><td>1</td><td>ENSONIG00000018181</td></tr>
	<tr><td>1</td><td>ENSONIG00000016307</td></tr>
	<tr><td>1</td><td>ENSONIG00000040474</td></tr>
	<tr><td>1</td><td>ENSONIG00000029818</td></tr>
	<tr><td>1</td><td>ENSONIG00000005310</td></tr>
	<tr><td>1</td><td>ENSONIG00000030886</td></tr>
	<tr><td>1</td><td>ENSONIG00000021199</td></tr>
	<tr><td>1</td><td>ENSONIG00000007138</td></tr>
	<tr><td>1</td><td>ENSONIG00000002867</td></tr>
	<tr><td>1</td><td>ENSONIG00000038261</td></tr>
	<tr><td>1</td><td>ENSONIG00000010015</td></tr>
	<tr><td>1</td><td>ENSONIG00000017936</td></tr>
	<tr><td>1</td><td>ENSONIG00000039471</td></tr>
	<tr><td>1</td><td>ENSONIG00000018265</td></tr>
	<tr><td>1</td><td>ENSONIG00000031275</td></tr>
	<tr><td>1</td><td>ENSONIG00000033303</td></tr>
	<tr><td>1</td><td>ENSONIG00000030635</td></tr>
	<tr><td>1</td><td>ENSONIG00000019905</td></tr>
	<tr><td>⋮</td><td>⋮</td></tr>
	<tr><td>2</td><td>ENSONIG00000041821</td></tr>
	<tr><td>2</td><td>ENSONIG00000003335</td></tr>
	<tr><td>2</td><td>ENSONIG00000029469</td></tr>
	<tr><td>2</td><td>ENSONIG00000029705</td></tr>
	<tr><td>2</td><td>ENSONIG00000016302</td></tr>
	<tr><td>2</td><td>ENSONIG00000010481</td></tr>
	<tr><td>2</td><td>ENSONIG00000013081</td></tr>
	<tr><td>2</td><td>ENSONIG00000021247</td></tr>
	<tr><td>2</td><td>ENSONIG00000035337</td></tr>
	<tr><td>2</td><td>ENSONIG00000010847</td></tr>
	<tr><td>2</td><td>ENSONIG00000021268</td></tr>
	<tr><td>2</td><td>ENSONIG00000017958</td></tr>
	<tr><td>2</td><td>ENSONIG00000032419</td></tr>
	<tr><td>2</td><td>ENSONIG00000037819</td></tr>
	<tr><td>2</td><td>ENSONIG00000000641</td></tr>
	<tr><td>2</td><td>ENSONIG00000042229</td></tr>
	<tr><td>2</td><td>ENSONIG00000014454</td></tr>
	<tr><td>2</td><td>ENSONIG00000033545</td></tr>
	<tr><td>2</td><td>ENSONIG00000011278</td></tr>
	<tr><td>2</td><td>ENSONIG00000028295</td></tr>
	<tr><td>2</td><td>ENSONIG00000020399</td></tr>
	<tr><td>2</td><td>ENSONIG00000007920</td></tr>
	<tr><td>2</td><td>ENSONIG00000035611</td></tr>
	<tr><td>2</td><td>ENSONIG00000003477</td></tr>
	<tr><td>2</td><td>ENSONIG00000041753</td></tr>
	<tr><td>2</td><td>ENSONIG00000033213</td></tr>
	<tr><td>2</td><td>ENSONIG00000013172</td></tr>
	<tr><td>2</td><td>ENSONIG00000014031</td></tr>
	<tr><td>2</td><td>ENSONIG00000007701</td></tr>
	<tr><td>2</td><td>ENSONIG00000016673</td></tr>
</tbody>
</table>




```R

```


```R
 library(dplyr)
```

TF enrichment


```R
#From biomart
TF_fam <- read.table("/home/xuanwenxin/xwx_data/re_jgg/GO_KEGG_database/Oreochromis_niloticus_TF.txt",
                     sep = "\t",header = TRUE)
```


```R
TF_fam_2 <- TF_fam[,c(3,4)] %>% distinct() %>% arrange(Family) %>% dplyr::select(term = Family,gene = Ensembl)
```


```R
library(clusterProfiler)
```


```R
lapply(split(cluster_gene$gene,f = cluster_gene$cluster) ,function(genes){
    genes2 <- genes[genes %in% TF_fam_2$gene]#挑选是TF的DEGs
    print(c("TF Num :",length(genes2)))
    res <- enricher(gene = genes2,TERM2GENE = TF_fam_2,pAdjustMethod = "BH",minGSSize = 3)
    return(res)
}) -> TF_enrich_res

names(TF_enrich_res) <- paste("cluster",names(split(cluster_gene$gene,f = cluster_gene$cluster)),sep = "")
```

    [1] "TF Num :" "362"     
    [1] "TF Num :" "73"      
    [1] "TF Num :" "187"     



```R
lapply(names(TF_enrich_res),function(x){
    x2 <- TF_enrich_res[[x]]@result
    x2$cluster <- x
    return(x2)
}) -> TF_enrich_res_list
TF_enrich_res_df <- do.call(rbind,TF_enrich_res_list)
```


```R
setwd("/home/xuanwenxin/xwx_data/re_jgg/output/enrichment")
```


```R
openxlsx::write.xlsx(TF_enrich_res_df,
                     "BulkRNA_TF_enrichment.xlsx")
```


```R
TF_enrich_res_df_filter <- TF_enrich_res_df %>% 
    filter(p.adjust < 0.05) %>% group_by(cluster) %>% 
    arrange(p.adjust,.by_group = TRUE) %>% ungroup
```


```R
openxlsx::write.xlsx(TF_enrich_res_df_filter,
                     "BulkRNA_TF_enrichment_filter.xlsx")
```


```R

```


```R
GO enrichment
```


```R
#GO_KEGG富集分析相关函数
enrichment_stat_go = function(bl,il) {
#  gos <- gene_go2
#  go_des <- go_des
  gos.bl <- gos[gos$gene %in% bl,]
  gos.il <- gos[gos$gene %in% il,]
  gos.bl.stat <- gos.bl %>% group_by(GO) %>% dplyr::summarise(GeneNumber=n()) %>% ungroup
  gos.il.stat <- gos.il %>% group_by(GO) %>% dplyr::summarise(GeneNumber=n()) %>% ungroup
  list.df <- gos.bl.stat %>% group_by(GO) %>% left_join(gos.il.stat,by="GO") %>% ungroup
  list.df[is.na(list.df)] <- 0
  names(list.df) <- c("GO","k","x")#k:bl x:il
  list.df$m = as.numeric(length(unique(sort(gos.il$gene))))
  list.df$n = as.numeric(length(unique(sort(gos.bl$gene))))-list.df$k
  list.df$Pvalue <- apply(list.df[,-1],1,FUN=function(x) {phyper(x[2]-1,x[1],x[4],x[3],lower.tail=F)})
  # if (id == "go"){
  #   list <- list[list$class == "biological_process" & list$k >= 10,]
  # }
  
  #list$p.adjust <- p.adjust(list$Pvalue,method="fdr",n=length(list$Pvalue))
  list.df <- list.df[order(list.df$Pvalue,decreasing=F),]
  list.df$fold = list.df$x*(list.df$n+list.df$k)/list.df$m/list.df$k ## should add fold enrichment value
  list.1 <- list.df %>% group_by(GO) %>% left_join(go_des,by="GO") %>% ungroup
  list.1 <- list.1[list.1$k>=10 & list.1$x>=3,]
  list.1$p.adjust <- p.adjust(list.1$Pvalue,method="fdr",n=length(list.1$Pvalue))
  list.1 <- list.1[list.1$Pvalue < 0.05,]
  list.1
}

enrichment_stat_kegg = function(bl,il) {
#  keggs <- keggs
#  kegg_des <- kegg_des
  
  keggs.bl <- keggs[keggs$gene %in% bl,]
  keggs.il <- keggs[keggs$gene %in% il,]
  keggs.bl.stat <- keggs.bl %>% group_by(kegg) %>% dplyr::summarise(GeneNumber=n()) %>% ungroup
  keggs.il.stat <- keggs.il %>% group_by(kegg) %>% dplyr::summarise(GeneNumber=n()) %>% ungroup
  list.df <- keggs.bl.stat %>% group_by(kegg) %>% left_join(keggs.il.stat,by="kegg") %>% ungroup
  list.df[is.na(list.df)]<-0
  names(list.df) <- c("kegg","k","x")
  list.df$m = as.numeric(length(unique(sort(keggs.il$gene))))
  list.df$n = as.numeric(length(unique(sort(keggs.bl$gene))))-list.df$k
  list.df$Pvalue <- apply(list.df[,-1],1,FUN=function(x) {phyper(x[2]-1,x[1],x[4],x[3],lower.tail=F)})
  # if (id == "kegg"){
  #   list <- list[list$class == "biological_process" & list$k >= 10,]
  # }
  
  #list.df$p.adjust <- p.adjust(list.df$Pvalue,method="fdr",n=length(list.df$Pvalue))
  list.df <- list.df[order(list.df$Pvalue,decreasing=F),]
  list.df$fold = list.df$x*(list.df$n+list.df$k)/list.df$m/list.df$k ## should add fold enrichment value
  list.1 <- list.df %>% group_by(kegg) %>% left_join(kegg_des,by="kegg") %>% ungroup
  list.1 <- list.1[list.1$k>=10 & list.1$x>=2,]
  list.1$p.adjust <- p.adjust(list.1$Pvalue,method="fdr",n=length(list.1$Pvalue))
  list.1 <- list.1[list.1$Pvalue < 0.05,]
  list.1
}
```


```R
load("/home/xuanwenxin/xwx_data/re_jgg/GO_KEGG_database/go_info.RData")
load("/home/xuanwenxin/xwx_data/re_jgg/GO_KEGG_database/kegg_info.RData")
```


```R
lapply(split(cluster_gene$gene,f = cluster_gene$cluster),function(x){
    GO_res <- enrichment_stat_go(bl = rownames(mat2),il = x)
    return(GO_res)
}) -> GO_res
```


```R
names(GO_res) <- paste("cluster",names(split(cluster_gene$gene,f = cluster_gene$cluster)),sep = "")
```


```R
openxlsx::write.xlsx(GO_res,
                     "/home/xuanwenxin/xwx_data/re_jgg/output/enrichment/Heatmap_cluster_gene_GO.xlsx")
```


```R
head(geneID_name)
geneID_name_uniq <- geneID_name %>% group_by(gene) %>% dplyr::slice(1) %>% ungroup
head(geneID_name_uniq)
```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>gene</th><th scope=col>Gene.name</th><th scope=col>eggnog_Preferred_name</th><th scope=col>eggnog_Description</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>ENSONIG00000000001</td><td>ppfibp1a</td><td>PPFIBP1 </td><td>PTPRF interacting protein, binding protein 1 (liprin beta 1)</td></tr>
	<tr><th scope=row>2</th><td>ENSONIG00000000002</td><td>syt10   </td><td>SYT10   </td><td>Synaptotagmin X                                             </td></tr>
	<tr><th scope=row>3</th><td>ENSONIG00000000002</td><td>syt10   </td><td>CYB5R3  </td><td>b5 reductase                                                </td></tr>
	<tr><th scope=row>4</th><td>ENSONIG00000000003</td><td>syt1a   </td><td>SYT1    </td><td>Synaptotagmin Ia                                            </td></tr>
	<tr><th scope=row>5</th><td>ENSONIG00000000005</td><td>ppp1r12a</td><td>PPP1R12A</td><td>Protein phosphatase 1 regulatory                            </td></tr>
	<tr><th scope=row>6</th><td>ENSONIG00000000006</td><td>rps16   </td><td>RPS16   </td><td>Belongs to the universal ribosomal protein uS9 family       </td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A tibble: 6 × 4</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>Gene.name</th><th scope=col>eggnog_Preferred_name</th><th scope=col>eggnog_Description</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>ENSONIG00000000001</td><td>ppfibp1a</td><td>PPFIBP1 </td><td>PTPRF interacting protein, binding protein 1 (liprin beta 1)</td></tr>
	<tr><td>ENSONIG00000000002</td><td>syt10   </td><td>SYT10   </td><td>Synaptotagmin X                                             </td></tr>
	<tr><td>ENSONIG00000000003</td><td>syt1a   </td><td>SYT1    </td><td>Synaptotagmin Ia                                            </td></tr>
	<tr><td>ENSONIG00000000005</td><td>ppp1r12a</td><td>PPP1R12A</td><td>Protein phosphatase 1 regulatory                            </td></tr>
	<tr><td>ENSONIG00000000006</td><td>rps16   </td><td>RPS16   </td><td>Belongs to the universal ribosomal protein uS9 family       </td></tr>
	<tr><td>ENSONIG00000000007</td><td>myf6    </td><td>MYF6    </td><td>Myogenic factor 6                                           </td></tr>
</tbody>
</table>




```R
library(dplyr)
lapply(split(cluster_gene$gene,f = cluster_gene$cluster),function(x){
    GO_res <- enrichment_stat_go(bl = rownames(mat2),il = x)
    
    gos2 <- gos %>% group_by(gene)  %>% left_join(geneID_name_uniq,by = "gene") %>% ungroup
    
    GO_gene <- gos2[gos2$gene %in% x,] %>% group_by(GO) %>% summarise(gene = paste(gene,collapse = "/"),
                                                                      Gene.name = paste(Gene.name,collapse = "/"),
                                                                      eggnog_Preferred_name = paste(eggnog_Preferred_name,collapse = "/")
                                                                     ) %>% ungroup
    
    GO_res2 <- GO_res %>% left_join(GO_gene,by = "GO")
    return(GO_res2)
}) -> GO_res
```


```R
names(GO_res) <- paste("cluster",names(split(cluster_gene$gene,f = cluster_gene$cluster)),sep = "")
```


```R
setwd("/home/xuanwenxin/xwx_data/re_jgg/output/enrichment")
```


```R
openxlsx::write.xlsx(GO_res,
                     "Heatmap_cluster_gene_GO.xlsx")
```

    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0016435/ENSONIG0 ... SONIG00000042622 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0014273/ENSONIG0 ... SONIG00000033068 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0003748/ENSONIG0 ... SONIG00000033068 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0012314/ENSONIG0 ... SONIG00000042622 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0019718/ENSONIG0 ... SONIG00000042622 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0015699/ENSONIG0 ... SONIG00000033068 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0011895/ENSONIG0 ... SONIG00000042622 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0016705/ENSONIG0 ... SONIG00000042622 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0014801/ENSONIG0 ... SONIG00000042622 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0010013/ENSONIG0 ... SONIG00000042622 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0001266/ENSONIG0 ... SONIG00000042622 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0042955/ENSONIG0 ... SONIG00000042622 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0008106/ENSONIG0 ... SONIG00000042622 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0020052/ENSONIG0 ... SONIG00000010473 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0038934/ENSONIG0 ... SONIG00000010473 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0005359/ENSONIG0 ... SONIG00000010473 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0008444/ENSONIG0 ... SONIG00000010473 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0008468/ENSONIG0 ... SONIG00000010473 is truncated. 
    Number of characters exeed the limit of 32767.”
    Warning message in wb$writeData(df = x, colNames = colNames, sheet = sheet, startCol = startCol, :
    “0000784/ENSONIG0 ... SONIG00000010473 is truncated. 
    Number of characters exeed the limit of 32767.”



```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```
