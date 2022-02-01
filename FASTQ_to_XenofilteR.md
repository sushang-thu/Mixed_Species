## FASTQ file splitting
Due to the line limit in XenofilteR processing, we need to split the raw-read fastq first.
I use the split2 function in the Seqkit pipeline.

```
## First check if you have already installed seqkit
conda list seqkit
## If not, install the seqkit
conda install seqkit
## you can use the helppage to learn the commands under seqkit
seqkit -h
seqkit split2 -h
```
We use `split2` because we are going to process pair-ended sequencing files.
```
## cd into the folder
seqkit split2 -s 50000000 -O fastq_split/ -1 BT1_R1_001.fastq.gz -2 BT1_R2_001.fastq.gz
## -s specifies splitting by sequence numbers
## -O specifies the output folder. If you dont assign a folder, seqkit will create a folder named $inputname.split in the input directory.
## -1 and -2 specify the paired reads.
## You will see the following info on the screen:
[INFO] flag -1/--read1 and -2/--read2 given, ignore: -
[INFO] split seqs from BT1_R1_001.fastq.gz and BT1_R2_001.fastq.gz
[INFO] split into 50000000 seqs per file
[INFO] write 50000000 sequences to file: fastq_split/BT1_R1_001.part_001.fastq.gz
[INFO] write 50000000 sequences to file: fastq_split/BT1_R2_001.part_001.fastq.gz
[INFO] write 50000000 sequences to file: fastq_split/BT1_R1_001.part_002.fastq.gz
[INFO] write 50000000 sequences to file: fastq_split/BT1_R2_001.part_002.fastq.gz
[INFO] write 14219233 sequences to file: fastq_split/BT1_R1_001.part_003.fastq.gz
[INFO] write 14219233 sequences to file: fastq_split/BT1_R2_001.part_003.fastq.gz
```
I am not sure about the exact maximal length that XenofilteR can process but 50000000 is OK.
You dont need to unzip the fastq gz files.
We can check the file info using the `stat` function in `seqkit`.
```
seqkit stat fastq_split/*
file                                      format  type    num_seqs        sum_len  min_len  avg_len  max_len
fastq_split/BT1_R1_001.part_001.fastq.gz  FASTQ   DNA   50,000,000  7,500,000,000      150      150      150
fastq_split/BT1_R1_001.part_002.fastq.gz  FASTQ   DNA   50,000,000  7,500,000,000      150      150      150
fastq_split/BT1_R1_001.part_003.fastq.gz  FASTQ   DNA   14,219,233  2,132,884,950      150      150      150
fastq_split/BT1_R2_001.part_001.fastq.gz  FASTQ   DNA   50,000,000  7,500,000,000      150      150      150
fastq_split/BT1_R2_001.part_002.fastq.gz  FASTQ   DNA   50,000,000  7,500,000,000      150      150      150
fastq_split/BT1_R2_001.part_003.fastq.gz  FASTQ   DNA   14,219,233  2,132,884,950      150      150      150

## Now we can see the BT1 sample R1 and R2 reads have been splitted into 3 files, separately. These files are from PE150 sequencing.
```

Next, we are going to apply commands to other files.
I have multiple pairs of fastq files.

BT1_R1_001.fastq.gz vs BT1_R2_001.fastq.gz
BT2_R1_001.fastq.gz vs BT2_R2_001.fastq.gz
...

To save coding time, I am going to use a `for` loop.
```
## Find the naming rules in your files
for i in *_R1_001.fastq.gz # This will exhaustively and non-repeatedly scan all the files ending with "_R1_001.fastq.gz"
do
var=$i # load the file name to a new object "var", you can name with any other word as long as it does not conflict with system settings. 
name=${var%%_*} # this code will remove anything after the first "_" in the string of object "var"
seqkit split2 -s 50000000 -O fastq_split/ -1 $name"_R1_001.fastq.gz" -2 $name"_R2_001.fastq.gz"
done
```
Now every pair of seq reads will be splitted one by one.

## Genomic alignment with hisat2
```
xenograft="/mnt/e/Rawdata/Yawei/Genewiz_2021"
cd $xenograft/fastq_split
## note that here I define an object called "xenograft" to store the directory of fastq files, so that I can quickly redirect to the folder by simply $xenograft
## create folders named "sam", "bam","sortbam" with the command "mkdir" under the folder $xenograft
```

Now I am going to start the alignment, onto human or mouse reference genomes.
```
for i in *_R1_001.part_001.fastq.gz
do
var=$i # load the file name to a new object "var", you can name with other words. 
name=${var%%_*} # remove anything after the first "_"
hisat2 -p 4 -t -x ~/human/gencode/GRCh38.p13.genome.fa -1 $name"_R1_001.part_001.fastq.gz" -2 $name"_R2_001.part_001.fastq.gz" -S ../sam/$name"_human.part_001.sam"
hisat2 -p 4 -t -x ~/mouse/gencode/GRCm39.genome.fa -1 $name"_R1_001.part_001.fastq.gz" -2 $name"_R2_001.part_001.fastq.gz" -S ../sam/$name"_mouse.part_001.sam" 
done

for i in *_R1_001.part_002.fastq.gz
do
var=$i # load the file name to a new object "var", you can name with other words. 
name=${var%%_*} # remove anything after the first "_"
hisat2 -p 4 -t -x ~/human/gencode/GRCh38.p13.genome.fa -1 $name"_R1_001.part_002.fastq.gz" -2 $name"_R2_001.part_002.fastq.gz" -S ../sam/$name"_human.part_002.sam"
hisat2 -p 4 -t -x ~/mouse/gencode/GRCm39.genome.fa -1 $name"_R1_001.part_002.fastq.gz" -2 $name"_R2_001.part_002.fastq.gz" -S ../sam/$name"_mouse.part_002.sam"
done

for i in *_R1_001.part_003.fastq.gz
do
var=$i # load the file name to a new object "var", you can name with other words. 
name=${var%%_*} # remove anything after the first "_"
hisat2 -p 4 -t -x ~/human/gencode/GRCh38.p13.genome.fa -1 $name"_R1_001.part_003.fastq.gz" -2 $name"_R2_001.part_003.fastq.gz" -S ../sam/$name"_human.part_003.sam"
hisat2 -p 4 -t -x ~/mouse/gencode/GRCm39.genome.fa -1 $name"_R1_001.part_003.fastq.gz" -2 $name"_R2_001.part_003.fastq.gz" -S ../sam/$name"_mouse.part_003.sam" 
done
```

### File convertion after alignment
```
cd ../sam
for i in *.sam
do
        samtools view -S $i -b > ../bam/$(basename $i .sam).bam
done

cd ../bam
for i in *.bam
do
        samtools sort $i -o ../sortbam/$(basename $i .bam).sort.bam
done
```

cd ../sortbam
mkdir human
mkdir mouse
Rscript --vanilla ~/Xenotest.R
# install the prerequisite packages first!
# BiocManager::install("Rsamtools") # don't forget the quotation marks!
# BiocManager::install(c("GenomicAlignments", "BiocParallel", "futile.logger")) # c() is a fixed grammar for array in R language.

library(XenofilteR)
bp.param <- SnowParam(workers = 1, type = "SOCK")
# sample.list creation
# learn how to create via https://bookdown.org/ndphillips/YaRrr/creating-matrices-and-dataframes.html

#setwd or cd into the folder that contains your sorted bam files.

sample <- c("BT1_human.part_001.sort.bam","BT1_mouse.part_001.sort.bam")
sample.list <- rbind(sample)
sample.list2 <- sample.list[,c(2,1)]
XenofilteR(sample.list, destination.folder = "./human", bp.param = bp.param, output.names = NULL, MM_threshold = 8)
XenofilteR(sample.list2, destination.folder = "./mouse", bp.param = bp.param, output.names = NULL, MM_threshold = 8)
## The developer suggested setting the MM-threshold to 8 (default 4) for samples with PE150

sample <- c("BT1_human.part_003.sort.bam","BT1_mouse.part_003.sort.bam")
sample.list <- rbind(sample)
sample.list2 <- sample.list[,c(2,1)]
XenofilteR(sample.list, destination.folder = "./human", bp.param = bp.param, output.names = NULL)
XenofilteR(sample.list2, destination.folder = "./mouse", bp.param = bp.param, output.names = NULL)
hisat2 -p 4 -t -x ~/mouse/gencode/GRCm39.genome.fa -1 $name"_1.fastq" -2 $name"_2.fastq" -S /mnt/e/bioinfo/sam/$name.sam

codes for analysing filtered mouse reads

cd /mnt/e/Rawdata/Yawei/Genewiz_2021/sortbam/mouse2/Filtered_bams/merged

for i in *_mouse.part_001.sort_Filtered.bam
do
var=$i # load the file name to a new object "var", you can name with other words. 
name=${var%%_*} # remove anything after the first "_"
samtools merge $name"_mouse_filtered.bam" $name*
done

featureCounts -T 8 -p -t exon -g gene_id \
-a ~/mouse/gencode/gencode.vM28.annotation.gtf \
-o 20220129mouseXenofilteR.txt \
AP1_mouse_filtered.bam AP2_mouse_filtered.bam AP3_mouse_filtered.bam \
BA2_mouse_filtered.bam BA3_mouse_filtered.bam BA4_mouse_filtered.bam \
BT1_mouse_filtered.bam BT2_mouse_filtered.bam BT3_mouse_filtered.bam \
OT3_mouse_filtered.bam OT5_mouse_filtered.bam OT6_mouse_filtered.bam \
SQ2_mouse_filtered.bam SQ4_mouse_filtered.bam SQ5_mouse_filtered.bam

## DESeq2 analyses

```
## Install DESeq2 in R if not already

## Load DESeq2
library(DESeq2)

## input
## To run a DESeq analysis, you need three objects:
## 1. cts: i.e. Counts of your gene expressions.
## 2. coldata: i.e. Grouping info
## 3. condition:

setwd("E:/Rawdata/Yawei/Genewiz_2021/sortbam/human/Filtered_bams/bamonly/merged")

## Load count matrix. Here I use the count matrix from FeatureCounts.
reads <- read.table("20220217humanXenofilteR.txt",header = T,quote = "\t",skip = 1)
reads2 <- read.csv("20220217humanXenofilteR.txt",header = T,sep = "",skip = 1)
## both are fine

## 
cts_00 <- reads[,c(1,7:15)]
colnames(cts_00)[2:10] <- c("BT1","BT2","BT3","OT3","OT5","OT6","SQ2","SQ4","SQ5")
coldata_xeno <- data.frame(c(rep("BT",3),rep("OT",3),rep("SQ",3)))
## Note that here we have 3 groups to run DESeq2 at the same time, later we can use the contrast function to specify which is the control group.
## if you only have two groups, untreated vs treated, usually the first group will be taken as the untreated control.
row.names(coldata_xeno) <- colnames(cts_00)[-1] # you dont need the first column, i.e. the gene id.
colnames(coldata_xeno) <- "condition"
head(coldata_xeno) # you will see the first 6 rows.
# If you want to see the whole dataframe, you can click the object on the right window or use 
View(coldata_xeno)
## The order of samples in coldata matrix must be identical with the order in count matrix.

rownames(cts_00) <- cts_00$Geneid
cts <- cts_00[,-1]
cts[1:6,1:3]# check the first 6 rows and 3 columns
head(cts)# by default, check the first 6 rows
## if the values in your count matrix are not all integers, then you need to round the values because DESeq2 only process integers.
## cts <- round(cts) # round the values of cts and load into an updated version of cts
## The count matrix from FeatureCounts has only integers so we dont need the round function here.

## create the DESeqDataSet object, i.e., dds.
dds <- DESeqDataSetFromMatrix(countData=cts,colData=coldata_xeno,design=~condition) # You can also name as dds_xeno or whatever.
# only count non-zero counts
dds <- dds[rowSums(counts(dds)) >= 10, ] # you may also use >=1
## If I dont filter the low counts, the total genes will be 61533 but after filtering there are only 32103.
dds <- DESeq(dds)
res <- results(dds)
res_BT2OT <- results(dds,contrast = c("condition","BT","OT")) ## Pay attention to the order of BT and OT, here the latter will be the control.
summary(res_BT2OT)
res_BT2OT
# If you simply want to know how many genes have adjusted p vlues less than 0.1 (default)
sum(res_BT2OT$padj < 0.1,na.rm = TRUE) ## You will get a total number of up and down genes from the summary.
## Or you can adjust the cut off to 0.05
res05 <- results(dds,contrast = c("condition","BT","OT"),alpha = 0.05)
summary(res05)

library(DESeq2)
res_BT2SQ <- results(dds,contrast = c("condition","BT","SQ")) 
res_OT2SQ <- results(dds,contrast = c("condition","OT","SQ")) 
write.csv(res_BT2SQ,"BTvsSQ.csv",row.names = TRUE)
write.csv(res_OT2SQ,"OTvsSQ.csv",row.names = TRUE)
BTSQ <- read.csv("BTvsSQ.csv")
OTSQ <- read.csv("OTvsSQ.csv")
BTSQ[,1] <- gsub("\\..*","",BTSQ[,1])
colnames(BTSQ)[1] <- "ensembl_gene_id"
BTSQ_annot <- merge(BTSQ,humanID,by = "ensembl_gene_id")
BTSQ_annot_pc <- subset(BTSQ_annot,gene_biotype == "protein_coding")

OTSQ[,1] <- gsub("\\..*","",OTSQ[,1])
colnames(OTSQ)[1] <- "ensembl_gene_id"
OTSQ_annot <- merge(OTSQ,humanID,by = "ensembl_gene_id")
OTSQ_annot_pc <- subset(OTSQ_annot,gene_biotype == "protein_coding")

BTSQ_DE <- subset(BTSQ_annot_pc,abs(log2FoldChange) >= log2(1.5) & padj<=0.05)
BTOT_DE <- subset(BTOT_annot_pc,abs(log2FoldChange) >= log2(1.5) & padj<=0.05)
OTSQ_DE <- subset(OTSQ_annot_pc,abs(log2FoldChange) >= log2(1.5) & padj<=0.05)
write.csv(BTOT_DE,"20220128BTOT_DE.csv",row.names = T)
write.csv(BTSQ_DE,"20220128BTSQ_DE.csv",row.names = T)
write.csv(OTSQ_DE,"20220128OTSQ_DE.csv",row.names = T)


install.packages("dplyr","tidyverse")
##separate(data = fc_BTOT,col = ensembl_gene_id,into = c("ensembl_gene_id","version"),sep = ".")
write.csv(res_BT2OT,"BTvsOT.csv",row.names = TRUE)
fc_BTOT <- read.csv("BTvsOT.csv")
colnames(fc_BTOT)[1] <- "ensembl_gene_id"
fc_BTOT[,8] <- gsub("\\..","",fc_BTOT[,1]) ## but those with version over 10 will be wrong.
fc_BTOT[,8] <- gsub("\\..*","",fc_BTOT[,1]) ## Now I add an * and it works
fc <- fc_BTOT[,c(8,2:7)]
colnames(fc)[1] <- "ensembl_gene_id"

## following needs verification

fc <- separate(data = fc_BTOT,col = "ensembl_gene_id",into = c("ensembl_gene_id","version"),sep = ".")
library(tidyverse)
humanID <- read.csv("humanID.csv")
##Note that ensembl gene ids in BTOT file have version info. We have to remove it.

gsub("\\.\\d","", "a.1") # d specifiy numbers
gsub("\\.\\d+","", "a.12") #more than 1 digits

# I used to use tidyr but I lost the successful codes.

BTOT_annot <- merge(fc,humanID,by = "ensembl_gene_id")
BTOT_annot_pc <- subset(BTOT_annot,gene_biotype == "protein_coding")
## Visualization

## a single gene

plotMA(res_BT2OT,ylim=c(-2,2))
idx <- identify(res05$baseMean,res05$log2FoldChange)
rownames(res05)[idx]
## This takes some time. Not recommended.
vsd <-vst(dds,blind = FALSE)
plotPCA(vsd,intgroup=c("condition"))

library(pheatmap)
install.packages("pheatmap")
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized = TRUE)),decreasing = TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])
pheatmap(assay(ntd)[select,],cluster_rows = FALSE,show_rownames = FALSE,cluster_cols = FALSE,annotation_col = df)
pheatmap(assay(vsd)[select,],cluster_rows = FALSE,show_rownames = FALSE,cluster_cols = FALSE,annotation_col = df)
```

