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
BT1_R1_001.fastq.gz 
BT1_R2_001.fastq.gz
BT2_R1_001.fastq.gz 
BT2_R2_001.fastq.gz
...

I am going to use a for loop.
```
for i in *_R1_001.fastq.gz # This will exhaustively and non-repeatedly scan all the files ending with "_R1_001.fastq.gz"
do
var=$i # load the file name to a new object "var", you can name with other words. 
name=${var%%_*} # remove anything after the first "_"
seqkit split2 -s 50000000 -O fastq_split/ -1 $name"_R1_001.fastq.gz" -2 $name"_R2_001.fastq.gz"
done
```

xenograft="/mnt/e/Rawdata/Yawei/Genewiz_2021"
cd $xenograft/fastq_split
hisat2 -p 4 -t -x ~/human/gencode/GRCh38.p13.genome.fa -1 BT1_R1_001.part_001.fastq.gz -2 BT1_R2_001.part_001.fastq.gz -S sam/BT1_human.part_001.sam 
hisat2 -p 4 -t -x ~/mouse/gencode/GRCm39.genome.fa -1 BT1_R1_001.part_001.fastq.gz -2 BT1_R2_001.part_001.fastq.gz -S sam/BT1_mouse.part_001.sam
hisat2 -p 4 -t -x ~/human/gencode/GRCh38.p13.genome.fa -1 BT1_R1_001.part_002.fastq.gz -2 BT1_R2_001.part_002.fastq.gz -S sam/BT1_human.part_002.sam 
hisat2 -p 4 -t -x ~/mouse/gencode/GRCm39.genome.fa -1 BT1_R1_001.part_002.fastq.gz -2 BT1_R2_001.part_002.fastq.gz -S sam/BT1_mouse.part_002.sam
hisat2 -p 4 -t -x ~/human/gencode/GRCh38.p13.genome.fa -1 BT1_R1_001.part_003.fastq.gz -2 BT1_R2_001.part_003.fastq.gz -S sam/BT1_human.part_003.sam 
hisat2 -p 4 -t -x ~/mouse/gencode/GRCm39.genome.fa -1 BT1_R1_001.part_003.fastq.gz -2 BT1_R2_001.part_003.fastq.gz -S sam/BT1_mouse.part_003.sam

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

cd ../sortbam
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
sample.list2 <- 
XenofilteR(sample.list, destination.folder = "./human", bp.param = bp.param, output.names = NULL)
XenofilteR(sample.list2, destination.folder = "./mouse", bp.param = bp.param, output.names = NULL)


hisat2 -p 4 -t -x ~/mouse/gencode/GRCm39.genome.fa -1 $name"_1.fastq" -2 $name"_2.fastq" -S /mnt/e/bioinfo/sam/$name.sam

