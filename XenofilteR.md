Here I will use XenofilteR to pickup human-specific reads and mouse-specific reads from two test files.
Please fork or clone the original XenofilteR package via https://github.com/PeeperLab/XenofilteR
This is an R package.
```
# install the prerequisite packages first!
# BiocManager::install("Rsamtools") # don't forget the quotation marks!
# BiocManager::install(c("GenomicAlignments", "BiocParallel", "futile.logger")) # c() is a fixed grammar for array in R language.
```
Load the library and set up environment.
```
library(XenofilteR)
bp.param <- SnowParam(workers = 1, type = "SOCK")
```

Download the test files or use your own files. Test files [link](https://github.com/PeeperLab/XenofilteR/tree/master/inst/extdata)
Use setwd or cd into the folder that contains your sorted bam files.

The XenofilteR requires you to provide the sorted bam files aligned to human or mouse reference genome separately from the same raw-read fastq.
Then you can assign the host and graft simply by putting the sorted bam aligned to graft first then that to the host into a sample.list.

Create a sample.list 
Here I create 2 list because I want to test the difference between host and graft settings.
`# learn how to create via https://bookdown.org/ndphillips/YaRrr/creating-matrices-and-dataframes.html`
```
sample <- c("Test_hg19_NRAS.bam","Test_mm10_NRAS.bam") # human as graft
sample.list <- rbind(sample)

sample2 <- c("Test_mm10_NRAS.bam","Test_hg19_NRAS.bam")
sample.list2 <- rbind(sample2)
```
Create 2 folders named human and mouse respectively. Then run the XenofilteR. Remember:
XenofilteR will filter out all reads that map to host only or map both to host and graft, so only graft-specific reads will be left after filteration.

```
XenofilteR(sample.list, destination.folder = "./human", bp.param = bp.param, output.names = NULL) 
# human as graft so human-reads will be left.
XenofilteR(sample.list2, destination.folder = "./mouse", bp.param = bp.param, output.names = NULL) 
# mouse as graft so mouse-reads will be left.
```
The results intepretation
Everytime you run XenofilteR, you will get a folder named `Filtered_bams` under your assigned output folder.
There will be 3 files in this folder: 
- yourinputname_Filtered.bam
- yourinputname_Filtered.bam.bai
- XenofilteR.log

For my test, I will get `Filtered_bams` folder under the folder `human` with the following 3 files:
- Test_hg19_NRAS_Filtered.bam
- Test_hg19_NRAS_Filtered.bam.bai
- XenofilteR.log

For the human as graft, XenofilteR "Filtered 59 read pairs out of 645  -  9.15 Percent";
For the mouse as graft, XenofilteR "Filtered 185 read pairs out of 233  -  79.4 Percent";
