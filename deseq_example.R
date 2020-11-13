#################################################
# LOAD REQUIRED PACKAGES

library("DESeq2")
library("ggplot2")
library("ggrepel")

#################################################
# PREP INPUT TABLE

# set your working dir
setwd("/Shares/CL_Shared/data/atma/1_QUARANTINE/9_NOV2020/testing_wiki_scripts/deseq_from_K27ac_chipseq/counts")

# read in the tab-delimited count table
countdata <- read.table("bamCountsWithinPeaks.tab", sep="\t", header = FALSE)

# check the table looks correct
head(countdata)

# create a new column which is the chr and coordinates combined
countdata$chr <- paste(countdata$V1,countdata$V2, sep=":")
countdata$coord <- paste(countdata$chr,countdata$V3, sep="-")
#countdata$coord

# set countdata$coord to be the rownames
rownames(countdata) <- countdata$coord
head(countdata)

# remove the first four columns (chr, start, stop, regionlabel)
# and the last two columns (coord, chr) since we don't need them anymore
# retaining only the counts
countdata <- countdata[, c(5, 6, 7, 8)] 
head(countdata)

# add bam names as column names
colnames(countdata) <- c("0hr_sw480_run1","0hr_sw480_run2","16hr_sw480_run1","16hr_sw280_run2")
head(countdata)

# convert table to matrix format
countdata <- as.matrix(countdata)
head(countdata)

# assign control vs treated samples
condition <- factor(c(rep("control", 2), rep("treated", 2)))

# create a "coldata" table containing the sample names with their appropriate condition (e.g. control versus cancer sample)
coldata <- data.frame(row.names=colnames(countdata), condition)
coldata

#################################################
# RUN DESEQ2

# construct a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
dds

# relevel to set the controls as the reference levels
dds$condition <- relevel(dds$condition, ref = "control")
dds

# run deseq2
dds <- DESeq(dds)
resultsNames(dds)

# call results
res_treated_vs_control <- results(dds, contrast=c("condition", "treated", "control"))
head(res_treated_vs_control)

# omit rows with counts of "N/A" 
res_treated_vs_control <- na.omit(res_treated_vs_control)
head(res_treated_vs_control)

# sort results by ascending adjusted pvalue
res_treated_vs_control <- res_treated_vs_control[order(res_treated_vs_control$padj), ]
head(res_treated_vs_control)

# report the number of rows with an adjusted pvalue less than 0.05
table(res_treated_vs_control$padj<0.05)

#################################################
# SAVE RESULTS TABLES

# save this output table to fiji
# table will be saved to the current working dir (set at the start of this script)
# if you want to save it elsewhere, include the path before the file name
write.table(res_treated_vs_control, file="res_treated_vs_control.tab", quote = FALSE, row.names = TRUE, col.names=NA, sep = "\t")

# also save a table of the raw counts and normalized counts
raw_counts = counts(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(raw_counts, file="raw_counts.tab", quote = FALSE, row.names = TRUE, col.names=NA, sep = "\t")
write.table(normalized_counts, file="normalized_counts.tab", quote = FALSE, row.names = TRUE, col.names=NA, sep = "\t")

# extract regions that are significantly different in treated samples compared to control
# separate into two files based on whether they are significantly up (treated > control) or down (control > treated)
sigUp = subset(res_treated_vs_control, padj<0.05 & log2FoldChange>1)
head(sigUp)
nrow(sigUp)
write.table(sigUp, file="sigUp.tab", quote = FALSE, row.names = TRUE, col.names=NA, sep = "\t")

sigDown = subset(res_treated_vs_control, padj<0.05 & log2FoldChange<(-1))
head(sigDown)
nrow(sigDown)
write.table(sigDown, file="sigDown.tab", quote = FALSE, row.names = TRUE, col.names=NA, sep = "\t")

#################################################
# PLOT STUFF

# make an MA plot
plotMA(res_treated_vs_control, ylim=c(-10,10))

# save a copy of the MA plot 
dev.copy(png,'treated_vs_control_MAplot.png')
dev.off()

# make a volcano plot using ggplot2
# first, make the results table a data frame
res_treated_vs_control <- as.data.frame(res_treated_vs_control)
res_treated_vs_control

ggplot(res_treated_vs_control, aes(log2FoldChange, -log10(padj)), colour="grey") +
  scale_color_discrete(name = 'Labels') +
  theme_bw() + 
  labs(y="-log10 adjusted pvalue", x = "log2 fold change") +
  # set all dots to be grey
  geom_point(data=res_treated_vs_control, colour = "grey") + 
  # if pvalue<0.05, change dot color to green
  geom_point(data=res_treated_vs_control[which(res_treated_vs_control $padj <0.05),], colour = "springgreen2") + 
  # if log2FC >1, change dot color to orange
  geom_point(data=res_treated_vs_control[which(abs(res_treated_vs_control $log2FoldChange)>1),], colour = "darkgoldenrod1") +
  # if both, change dot color to blue
  geom_point(data=res_treated_vs_control[which(abs(res_treated_vs_control $log2FoldChange)>1 & res_treated_vs_control$padj<0.05),], colour = "royalblue1") +
  # add text labels to the most significant regions
  geom_text_repel(data =res_treated_vs_control[which(res_treated_vs_control $padj <0.000005),], mapping = aes(log2FoldChange, -log10(padj), label = rownames(res_treated_vs_control[which(res_treated_vs_control $padj <0.000005),])),size = 4,force = 1)

# save a copy of the volcano plot
dev.copy(png, res=200, height = 1000, width = 1000, pointsize=4, 'treated_vs_control_volcano.png')
dev.off()
