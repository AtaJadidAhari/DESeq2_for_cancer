# import needed libs
library(DESeq2)
library(ggplot2)
library(pheatmap)

#load count matrix
expressions <- read.table("expression_counts.txt",
                         header=TRUE, row.names=1)
#change row and col to match DSEQ
exp_transpose = t(expressions)
write.csv(exp_transpose, "exp_transposed.csv")



#create annotation file for DESEQ
name <- row.names(expressions)
condition <- c("disease", "disease", "disease", 'control', 'control', 'control')
annot <- data.frame(name, condition)


#create DESeqdataset
dds <- DESeqDataSetFromMatrix(countData = exp_transpose,
                              colData = annot,
                              design = ~condition)



#remove lowly expressed genes
keep <- rowSums(counts(dds)) >= 10
high_dds <- dds[keep, ]

#run DESeq
ddsDE <- DESeq(high_dds)



#export normalized read counts
norm_counts <- counts(ddsDE, normalized=T)
write.csv(norm_counts, "normal.csv")



# results
res <- results(ddsDE, alpha = 0.05) #adjusted p-val
res_ordered <- res[order(res$padj), ]
write.csv(res, "results.csv")
summary(res)
# from LFC (log fold change) we get that:

resultsNames(ddsDE)
#  "Type_disease_vs_control"
# 38 up regulated
# 35 down regulated



#check more

normcoutns <- read.csv("normal.csv", row.names = 1)
deSeqRes <- read.csv("results.csv", row.names = 1)

# plotMA
plotMA(ddsDE, ylim = c(-6, 6))
deSeqRes$sig <- ifelse(deSeqRes$padj <= 0.05, "yes", "no")
ggplot(deSeqRes, aes(x=log10(baseMean), y=log2FoldChange, color=sig)) + geom_point()

# volcano plot
ggplot(deSeqRes, aes(x=log2FoldChange, y = -log10(padj), color=sig)) + geom_point()



# select significant genes
significants <- subset(deSeqRes, padj <= 0.05)

# select up and down regulated ones
up_regulated <- row.names(subset(significants, log2FoldChange > 0))
down_regulated <- row.names(subset(significants, log2FoldChange < 0))

# load cosmic data
cosmic_data <-read.csv("Census_allSat_Jan81918072022.csv", row.names = 1)
genes_in_cosmic <- row.names(cosmic_data)

# intersect my significant genes with cosmic data
up_intersect = intersect(up_regulated, genes_in_cosmic)
down_intersect = intersect(down_regulated, genes_in_cosmic)
