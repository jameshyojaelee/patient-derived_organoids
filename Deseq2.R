library(tximport)


setwd("D:/patient-derived_organoids/Deseq2")
dir <- "D:/patient-derived_organoids/quant_kallisto"

library(ensembldb)
library(EnsDb.Hsapiens.v86)

EnsDb <- EnsDb.Hsapiens.v86
EnsDb

k <- keys(EnsDb, keytype = "TXNAME") #or GENENAME
tx2gene <- select(EnsDb, k, "GENENAME", "TXNAME") # close library(dplyr) before running this line
head(tx2gene)


samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
files <- file.path(dir, samples$run, "abundance.h5")
names(files) <- c("hF2_TAK",
                  "hF2_Untreated",
                  "hF3_TAK",
                  "hF3_Untreated",
                  "hF23_TAK",
                  "hF23_Untreated",
                  "hF44_TAK",
                  "hF44_Untreated",
                  "hM1e_TAK",
                  "hM1e_Untreated",
                  "hM19a_TAK",
                  "hM19a_Untreated")

txi.kallisto <- tximport(files, type="kallisto", tx2gene = tx2gene, txOut=TRUE) #transcript level
head(txi.kallisto$counts)

txi.kallisto2 <- tximport(files, type="kallisto", tx2gene = tx2gene, ignoreTxVersion=TRUE) #gene level
head(txi.kallisto2$counts)


################################################################################

library(DESeq2)
ddsTxi <- DESeqDataSetFromTximport(txi.kallisto2,
                                   colData = samples,
                                   design = ~ condition)
head(counts(ddsTxi))
ddst <- estimateSizeFactors(ddsTxi)
sizeFactors(ddst)

raw_counts <- counts(ddst, normalized=FALSE)
write.table(raw_counts, file="raw_counts.txt", sep="\t", quote=F, col.names=NA)


normalized_counts <- counts(ddst, normalized=TRUE)

normalized_counts <- as.data.frame(normalized_counts)
library(dplyr)
normalized_counts <- tibble::rownames_to_column(normalized_counts, "NAME") #convert the matrix names to a column
#normalized_counts$TXID <- gsub("\\..*","",normalized_counts$TXID) #delete anything after dot


write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)

tpm <- 1e6*t(t(normalized_counts[c(1:nrow(normalized_counts)),])/colSums(normalized_counts))

write.table(tpm, file="tpm.txt", sep="\t", quote=F, col.names=NA)

head(tpm) # take a look at the TPM values
colSums(tpm) # make sure that TPM values add up to the same number for each sample
log2.tpm <- log2(tpm+1)


#Run DeSeq2 and save results
dds <- DESeq(ddsTxi)


res <- results(dds)
res.df <- as.data.frame(res)
head(res.df)

#poltMA
plotMA(res)

res.df[is.na(res.df[,2]),2] <- 0 # replace null in log2 fold change by 0
res.df[is.na(res.df[,6]),6] <- 1 # replace null in adjusted p-values by 1
DE.bools <- abs(res.df[,2]) > 1 & res.df[,6] < 0.05 # filter for fold change greater than 2 and adjusted p-value less than 0.05

heatmap(log2.tpm[DE.bools,])


################# PCA analysis#########################
library("pcaExplorer")

pcaExplorer(dds=dds)
