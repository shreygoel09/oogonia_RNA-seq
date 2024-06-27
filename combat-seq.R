library(sva)

workdir <- "/Users/Shrey/Shrey_Work/Duke/Research/Chatterjee_Lab/Oogonia_Shrey/Revisions"
setwd(workdir)


# read counts
raw_counts <- read.table('raw_counts.csv', header = TRUE, sep = ",", check.names=TRUE)
row.names(raw_counts) <- raw_counts$Geneid
raw_counts <- raw_counts[, -which(names(raw_counts) == "Geneid")]

# run ComBat-seq
raw_counts <- matrix(rnbinom(2729205, size=10, prob=0.1), nrow=60649, ncol=45)
batch <- c(rep('hiPSC', 2), rep('hPGCLC', 2), rep('NANOS3', 12), rep('D3', 3), rep('Yatsenko', 3), rep('Yamashiro', 12), rep('Zhang', 11))
#group <- rep(c('hiPSC', 'hPGCLC', 'Oogonia-Like Cells'), times=c(2, 12, 3))
counts <- ComBat_seq(raw_counts, batch=batch, group=NULL)

write.csv(raw_counts, file = "combat-seq_counts.csv", row.names = FALSE)




pranam_counts <- read.table('pranam_genes_logtpm_mat.csv', header = TRUE, sep = ',', check.names=TRUE)
rownames(pranam_counts) <- pranam_counts$X
pranam_counts <- subset(pranam_counts, select = -X)

pranam_counts <- matrix(rnbinom(1105, size=10, prob=0.1), nrow=65, ncol=17)
batch <- c(rep('hiPSC', 2), rep('NANOS3', 12), rep('D3', 3))
# hiPSC <- rep(c(1, 0), times=c(11, 54))
# PGC <- rep(c(0, 1, 0), times=c(11, 14, 40))
# Oogonia <- rep(c(0, 1, 0), times=c(25, 24, 16))
# Oocyte <- rep(c(0, 1), times=c(49, 16))
# covar_mat <- cbind(hiPSC, PGC, Oogonia, Oocyte)
# pranam_counts <- ComBat_seq(pranam_counts, batch=batch, group=NULL, covar_mod=covar_mat)
adjusted_pranam_counts <- ComBat_seq(pranam_counts, batch=batch, group=NULL)

#group <- rep(c('hiPSC', 'PGC', 'Oogonia', 'Oocyte'), times=c(11, 14, 24, 16))

