library(tidyverse)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(stringr)
library(pcaExplorer)
library(devtools)
library(ggfortify)
library(easyGgplot2)
library(factoextra)
library(cluster)
library(RColorBrewer)
library(biomaRt)


workdir <- "/Users/Shrey/Shrey_Work/Duke/Research/Chatterjee_Lab/Oogonia_Shrey/Revisions/may_2024"
setwd(workdir)


# read counts
logtpm <- read.table('log2(norm_counts+1)_2024-02-10.csv', header=TRUE, sep=',', check.names=TRUE)
logtpm <- logtpm[, !grepl("^X$", colnames(logtpm))]
logtpm$Geneid <- sub("\\..*", "", logtpm$Geneid)

# read genes
new_genes <- read.table('new_genes_list_jan2024.csv', header=TRUE, sep=",", check.names =TRUE)


# convert Ensembl gene ID to gene symbol
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- logtpm$Geneid
gene_names <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = ensembl_ids,
                    mart = ensembl)

logtpm <- merge(logtpm, gene_names, by.x = "Geneid", by.y = "ensembl_gene_id", all.x = TRUE)
logtpm$Geneid <- logtpm$external_gene_name
logtpm <- logtpm[, !(names(logtpm) == "external_gene_name")]
colnames(logtpm)[colnames(logtpm) == "Geneid"] <- "Genes"


# remove other samples
logtpm <- na.omit(logtpm)
logtpm <- logtpm[!is.na(logtpm$Genes) & logtpm$Genes != "", ]

write.csv(logtpm, "log2(norm_counts+1)_feb2024_gene_names.csv", row.names = FALSE)

new_genes <- na.omit(new_genes)


# subset based on gene symbols
new_genes_logtpm <- subset(logtpm, Genes %in% new_genes$Genes)


new_genes <- subset(new_genes, Genes %in% unique(new_genes_logtpm$Genes))
new_genes_logtpm <- new_genes_logtpm[match(new_genes$Genes, new_genes_logtpm$Genes), ]
new_genes_logtpm[, "Genes"] <- new_genes_logtpm$Genes
new_genes_logtpm <- distinct(new_genes_logtpm)
rownames(new_genes_logtpm) <- new_genes_logtpm$Genes
new_genes_logtpm <- subset(new_genes_logtpm,select=-c(Genes))
new_genes_logtpm_mat <- data.matrix(new_genes_logtpm)


# NEW GENES HEATMAP
mat <- new_genes_logtpm_mat
new_genes_column <- mat[, 1]
#mat <- scale(mat)
new_genes <- new_genes[new_genes$Genes %in% rownames(mat), ]

# If you want to include the first column in the scaled matrix, you can bind it back
#scaled_matrix <- cbind(first_column, scaled_matrix)


f1 <- colorRamp2(seq(-6, 6, length = 3), 
                 c("blue", "white", "red"), 
                 space = "RGB")

f2 <- colorRamp2(seq(min(mat[which(!is.na(mat))]), max(mat[which(!is.na(mat))]), length = 3), 
                 c("blue", "white", "red"), 
                 space = "RGB")

f3 <- colorRamp2(c(min(mat[which(!is.na(mat))]), 
                   2*min(mat[which(!is.na(mat))])/3,
                   min(mat[which(!is.na(mat))])/3,
                   0, 
                   max(mat[which(!is.na(mat))])), c("#08306B","#2171B5","#DEEBF7", "white", "red"))

f4 <- colorRamp2(c(min(mat[which(!is.na(mat))]), 
                   2*min(mat[which(!is.na(mat))])/3,
                   min(mat[which(!is.na(mat))])/3,
                   0,
                   max(mat[which(!is.na(mat))])/3,
                   2*max(mat[which(!is.na(mat))])/3,
                   max(mat[which(!is.na(mat))])), c("#08306B","#2171B5","#DEEBF7", "white","#FB7353","#EF3B2C","#A50F15"))
colmap <- f4

colors_stage <- c("hiPSC" = "#000000",
                  "hPGCLC" = "grey18",
                  "RA-FGC" = "grey30",
                  "Oogonia" = "snow4",
                  "Maturing_Oocyte" = "snow3",
                  "Meiotic_Oocyte" = "snow2"
)

font_colors <- c("hiPSC" = "white",
                 "hPGCLC" = 'white',
                 "RA-FGC" = "white",
                 "Oogonia" = "white",
                 "Maturing_Oocyte" = "black",
                 "Meiotic_Oocyte" = "black"
)


new_genes_heat <- Heatmap(mat,
                           name = 'Log2FC',
                           cluster_rows = F, 
                           cluster_columns = F,
                           col = colmap,
                           row_names_gp = gpar(fontsize = 7.8),
                           show_row_names=TRUE,
                           column_names_rot = 90,
                           #row_names_gp = NULL,
                           column_names_gp = gpar(fontsize = 9),
                           #row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 10)),
                           column_names_max_height = max_text_height(colnames(mat), gp = gpar(fontsize = 4)),
                           heatmap_width = unit(9, "in"),
                           heatmap_height = unit(7, "in"),
                           heatmap_legend_param = list(title = gt_render("<span style='color:black'>**log2(Normalized Counts + 1)**</span>"),
                                                       title_gp = gpar(fontsize=16),
                                                       labels_gp = gpar(fontsize=16)),
                           row_split=factor(c(new_genes$Partition),levels=c("TF_Control","PGCLC","Oogonia","Oocyte","Meiosis")),
                           row_title_rot = 0,
                           column_split= NULL,
                           row_gap = unit(3, "mm"),
                           column_gap = unit(15, "mm"),
                           row_title_gp = gpar(
                             fill = colors_stage,
                             col = font_colors,
                             border = c(rep("white",6))),
                           #row_gap_color = "black"
                          
)
new_genes_heat



png(file="log2(norm_counts+1)_heat.png", width=13, height=12, units='in', res=300)
new_genes_heat
dev.off()



