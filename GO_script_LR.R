#Project "Molecular mechanisms and energetic costs of recovery from freezing in the polyextremophile midge Belgica antarctica"

#Compiled and annotated codes for part 3b: transcriptomics, go

#expects files:
#"up_FC1.txt" - (and all the other lists of DEGs separated by direction)
#gene_length_processed.txt (generated in the command line. see extended methods if needed)
#gene2go_combined.txt (generated in the command line. see extended methods if needed). needs gene_id, GO_terms and length columns

#run this code to get the list of DEG by gene_id and L2FC (e.g., up_FC1.txt):
degs <- read.csv("DEGs_full_Control1_vs_Frozen1.csv")

#extract genes with positive log2FoldChange (up-regulated)
up_genes <- degs[degs$log2FoldChange > 0, ]
down_genes <- degs[degs$log2FoldChange < 0, ]
#save to .txt files
write.table(up_genes, "up_FC1.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(down_genes, "down_FC1.txt", sep = "\t", row.names = FALSE, quote = FALSE)

setwd("_")

# install and load required package 
# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("goseq", force=TRUE)
library(goseq)

# load gene2go data
# format expected: 
# column 1: gene_id
# column 2: GO_terms
# column 3: length 

gene2go_data <- read.delim("gene_length_processed.txt", header = TRUE, stringsAsFactors = FALSE) #this only has length information for genes that have go terms associated to them
head(gene2go_data)

# load significant genes from file 
significant_genes <- read.delim("up_FC1.txt", header = TRUE, stringsAsFactors = FALSE)
gene_ids <- significant_genes$gene_id

# define the background gene set
# the background should be ALL genes tested in your original DESeq2 analysis
background_genes <- unique(gene2go_data$gene_id)

# create the significance vector
# 1=significant, 0=not significant
gene_vector <- as.integer(background_genes %in% gene_ids)
names(gene_vector) <- background_genes

# create a named vector of gene lengths from gene2go file
gene_lengths <- setNames(gene2go_data$length, gene2go_data$gene_id)

# subset to include only genes present in background
gene_lengths_vector <- gene_lengths[names(gene_vector)]

# create Probability Weighting Function (PWF)
pwf <- nullp(gene_vector, bias.data = gene_lengths_vector, plot.fit = TRUE)

# GO annotation mapping, create the gene2cat list required by goseq
go_mapping <- strsplit(gene2go_data$GO_terms, ",")
names(go_mapping) <- gene2go_data$gene_id

# run GO enrichment analysis
go_results <- goseq(pwf, gene2cat = go_mapping, method = "Wallenius")

# adjust for multiple testing
go_results$padj <- p.adjust(go_results$over_represented_pvalue, method = "BH")

# filter significant results
significant_go_terms <- go_results[go_results$padj < 0.05, ]
head(significant_go_terms)

# sort by adjusted p-value
significant_go_terms <- significant_go_terms[order(significant_go_terms$padj), ]
head(significant_go_terms)

# save results
#write.csv(go_results, file = "up_FC1_results.csv", row.names = FALSE)
#write.csv(significant_go_terms, file = "up_FC1_significant_GO_terms.csv", row.names = FALSE)

##description of output columns
#category: The GO term identifier
#over_represented_pvalue: raw p-value testing if the GO term is enriched in significant gene list
#under_represented_pvalue: raw p-value testing if the GO term is not enriched
#numDEInCat: number of significant genes in list that are annotated with the respective GO term
#numInCat: total number of genes in the background set that are annotated with the GO term
#term: description of the GO term
#ontology: which GO category: BP (Biological Process), MF (Molecular Function), or CC (Cellular Component)
#padj: adjusted p-value for over-representation (Benjamini-Hochberg FDR correction) 

# summary
nrow(significant_go_terms) #number of significant GO terms
length(gene_ids) #number of input significant genes
length(background_genes) #number of genes in background se
length(gene_vector) # number of genes used in analysis

#now you may repeat the script for the timepoint/direction/contrast of interest.

#
RStudio.Version()$citation
sessionInfo()
#
