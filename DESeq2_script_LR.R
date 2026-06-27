#Project "Molecular mechanisms and energetic costs of recovery from freezing in the polyextremophile midge Belgica antarctica"

#Compiled and annotated codes for part 3a: transcriptomics, deseq2

#Set directory
setwd("__")

#expects files:
#count_recovery.csv (count matrix from stringtie/prepde.py)
#galaxy_eggnog_arthropoda_annotations_processed.tsv (annotation from eggnog/galaxy)
#phen_recovery.txt (phenotype file, consult deseq2 manual if needed)

#Load required packages
library("DESeq2")
library("ggplot2")
library("ashr")
library("pheatmap")
library("RColorBrewer")
library("apeglm")
library("EnhancedVolcano")
library("scales")


# load data ####
cts <- as.matrix(read.csv("count_recovery.csv", row.names="gene_id"))
head(cts,2)
coldata <- read.csv("phen_recovery.txt", sep="\t", row.names=1)
coldata$condition <- factor(coldata$condition) #need to modify conditions as factors
head(coldata,2)

# using featureCounts to construct the deseq dataset
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
head(dds)

#chunk below enforces the treatment names to be treated as factors, which in this case bypass the alphabetical order of C0, C1, C15, C3, C7 -> C0, C1, C3, C7, C15
dds$condition <- factor(dds$condition, levels = c("Control0","Control1", "Control3", "Control7", "Control15", 
                                                  "Frozen0", "Frozen1", "Frozen3", "Frozen7", "Frozen15"))

#Pre-filtering
dds <- dds[ rowSums(counts(dds)) > 1, ] #filter genes with 0 or 1 reads
dds <- DESeq(dds)

# for loop to iterate trough all conditions. generates results for pairwise comparisons of interest (stored in the list below)
conditions <- levels(coldata$condition)
comparisons <- list(
  c("Frozen0", "Control0"),
  c("Frozen1", "Control1"),
  c("Frozen3", "Control3"),
  c("Frozen7", "Control7"),
  #  c("Control1", "Control0"),
  c("Frozen15", "Control15")
)

# create lists to store all results
results_list <- list()
deg_lists_raw <- list()
deg_lists <- list()

# load annotation and subset fo merging with final datasets later
anno_arthropoda <- read.delim("galaxy_eggnog_arthropoda_annotations_processed.tsv", header = TRUE, stringsAsFactors = TRUE)
colnames(anno_arthropoda)

# in this case, i'm subseting specific columns to merge only columns of interest with deseq2 output
anno_arthropoda <- anno_arthropoda[, c("gene_id", "seed_ortholog", "Description", 
                                       "Preferred_name", "COG_category")]

for (comp in comparisons) {
  contrast_name <- paste(comp[1], "vs", comp[2], sep = "_")
  res <- results(dds, contrast = c("condition", comp[1], comp[2]))
  
  # apply LFC shrinkage for visualization
  res_shrink <- lfcShrink(dds, 
                          contrast = c("condition", comp[1], comp[2]), 
                          type = "ashr")
  
  # order results table by the smallest p value if needed:
  # res_ordered <- res_shrink[order(res_shrink$pvalue),] 
  
  # store both regular and shrunken results
  results_list[[contrast_name]] <- list(
    raw = res,
    shrunk = res_shrink
  )
  
  deg_lists_raw[[contrast_name]] <- as.data.frame(res)
  deg_lists[[contrast_name]] <- as.data.frame(res_shrink)
  
  # making a subset to report only the significant DEGs
  res_shrink_df <- subset(res_shrink, padj < 0.05)
  res_shrink_df <- as.data.frame(res_shrink_df)
  res_shrink_df$gene_id <- rownames(res_shrink_df)
  
  # merge annotations with final datasets
  res_annotated <- merge(res_shrink_df, anno_arthropoda, by = "gene_id", all.x = TRUE)
  
  # this codes outputs the FULL DEG table (not filtered by only significant DEGs)
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_annotated_all <- merge(res_df, anno_arthropoda, by = "gene_id", all.x = TRUE)
  
  # save FULL DEG table
 #write.csv(res_annotated_all, paste0("DEGs_full_", contrast_name, ".csv"), row.names = FALSE)

    # volcano plot per comparison
  #fist prep data to cap FC and set out-of-bound flags
#  res_shrink_plot <- as.data.frame(res_shrink) %>%
#    mutate(
#      FC_plot   = case_when(
#        log2FoldChange >  1 ~  1,
#        log2FoldChange < -1 ~ -1,
#        TRUE                ~ log2FoldChange
#      ),
#      oob_right = log2FoldChange >  1,
#      oob_left  = log2FoldChange < -1
#    )
  
  
  volcano_plot <- EnhancedVolcano(res_shrink,
                                  #lab = rownames(res_shrink),
                                  lab = rep("", nrow(res_shrink)),
                                  x = 'log2FoldChange',
                                  y = 'pvalue',
                                  #title = contrast_name,
                                  title = NULL,
                                  subtitle = NULL,
                                  caption = NULL,
                                  pCutoff = 0.05,
                                  FCcutoff = 0,
                                  pointSize = 0.5,
                                  xlim = c(-1, 1),
                                  ylim = c(0, 10),
                                  labSize = 4.0,
                                  col=c('grey40', 'grey40', 'grey40', '#54B9E3'),
                                  colAlpha = 1,
                                  gridlines.major = FALSE,
                                  gridlines.minor = FALSE,
                                  legendLabels = c("NS",
                                                   "log2FC",
                                                   "p<0.05",
                                                   "-1<log2FC>1")) + 
    theme(axis.title = element_text(size = 12)) + 
    theme(axis.text.x = element_text(size = 12)) + 
    theme(axis.text.y = element_text(size = 12)) + 
    theme(legend.position = "none") #+ coord_flip()
  volcano_plot

#volcano_plot <- volcano_plot +
#    geom_text(data = filter(res_shrink_plot, oob_right),
#              aes(x = 1, y = -log10(pvalue), label = "▶"),
#              size = 4, color = "#54B9E3",
#              inherit.aes = FALSE) +
#    geom_text(data = filter(res_shrink_plot, oob_left),
#              aes(x = -1, y = -log10(pvalue), label = "◀"),
#              size = 4, color = "#54B9E3",
#              inherit.aes = FALSE)
#volcano_plot
  
#save if needed
 ggsave(paste0("volcano_", contrast_name, "_oob.pdf"), 
         plot = volcano_plot,
         device = cairo_pdf, family = "Arial",
         width = 2.5, height = 3)
  
} 


#this should generate lists of degs (including significant and non significant) that can be used for kegg and go analyses
#should also generate volcano plots, which, together, make up figure 2 in the paper

#
RStudio.Version()$citation
sessionInfo()
#




