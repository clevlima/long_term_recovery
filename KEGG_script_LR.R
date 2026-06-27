#Project "Molecular mechanisms and energetic costs of recovery from freezing in the polyextremophile midge Belgica antarctica"

#Compiled and annotated codes for part 3c: transcriptomics, kegg

#Set directory
setwd("__")

#expects files:
#count_recovery.csv (count matrix, same as used in deseq2)
#gene2kegg_arthropoda.tsv (generated in the command line, see extended methods if needed)


#Required packages
library(gage)
library(pathview)
library(ggplot2)

#load the count data
cnts <- as.matrix(read.table(file = "count_recovery.csv", header = TRUE, row.names = 1, sep = ","))

## setting the kegg species list to generic whatever
kg.ko=kegg.gsets("ko")
kegg.gs=kg.ko$kg.sets[kg.ko$sigmet.idx]

# Loading the reference library
YRWA <- read.delim("gene2kegg_arthropoda.tsv", header = FALSE, stringsAsFactors = TRUE)
head(YRWA)
str(YRWA)

# parsing in the library and adding up all the KEGG term counts
D2 <- as.matrix(YRWA)
D2 <- D2[, c(1,2)]
cnts=data.frame(cnts)
cnts.data <- mol.sum(cnts, id.map=D2) 

# removing rows where total counts = 0
cnts=cnts.data
sel.rn=rowSums(cnts) !=0
cnts=cnts[sel.rn,]

# normalizing
libsizes=colSums(cnts)
size.factor=libsizes/exp(mean(log(libsizes)))
cnts.norm=t(t(cnts)/size.factor)
range(cnts.norm)
cnts.norm=log2(cnts.norm+8) #log2 + 8 is recommended in the vignette
range(cnts.norm)

# if you want to just output a file of normalized counts
#write.table(cnts.norm, "normalized_counts_eggnog.txt", sep="\t")
cnts.norm
#day 0
#significant pathways
cnts.kegg.p_0 <- gage(cnts.norm, gsets = kegg.gs, ref = 1:4, samp = 21:24, compare ="unpaired")
cnts.d_0 = cnts.norm[, 21:24]-rowMeans(cnts.norm[, 1:4])

# vector of colors (e.g., 100 shades from grey40 to white to purple) for heatmap
heatmap_colors <- colorRampPalette(c("grey40", "white", "#54B9E3"))(100)

# this definse fixed breaks from -4 to +4, which ensures that the color gradient is consistent across all heatmaps
my_breaks <- seq(-4, 4, length.out = 101)  # 101 breaks for 100 colors

#making heatmap with these pathways and saving it in current directory
kegg.sig_0 <- sigGeneSet(cnts.kegg.p_0, 
                         cutoff = 0.05, 
                         outname = "Heatmap_FC0", 
                         pdf.size = c(8.5, 11), 
                         margins = c(10, 30), #These margins worked ok for my on Adobe Acrobat. Change the cutoff if desired
                         cexRow = 1.5,    # row label size
                         cexCol = 1.5,    # column label size
                         cex.main = 1.5,  # main title size
                         breaks = my_breaks,
                         col = heatmap_colors)

kegg.sig_0

write.table(rbind(kegg.sig_0$greater, kegg.sig_0$less), #returns a table with both significant up- and downregulated pathways
            file = "DE_FC0.txt", sep = "\t")


#day 1
#significant pathways
cnts.kegg.p_1 <- gage(cnts.norm, gsets = kegg.gs, ref = 5:8, samp = 25:28, compare ="unpaired")
cnts.d_1 = cnts.norm[, 25:28]-rowMeans(cnts.norm[, 5:8])

kegg.sig_1 <- sigGeneSet(cnts.kegg.p_1, 
                         cutoff = 0.05, 
                         outname = "Heatmap_FC1", 
                         pdf.size = c(8.5, 11), 
                         margins = c(10, 30), 
                         cexRow = 1,    
                         cexCol = 1.25, 
                         cex.main = 1.5,  
                         breaks = my_breaks,
                         col = heatmap_colors)
kegg.sig_1

write.table(rbind(kegg.sig_1$greater, kegg.sig_1$less), 
            file = "DE_FC1.txt", sep = "\t")


#day 3
cnts.kegg.p_3 <- gage(cnts.norm, gsets = kegg.gs, ref = 9:12, samp = 29:32, compare ="unpaired")
cnts.d_3= cnts.norm[, 29:32]-rowMeans(cnts.norm[, 9:12])

kegg.sig_3<-sigGeneSet(cnts.kegg.p_3, 
                       cutoff = 0.05, 
                       outname = "Heatmap_FC3", 
                       pdf.size = c(8.5, 11), 
                       margins = c(10, 30),
                       cexRow = 1.5,    
                       cexCol = 1.5,    
                       cex.main = 1.5, 
                       col = heatmap_colors,
                       breaks = my_breaks)
kegg.sig_3

write.table(rbind(kegg.sig_3$greater, kegg.sig_3$less), #returns a table with both significant up- and downregulated pathways
            file = "DE_FC3.txt", sep = "\t")


#day 7
cnts.kegg.p_7 <- gage(cnts.norm, gsets = kegg.gs, ref = 13:16, samp = 33:36, compare ="unpaired")
cnts.d_7 = cnts.norm[, 33:36]-rowMeans(cnts.norm[, 13:16])

kegg.sig_7<-sigGeneSet(cnts.kegg.p_7, 
                       cutoff = 0.05, 
                       outname = "Heatmap_FC7", 
                       pdf.size = c(8.5, 11), 
                       margins = c(10, 30),
                       cexRow = 1.5,    
                       cexCol = 1.5,    
                       cex.main = 1.5,
                       breaks = my_breaks,
                       col = heatmap_colors)
kegg.sig_7

write.table(rbind(kegg.sig_7$greater, kegg.sig_7$less),
            file = "DE_FC7.txt", sep = "\t")



#day 15
cnts.kegg.p_15 <- gage(cnts.norm, gsets = kegg.gs, ref = 17:20, samp = 37:40, compare ="unpaired")
cnts.d_15 = cnts.norm[, 37:40]-rowMeans(cnts.norm[, 17:20])


#making heatmap with these pathways and saving it in current directory
kegg.sig_15<-sigGeneSet(cnts.kegg.p_15, 
                        cutoff = 0.05, 
                        outname = "Heatmap_FC15", 
                        pdf.size = c(8.5, 11), 
                        margins = c(10, 30),
                        cexRow = 1.5,    
                        cexCol = 1.5,    
                        cex.main = 1.5,
                        breaks = my_breaks,
                        col = heatmap_colors)
kegg.sig_15

write.table(rbind(kegg.sig_15$greater, kegg.sig_15$less),
            file = "DE_FC15.txt", sep = "\t")


#download pathways of interest individually this puts out an overlay of the differentially expressed genes from data (needs to change the pathway.id and don't forget to change the output name!!)

#this first gives a file with expression value separated by replicate
pv.out <- pathview(gene.data = cnts.d_0, pathway.id = "03050", 
                   species = "ko", out.suffix = "FC0_proteasome", 
                   kegg.native = T, 
                   low=list(gene="grey40"), high=list(gene="#54B9E3")) 


pv.out <- pathview(gene.data = cnts.d_1, pathway.id = "04810", 
                   species = "ko", out.suffix = "FC1_cytoskeleton", 
                   kegg.native = T, 
                   low=list(gene="grey40"), high=list(gene="#54B9E3")) 


#this code below does the same as above but without the replicate separation 
#(i think it takes the mean expression values of all replicates)
pv.out <- pathview(
  gene.data   = cnts.d_0[, 1],
  pathway.id  = "03050",
  species     = "ko",
  out.suffix  = "FC0_proteasome",
  kegg.native = TRUE,
  low         = list(gene = "grey40"),
  mid         = list(gene = "white"),
  high        = list(gene = "#54B9E3")
)

pv.out <- pathview(
  gene.data   = cnts.d_1[, 1],
  pathway.id  = "04810",
  species     = "ko",
  out.suffix  = "FC1_cytoskeleton",
  kegg.native = TRUE,
  low         = list(gene = "grey40"),
  mid         = list(gene = "white"),
  high        = list(gene = "#54B9E3")
)


#these should generate S2 through S5.For figs S3 and S4, I used the pathview graphs with replicates.

#
sessionInfo()
# 
