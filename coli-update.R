#---------Bulk RNAseq Analysis

#Set up the working directory----
setwd("Desktop/Decode_Workshop/Bulk_RNAseq/")

#Loading required packages----
#install.packages("R.utils")
library(R.utils)

#install.packages("data.table")
library(data.table)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("RUVSeq")
library(RUVseq)

#(if the above installation of RUVseq didn't work, try this)
#source("http://bioconductor.org/biocLite.R")
#biocLite("RUVSeq")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
library(DESeq2)

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

#install.packages("pheatmap")
library(pheatmap)

#install.packages("RColorBrewer")
library(RColorBrewer)

#install.packages("ggplot2")
library(ggplot2)


#Working on actual data----
f1 <- read.csv("merged_counts.txt", sep="\t", row.names = 1)

#Normal stage count----Ciprofloxacin
count_of_T1=length(grep(x = colnames(f1), pattern = "^Control."))
count_of_T1

#Drug stage count----
count_of_T2=length(grep(x = colnames(f1), pattern = "^Colistin."))
count_of_T2

#Filtering----
filter <- apply(f1, 1, function(x) length(x[x>0])>=1)
filtered <- f1[filter,]

genes <- rownames(filtered)
genes

# Set the conditions
t1<-rep(c("Control"),each=count_of_T1)
t2<-rep(c("Colistin"),each=count_of_T2)
x<-c(t1,t2)
x<-as.factor(x)
x
set <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(x,row.names=colnames(filtered)))
set

#Setting Color theme
colors <- brewer.pal(3, "Set2")

#Plotting basic graphs
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

# Normalization with help of RuvSeq for removing unwanted variations
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

differences <- makeGroups(x)
differences
set3 <- RUVs(set, genes, k=1, differences)
pData(set3)

# Assuming 'x' is a column in pData(set3)
pData(set3)$x <- relevel(pData(set3)$x, ref = "Control")

# Checking factor levels after releveling
levels(pData(set3)$x)

# Ensuring pData and countData dimensions match
all(rownames(pData(set3)) == colnames(counts(set3)))


#Computing Differential expression
dds <- DESeqDataSetFromMatrix(countData = counts(set3),colData = pData(set3), design= ~ W_1 + x)

dds <- DESeq(dds)

res <- results(dds)
head(res)

write.table(res, file = "coli.2mg.15min-update.csv", sep = ",", col.names = NA, qmethod = "double")


#Volcano plot for all the degs with Enchanced Valcano Package for downregulated and upregulated genes

# Load EnhancedVolcano package
library(EnhancedVolcano)

# Create a more colorful volcano plot 
res$regulation <- ifelse(res$padj < 0.05 & res$log2FoldChange > 1, "Upregulated",
                         ifelse(res$padj < 0.05 & res$log2FoldChange < -1, "Downregulated", "Not Significant"))

volcano_plot <- EnhancedVolcano(res,
                                lab = NA,                          
                                x = 'log2FoldChange',              
                                y = 'padj',                        
                                pCutoff = 0.05,                    
                                FCcutoff = 1.0,                    
                                pointSize = 2.0,                   
                                colAlpha = 1,                      
                                legendPosition = 'right',          
                                legendLabSize = 10,                
                                legendIconSize = 3.0,              
                                drawConnectors = FALSE,            
                                title = "Volcano Plot of Colistin.2mg.15min DEGs",
                                caption = "Thresholds: |Log2FC| > 1, Adjusted P-value < 0.05",
                                border = 'full',
                                borderWidth = 1,
                                borderColour = 'black') +
  geom_point(aes(color = regulation), size = 2, alpha = 1) + 
  scale_color_manual(values = c("Upregulated" = "forestgreen", 
                                "Downregulated" = "red2", 
                                "Not Significant" = "grey70")) +
  theme(axis.title.x = element_text(size = 12),     
        axis.title.y = element_text(size = 12),     
        axis.text.x = element_text(size = 12),       # X-axis tick labels size
        axis.text.y = element_text(size = 12),       # Y-axis tick labels size
        legend.text = element_text(size = 10),       
        plot.title = element_text(size = 14),   
        plot.caption = element_text(size = 10, hjust = 0.5),
        panel.grid.major = element_blank(),           
        panel.grid.minor = element_blank(),           
        plot.background = element_rect(fill = "white")) 

# Print the plot
print(volcano_plot)

# Save the volcano plot as a high-quality TIFF image for publication
ggsave("volcano_plot_high_quality.png", 
       plot = volcano_plot,         # The plot object to save
       device = "png",             # Save as TIFF
       width = 6,                   # Width of the image in inches
       height = 4,                  # Height of the image in inches
       units = "in",                # Measurement unit for width and height
       dpi = 300)                   # Resolution in DPI (300 DPI for high quality)



2. # Example of selected genes
selected_genes <- c("J696_03586", "J696_01953", "J696_00992","J696_01562","J696_01285","J696_03105","J696_03951","J696_02654","J696_02935")  # Replace with actual gene IDs

# Subset the data to include only the selected genes
f1_subset <- f1[rownames(f1) %in% selected_genes, ]
pheatmap(f1_subset, 
         scale = "row",       # Scale expression data by row
         cluster_rows = TRUE,  # Cluster genes (rows)
         cluster_cols = TRUE,  # Cluster samples (columns)
         show_rownames = TRUE, # Show gene names
         show_colnames = TRUE) # Show sample names
pheatmap(f1_subset, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 8, 
         fontsize_col = 8, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row")
# Save the heatmap as a high-quality TIFF image
tiff("heatmap_high_quality.tiff", width = 8, height = 6, units = "in", res = 300)
pheatmap(f1_subset, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 8, 
         fontsize_col = 8, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row")
dev.off()  # Close the device to finalize saving the file


3. Plot the downregulated and upregulated gene using ggplot
# Load necessary libraries
library(ggplot2)
filtered_res$regulation <- ifelse(filtered_res$padj < 0.05 & abs(filtered_res$log2FoldChange) > 1, 
                                  ifelse(filtered_res$log2FoldChange > 1, "Upregulated", "Downregulated"), 
                                  "Non-significant")
library(ggplot2)

ggplot(filtered_res, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Upregulated" = "forestgreen", 
                                "Downregulated" = "red2", 
                                "Non-significant" = "grey70")) +  # Color for non-significant genes
  labs(title = "Upregulated, Downregulated, and Non-significant Genes",
       x = bquote(~Log[2]~ 'Fold Change'),
       y = bquote(~-Log[10]~ 'Adjusted P-value')) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),       # X-axis tick labels size
        axis.text.y = element_text(size = 12),       # Y-axis tick labels size
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),           
        panel.grid.minor = element_blank(),           
        plot.title = element_text(size = 14),
        axis.line = element_line(color = "black"))

