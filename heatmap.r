# Installing Dependencies
if (!require("gplots")){
  install.packages("gplots",dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")){
  install.packages("RColorBrewer",dependencies = TRUE)
  library(RColorBrewer)
}
if(!require("Heatplus")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("Heatplus")  # annHeatmap or annHeatmap2
  library(Heatplus)
}

# MAPPING
#read gene_presence_absence.Rtab file
file_gene = file.choose()
filename_gene = basename(file_gene)

gene_data <- read.table(filename_gene, sep="\t", header = TRUE)
rownames(gene_data) <- gene_data[, 1]
colnames(gene_data) <- gsub('_02062016', "",colnames(gene_data),fixed=TRUE)
gene_data <- gene_data[, -1]

# GROUPPING
#read grouping.txt file
file_groups = file.choose()
filename_groups = basename(file_groups)

groups <- read.table(filename_groups, sep="\t", header = FALSE)
groups<- groups[,c(-4,-5)]
groups[,3]<-gsub('forstgreen','forestgreen', groups[,3],fixed = TRUE)
groups <- groups[match(colnames(gene_data), as.character(groups[,1])) ,]
stopifnot(all(as.character(groups[,1]) == colnames(gene_data)))

# Heatmap colour pallete
my_palette <- colorRampPalette(c("white", "blue", "darkblue"))(n = 299)


##### 4000 GENES WITH HIGHER VARIANCE #####

# variance - select 4000 genes with higher variance
variance_rnaSeq <- apply(gene_data, 1, var)
variance_rnaSeq <- sort(variance_rnaSeq, decreasing = TRUE)
variance_rnaSeq <- variance_rnaSeq[1:4000]

#Fixing col names
data_variance=as.matrix(gene_data[match(names(variance_rnaSeq), rownames(gene_data)),])
col_name=colnames(data_variance)
row_name=rownames(data_variance)

# creates a 5 x 5 inch image
png("clustergram_4000variable.png",        
    width = 13*300,        # 5 x 300 pixels
    height = 9*300,
    res = 500,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap.2(data_variance, 
          scale='none',
          notecol="black", 
          labCol=col_name,
          density.info="none",  
          trace="none",         
          margins =c(12,9),
          ColSideColors = as.character(groups[,3]),
          col=my_palette,
          labRow='', 
          dendrogram="both")
dev.off()

##### ALL DATA #####

#Fixing col names
data_all=as.matrix(gene_data)
col_name=colnames(data_all)
row_name=rownames(data_all)

# creates a 5 x 5 inch image
png("clustergram_allData.png",        
    width = 13*300,        # 5 x 300 pixels
    height = 9*300,
    res = 500,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap.2(data_all, 
          scale='none',
          notecol="black", 
          labCol=col_name,
          density.info="none",  
          trace="none",         
          margins =c(12,9),
          ColSideColors = as.character(groups[,3]),
          col=my_palette,
          labRow='', 
          dendrogram="both")
dev.off()