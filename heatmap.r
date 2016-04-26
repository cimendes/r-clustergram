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
if(!require("vegan")){
  install.packages("vegan", dependencies=TRUE)
  library(vegan)
}
if (!require("ape")){
  install.packages("ape",dependencies = TRUE)
  library(ape)
}


# MAPPING

# Unbiased
rpkm <- read.table("gene_presence_absence.Rtab", sep="\t", header = TRUE)
rownames(rpkm) <- rpkm[, 1]
colnames(rpkm) <- gsub('_02062016', "",colnames(rpkm),fixed=TRUE)
rpkm <- rpkm[, -1]

#grouping
groups <- read.table("grouping.txt", sep="\t", header = FALSE)
groups<- groups[,c(-4,-5)]
groups[,3]<-gsub('forstgreen','forestgreen', groups[,3],fixed = TRUE)
groups <- groups[match(colnames(rpkm), as.character(groups[,1])) ,]
stopifnot(all(as.character(groups[,1]) == colnames(rpkm)))
#rpkm<-cbind(rpkm,groups[,3])

# Variance
variance_rnaSeq <- apply(rpkm, 1, var)
variance_rnaSeq <- sort(variance_rnaSeq, decreasing = TRUE)
variance_rnaSeq <- variance_rnaSeq[1:4000]

#Dendogram from Mega
#cPhylo<-read.tree(file="ML_Bt500.nwk")
#dg=chronos(cPhylo)

# Heatmap
colour <- colorRampPalette(c("blue","white"),space="rgb")(75)
my_palette <- colorRampPalette(c("white", "blue", "darkblue"))(n = 299)
scaleyellowred <- colorRampPalette(c("red",'pink', "lightyellow"), space = "rgb")(100)

#Fixing col names
data=as.matrix(rpkm[match(names(variance_rnaSeq), rownames(rpkm)),])
col_name=colnames(data)
#col_name=gsub('_02062016', "",col_name,fixed=TRUE)
row_name=rownames(data)


#dendogram-vegan - NOT WORTH!
#calculate the Bray-Curtis dissimilarity matrix on the full dataset:
#data.dist <- vegdist(data, method = "bray")
# Do average linkage hierarchical clustering. Other options are 'average', 'complete' or 'single'. You'll need to choose the one that best fits the needs of your situation and your data.
#row.clus <- hclust(data.dist, "average")

# creates a 5 x 5 inch image
png("cluster_4000_blue.png",        
    width = 13*300,        # 5 x 300 pixels
    height = 9*300,
    res = 500,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap.2(data, 
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