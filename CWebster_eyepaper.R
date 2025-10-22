#all R script for publication - seasonal and daily transcriptome of mammalian eye
# Cara Webster - September 2025

#RNA read quality checking, trimming, read alignment with genome, and gene expression counts were all completed outside of R using commandline
#htseq-count data was generated for each biological replicate to then be inputted into R
#any script with seasonal data listed was also used with daily data, just with substituted naming conventions

###DESEQ2
library(DESeq2)
library(sva)
directory <- "C:/Documents/seasonal_counts"
sampleFiles <- grep("1.STAR",list.files(directory),value=TRUE)
metaData <- read.csv("seasonal_eye2.csv", fileEncoding = 'UTF-8-BOM')
metaData1 <- metaData[order(metaData$file_name),]

ID <- dplyr::pull(metaData1, ID)
season <- dplyr::pull(metaData1, season)
file_name <- dplyr::pull(metaData1, file_name)
batch <- dplyr::pull(metaData1, batch)

sampleTable <- data.frame(sampleName=sampleFiles, fileName=sampleFiles, season=season, ID=ID, batch=batch)
sampleTable$season <-factor(sampleTable$season)
sampleTable$batch <- factor(sampleTable$batch)
#for daily dataset, sex was also made into factor
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory = directory, design= ~ batch + season)
#for daily time point data, the following ddsHTSeq object was created
#ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory = directory, design= ~ sex + time)
countdata <- assay(ddsHTSeq)
countdata <- countdata[rowSums(countdata) > 0, ]
mod <- model.matrix(~ batch + season, data = colData(ddsHTSeq))
mod0 <- model.matrix(~ batch, data = colData(ddsHTSeq))
svobj <- svaseq(countdata, mod, mod0)

for (i in 1:ncol(svobj$sv)) {
  colData(ddsHTSeq)[[paste0("SV", i)]] <- svobj$sv[, i]
}

design(ddsHTSeq) <- ~ batch + SV1 + SV2 + season
#design(ddsHTSeq) <- ~ SV1 + SV2 + sex + time
dds <- DESeq(ddsHTSeq, test = "LRT", full = ~ batch + SV1 + SV2 + season, reduced = ~ batch + SV1 + SV2)
#final dds model used for daily time point
#dds <- DESeq(ddsHTSeq, test = "LRT", full = ~ SV1 + SV2 + sex + season, reduced = ~ SV1 + SV2 + sex)
dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 5 ) >= 3
eye_dds <- dds[idx,]
eye_res <- results(eye_dds)

###WGCNA - original code modified from fuzzyatelin.github.io
library(WGCNA)
normalized_counts <- counts(eye_dds, normalized=TRUE)
vst_data <- DESeq2::vst(eye_dds)
wgcna_input <- assay(vst_data)
wgcna_input_transposed <- t(wgcna_input)
spt <- pickSoftThreshold(wgcna_input_transposed)
#soft threshold of 5 used for seasonal, 6 for daily
softPower <- 5
adjacency <- adjacency(wgcna_input_transposed, power = softPower)

TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1-TOM
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")
Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
ModuleColors <- labels2colors(Modules)
MElist <- moduleEigengenes(wgcna_input_transposed, colors = ModuleColors)
MEs <- MElist$eigengenes
ME.dissimilarity = 1-cor(MElist$eigengenes, use = "complete")
METree = hclust(as.dist(ME.dissimilarity), method = "average")
merge <- mergeCloseModules(wgcna_input_transposed, ModuleColors,
                           cutHeight = 0.25)
mergedColors = merge$colors
mergedMEs = merge$newMEs

samples <- rownames(wgcna_input_transposed)
traitRows <- match(samples, metaData$file_name)
metaData_ordered <- metaData[match(samples, metaData$file_name), ]
stopifnot(all(metaData_ordered$file_name == samples))
datTraits <- metaData_ordered[, -1]
rownames(datTraits) <- metaData_ordered$file_name
datTraits <- binarizeCategoricalColumns(datTraits)
datTraits$season.fall.vs.all <- c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                  0, 0, 0, 0)

datTraits <- datTraits[, -c(1, 5)]
names(datTraits)[1] <- "spring"
names(datTraits)[2] <- "summer"
names(datTraits)[3] <- "winter"
names(datTraits)[4] <- "fall"

nGenes = ncol(wgcna_input_transposed)
nSamples = nrow(wgcna_input_transposed)
module.trait.correlation = cor(mergedMEs, datTraits, use = "p")
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples)

textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
                   signif(module.trait.Pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.correlation)
boldEntries <- abs(module.trait.correlation) > 0.5 | module.trait.Pvalue < 0.05
textMatrixCustom <- ifelse(boldEntries,
                           paste0("*", signif(module.trait.correlation, 2), "\n(",
                                  signif(module.trait.Pvalue, 1), ")"),
                           paste(signif(module.trait.correlation, 2), "\n(",
                                 signif(module.trait.Pvalue, 1), ")"))
dim(textMatrixCustom) <- dim(module.trait.correlation)
labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(datTraits),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrixCustom,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = "")

#define variable from datTrait, for now I am going to use winter as an example
winter = as.data.frame(datTraits$winter)
names(winter) = "winter"
modNames = substring(names(mergedMEs), 3)

geneModuleMembership = as.data.frame(cor(wgcna_input_transposed, mergedMEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

geneTraitSignficance = as.data.frame(cor(wgcna_input_transposed, winter, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignficance), nSamples))
names(geneTraitSignficance) = paste("GS.", names(winter), sep = "")
names(GSPvalue) = paste("p.GS.", names(winter), sep = "")

###mFuzz
library(Mfuzz)
filtered <- dds[genefilter::rowVars(assay(eye_dds)) > 0.5, ]
vst_data <- vst(filtered, blind=TRUE)
grouped <- rowsum(t(assay(vst_data)), group = eye_dds$season)
replicates_per_group <- as.numeric(table(vst_data$season))
averaged <- t(grouped / replicates_per_group)
season_order <- c("fall", "winter", "spring", "summer")

eset <- new("ExpressionSet", exprs=averaged)
eset_standard <- standardise(eset)
cl <- mfuzz(eset_standard, c=10, m=2)
mfuzz.plot2(eset_standard, cl=cl, mfrow = c(2, 2), time.labels = c("fall", "winter", "spring", "summer"),
            centre = TRUE, centre.col = "black")
cl$size

acore_list <- acore(eset_standard, cl=cl, min.acore = 0.2)
#input_genes <- list of significant gene names
get_gene_scores <- function(input_genes, acore_list) {
  results_score <- data.frame(gene = input_genes, membership_score = I(vector("list", length(input_genes))), 
                              clusters = I(vector("list", length(input_genes)))
  )
  for (i in seq_along(input_genes)) {
    gene <- input_genes[i]
    score <- c()
    gene_clusters <- c()
    
    for (j in seq_along(acore_list)) {
      df <- acore_list[[j]]
      
      matches <- df[df$NAME == gene, ]
      
      if (nrow(matches) > 0) {
        score <- c(score, as.numeric(matches$MEM.SHIP))
        gene_clusters <- c(gene_clusters, j)
      }
    }
    results_score$membership_score[[i]] <- if (length(score) > 0) score else NA
    results_score$clusters[[i]] <- if (length(gene_clusters) > 0) gene_clusters else NA
  }
  return(results_score)
}
results_score <- get_gene_scores(input_genes, acore_list)

###gprofiler2 GO analysis
#custom_id <- upload_GMT_file("custom_filtered.gmt")
#gp__M4Ac_Y300_5Eg <- new custom ID based on tadbra genome
res <- gost(query = sig_season, significant = FALSE, organism = 'gp__M4Ac_Y300_5Eg')

eye_res <- as.data.frame(res$result)
