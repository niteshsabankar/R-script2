library("flashClust")
library("WGCNA")
library("DESeq2")

enableWGCNAThreads()

outputDir <- "blastp_90_merged_3/"
merged = T
dir.create(outputDir)
### 1. Data pre-processing ### 
# Read expression matrices
datExprA1 <- read.table("Data/geph-original", header=TRUE, row.names = 1, sep="\t")
datExprA2 <- read.table("Data/ehux-original", header=TRUE, row.names = 1, sep="\t")
datExprA1 <- datExprA1[, -c(1, 4, 7, 10)]     # Remove sample 1, 4, 7, 10
datExprA2 <- datExprA2[, -c(1, 4, 7, 10)]

isexpr <- rowSums(datExprA1 > 10) >= 2     # Remove row with low expressions
datExprA1 <- datExprA1[isexpr, ]
isexpr <- rowSums(datExprA2 > 10) >= 2
datExprA2 <- datExprA2[isexpr, ]

# Read match file
#csv  <- read.csv("Data/Ehux_JGI_blastn_GCA_95.csv", header=TRUE)
csv  <- read.csv("Data/Ehux_JGI_blastp_GCA_90.csv", header=TRUE)
csv <- csv[!duplicated(csv[, "geph_ID"]), ]

new <- data.frame(csv, datExprA1[match(csv$geph_ID, row.names(datExprA1)), ])   # add all columns of geph to csv and create 'new' matrix
new <- data.frame(new, datExprA2[match(csv$ehux_ID, row.names(datExprA2)), ])

rownames(new) <- new$ehux_ID

new$ehux_ID <- NULL       # Remove unneccesary columns
new$geph_ID <- NULL
new$match <- NULL

new <- na.omit(new)
datExprA1 <- new[, c(1:8)]
datExprA2 <- new[, c(9:16)]

#Normalization
samples <- data.frame(samples = colnames(datExprA1))
ds <- DESeqDataSetFromMatrix(countData=datExprA1, colData=samples, design=~samples)  # Creating a DESeqDataSet object
colnames(ds) <- colnames(datExprA1)
dds <- estimateSizeFactors(ds)
log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1)
rs <- rowSums(counts(dds))
datExprA1 <- log.norm.counts[rs > 0,]

samples <- data.frame(samples = colnames(datExprA2))
ds <- DESeqDataSetFromMatrix(countData=datExprA2, colData=samples, design=~samples)  # Creating a DESeqDataSet object
colnames(ds) <- colnames(datExprA2)
dds <- estimateSizeFactors(ds)
log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1)
rs <- rowSums(counts(dds))
datExprA2 <- log.norm.counts[rs > 0,]

# Write Normalized data sets to files.
# write.csv(datExprA1, file="Data/geph-normalized")
# write.csv(datExprA2, file="Data/ehux-normalized")

### 2. Plot overall data correlations ###
softPower = 18
rankExprA1= rank(rowMeans(datExprA1))
rankExprA2= rank(rowMeans(datExprA2))

pdf(file=paste0(outputDir, "generalNetworkProperties.pdf"), height=10, width=12)
par(mar=c(5,10,4,10)+.1)
verboseScatterplot(rankExprA1,rankExprA2, xlab="Ranked Expression (G. Oceanica)", ylab="Ranked Expression (E. huxleyi)", pch = 1)
#verboseScatterplot(rankConnA1,rankConnA2, xlab="Ranked Connectivity (G. Oceanica)",ylab="Ranked Connectivity (E. huxleyi)") 
dev.off()

### 3. Put expression data into a multi-set format ####
nSets = 2
setLabels = c("G. oceanica", "E. huxleyi")
shortLabels = c("Geph", "Ehux")

# Form multi-set expression data:
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(t(datExprA1)))
multiExpr[[2]] = list(data = as.data.frame(t(datExprA2)))

# Define data set dimensions
exprSize = checkSets(multiExpr)
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;

### 4. Construct networks and Compute distances ###
# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate adjacencies in each individual data set
for (set in 1:nSets)
  adjacencies[set, , ] = adjacency(multiExpr[[set]]$data, power = softPower, type = "signed")
  #adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower;
# Initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate TOMs in each individual data set
for (set in 1:nSets)
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ], TOMType = "signed");


#Scaling of Topological Overlap Matrices to make them comparable across sets
scaleP = 0.95
set.seed(12345)
# Sample sufficiently large number of TOM entries
tomSamples = as.integer(1/(1-scaleP) * 1000);
scaleSample = sample(nGenes*(nGenes-1)/2, size = tomSamples)

TOMScalingSamples = list();
scaleQuant = rep(1, nSets)
scalePowers = rep(1, nSets)
for (set in 1:nSets)
{
  # Select the sampled TOM entries
  TOMScalingSamples[[set]] =  as.dist(TOM[set, , ])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] = quantile(TOMScalingSamples[[set]], probs = scaleP, type = 8);
  
  if(set > 1) 
  {
    scalePowers[set] =  log(scaleQuant[1])/log(scaleQuant[set]); 
    TOM[set, ,] = TOM[set, ,]^scalePowers[set]
  }
}

scaledTOMSamples = list();
for (set in 1:nSets)
  scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]

pdf(file = paste0(outputDir, "TOMScaling-QQPlot.pdf"), wi = 6, he = 6);
sizeGrWindow(6, 6)
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6,
                    xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]), main = "Q-Q plot of TOM", pch = 20)
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()

### Create modules ####
consensusTOM = pmin(TOM[1, , ], TOM[2, , ])
consTree = flashClust(as.dist(1-consensusTOM), method = "ward");

pdf(file = paste0(outputDir, "dendrogram.pdf"),height=10,width=14)
plot(consTree,xlab="",sub="",main="consensus Tree", labels=FALSE,hang=0.04);
dev.off() 

unmergedLabels = cutreeDynamic(dendro = consTree,
                             distM = 1-consensusTOM,
                             deepSplit = 2,
                             cutHeight = 30,
                             minClusterSize = 30,
                             pamRespectsDendro = FALSE);
unmergedColors = labels2colors(unmergedLabels)
print(table(unmergedLabels))

print(length(table(unmergedColors)))
pdf(file=paste0(outputDir, "Module_choices.pdf"), height=10,width=14);
plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

### Merging of modules ####
# Calculate module eigengenes
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs);
# Cluster consensus modules
consMETree = hclust(as.dist(consMEDiss), method = "average");
# Plot the result
sizeGrWindow(7,6)
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h=0.15, col = "red")
merge = mergeCloseModules(multiExpr, unmergedColors, cutHeight = 0.15, verbose = 3)
#moduleLabels = merge$colors
#moduleColors = labels2colors(moduleLabels)
mergedColors = merge$colors
consMEs = merge$newMEs;
pdf(file=paste0(outputDir, "merged_dendrogram.pdf"), width=6, height=6)
#sizeGrWindow(9,6)
plotDendroAndColors(consTree, cbind(unmergedColors, mergedColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
print(table(mergedColors))
print(length(table(mergedColors)))

#### Plot module conservation in individual trees ####
moduleColors = unmergedColors # select unmerged or merged colors
moduleMEs = unmergedMEs       # select unmergedMEs or consMEs 
if(merged) {
  moduleColors = mergedColors
  moduleMEs = consMEs
}
geneTreeA1 = flashClust(as.dist(1 - TOM[1, , ]), method = "ward")
geneTreeA2 = flashClust(as.dist(1 - TOM[2, , ]), method = "ward")
pdf(file = paste0(outputDir, "Conserved_modules.pdf"),height=8,width=12)
plotDendroAndColors(geneTreeA1, moduleColors, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors (G. Oceanica)")
plotDendroAndColors(geneTreeA2, moduleColors, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors (E. huxleyi)")
dev.off() 
dev.off() 

#multiData = list(A1=list(data=t(datExprA1)), A2=list(data=t(datExprA2)))
names(multiExpr) = c("A1", "A2")
multiColor = list(A1 = moduleColors, A2 = moduleColors)

mp = list()
for (set in 1:nSets) 
{
  mp[[set]] = modulePreservation(multiExpr, multiColor, referenceNetworks=set, verbose=3, networkType="signed", nPermutations=30)
}
stats = mp[[1]]$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
stats[order(-stats[,2]),c(1:2)] 
stats2 = mp[[2]]$preservation$Z$ref.A2$inColumnsAlsoPresentIn.A1
stats2[order(-stats2[,2]),c(1:2)] 

#### Analyze top preserved modules ###
noOfModules = 14	                                  # Number of top modules to be considered for
                                                    # module-trait correlation and Gene Enrichment Analysis
topModules <- rownames(head(stats[order(-stats[,2]),c(1:2)], noOfModules))
topModules = topModules[topModules  != "gold"]

topME = list()
for (set in 1:nSets)
{
  topME[[set]] = (moduleMEs[[set]]$data)[, paste0("ME", topModules)]
}

# Compute correlations of top Module eigengenes
noOfModules = length(topModules)
ME_cor = cor(topME[[1]], topME[[2]], method = "pearson")
MEPvalue = corPvalueStudent(ME_cor, noOfModules)
textME = paste(signif(ME_cor, 2), "\n(", signif(MEPvalue, 1), ")", sep = "")
pdf(paste0(outputDir, "MECorrelation.pdf"), height = 8, width = 10)
#par(mar = c(8,9,4,1)+.1)
#frame()
par(mfrow = c(1,1), mar=c(5,6,3,1))
labeledHeatmap(Matrix = ME_cor,
               xLabels = topModules,
               xSymbols = topModules,
               yLabels = topModules,
               ySymbols = topModules,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textME,
               setStdMargins = FALSE, 
               cex.text = 1, 
               zlim = c(-1,1), 
               main = paste("Module Eigengenes relationships. x = E. huxleyi, y = G. Oceanica"))
dev.off()

#Plot correlation heatmap between modules and traits
datTraits = read.csv("Data/datTraits.csv", row.names = 1)
for (set in 1:nSets) {
  ME_trait_cor = cor(topME[[set]], datTraits)
  ME_trait_pvalue = corPvalueStudent(ME_trait_cor, 8)
  text_ME_trait = paste(signif(ME_trait_cor, 2), "\n(", signif(ME_trait_pvalue, 1), ")", sep = "")
  pdf(paste0(outputDir, "ME_Trait_", shortLabels[[set]], ".pdf"), height=6, width = 10)
  par(mfrow = c(1,1), mar=c(5,6,3,1))
  labeledHeatmap(Matrix = ME_trait_cor,
                 xLabels = names(datTraits),
                 yLabels = topModules,
                 ySymbols = topModules,
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = text_ME_trait,
                 setStdMargins = FALSE, 
                 cex.text = 1, 
                 zlim = c(-1,1), 
                 main = paste(setLabels[[set]], "Module-trait relationships"))
  dev.off()
}


A1dir <- paste0(outputDir, "A1/")
A2dir <- paste0(outputDir, "A2/")
dir.create(A1dir)
dir.create(A2dir)

#### Output genes in modules ####
write.table(rownames(t(multiExpr[[1]]$data)), file = paste0(A1dir, "all.ids"), quote = FALSE, col.names = FALSE, row.names = FALSE)
# write ids in top modules
for (mod in topModules) 
{
  write.table(rownames(t(multiExpr[[1]]$data)[moduleColors == mod, ]), 
              file=paste0(A1dir, mod, ".ids"), 
              quote = FALSE, 
              col.names = FALSE,
              row.names = FALSE)
}

write.table(rownames(t(multiExpr[[2]]$data)), file = paste0(A2dir, "all.ids"), quote = FALSE, col.names = FALSE, row.names = FALSE)

for (mod in topModules) 
{
  write.table(rownames(t(multiExpr[[2]]$data)[moduleColors == mod, ]), 
              file=paste0(A2dir, mod, ".ids"), 
              quote = FALSE, 
              col.names = FALSE,
              row.names = FALSE)
}
