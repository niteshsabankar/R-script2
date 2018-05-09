options(warn=-1)							# Ignore warnings

library("flashClust")
library("WGCNA")
library("DESeq2")

datExprA1 <- read.table("Data/geph-original", header = TRUE, row.names = 1, sep = "\t")		# Read G. Oceanica RNA-seq file
datExprA2 <- read.table("Data/ehux-original", header = TRUE, row.names = 1, sep = "\t")		# Read E. huxleyi RNA-seq file
datTraits <- read.csv("Data/datTraits.csv", row.names = 1)					# Read traits file
csv  <- read.csv("Data/Ehux_JGI_blastp_GCA_90.csv", header=TRUE)				# Read gene match file

csv <- csv[!duplicated(csv[, "geph_ID"]), ]													# Remove duplicates from matched genes

datExprA1 <- datExprA1[, -c(1, 4, 7, 10)]     # Remove sample 1, 4, 7, 10
datExprA2 <- datExprA2[, -c(1, 4, 7, 10)]

isexpr <- rowSums(datExprA1 > 10) >= 2     # Remove low count rows
datExprA1 <- datExprA1[isexpr, ]
isexpr <- rowSums(datExprA2 > 10) >= 2
datExprA2 <- datExprA2[isexpr, ]

new <- data.frame(csv, datExprA1[match(csv$geph_ID, row.names(datExprA1)), ])   # Add all columns of geph to
new <- data.frame(new, datExprA2[match(csv$ehux_ID, row.names(datExprA2)), ])	# csv and create 'new' matrix

rownames(new) <- new$ehux_ID

new$ehux_ID <- NULL       # Remove unneccesary columns
new$geph_ID <- NULL
new$match <- NULL

new <- na.omit(new)
dim(new)

datExprA1 <- new[, c(1:8)]
datExprA2 <- new[, c(9:16)]

#Normalization
samples <- data.frame(samples = colnames(datExprA1))
ds <- DESeqDataSetFromMatrix(countData = datExprA1, colData = samples, design =~ samples)  # Creating a DESeqDataSet object
colnames(ds) <- colnames(datExprA1)
dds <- estimateSizeFactors(ds)
log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1)
rs <- rowSums(counts(dds))
datExprA1 <- log.norm.counts[rs > 0, ]

samples <- data.frame(samples = colnames(datExprA2))
ds <- DESeqDataSetFromMatrix(countData = datExprA2, colData = samples, design =~ samples)  # Creating a DESeqDataSet object
colnames(ds) <- colnames(datExprA2)
dds <- estimateSizeFactors(ds)
log.norm.counts <- log2(counts(dds, normalized = TRUE) + 1)
rs <- rowSums(counts(dds))
datExprA2 <- log.norm.counts[rs > 0, ]

#-------------------------------------------------------------------------------------
# Generate graphs for choosing softpower

# For Geph
powers = c(c(1:8), seq(from = 10, to = 20, by = 4), seq(from = 22, to = 90, by = 5))
sft = pickSoftThreshold(t(datExprA1), powerVector = powers, verbose = 5)

pdf("GephSft.pdf",height = 5, width = 9)
par(mfrow = c(1,2))
cex1 = 0.7

plot(sft$fitIndices[,1], 
	-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	xlab = "Soft Threshold (power)",
	ylab = "Scale Free Topology Model Fit,signed R^2",
	type = "n",
	main = paste("Scale independence"))

text(sft$fitIndices[,1],
	-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	labels = powers,
	cex = cex1,
	col = "red")

abline(h = 0.90, col = "red")

plot(sft$fitIndices[,1],
	sft$fitIndices[,5],
	xlab = "Soft Threshold (power)",
	ylab = "Mean Connectivity",
	type = "n",
	main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")
dev.off()

#--------------------------------------
# For Ehux
powers = c(c(1:22))

sft = pickSoftThreshold(t(datExprA2), powerVector = powers, verbose = 5)

pdf("EhuxSft.pdf", height = 5, width = 9)
par(mfrow = c(1,2))
cex1 = 0.7

plot(sft$fitIndices[,1],
	-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	xlab = "Soft Threshold (power)",
	ylab = "Scale Free Topology Model Fit,signed R^2",
	type = "n",
	main = paste("Scale independence"))

text(sft$fitIndices[,1],
	-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	labels = powers,
	cex = cex1,
	col = "red")

abline(h = 0.90, col = "red")

plot(sft$fitIndices[,1], 
	sft$fitIndices[,5],
	xlab = "Soft Threshold (power)",
	ylab = "Mean Connectivity", 
	type = "n",
	main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")
dev.off()

#-------------------------------------------------------------------------------------
# Generate ranked expresion graph which shows the correlation between two data sets

softPower = 18		# change this variable to change the softpower

rankExprA1 = rank(rowMeans(datExprA1))
rankExprA2 = rank(rowMeans(datExprA2))

pdf("generalNetworkProperties.pdf", height = 10, width = 12)
par(mar = c(5,10,4,10)+.1)
verboseScatterplot(rankExprA1,rankExprA2, 
		xlab = "Ranked Expression (G. Oceanica)", 
		ylab = "Ranked Expression (E. huxleyi)",
		pch = 1)
dev.off()

adjacencyA1 = adjacency(t(datExprA1), power = softPower,type = "signed")
diag(adjacencyA1) = 0
TOMA1 = TOMsimilarity(adjacencyA1, TOMType = "signed")
geneTreeA1 = flashClust(as.dist(1 - TOMA1), method = "ward")

adjacencyA2 = adjacency(t(datExprA2), power = softPower, type = "signed")
diag(adjacencyA2) = 0
TOMA2 = TOMsimilarity(adjacencyA2, TOMType = "signed")
geneTreeA2 = flashClust(as.dist(1 - TOMA2), method = "ward")

#-------------------------------------------------------------------------------------
#Scaling of Topological Overlap Matrices to make them comparable across sets

setLabels = c("G. Oceanica", "E. huxleyi")
nGenes = length(rownames(datExprA1))
nSets = 2
scaleP = 0.95
set.seed(12345)
nSamples = as.integer(1/(1-scaleP) * 1000)
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list()
scaleQuant = rep(1, nSets)
scalePowers = rep(1, nSets)

TOMScalingSamples[[1]] = as.dist(TOMA1)[scaleSample]
TOMScalingSamples[[2]] = as.dist(TOMA2)[scaleSample]
scaleQuant[1] = quantile(TOMScalingSamples[[1]], probs = scaleP, type = 8)
scaleQuant[2] = quantile(TOMScalingSamples[[2]], probs = scaleP, type = 8)

scalePowers[2] = log(scaleQuant[1])/log(scaleQuant[2])
TOMA2 = TOMA2^scalePowers[2]

scaledTOMSamples = list()
for (set in 1:2)
scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]

pdf(file = "TOMScaling-QQPlot.pdf", wi = 6, he = 6)
qqUnscaled = qqplot(TOMScalingSamples[[1]],
		TOMScalingSamples[[2]], 
		plot.it = TRUE, 
		cex = 0.6,
		xlab = paste("TOM in", setLabels[1]), 
		ylab = paste("TOM in", setLabels[2]), 
		main = "Q-Q plot of TOM", pch = 20)

qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)

points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)

abline(a = 0, b = 1, col = "blue")

legend("topleft", 
		legend = c("Unscaled TOM", "Scaled TOM"), 
		pch = 20, 
		col = c("black", "red"))
dev.off()

#-------------------------------------------------------------------------------------
# Generate consensus TOM using which modules are decided. Change 'cutHeight' and 'deepSplit'
# option in cutreeDynamic function to change number of modules and size of modules

consensusTOM = pmin(TOMA1, TOMA2)

consTree = flashClust(as.dist(1-consensusTOM), method = "ward")

pdf("dendrogram.pdf", height = 10, width = 14)
plot(consTree, xlab = "", sub = "", main = "consensus Tree", labels = FALSE, hang = 0.04)
dev.off() 

moduleLabels = cutreeDynamic(dendro = consTree,
			distM = 1-consensusTOM,
			deepSplit = 2,
			cutHeight = 30,
			minClusterSize = 30,
			pamRespectsDendro = FALSE)

moduleColors = labels2colors(moduleLabels)
table(moduleColors)

pdf("Module_choices.pdf", height = 10, width = 14)
plotDendroAndColors(consTree, 
		moduleColors, 
		"Dynamic Tree Cut",
		dendroLabels = FALSE, 
		hang = 0.03,
		addGuide = TRUE, 
		guideHang = 0.05)
dev.off()

#-------------------------------------------------------------------------------------
# Impose modules generated using consensud tree on two datasets and plot the graph

pdf("Final_modules.pdf", height = 8, width = 12)
plotDendroAndColors(geneTreeA1, 
		moduleColors, 
		"Modules", 
		dendroLabels = FALSE, 
		hang = 0.03, 
		addGuide = TRUE,
		guideHang = 0.05, 
		main = "Gene dendrogram and module colors (G. Oceanica)")

plotDendroAndColors(geneTreeA2, 
		moduleColors, 
		"Modules", 
		dendroLabels = FALSE, 
		hang = 0.03, 
		addGuide = TRUE,
		guideHang = 0.05, 
		main = "Gene dendrogram and module colors (E. huxleyi)")
dev.off() 

#-------------------------------------------------------------------------------------
# Check module preservation between two datasets. This will calculate Z-score which 
# shows to what extend gene expressions in each module are same in both data sets

multiData = list(A1 = list(data = t(datExprA1)), A2 = list(data = t(datExprA2)))
multiColor = list(A1 = moduleColors)
mp=modulePreservation(multiData,
		multiColor, 
		referenceNetworks = 1, 
		verbose = 3, 
		networkType = "signed", 
		nPermutations = 30)

stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
stats[order(-stats[,2]),c(1:2)]

write.csv(stats[order(-stats[,2]),c(1:2)], file = "PreservationTable.csv")

#-------------------------------------------------------------------------------------
# Decide how many modules to be considered for further analysis. Modules starting from
# top in module preservation list are considered for analysis as they are sorted according 
# to their Z-score  

noOfModules = 12	# Number of top modules to be considered for
					# module-trait correlation and Gene Enrichment Analysis

topModules <- rownames(head(stats[order(-stats[,2]),c(1:2)], noOfModules))

MEtopModules=character()
index = 1
for (i in topModules) {
	MEtopModules[index] <-  paste("ME", i, sep = "")		# Add string 'ME' to before each module 
	index <- index + 1						# for module eigengene calculations and graph
}

#-------------------------------------------------------------------------------------
# Calculate module eigen of each module in two data sets

PCs1A = moduleEigengenes(t(datExprA1), colors = moduleColors)
ME_1A = PCs1A$eigengenes

PCs2A = moduleEigengenes(t(datExprA2), colors = moduleColors)
ME_2A = PCs2A$eigengenes

# Keep module Eigengenes of only 10 highly preserved modules
ME_1A_preserved <- ME_1A[, MEtopModules]		# Considering only top modules for further analysis
ME_2A_preserved <- ME_2A[, MEtopModules]		# according to noOfModules chosen by the user.

#-------------------------------------------------------------------------------------
# Generate correlation heatmap between module eigengenes of two datasets.

ME_cor <- cor(ME_1A_preserved, ME_2A_preserved, use = "p")
MEPvalueA1 = corPvalueStudent(ME_cor, noOfModules)
textME = paste(signif(ME_cor, 2), "\n(", signif(MEPvalueA1, 1), ")", sep = "")
pdf("MECorrelation.pdf", height = 8, width = 10)
par(mar = c(8,9,4,1)+.1)

labeledHeatmap(Matrix = ME_cor,
            xLabels = names(ME_2A_preserved),
            xSymbols = names(ME_2A_preserved),
            yLabels = names(ME_1A_preserved),
            ySymbols = names(ME_1A_preserved),
            colorLabels = FALSE,
            colors = blueWhiteRed(50),
            textMatrix = textME,
            setStdMargins = FALSE, 
            cex.text = 1, 
            zlim = c(-1,1), 
            main = paste("Module Eigengenes relationships. x = E. huxleyi, y = G. Oceanica"))
dev.off()

#-------------------------------------------------------------------------------------
# Generate Correlation heatmaps between module eigengene and traits for two data sets

# For Geph
moduleTraitCorA1 = cor(ME_1A_preserved, datTraits, use = "p")
moduleTraitPvalueA1 = corPvalueStudent(moduleTraitCorA1, 8)

#Print correlation heatmap between modules and traits
textMatrixA1 = paste(signif(moduleTraitCorA1, 2), "\n(", signif(moduleTraitPvalueA1, 1), ")", sep = "")
dim(textMatrixA1) = dim(moduleTraitCorA1)
pdf("heatmapA1.pdf", height = 6, width = 10)
par(mar = c(5,9,4,1)+.1)

labeledHeatmap(Matrix = moduleTraitCorA1,
            xLabels = names(datTraits),
            yLabels = names(ME_1A_preserved),
            ySymbols = names(ME_1A_preserved),
            colorLabels = FALSE,
            colors = blueWhiteRed(50),
            textMatrix = textMatrixA1,
            setStdMargins = FALSE, 
            cex.text = 1, 
            zlim = c(-1,1), 
            main = paste("G. Oceanica Module-trait relationships"))
dev.off()

#--------------------------------------
# For Ehux
moduleTraitCorA2 = cor(ME_2A_preserved, datTraits, use = "p")
moduleTraitPvalueA2 = corPvalueStudent(moduleTraitCorA2, 8)

#Print correlation heatmap between modules and traits
textMatrixA2 = paste(signif(moduleTraitCorA2, 2), "\n(", signif(moduleTraitPvalueA2, 1), ")", sep = "")
dim(textMatrixA2) = dim(moduleTraitCorA2)

#display the corelation values with a heatmap plot

pdf("heatmapA2.pdf", height = 6, width = 10)
par(mar = c(5,9,4,1)+.1)

labeledHeatmap(Matrix = moduleTraitCorA2,
            xLabels = names(datTraits),
            yLabels = names(ME_2A_preserved),
            ySymbols = names(ME_2A_preserved),
            colorLabels = FALSE,
            colors = blueWhiteRed(50),
            textMatrix = textMatrixA2,
            setStdMargins = FALSE, 
	        cex.text = 1,
            zlim = c(-1,1),
            main = paste("E. huxleyi Module-trait relationships"))
dev.off()

#-------------------------------------------------------------------------------------
# Combined HeatMap for module eigen genes and traits for two data sets to assess similarity
# between two data set's ME and traits correlation

# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCorA1), ncol(moduleTraitCorA1))
consensusPvalue = matrix(NA, nrow(moduleTraitCorA1), ncol(moduleTraitCorA1))

# Find consensus negative correlations
negative = moduleTraitCorA1 < 0 & moduleTraitCorA2 < 0
consensusCor[negative] = pmax(moduleTraitCorA1[negative], moduleTraitCorA2[negative])
consensusPvalue[negative] = pmax(moduleTraitPvalueA1[negative], moduleTraitPvalueA2[negative])

# Find consensus positive correlations
positive = moduleTraitCorA1 > 0 & moduleTraitCorA2 > 0
consensusCor[positive] = pmin(moduleTraitCorA1[positive], moduleTraitCorA2[positive])
consensusPvalue[positive] = pmax(moduleTraitPvalueA1[positive], moduleTraitPvalueA2[positive])

textMatrix = paste(signif(consensusCor, 2), "\n(", signif(consensusPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCorA2)

pdf("heatmapA1-A2.pdf", height = 6,width = 10)
par(mar = c(5,9,4,1)+.1)

labeledHeatmap(Matrix = consensusCor,
            xLabels = names(datTraits),
            yLabels = names(ME_1A_preserved),
            ySymbols = names(ME_1A_preserved),
            colorLabels = FALSE,
            colors = blueWhiteRed(50),
            textMatrix = textMatrix,
            setStdMargins = FALSE, 
			cex.text = 1,
       	    zlim = c(-1,1),
	    main = paste("G. Oceanica - E. huxleyi Module-trait relationships"))

dev.off()
#-------------------------------------------------------------------------------------
# Generate gene list for each module and store it in text files for two data sets
# Requires Directories A1 and A2 in project directory for saving text files
# A1 for storing module gene text file for dataset-1 and A2 for dataset-2.

modulesA1 <- as.data.frame(datExprA1)
modulesA1$module <- moduleColors

modulesA2 <- as.data.frame(datExprA2)
modulesA2$module <- moduleColors

for (i in topModules) {
	i <- noquote(i)
	A1path <- paste("A1/", paste(i, "A1", sep = ""), sep = "")
	A2path <- paste("A2/", paste(i, "A2", sep = ""), sep = "")

	write.table(rownames(modulesA1[grep(paste("^", i, "$", sep = ""), modulesA1$module),]), 
		file = A1path, 
		col.names = FALSE, 
		row.names = FALSE, 
		quote = FALSE)
	write.table(rownames(modulesA2[grep(paste("^", i, "$", sep = ""), modulesA2$module),]), 
		file = A2path, 
		col.names = FALSE, 
		row.names = FALSE, 
		quote = FALSE)
}