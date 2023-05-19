# Ncounter_normalization

```r
setwd("C:/1Philekorea/Nanostring/nCounter/WorldKimchi Research Institute/48 RCC file")

############## Required R packages ############
library(ggplot2)
library(NanoStringNorm)
library(data.table)
library(NormqPCR)
library(ComplexHeatmap)
library(dplyr)
library(matrixStats)
library(plotrix)
library(scales)
library()

install.packages("rlang")
install.packages("NanoStringNorm", force = TRUE)
install.packages("gdata")
install.packages("C:/1Philekorea/Nanostring/nCounter/WorldKimchi Research Institute/NanoStringNorm", repos = NULL, type = "source")
install.packages("data.table")
remove.packages("vctrs")
remove.packages("remotes")
remove.packages("remotes", lib = .libPaths())
install.packages("remotes")
remotes::install_version("vctrs", version = "0.6.0")


Libraries <- c('matrixStats', 'ruv', 'NanoStringNorm', 'data.table','NormqPCR',
               'ComplexHeatmap','dplyr', 'plotrix', 'tuple', 'scales', 'ggplot2')
lapply(Libraries, require, character.only = TRUE)


#####  Reading sample annotation and raw counts #####
#***** Reading Nanostring gene expression raw data
Nano_ExpressionMatrix <- read.delim('WorldKimchirawcounts.txt',
                                    stringsAsFactors = FALSE, header = TRUE, as.is = TRUE)
dim(Nano_ExpressionMatrix) # 614 165
# all Endogenous genes (600) and Nanostring negative and positive spike-in controls (14)

#***** Reading sample and clinical information
Nano_SampleInfo <- read.delim('WorldKimchisampleinform.txt',
                              stringsAsFactors = FALSE, header = TRUE, as.is = TRUE)
dim(Nano_SampleInfo) # 48 7
table(Nano_SampleInfo$Tissues)

#***** Expression matrix
Nano_RawCounts <- as.matrix(Nano_ExpressionMatrix[ , 4:ncol(Nano_ExpressionMatrix)])
dim(Nano_RawCounts) # 54 48
sum(Nano_RawCounts == 0) # 0
row.names(Nano_RawCounts) <- Nano_ExpressionMatrix$Name
Nano_SampleInfo$SampleNames <- colnames(Nano_RawCounts)


######  Exploratory  data analysis – raw counts ##########
#***** Colors for each cartridges
Color_Batches <- c('purple','orange','darkred','blue','chartreuse',
                   'darkgoldenrod4','tan2','darkgreen','red3','darkmagenta',
                   'deeppink','violet','navy','red','dodgerblue')


#***** box plot of raw data - Endogenous genes only
RawCounts_log <- log2(Nano_RawCounts[1:48 , ]) # Excluding 13 housekeeping genes and Nanostring spike-ins
par(mar = c(6.5,6.5,2.3,0), mgp = c(3.7 , 1 , 0))
boxplot(RawCounts_log, las = 1, cex.axis = 2, ylab = '' , xlab = '', cex.lab = 4,
        xaxt = 'n', yaxt = 'n', main = 'Unnormalized counts', cex.main = 3.5,
        outline = FALSE, names = FALSE, frame = FALSE,
        whisklty = 3, whisklwd = 1.5, staplelty = 1, notch = TRUE, boxlwd = 2,
        staplelwd = 0 , boxcol = Color_Batches[factor(Nano_SampleInfo$Cartridges)],
        border = Color_Batches[factor(Nano_SampleInfo$Cartridges)] , col='gray87')
box(lwd = 7, bty = 'l')
axis(1, cex.axis=1, at = c(1, seq(20,166,20)), cex.axis = 2.5, lwd.ticks = 4, mgp = c(3.5,1.6,0))
axis(2, at = c(0, seq(3,15,3)), mgp = c(3.5,.9,0), lwd.ticks = 4, las = 1, cex.axis=3)
mtext(expression(paste('Samples', '(', 'n'[samples], '=', '162', ')')), 1, line = 4.5, cex = 2.5)
mtext(expression(paste(Log[2],' (raw counts)')), 2, line = 3.5, cex = 2.7)

#***** RLE plot - Figure 1 A - Unormalized
par(mar = c(6.5,6.5,2.3,0))
boxplot(RawCounts_log - rowMedians(RawCounts_log),
        main = '', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', ylim = c(-4.1,4.1),
        outline = FALSE, names = FALSE, frame = FALSE, whisklty = 3, whisklwd = 1.5, staplelty = 1, notch = TRUE, boxlwd = 2,
        staplelwd = 0 , boxcol = Color_Batches[factor(Nano_SampleInfo$Cartridges)],
        border = Color_Batches[factor(Nano_SampleInfo$Cartridges)], col = 'gray87')
box(lwd = 7, bty = 'l')
title('Unnormalized counts', line = -2, cex.main = 3.5)
Median_RawData <- apply(RawCounts_log - rowMedians(RawCounts_log), 2, median)
points(c(1:ncol(RawCounts_log)), Median_RawData, col = Color_Batches[factor(Nano_SampleInfo$Cartridges)], pch = 19, cex = 1.2)
axis(2, mgp = c(3.5, .9 ,0), lwd.ticks=6, las=1, cex.axis=3)
mtext('RLE', 2, line = 3.5, cex = 3.5)
abline(h = 0, col = 'black', lwd = 5, lty = 2)
par(lwd = 3)
axis.break(2, -4.2, style = 'zigzag', brw = .02)
legend(160, 4.1, legend = c(1,2,3,'.','.','.', 13,14,15),
       col = c(Color_Batches[1:3], rep('white', 3), Color_Batches[13:15]),
       pch = 19, bty = 'n', cex = 1.4)
text(x = 162, y = 4.2 ,labels  = 'Cartridges', cex = 1.5)
rm(Median_RawData)

#***** Average plot
### Average of Nanostring positive spike-ins controls
Mean_NegativeControlProbes <- apply(Nano_RawCounts[41:48, ], 2, mean)
### Average of Nanostring positive spike-ins controls
Mean_PositiveControlProbes <- apply(Nano_RawCounts[49:54, ], 2, mean)
### Average of housekeeping genes
Mean_HousekeepingGenes <- apply(Nano_RawCounts[38:40, ], 2, mean)
### Library size
LibrarySize <- colSums(Nano_RawCounts [ 1:37, ])

#### Average plots -  Supplementary Figure 2
LibrarySize <- colSums(Nano_RawCounts[1:37, ])
par(mar = c(6, 7, 0, 0))
plot(log2(LibrarySize), ylim = c(0, 22), bty = 'l', typ = 'n', ylab = '', xlab = '', xaxt = 'n', yaxt = 'n')
axis(1, cex.axis = 1, at = c(1, seq(20, 166, 20)), cex.axis = 2, lwd.ticks = 4, mgp = c(3.5, 1.4, 0))
axis(2, cex.axis = 1, cex.axis = 2, lwd.ticks = 4, mgp = c(3.5, 1, 0), las = 1)
mtext(expression(paste(Log [2], ' (raw counts)')), 2, line = 3.5, cex = 3)
mtext(expression(paste('Samples', ' (', 'n'[samples], '=', '48', ')')), 1, line = 4.5, cex = 2.5)
X <- c(0, 6, 18, 30, 42, 53, 65, 77, 89, 101, 113, 125, 137, 143, 155, 167)
GradiantColors <- paste0('gray', seq(90, 20, by = -5))
for (i in 1:15) rect(X[i], -5, X[i + 1], 22, col = GradiantColors[i], lty = 0)
points(log2(LibrarySize), col = alpha('darkgoldenrod1', .8), pch = 15, cex = 1.8, lwd = 1.5)
lines(smooth.spline(c(1:48), log2(LibrarySize), df = 40), col = 'darkgoldenrod1', lwd = 4)
points(log2(Mean_HousekeepingGenes) + 3, cex = 2, col = alpha('red', .6), pch = 19)
lines(smooth.spline(c(1:48), log2(Mean_HousekeepingGenes) + 3, df = 40), col = 'red', lwd = 4)
points(log2(Mean_PositiveControlProbes) - 1, col = alpha('cyan2', .8), pch = 18, cex = 2.5)
lines(smooth.spline(c(1:48), log2(Mean_PositiveControlProbes) - 1, df = 40), col = 'cyan1', lwd = 4)
points(log2(Mean_NegativeControlProbes), col = alpha('green2', .6), pch = 17, cex = 1.8)
lines(smooth.spline(c(1:48), log2(Mean_NegativeControlProbes), df = 40), col = 'green2', lwd = 4)
box(lwd = 5, bty = 'l')
rm(X, GradiantColors, Mean_NegativeControlProbes,
   Mean_PositiveControlProbes, Mean_HousekeepingGenes,
   LibrarySize)
# Add legend with reduced size and adjusted position
legend('topright', legend = c('Library Size', 'Housekeeping Genes', 'Positive Control Probes', 'Negative Control Probes'), 
       col = c('darkgoldenrod1', 'red', 'cyan2', 'green2'), 
       pch = c(15, 19, 18, 17), 
       cex = 1.2,  # Adjust the value to reduce the size of the legend text
       lwd = c(1.5, 4, 4, 4),
       bty = 'n', 
       inset = c(0, 0.02))  # Adjust the values to adjust the position of the legend

#########  Nanostring normalization - uisng NanoStringNorm R package #######

#***** Housekeeping genes - 3 genes
Selected_HKgenes <- Nano_ExpressionMatrix$Name[ 38:40]
Nano_ExpressionMatrix$Code.Class[Nano_ExpressionMatrix$Name %in% Selected_HKgenes] <- "Housekeeping"
table(Nano_ExpressionMatrix$Code.Class)
# Endogenous  Housekeeping  Negative     Positive
# 38           3            8            6
rm(Selected_HKgenes)

#***** Nanostring normalization
### These options recommended by Nanostring and are commonly used:
# CodeCount = geo.mean
# Background = mean.2sd
# SampleContent = housekeeping.geo.mean
dim(Nano_ExpressionMatrix) ##  54 51
Nanostring_normalized <- NanoStringNorm(x = Nano_ExpressionMatrix,
                                        CodeCount = 'geo.mean',
                                        Background = "mean.2sd",
                                        SampleContent = 'housekeeping.geo.mean',
                                        round.values = FALSE,
                                        take.log = FALSE ,
                                        return.matrix.of.endogenous.probes = TRUE)

sum(is.na(Nanostring_normalized)) ## 0
NanostringNormalized <- as.matrix(log2(Nanostring_normalized + 1))
all(colnames(NanostringNormalized) == Nano_SampleInfo$SampleNames) # TRUE

#***** RLE plots - Figure 1 A - nCounter normalization
par(mar = c(6.5,6.5,2.3,0))
boxplot(NanostringNormalized - rowMedians(NanostringNormalized), main = '', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', ylim = c(-4.1,4.1),
        outline = FALSE, names = FALSE, frame = FALSE, whisklty = 3, whisklwd = 1.5, staplelty = 1, notch = FALSE, boxlwd = 2,
        staplelwd = 0 , boxcol = Color_Batches[factor(Nano_SampleInfo$Cartridges)],
        border = Color_Batches[factor(Nano_SampleInfo$Cartridges)], col = 'gray87')
box(lwd = 7, bty = 'l')
title('nCounter normalized', line = -2, cex.main = 3.5)
Median_Nano <- apply(NanostringNormalized - rowMedians(NanostringNormalized), 2, median)
points(c(1:ncol(NanostringNormalized)), Median_Nano, col = Color_Batches[factor(Nano_SampleInfo$Cartridges)], pch = 19, cex = 1.2)
axis(2, mgp = c(3.5, .9 ,0), lwd.ticks=6, las=1, cex.axis=3)
mtext('RLE', 2, line = 3.5, cex = 3.5)
abline(h = 0, col = 'black', lwd = 5, lty = 2)
par(lwd = 3)
axis.break(2, -4.2, style = 'zigzag', brw = .02)
rm(Median_Nano)

#***** Creating replicate matrix
length(Nano_SampleInfo$Patient.barcodes) # 162
length(unique(Nano_SampleInfo$Patient.barcodes)) # 135

ReplicateMatrix <- ruv::replicate.matrix(Nano_SampleInfo$Patient.barcodes)
dim(ReplicateMatrix) # 162 133

### Making sure that every row has got only one number
par(mfrow = c(2,1))
barplot(colSums(ReplicateMatrix))
barplot(rowSums(ReplicateMatrix))
par(mfrow = c(1,1))

#***** Finding control genes
### Step 1: Using all genes as a set of negative control genes
# Performing RUVIII
dataRUV <- t(log2(Nano_RawCounts[1:38, ]))
RUVcorrected <- ruv::RUVIII(Y = dataRUV, M = ReplicateMatrix, ctl = c(1:38))
RUVcorrected <- t(RUVcorrected)

### Step 2: Selecting the most stable genes from step 1
LowVarGenes <- apply(RUVcorrected, 1, var)
ControlGenes <- which(LowVarGenes < .5)
length(ControlGenes) # 492 genes


### Step3: Performing RUV-III using the ControlGenes set
RUVcorrected <- RUVIII(Y = dataRUV, M = ReplicateMatrix, ctl = ControlGenes, k = 6)
RUVcorrected <- t(RUVcorrected)
dim(RUVcorrected)


#***** RLE plots - Figure 1 A, RUV-III normalization



boxplot(RUVcorrected - rowMedians(RUVcorrected), main = '', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', ylim = c(-4.1, 4.1),
        outline = FALSE, names = FALSE, frame = FALSE, whisklty = 3, whisklwd = 1.5, staplelty = 1, notch = FALSE, boxlwd = 2,
        staplelwd = 0, boxcol = Color_Batches[factor(Nano_SampleInfo$Cartridges)],
        border = Color_Batches[factor(Nano_SampleInfo$Cartridges)], col = 'gray87')



par(mar = c(6.5,6.5,2.3,0))
boxplot(RUVcorrected - rowMedians(RUVcorrected), main = '', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', ylim = c(-4.1,4.1),
        outline = FALSE, names = FALSE, frame = FALSE, whisklty = 3, whisklwd = 1.5, staplelty = 1, notch = FALSE, boxlwd = 2,
        staplelwd = 0 , boxcol = Color_Batches[factor(Nano_SampleInfo$Cartridges)],
        border = Color_Batches[factor(Nano_SampleInfo$Cartridges)], col = 'gray87')
box(lwd = 7, bty = 'l')
title('RUV-III normalized', line = -2, cex.main = 3.5)
Median_RUV <- apply(RUVcorrected - rowMedians(RUVcorrected), 2, median)
points(c(1:ncol(RUVcorrected)), Median_RUV, col = Color_Batches[factor(Nano_SampleInfo$Cartridges)], pch = 19, cex = 1.2)
axis(1, cex.axis = 1, at = c(1, seq(20,166,20)), cex.axis = 2.5, lwd.ticks = 6,  mgp = c(3.5,1.6,0))
axis(2, mgp = c(3.5,.9,0), lwd.ticks = 6, las = 1, cex.axis = 3)
mtext(expression(paste('Samples','(', 'n'[samples], '=', '162', ')')),1 ,line = 5, cex = 3)
mtext('RLE', 2, line = 3.5, cex = 3.5)
abline(h = 0, col = 'black', lwd = 5, lty = 2)
```r



