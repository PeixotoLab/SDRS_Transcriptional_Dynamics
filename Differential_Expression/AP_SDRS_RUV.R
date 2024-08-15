### RUV.R
## Author: Alexander Popescu (June 2023)

## Adapted from Muheim et al., 2023; Ingiosi et al., 2019

## Pipeline for filtering, RUVs normalization (RUVSeq), and differential expression analysis (edgeR) of the matrix of gene counts from Tximeta


### Packages ----

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(limma)  # 3.56.2
library(statmod)  # 1.5.0
library(edgeR)  # 3.42.4
library(EDASeq)  # 2.34.0
library(RUVSeq)  # 1.34.0
library(RColorBrewer)  # 1.1-3
library("org.Mm.eg.db")  # 3.17.0
library(dplyr)  # 1.1.2
library(readxl)  # 1.4.2
library(here)  # 1.0.1

here::i_am("AP_SDRS_RUV.R")

sink('sessionInfo_RUV.txt')
sessionInfo()
sink()


### Data ----

# Load counts matrix (text file) created at the end of Tximeta
sdrs <- read.table(file = here("Data", "gse_HCSDRS_WT_SleepIntegration_salmon.txt"),
                   row.names = 1, header = TRUE)

head(sdrs)
colnames(sdrs) # HC, RS, SD

# Remove gene ID version numbers (after ".")
rownames(sdrs) <- lapply(rownames(sdrs), sub, pattern = "\\.\\d+$", replacement = "")

# If needed, rename columns to make plots more readable and reorder them #

head(sdrs)
colnames(sdrs)
dim(sdrs) # 54307 x 42

write.table(sdrs, file = here("Data", "gse_HCSDRS_WT_SleepIntegration_salmon_trunc.txt"),
            sep = "\t")


### Filter ----

## We want to discard genes that are not expressed or have very low read counts ##

# We made a function that states, firstly, the minimum number of reads we want and, secondly, the minimum number of samples that need to have the minimum number of reads

x <- factor(rep(c("HC3", "HC5", "HC6", "HC7", "HC12", "SD3", "SD5", "SD6", "RS2", "RS6"),
                c(3, 5, 5, 5, 5, 3, 5, 3, 5, 3)),
            levels = c("HC3", "HC5", "HC6", "HC7", "HC12", "SD3", "SD5", "SD6", "RS2", "RS6"))

names(x) <- colnames(sdrs)

## We filter out genes present less than 10x across more than 5 samples ##

# These numbers will depend on your dataset and the number of conditions and replicates
# Check how many positive controls remain in your dataset with different parameters

filter <- apply(sdrs, 1, function(x) length(x[which(x > 10)]) > 5)
filtered <- as.matrix(sdrs)[filter, ]

head(filtered)
dim(filtered) # 18872 x 42

# Save list of expressed genes to use as a background for functional analysis
write.table(rownames(filtered),
            file = here("Data", "gse_HCSDRS_WT_SleepIntegration_salmon_background.txt"),
            sep = "\t")


### Positive Controls ----

## To assess the reproducibility of results from differential expression and compare the relative performances of UQ + RUV vs UQ only, we use:

# A set of positive controls assembled independent of brain region, technology, and lab after 5-6 hours of SD
# Genes previously reported to be differentially expressed (adj. p-value < 0.01) at either 2 or 6 hours of RS following 5–6 hours of SD in the prefrontal cortex

# These genes are from Additional File 2 of Gerstner et al., 2016, with Ensembl IDs and annotation updated to Ensembl 105
# See `Controls.R` script for additional details

file <- here("Data", "Supplemental Table S1.xlsx")

pos.sd.ctrls <- read_excel(file, sheet = 2)
pos.sd.ctrls <- pos.sd.ctrls$`ENSEMBL ID` # 677
sd.pos <- intersect(pos.sd.ctrls, rownames(filtered)) # 670

pos.rs2.ctrls <- read_excel(file, sheet = 3)
pos.rs2.ctrls <- pos.rs2.ctrls$`ENSEMBL ID` # 176
rs2.pos <- intersect(pos.rs2.ctrls, rownames(filtered)) # 170

pos.rs6.ctrls <- read_excel(file, sheet = 4)
pos.rs6.ctrls <- pos.rs6.ctrls$`ENSEMBL ID` # 45
rs6.pos <- intersect(pos.rs6.ctrls, rownames(filtered)) # 40


### UQ Normalization ----

## Upper Quartile normalization (EDASeq) is used to correct for library size and gives us a better recovery of positive controls than TMM (edgeR) does

uq <- betweenLaneNormalization(filtered, which = "upper")
dim(uq) # 18872 x 42

## PCA and RLE Plots ##

pchvec1 <- ifelse(grepl("HC7|RS2", colnames(uq)), 15,
               ifelse(grepl("HC5|SD5", colnames(uq)), 17, 19))

pchvec2 <- ifelse(grepl("HC7|RS2", levels(x)), 15,
                  ifelse(grepl("HC5|SD5", levels(x)), 17, 19))

colors <- factor(rep(c("#AAAAAA", "#888888", "#666666", "#444444", "#222222",
                       "#FF6666", "#FF3333", "#CC0000",
                       "#6495ED", "#336699"),
                     c(3, 5, 5, 5, 5, 3, 5, 3, 5, 3)),
                 levels = c("#AAAAAA", "#888888", "#666666", "#444444", "#222222",
                            "#FF6666", "#FF3333", "#CC0000",
                            "#6495ED", "#336699")) 

# RLE plots reveal confounders when the mean and the variance are not similar
plotRLE(uq, col = as.character(colors), outline = FALSE, las = 3, ylim = c(-0.8, 0.8),
        ylab = "Relative Log Expression", main = "Upper Quartile")

# PCA plots show principal components by which variance can be explained
plotPCA(uq, labels = FALSE, pch = pchvec1, col = as.character(colors),
        main = "Upper Quartile", cex = 1.25, xlim = c(-0.3, 0.3), ylim = c(-0.3, 0.3))

legend("topright",
       legend = levels(x),
       col = levels(colors),
       pch = pchvec2,
       bty = 'n',
       text.col = "black", 
       horiz = F,
       inset = c(0.025, 0.025))


### UQ DE Analysis ----

## We get differential expression results after UQ normalization only as a benchmark by which improvement with RUV can be assessed, not for the final gene lists

# A generalized linear model (GLM) is used to model the counts
# Quasi-likelihood F-test (QLF) is used to estimate the parameters in the GLM

# SD5 and SD6 are grouped together in the contrasts, as are HC5 and HC6, to allow for the removal of the lab effect (see RUV) and a single estimate of the effects of long SD

design <- model.matrix(~x - 1)
design

y <- DGEList(counts = filtered, group = x)
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateDisp(y, design, verbose = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

qlf1 <- glmQLFTest(fit, contrast = c(-1, 0, 0, 0, 0, 1, 0, 0, 0, 0)) # SD3 vs HC3
qlf2 <- glmQLFTest(fit, contrast = c(0, -1, -1, 0, 0, 0, 1, 1, 0, 0)) # SD5-6 vs HC5-6
qlf3 <- glmQLFTest(fit, contrast = c(0, 0, 0, -1, 0, 0, 0, 0, 1, 0)) # RS2 vs HC7
qlf4 <- glmQLFTest(fit, contrast = c(0, 0, 0, 0, -1, 0, 0, 0, 0, 1)) # RS6 vs HC12

summary(decideTests(qlf1))
summary(decideTests(qlf2))
summary(decideTests(qlf3))
summary(decideTests(qlf4))

topUQSD1 <- topTags(qlf1, n = Inf)$table
topUQSD2 <- topTags(qlf2, n = Inf)$table
topUQSD3 <- topTags(qlf3, n = Inf)$table
topUQSD4 <- topTags(qlf4, n = Inf)$table

## Number of Differentially Expressed Genes (DEGs) ##

# In addition to the FDR cutoff, we applied an expression cutoff (log2 CPM < 0)

topUQSD1.cut <- topUQSD1[topUQSD1$FDR < 0.05 & topUQSD1$logCPM > 0, ]
topUQSD2.cut <- topUQSD2[topUQSD2$FDR < 0.05 & topUQSD2$logCPM > 0, ]
topUQSD3.cut <- topUQSD3[topUQSD3$FDR < 0.05 & topUQSD3$logCPM > 0, ]
topUQSD4.cut <- topUQSD4[topUQSD4$FDR < 0.05 & topUQSD4$logCPM > 0, ]

nrow(topUQSD1.cut)  # 2080
sum(topUQSD1.cut$logFC > 0)  # 852
sum(topUQSD1.cut$logFC < 0)  # 1228

nrow(topUQSD2.cut)  # 7041
sum(topUQSD2.cut$logFC > 0)  # 2893
sum(topUQSD2.cut$logFC < 0)  # 4148

nrow(topUQSD3.cut)  # 3288
sum(topUQSD3.cut$logFC > 0)  # 1443
sum(topUQSD3.cut$logFC < 0)  # 1845

nrow(topUQSD4.cut)  # 1511
sum(topUQSD4.cut$logFC > 0)  # 877
sum(topUQSD4.cut$logFC < 0)  # 634

## Recovery of Positive Controls ##

sum(rownames(topUQSD2.cut) %in% sd.pos)
(sum(rownames(topUQSD2.cut) %in% sd.pos)/(length(sd.pos))) * 100
# 582 (86.9%)

sum(rownames(topUQSD3.cut) %in% rs2.pos)
(sum(rownames(topUQSD3.cut) %in% rs2.pos)/(length(rs2.pos))) * 100
# 135 (79.4%)

sum(rownames(topUQSD4.cut) %in% rs6.pos)
(sum(rownames(topUQSD4.cut) %in% rs6.pos)/(length(rs6.pos))) * 100
# 19 (47.5%)


### UQ Volcano Plots ----

x1 <- expression(log[2] ~ fold ~ change)
y1 <- expression(-log[10] ~ pvalue)

# Volcano Plot 1: SD3 vs HC3
plot(topUQSD1[, 1], -log10(topUQSD1$PValue), pch = 20, col = "#DDDDDD", cex = 0.6,
     main = "Upper Quartile SD3 vs HC3", ylab = y1, xlab = x1,
     ylim = c(0, 25), xlim = c(-6, 6))

# de
de <- rownames(topUQSD1.cut)
points(topUQSD1[de, 1], -log10(topUQSD1[de, "PValue"]), pch = 20, col = "#B8B8B8",
       cex = 0.6, lwd = 2)

# Volcano Plot 2: SD5-6 vs HC5-6
plot(topUQSD2[, 1], -log10(topUQSD2$PValue), pch = 20, col = "#DDDDDD", cex = 0.6,
     main = "Upper Quartile SD5-6 vs HC5-6", ylab = y1, xlab = x1,
     ylim = c(0, 25), xlim = c(-6, 6))

# de
de <- rownames(topUQSD2.cut)
points(topUQSD2[de, 1], -log10(topUQSD2[de, "PValue"]), pch = 20, col = "#B8B8B8",
       cex = 0.6, lwd = 2)

# positive controls
uq.sd.pos <- intersect(de, sd.pos)
points(topUQSD2[uq.sd.pos, 1], -log10(topUQSD2[uq.sd.pos, "PValue"]), pch = 20,
       col = "#666666", cex = 0.6, lwd = 2)

# Volcano Plot 3: RS2 vs HC7
plot(topUQSD3[, 1], -log10(topUQSD3$PValue), pch = 20, col = "#DDDDDD", cex = 0.6,
     main = "Upper Quartile RS2 vs HC7", ylab = y1, xlab = x1,
     ylim = c(0, 25), xlim = c(-6, 6))

# de
de <- rownames(topUQSD3.cut)
points(topUQSD3[de, 1], -log10(topUQSD3[de, "PValue"]), pch = 20, col = "#B8B8B8",
       cex = 0.6, lwd = 2)

# positive controls
uq.rs2.pos <- intersect(de, rs2.pos)
points(topUQSD3[uq.rs2.pos, 1], -log10(topUQSD3[uq.rs2.pos, "PValue"]), pch = 20,
       col = "#666666", cex = 0.6, lwd = 2)

# Volcano Plot 4: RS6 vs HC12
plot(topUQSD4[, 1], -log10(topUQSD4$PValue), pch = 20, col = "#DDDDDD", cex = 0.6,
     main = "Upper Quartile RS6 vs HC12", ylab = y1, xlab = x1,
     ylim = c(0, 25), xlim = c(-6, 6))

# de
de <- rownames(topUQSD4.cut)
points(topUQSD4[de, 1], -log10(topUQSD4[de, "PValue"]), pch = 20, col = "#B8B8B8",
       cex = 0.6, lwd = 2)

# positive controls
uq.rs6.pos <- intersect(de, rs6.pos)
points(topUQSD4[uq.rs6.pos, 1], -log10(topUQSD4[uq.rs6.pos, "PValue"]), pch = 20,
       col = "#666666", cex = 0.6, lwd = 2)


### RUV Normalization (k = 3) ----

## Remove unwanted variation (RUV) normalization (Risso et al., 2014) is used to estimate factors of unwanted variation for the matrix of gene counts following filtering and UQ normalization

## RUVs uses technical replicates or negative control samples (for which the covariates of interest are constant) to estimate the factors of unwanted variation W

# Each row in the 'groups' matrix represents one replicate group and each column represents a replicate in the group, hence the padding with -1's

# Grouping HC5/HC6 together and SD5/SD6 together is a reasonable biological assumption and allows us to remove the "lab effect" (differences in sleep deprivation procedures across labs etc.) since RUVs can only account for variation that occurs within replicate groups

groups <- matrix(data = c(1:3, rep(-1, 7), 4:13, 14:18, rep(-1, 5), 19:23, rep(-1, 5),
                          24:26, rep(-1, 7), 27:34, rep(-1, 2),
                          35:39, rep(-1, 5), 40:42, rep(-1, 7)),
                 nrow = 8, byrow = TRUE)

groups

## Although RUVs performs well even when all of the genes in the filtered matrix are used, it is preferable to use a list of negative control or housekeeping genes

# Our list of genes less likely to be affected by 5-6 hours of sleep deprivation in the prefrontal cortex comes from Gerstner et al., 2016, see 'Controls.R' and Manuscript for additional details

neg.sd.ctrls <- read_excel(file, sheet = 1)
neg.sd.ctrls <- neg.sd.ctrls$`ENSEMBL ID` # 5494

neg <- intersect(neg.sd.ctrls, rownames(uq)) # 3422
# if use you use all of the genes in the filtered matrix: neg <- rownames(filtered)

## RUVs normalization with k = 3 (k is the number of factors of unwanted variation) maximized the number of differentially expressed genes and the number of positive controls we recover

s <- RUVs(x = uq, cIdx = neg, scIdx = groups, k = 3)

## RLE and PCA plots ##

plotRLE(s$normalizedCounts, col = as.character(colors), outline = FALSE, las = 3,
        ylim = c(-0.6, 0.6), ylab = "Relative Log Expression", main = "RUV (k = 3)")

plotPCA(s$normalizedCounts, labels = FALSE, pch = pchvec1, col = as.character(colors),
        main = "RUV (k = 3)", cex = 1.25, xlim = c(-0.3, 0.3), ylim = c(-0.4, 0.4))

legend("topright",
       legend = levels(x),
       col = levels(colors),
       pch = pchvec2,
       bty = 'n',
       text.col = "black", 
       horiz = F,
       inset = c(0.01, 0.025))


### RUV DE Analysis ----

# A mixed model (GLM plus unwanted factors estimated by RUV considered as additional covariates) is used to model the counts

design <- model.matrix(~x + s$W - 1)
design

y <- DGEList(counts = filtered, group = x)
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateDisp(y, design, verbose = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
fit$design

# For each k, add a 0 to the end of contrast matrix (check fit$design matrix to see why)

qlf1 <- glmQLFTest(fit, contrast = c(-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0))
qlf2 <- glmQLFTest(fit, contrast = c(0, -1, -1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0))
qlf3 <- glmQLFTest(fit, contrast = c(0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0))
qlf4 <- glmQLFTest(fit, contrast = c(0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0))

summary(decideTests(qlf1))
# -1*xHC3 1*xSD3
# Down             1542
# NotSig          16354
# Up                976

summary(decideTests(qlf2))
# -1*xHC5 -1*xHC6 1*xSD5 1*xSD6
# Down                            5164
# NotSig                         10320
# Up                              3388

summary(decideTests(qlf3))
# -1*xHC7 1*xRS2
# Down             2301
# NotSig          14668
# Up               1903

summary(decideTests(qlf4))
# -1*xHC12 1*xRS6
# Down               953
# NotSig           16740
# Up                1179

topSD1 <- topTags(qlf1, n = Inf)$table
topSD2 <- topTags(qlf2, n = Inf)$table
topSD3 <- topTags(qlf3, n = Inf)$table
topSD4 <- topTags(qlf4, n = Inf)$table

## Number of DEGs (FDR < 0.05 and log2 CPM > 0 cutoffs) ##

topSD1.cut <- topSD1[topSD1$FDR < 0.05 & topSD1$logCPM > 0, ]
topSD2.cut <- topSD2[topSD2$FDR < 0.05 & topSD2$logCPM > 0, ]
topSD3.cut <- topSD3[topSD3$FDR < 0.05 & topSD3$logCPM > 0, ]
topSD4.cut <- topSD4[topSD4$FDR < 0.05 & topSD4$logCPM > 0, ]

nrow(topSD1.cut)  # 2315
sum(topSD1.cut$logFC > 0)  # 942
sum(topSD1.cut$logFC < 0)  # 1373

nrow(topSD2.cut)  # 7493
sum(topSD2.cut$logFC > 0)  # 3103
sum(topSD2.cut$logFC < 0)  # 4390

nrow(topSD3.cut)  # 3908
sum(topSD3.cut$logFC > 0)  # 1785
sum(topSD3.cut$logFC < 0)  # 2123

nrow(topSD4.cut)  # 1989
sum(topSD4.cut$logFC > 0)  # 1104
sum(topSD4.cut$logFC < 0)  # 885

## Recovery of Controls ##

sum(rownames(topSD2.cut) %in% sd.pos)
(sum(rownames(topSD2.cut) %in% sd.pos)/length(sd.pos)) * 100
# 593 (88.5%)

sum(rownames(topSD2.cut) %in% neg)
(sum(rownames(topSD2.cut) %in% neg)/length(neg)) * 100
# 1093 (31.9%)

sum(rownames(topSD3.cut) %in% rs2.pos)
(sum(rownames(topSD3.cut) %in% rs2.pos)/length(rs2.pos)) * 100
# 138 (81.2%)

sum(rownames(topSD4.cut) %in% rs6.pos)
(sum(rownames(topSD4.cut) %in% rs6.pos)/length(rs6.pos)) * 100
# 21 (52.5%)

## Full Lists ##

write.table(x = topSD1, file = here("Data", "SDRS_Full_Gene_SD3_HC3_k=3.txt"), sep = "\t")
write.table(x = topSD2, file = here("Data", "SDRS_Full_Gene_SD56_HC56_k=3.txt"), sep = "\t")
write.table(x = topSD3, file = here("Data", "SDRS_Full_Gene_RS2_HC7_k=3.txt"), sep = "\t")
write.table(x = topSD4, file = here("Data", "SDRS_Full_Gene_RS6_HC12_k=3.txt"), sep = "\t")

## DEG Lists ##

write.table(x = topSD1.cut, file = here("Data", "SDRS_0.05_Gene_SD3_HC3_k=3.txt"), sep = "\t")
write.table(x = topSD2.cut, file = here("Data", "SDRS_0.05_Gene_SD56_HC56_k=3.txt"), sep = "\t")
write.table(x = topSD3.cut, file = here("Data", "SDRS_0.05_Gene_RS2_HC7_k=3.txt"), sep = "\t")
write.table(x = topSD4.cut, file = here("Data", "SDRS_0.05_Gene_RS6_HC12_k=3.txt"), sep = "\t")

topSD1.cut.up <- topSD1.cut[topSD1.cut$logFC > 0, ]
topSD1.cut.down <- topSD1.cut[topSD1.cut$logFC < 0, ]
topSD2.cut.up <- topSD2.cut[topSD2.cut$logFC > 0, ]
topSD2.cut.down <- topSD2.cut[topSD2.cut$logFC < 0, ]
topSD3.cut.up <- topSD3.cut[topSD3.cut$logFC > 0, ]
topSD3.cut.down <- topSD3.cut[topSD3.cut$logFC < 0, ]
topSD4.cut.up <- topSD4.cut[topSD4.cut$logFC > 0, ]
topSD4.cut.down <- topSD4.cut[topSD4.cut$logFC < 0, ]

## DEG Lists (Up & Down) ##

write.table(x = topSD1.cut.up, file = here("Data", "SDRS_0.05_Gene_SD3_HC3_k=3_UP.txt"),
            sep = "\t")
write.table(x = topSD1.cut.down, file = here("Data", "SDRS_0.05_Gene_SD3_HC3_k=3_DOWN.txt"),
            sep = "\t")
write.table(x = topSD2.cut.up, file = here("Data", "SDRS_0.05_Gene_SD56_HC56_k=3_UP.txt"),
            sep = "\t")
write.table(x = topSD2.cut.down, file = here("Data", "SDRS_0.05_Gene_SD56_HC56_k=3_DOWN.txt"),
            sep = "\t")
write.table(x = topSD3.cut.up, file = here("Data", "SDRS_0.05_Gene_RS2_HC7_k=3_UP.txt"),
            sep = "\t")
write.table(x = topSD3.cut.down, file = here("Data", "SDRS_0.05_Gene_RS2_HC7_k=3_DOWN.txt"),
            sep = "\t")
write.table(x = topSD4.cut.up, file = here("Data", "SDRS_0.05_Gene_RS6_HC12_k=3_UP.txt"),
            sep = "\t")
write.table(x = topSD4.cut.down, file = here("Data", "SDRS_0.05_Gene_RS6_HC12_k=3_DOWN.txt"),
            sep = "\t")


### RUV Volcano Plots ----

## Histograms of P-value Distributions ##

hist(topSD1$PValue, main = "RUV k = 3", xlab = "p-value", breaks = 100, ylim = c(0, 3000))
hist(topSD2$PValue, main = "RUV k = 3", xlab = "p-value", breaks = 100, ylim = c(0, 8000))
hist(topSD3$PValue, main = "RUV k = 3", xlab = "p-value", breaks = 100, ylim = c(0, 4500))
hist(topSD4$PValue, main = "RUV k = 3", xlab = "p-value", breaks = 100, ylim = c(0, 3000))

## Volcano Plots ##

# Volcano Plot 1: SD3 vs HC3
plot(topSD1[, 1], -log10(topSD1$PValue), pch = 20, col = "#DDDDDD", cex = 0.6,
     main = "RUV SD3 vs HC3", ylab = y1, xlab = x1,
     ylim = c(0, 25), xlim = c(-6, 6))

# de
de <- rownames(topSD1.cut)
points(topSD1[de, 1], -log10(topSD1[de, "PValue"]), pch = 20, col = "#B8B8B8",
       cex = 0.6, lwd = 2)

# Volcano Plot 2: SD5-6 vs HC5-6
plot(topSD2[, 1], -log10(topSD2$PValue), pch = 20, col = "#DDDDDD", cex = 0.6,
     main = "RUV SD5-6 vs HC5-6", ylab = y1, xlab = x1,
     ylim = c(0, 25), xlim = c(-6, 6))

# de
de <- rownames(topSD2.cut)
points(topSD2[de, 1], -log10(topSD2[de, "PValue"]), pch = 20, col = "#B8B8B8",
       cex = 0.6, lwd = 2)

# positive controls
ruv.sd.pos <- intersect(de, sd.pos)
points(topSD2[ruv.sd.pos, 1], -log10(topSD2[ruv.sd.pos, "PValue"]), pch = 20,
       col = "#666666", cex = 0.6, lwd = 2)

# Volcano Plot 3: RS2 vs HC7
plot(topSD3[, 1], -log10(topSD3$PValue), pch = 20, col = "#DDDDDD", cex = 0.6,
     main = "RUV RS2 vs HC7", ylab = y1, xlab = x1,
     ylim = c(0, 25), xlim = c(-6, 6))

# de
de <- rownames(topSD3.cut)
points(topSD3[de, 1], -log10(topSD3[de, "PValue"]), pch = 20, col = "#B8B8B8",
       cex = 0.6, lwd = 2)

# positive controls
ruv.rs2.pos <- intersect(de, rs2.pos)
points(topSD3[ruv.rs2.pos, 1], -log10(topSD3[ruv.rs2.pos, "PValue"]), pch = 20,
       col = "#666666", cex = 0.6, lwd = 2)

# Volcano Plot 4: RS6 vs HC12
plot(topSD4[, 1], -log10(topSD4$PValue), pch = 20, col = "#DDDDDD", cex = 0.6,
     main = "RUV RS6 vs HC12", ylab = y1, xlab = x1,
     ylim = c(0, 25), xlim = c(-6, 6))

# de
de <- rownames(topSD4.cut)
points(topSD4[de, 1], -log10(topSD4[de, "PValue"]), pch = 20, col = "#B8B8B8",
       cex = 0.6, lwd = 2)

# positive controls
ruv.rs6.pos <- intersect(de, rs6.pos)
points(topSD4[ruv.rs6.pos, 1], -log10(topSD4[ruv.rs6.pos, "PValue"]), pch = 20, col =
         "#666666", cex = 0.6, lwd = 2)

