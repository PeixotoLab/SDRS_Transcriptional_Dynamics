### Controls.R
## Author: Alexander Popescu (June 2023)

## Adapted from the Bioconductor vignette for biomaRt "Accessing Ensembl annotation with biomaRt" (Smith, Durinck, and Huber)

## For control gene lists from Gerstner et al. (2016), map Affymetrix probe IDs to Ensembl IDs (Ensembl release 105)

## We use Perl (see Annotate.pl) for annotation of Ensembl IDs but this can also be done with biomaRt


### Packages ----

library("biomaRt")  # 2.56.0
library(readxl)  # 1.4.2
library(here)  # 1.0.1

here::i_am("AP_SDRS_Controls.R")
path <- "Annotate_Files"

sink('sessionInfo_MapControls.txt')
sessionInfo()
sink()


### Determine Mart and Platform ----

listEnsembl() # genes

## Use version corresponding to GRCm39 and Gencode M28 (version 105) used to run Salmon rather than the most recent version (version 109), see www.gencodegenes.org/mouse/releases.html

# with useMart() you must specify a host, which incorporates the version
mouse.anno <- useMart("ensembl", "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org")

# with useEnsembl() you must specify the version if it is not the most recent one
mouse.anno <- useEnsembl(biomart = "genes",
                         dataset = "mmusculus_gene_ensembl",
                         version = 105)

attributes <- listAttributes(mouse.anno)
dim(attributes) # 2983 x 3
filters <- listFilters(mouse.anno)
dim(filters) # 396 x 2

## Search for specific microarray platform GPL from GEO
# From GSE78215: [MoGene-2_1-st] Affymetrix Mouse Gene 2.1 ST Array [transcript (gene) version]
# Find matching platform name in mouse.anno@attributes or mouse.anno@filters

searchAttributes(mouse.anno, "mogene")
platformName <- attributes[103, ]
platform <- platformName$name # "affy_mogene_2_1_st_v1"

# For each control list, we use the platform as a filter in getBM()


### SD Positive Controls ----

file <- here(path, "BMC_Genomics_AF2_SD_PosCtrls.xlsx")
sd.pos.ctrls <- read_xlsx(file)
dim(sd.pos.ctrls) # 679 x 10

# see what probe ID types were used (can query with either)
colnames(sd.pos.ctrls)[2:3]

affyids <- sd.pos.ctrls$"Mouse Gene 2.1 ST probeset ID" # vector of probe IDs
length(affyids) # 679

## Use bioMaRt to map probe IDs to Ensembl IDs

# modify attributes based on desired query output
affyToEnsembl <- getBM(attributes = c(platform, "ensembl_gene_id", "external_gene_name"),
                       filters = platform,
                       values = affyids,
                       mart = mouse.anno)

dim(affyToEnsembl) # 718 x 3

# discard probes that map to multiple unique Ensembl IDs
duplicates <- duplicated(affyToEnsembl$affy_mogene_2_1_st_v1)
affyToEnsembl <- affyToEnsembl[!duplicates, ]
dim(affyToEnsembl) # 677 x 3

# 668 Ensembl IDs were unchanged with the mapping
length(intersect(affyToEnsembl$ensembl_gene_id, sd.pos.ctrls$`Ensembl ID`))

# export the mapped list of Ensembl IDs as a text file for annotation with Perl
write.table(affyToEnsembl$ensembl_gene_id,
            here(path, "SD_Pos_Ctrls_AF2.txt"),
            sep = "\t")


### SD Negative Controls ----

outfile <- "AdditionalFile4_BMC_Full.xlsx"

# file used for negative controls (genes less likely to be affected by the treatment, SD or RS) is available on GitHub
if (!file.exists(outfile)) {
  download.file("https://github.com/drighelli/peixoto/blob/master/data/controls/Additional%20File%204%20full%20list%20of%20BMC%20genomics%20SD&RS2.xlsx?raw=true",
                destfile = outfile) }

# Sheet 1: SD, Sheet 2: RS1, Sheet 3: RS2, Sheet 4: RS3, Sheet 5: RS6

sd.ctrls <- read_excel(outfile, sheet = 1)

# keep only genes with an adj p-value > 0.9
sd.ctrls <- sd.ctrls[order(sd.ctrls$adj.P.Val),]
sd.neg.ctrls <- sd.ctrls[sd.ctrls$adj.P.Val > 0.9, ]

# check the number of negative controls before remapping
sd.neg.ctrls.old <- sd.neg.ctrls$"Ensembl ID"
sd.neg.ctrls.old <- na.omit(sd.neg.ctrls.old)
sd.neg.ctrls.old <- unique(sd.neg.ctrls.old)
length(sd.neg.ctrls.old) # 4046

affyids <- sd.neg.ctrls$"Affymetrix Probeset ID"
length(affyids) # 7874

# use bioMaRt to map probe IDs to Ensembl IDs
affyToEnsembl <- getBM(attributes = c(platform, "ensembl_gene_id", "external_gene_name"),
                       filters = platform,
                       values = affyids,
                       mart = mouse.anno)

dim(affyToEnsembl) # 8139 x 3

# discard probes that map to multiple unique Ensembl IDs
duplicates <- duplicated(affyToEnsembl$affy_mogene_2_1_st_v1)
affyToEnsembl <- affyToEnsembl[!duplicates, ]
dim(affyToEnsembl) # 6043 x 3

# 3890 Ensembl IDs were unchanged with the mapping
length(intersect(affyToEnsembl$ensembl_gene_id, sd.neg.ctrls.old))

# subset Ensembl IDs
sd.neg.ctrls <- affyToEnsembl$ensembl_gene_id

# remove NAs and duplicates
sd.neg.ctrls <- na.omit(sd.neg.ctrls)
sd.neg.ctrls <- unique(sd.neg.ctrls)
length(sd.neg.ctrls) # 5559

write.table(sd.neg.ctrls,
            here(path, "SD_Neg_Ctrls_AF4.txt"),
            sep = "\t")


### RS2 Positive Controls ----

rs2.full <- read_excel(outfile, sheet = 3)

# keep only genes with an adj p-value < 0.01
rs2.full <- rs2.full[order(rs2.full$adj.P.Val), ]
rs2.pos.ctrls <- rs2.full[rs2.full$adj.P.Val < 0.01, ]

affyids <- rs2.pos.ctrls$"Affymetrix Probeset ID"
length(affyids) # 186

# use bioMaRt to map probe IDs to Ensembl IDs
affyToEnsembl <- getBM(attributes = c(platform, "ensembl_gene_id", "external_gene_name"),
                       filters = platform,
                       values = affyids,
                       mart = mouse.anno)

dim(affyToEnsembl) # 187 x 3

# discard probes that map to multiple unique Ensembl IDs
duplicates <- duplicated(affyToEnsembl$affy_mogene_2_1_st_v1)
affyToEnsembl <- affyToEnsembl[!duplicates, ]
dim(affyToEnsembl) # 181 x 3

# check the number of positive controls before remapping
rs2.pos.ctrls.old <- rs2.pos.ctrls$"Ensembl ID"
rs2.pos.ctrls.old <- na.omit(rs2.pos.ctrls.old)
rs2.pos.ctrls.old <- unique(rs2.pos.ctrls.old)
length(rs2.pos.ctrls.old) # 163

# 159 Ensembl IDs were unchanged with the mapping
length(intersect(affyToEnsembl$ensembl_gene_id, rs2.pos.ctrls.old))

# subset Ensembl IDs
rs2.pos.ctrls <- affyToEnsembl$ensembl_gene_id

# remove NAs and duplicates
rs2.pos.ctrls <- na.omit(rs2.pos.ctrls)
rs2.pos.ctrls <- unique(rs2.pos.ctrls)
length(rs2.pos.ctrls) # 176

write.table(rs2.pos.ctrls,
             here(path, "RS2_Pos_Ctrls_AF4.txt"),
             sep = "\t")


### RS6 Positive Controls ----

rs6.full <- read_excel(outfile, sheet = 5)

# keep only genes with an adj p-value < 0.01
rs6.full <- rs6.full[order(rs6.full$adj.P.Val), ]
rs6.pos.ctrls <- rs6.full[rs6.full$adj.P.Val < 0.01, ]

affyids <- rs6.pos.ctrls$"Affymetrix Probeset ID"
length(affyids) # 49

# use bioMaRt to map probe IDs to Ensembl IDs
affyToEnsembl <- getBM(attributes = c(platform, "ensembl_gene_id", "external_gene_name"),
                       filters = platform,
                       values = affyids,
                       mart = mouse.anno)

dim(affyToEnsembl) # 52 x 3

# discard probes that map to multiple unique Ensembl IDs
duplicates <- duplicated(affyToEnsembl$affy_mogene_2_1_st_v1)
affyToEnsembl <- affyToEnsembl[!duplicates, ]
dim(affyToEnsembl) # 46 x 3

# check the number of positive controls before remapping
rs6.pos.ctrls.old <- rs6.pos.ctrls$"Ensembl ID"
rs6.pos.ctrls.old <- na.omit(rs6.pos.ctrls.old)
rs6.pos.ctrls.old <- unique(rs6.pos.ctrls.old)
length(rs6.pos.ctrls.old) # 36

# 36 Ensembl IDs were unchanged with the mapping
length(intersect(affyToEnsembl$ensembl_gene_id, rs6.pos.ctrls.old))

# subset Ensembl IDs
rs6.pos.ctrls <- affyToEnsembl$ensembl_gene_id

# remove NAs and duplicates
rs6.pos.ctrls <- na.omit(rs6.pos.ctrls)
rs6.pos.ctrls <- unique(rs6.pos.ctrls)
length(rs6.pos.ctrls) # 45

write.table(rs6.pos.ctrls,
            here(path, "RS6_Pos_Ctrls_AF4.txt"),
            sep = "\t")


### RS2 Negative Controls ----

# keep only genes with an adj p-value > 0.9
rs2.neg.ctrls <- rs2.full[rs2.full$adj.P.Val > 0.9, ]

affyids <- rs2.neg.ctrls$"Affymetrix Probeset ID"
length(affyids) # 12052

# use bioMaRt to map probe IDs to Ensembl IDs
affyToEnsembl <- getBM(attributes = c(platform, "ensembl_gene_id", "external_gene_name"),
                       filters = platform,
                       values = affyids,
                       mart = mouse.anno)

dim(affyToEnsembl) # 12302 x 3

# discard probes that map to multiple unique Ensembl IDs
duplicates <- duplicated(affyToEnsembl$affy_mogene_2_1_st_v1)
affyToEnsembl <- affyToEnsembl[!duplicates, ]
dim(affyToEnsembl) # 9453 x 3

# check the number of negative controls before remapping
rs2.neg.ctrls.old <- rs2.neg.ctrls$"Ensembl ID"
rs2.neg.ctrls.old <- na.omit(rs2.neg.ctrls.old)
rs2.neg.ctrls.old <- unique(rs2.neg.ctrls.old)
length(rs2.neg.ctrls.old) # 6456

# 6202 Ensembl IDs were unchanged with the mapping
length(intersect(affyToEnsembl$ensembl_gene_id, rs2.neg.ctrls.old))

# subset Ensembl IDs
rs2.neg.ctrls <- affyToEnsembl$ensembl_gene_id

# remove NAs and duplicates
rs2.neg.ctrls <- na.omit(rs2.neg.ctrls)
rs2.neg.ctrls <- unique(rs2.neg.ctrls)
length(rs2.neg.ctrls) # 8593

write.table(rs2.neg.ctrls,
            here(path, "RS2_Neg_Ctrls_AF4.txt"),
            sep = "\t")


### RS6 Negative Controls ----

# keep only genes with an adj p-value > 0.9
rs6.neg.ctrls <- rs6.full[rs6.full$adj.P.Val > 0.9, ]

affyids <- rs6.neg.ctrls$"Affymetrix Probeset ID"
length(affyids) # 8290

# use bioMaRt to map probe IDs to Ensembl IDs
affyToEnsembl <- getBM(attributes = c(platform, "ensembl_gene_id", "external_gene_name"),
                       filters = platform,
                       values = affyids,
                       mart = mouse.anno)

dim(affyToEnsembl) # 9066 x 3

# discard probes that map to multiple unique Ensembl IDs
duplicates <- duplicated(affyToEnsembl$affy_mogene_2_1_st_v1)
affyToEnsembl <- affyToEnsembl[!duplicates, ]
dim(affyToEnsembl) # 6546 x 3

# check the number of negative controls before remapping
rs6.neg.ctrls.old <- rs6.neg.ctrls$"Ensembl ID"
rs6.neg.ctrls.old <- na.omit(rs6.neg.ctrls.old)
rs6.neg.ctrls.old <- unique(rs6.neg.ctrls.old)
length(rs6.neg.ctrls.old) # 4494

# 4337 Ensembl IDs were unchanged with the mapping
length(intersect(affyToEnsembl$ensembl_gene_id, rs6.neg.ctrls.old))

# subset Ensembl IDs
rs6.neg.ctrls <- affyToEnsembl$ensembl_gene_id

# remove NAs and duplicates
rs6.neg.ctrls <- na.omit(rs6.neg.ctrls)
rs6.neg.ctrls <- unique(rs6.neg.ctrls)
length(rs6.neg.ctrls) # 6011

write.table(rs6.neg.ctrls,
            here(path, "RS6_Neg_Ctrls_AF4.txt"),
            sep = "\t")
