### Tximeta.R

## Author: Stephanie Hicks, with modifications by Christine Muheim (March 2022) and Alexander Popescu (June 2023)

## Use Tximeta to read in quant files from Salmon quant into Summarized Experiment files and output matrix of gene counts for downstream analysis


### Packages ----

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(tximeta)  # 1.18.0
library(SummarizedExperiment)  # 1.30.1
library(org.Mm.eg.db)  # 3.17.0
library(here)  # 1.0.1

here::i_am("AP_SDRS_Tximeta.R")

sink('sessionInfo_Tximeta.txt')
sessionInfo()
sink()


### Create Linked Transcriptome ----

# update paths for index and quant files as needed
index_dir <- here("Salmon/Index", "gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys")
fasta_path <- here("Salmon/Index", "gencode.vM28.transcripts.fa.gz")
gtf_path <- here("Salmon/Index", "gencode.vM28.annotation.gtf.gz")
json_file <- here("Salmon/Index", paste0(basename(index_dir), ".json"))

# index is added to the cache automatically
makeLinkedTxome(indexDir = index_dir, 
                source = "GENCODE", organism = "Mus musculus", 
                release = "M28", genome = "GRCm39", 
                fasta = fasta_path,
                gtf = gtf_path,
                write = TRUE, jsonFile = json_file)


### Build Sample Table ----

all_files <- list.files(here("Salmon/Results"))

# modify common identifier as needed
all_files <- stringr::str_subset(all_files, "WT")
file_paths <- here("Salmon/Results", all_files, "quant.sf")

coldata <- data.frame(files = file_paths,
                      names = stringr::str_sub(all_files, start = 5, end = -7),
                      stringsAsFactors = FALSE)

coldata
all(file.exists(coldata$files)) # TRUE


### Run Tximeta ----

## Given the coldata, Tximeta can now import the quant files into a Summarized Experiment object

## Note that Tximeta must be run separately for the gene level and the transcript level (not analyzed here)

# For the gene level, drop the inferential replicates from Salmon
# Set 'useHub = FALSE' to avoid downloading and using the reference transcriptome from AnnotationHub
# We do not specify 'countsFromAbundance = scaledTPM' because that is redundant with UQ normalization (see RUV.R)

se <- tximeta(coldata, type = "salmon", txOut = TRUE, useHub = FALSE)
seG <- tximeta(coldata, type = "salmon", txOut = TRUE, dropInfReps = TRUE, useHub = FALSE)

# add gene IDs to se object
se <- addIds(se, "SYMBOL")
mcols(se)

colData(se)
assayNames(se)
rowRanges(se)

# summarize to gene-level counts
gse <- summarizeToGene(seG)
gse <- addIds(gse, "SYMBOL", gene = TRUE)

colData(gse)
assayNames(gse)
rowRanges(gse)


### Counts Matrix ----

## Use the Summarized Experiment objects from the previous step

# transcript-level
transcriptCounts <- assays(se)[["counts"]]
write.table(transcriptCounts,
            file = here("Data", "se_HCSDRS_WT_SleepIntegration_salmon.txt"),
            sep = "\t")

# gene-level
geneCounts <- assays(gse)[["counts"]]
write.table(geneCounts,
            file = here("Data", "gse_HCSDRS_WT_SleepIntegration_salmon.txt"),
            sep = "\t")
