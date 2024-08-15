### BubblePlot.R
## Author: Alexander Popescu (January 2024) with modifications by Caitlin Ottaway (July 2024)

## Adapted from Katie Ford's code for the scRNA-seq pipeline

## Use ggplot2 to visualize EASE score (Fischer Exact p-value) and fold enrichment of all clustered and unclustered terms from functional enrichment analysis of intersected lists


### Packages ----

library(ggplot2) # v 3.5.1
library(shades) # v 1.4.0
library(grid) # v 4.4.0
library(psych) # v 2.4.3
library(dplyr) # v 1.1.4
library(magrittr) # v 2.0.3
library(stringr) # v 1.5.1
library(readxl) # v 1.4.3

# sink('sessionInfo_BubblePlot.txt')
# sessionInfo()
# sink()


### Import Results from DAVID ----

SD56 <- read_excel("Data/Dec_2023_Functional_Annotation_Output.xlsx", sheet = 1)
SD3_SD56 <- read_excel("Data/Dec_2023_Functional_Annotation_Output.xlsx", sheet = 2)
FR_Union <- read_excel("Data/Dec_2023_Functional_Annotation_Output.xlsx", sheet = 3)

SD56_RS2 <- read_excel("Data/Dec_2023_Functional_Annotation_Output.xlsx", sheet = 4)
SD3_SD56_RS2 <- read_excel("Data/Dec_2023_Functional_Annotation_Output.xlsx", sheet = 5)
SR_Union <- read_excel("Data/Dec_2023_Functional_Annotation_Output.xlsx", sheet = 6)

SD56_RS2_RS6 <- read_excel("Data/Dec_2023_Functional_Annotation_Output.xlsx", sheet = 7)
SD3_SD56_RS2_RS6 <- read_excel("Data/Dec_2023_Functional_Annotation_Output.xlsx", sheet = 8)
NR_Union <- read_excel("Data/Dec_2023_Functional_Annotation_Output.xlsx", sheet = 9)

RS2 <- read_excel("Data/Dec_2023_Functional_Annotation_Output.xlsx", sheet = 10)
RS6 <- read_excel("Data/Dec_2023_Functional_Annotation_Output.xlsx", sheet = 11)


### Separate by Clustered/Unclustered and Upregulated/Downregulated ----

## Separate clustered (not always present) vs unclustered terms and separate by upregulated vs downregulated genes

## Remove headers, change the column names, and, for clustered terms, remove the line with cluster number and enrichment

SD3_SD56_Up_Clustered <- SD3_SD56[1:22, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

SD3_SD56_Up_Unclustered <- SD3_SD56[25:29, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

SD3_SD56_Down_Unclustered <- SD3_SD56[33:36, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

FR_Union_Up_Clustered <- FR_Union[1:25, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

FR_Union_Up_Unclustered <- FR_Union[28:36, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

FR_Union_Down_Clustered <- FR_Union[39:49, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

FR_Union_Down_Unclustered <- FR_Union[52:69, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

SD56_Up_Clustered <- SD56[1:20, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

SD56_Up_Unclustered <- SD56[23:31, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

SD56_Down_Clustered <- SD56[34:51, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

SD56_Down_Unclustered <- SD56[54:75, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))


SD3_SD56_RS2_Up_Clustered <- SD3_SD56_RS2[1:13, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

SD3_SD56_RS2_Up_Unclustered <- SD3_SD56_RS2[16:18, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

SD3_SD56_RS2_Down_Unclustered <- SD3_SD56_RS2[22:25, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

SR_Union_Up_Clustered <- SR_Union[1:18, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

SR_Union_Up_Unclustered <- SR_Union[21:25, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

SR_Union_Down_Clustered <- SR_Union[28:51, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

SR_Union_Down_Unclustered <- SR_Union[54:66, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

SD56_RS2_Up_Clustered <- SD56_RS2[1:15, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

SD56_RS2_Up_Unclustered <- SD56_RS2[18:22, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

SD56_RS2_Down_Clustered <- SD56_RS2[25:58, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

SD56_RS2_Down_Unclustered <- SD56_RS2[61:69, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))


SD3_SD56_RS2_RS6_Up_Clustered <- SD3_SD56_RS2_RS6[1:15, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

SD3_SD56_RS2_RS6_Up_Unclustered <- SD3_SD56_RS2_RS6[18:20, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

SD3_SD56_RS2_RS6_Down_Clustered <- SD3_SD56_RS2_RS6[23:28, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

SD3_SD56_RS2_RS6_Down_Unclustered <- SD3_SD56_RS2_RS6[31:32, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

NR_Union_Up_Clustered <- NR_Union[1:23, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

NR_Union_Up_Unclustered <- NR_Union[26:29, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

NR_Union_Down_Clustered <- NR_Union[32:36, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

NR_Union_Down_Unclustered <- NR_Union[39:44, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

SD56_RS2_RS6_Up_Unclustered <- SD56_RS2_RS6[2:8, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

SD56_RS2_RS6_Down_Clustered <- SD56_RS2_RS6[11:18, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

SD56_RS2_RS6_Down_Unclustered <- SD56_RS2_RS6[21:24, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))


RS2_Up_Clustered <- RS2[1:17, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

RS2_Up_Unclustered <- RS2[20:28, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

RS2_Down_Clustered <- RS2[31:35, ] %>%
  na.omit() %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR")) %>%
  subset(Category != "Category")

RS2_Down_Unclustered <- RS2[38:43, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

RS6_Up_Unclustered <- RS6[2:16, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))

RS6_Down_Unclustered <- RS6[20:21, ] %>%
  set_colnames(c("Category", "Term", "Gene_Count", "%", "p_Value", "Genes", "Gene_Names",
                 "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", "Bonferroni",
                 "Benjamini", "FDR"))


### Combine Clustered/Unclustered ----

## Combine clustered and unclustered terms for plotting and add columns to track order and clustering for terms

SD3_SD56_Up <- rbind(SD3_SD56_Up_Clustered, SD3_SD56_Up_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Up", "Cluster 2 Up", "Cluster 3 Up", "Unclustered Up"),
                     times = c(4, 4, 6, 5)))

SD3_SD56_Down <- SD3_SD56_Down_Unclustered %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep("Unclustered Down", 4))

FR_Union_Up <- rbind(FR_Union_Up_Clustered, FR_Union_Up_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Up", "Cluster 2 Up", "Cluster 3 Up", "Unclustered Up"),
                     times = c(5, 3, 9, 9)))

FR_Union_Down <- rbind(FR_Union_Down_Clustered, FR_Union_Down_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Down", "Cluster 2 Down", "Unclustered Down"),
                     times = c(3, 3, 18)))

SD56_Up <- rbind(SD56_Up_Clustered, SD56_Up_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Up", "Cluster 2 Up", "Cluster 3 Up", "Unclustered Up"),
                     times = c(5, 3, 4, 9)))

SD56_Down <- rbind(SD56_Down_Clustered, SD56_Down_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Down", "Cluster 2 Down", "Cluster 3 Down",
                       "Unclustered Down"),
                     times = c(3, 4, 3, 22)))

SD3_SD56_RS2_Up <- rbind(SD3_SD56_RS2_Up_Clustered, SD3_SD56_RS2_Up_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Up", "Cluster 2 Up", "Unclustered Up"),
                     times = c(4, 4, 3)))

SD3_SD56_RS2_Down <- SD3_SD56_RS2_Down_Unclustered %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep("Unclustered Down", 4))

SR_Union_Up <- rbind(SR_Union_Up_Clustered, SR_Union_Up_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Up", "Cluster 2 Up", "Cluster 3 Up", "Unclustered Up"),
                     times = c(3, 4, 3, 5)))

SR_Union_Down <- rbind(SR_Union_Down_Clustered, SR_Union_Down_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Down", "Cluster 2 Down", "Cluster 3 Down",
                       "Unclustered Down"),
                     times = c(10, 3, 3, 13)))

SD56_RS2_Up <- rbind(SD56_RS2_Up_Clustered, SD56_RS2_Up_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Up", "Cluster 2 Up", "Unclustered Up"),
                     times = c(6, 4, 5)))

SD56_RS2_Down <- rbind(SD56_RS2_Down_Clustered, SD56_RS2_Down_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Down", "Cluster 2 Down", "Cluster 3 Down", "Cluster 4 Down",
                       "Unclustered Down"),
                     times = c(11, 3, 6, 3, 9)))

SD3_SD56_RS2_RS6_Up <- rbind(SD3_SD56_RS2_RS6_Up_Clustered,
                             SD3_SD56_RS2_RS6_Up_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Up", "Cluster 2 Up", "Unclustered Up"),
                     times = c(3, 7, 3)))

SD3_SD56_RS2_RS6_Down <- rbind(SD3_SD56_RS2_RS6_Down_Clustered,
                               SD3_SD56_RS2_RS6_Down_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Down", "Unclustered Down"),
                     times = c(4, 2)))

NR_Union_Up <- rbind(NR_Union_Up_Clustered, NR_Union_Up_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Up", "Cluster 2 Up", "Unclustered Up"),
                     times = c(6, 12, 4)))

NR_Union_Down <- rbind(NR_Union_Down_Clustered, NR_Union_Down_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Down", "Unclustered Down"),
                     times = c(3, 6)))

SD56_RS2_RS6_Up <- SD56_RS2_RS6_Up_Unclustered %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep("Unclustered Up", 7))

SD56_RS2_RS6_Down <- rbind(SD56_RS2_RS6_Down_Clustered, SD56_RS2_RS6_Down_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Down", "Unclustered Down"),
                     times = c(6, 4)))

RS2_Up <- rbind(RS2_Up_Clustered, RS2_Up_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Up", "Cluster 2 Up", "Unclustered Up"),
                     times = c(4, 8, 9)))

RS2_Down <- rbind(RS2_Down_Clustered, RS2_Down_Unclustered) %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep(c("Cluster 1 Down", "Unclustered Down"),
                     times = c(3, 6)))

RS6_Up <- RS6_Up_Unclustered %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep("Unclustered Up", 15))

RS6_Down <- RS6_Down_Unclustered %>%
  mutate(Order = row_number()) %>%
  mutate(Group = rep("Unclustered Down", 2))


### Determine Unique Clusters ----

## Steps: subset cluster, select genes, remove commas/quotes/whitespace and duplicates

SD3_SD56_Up_Cluster1_Genes <- subset(SD3_SD56_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD3_SD56_Up_Cluster2_Genes <- subset(SD3_SD56_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD3_SD56_Up_Cluster3_Genes <- subset(SD3_SD56_Up, Group == "Cluster 3 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

FR_Union_Up_Cluster1_Genes <- subset(FR_Union_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

FR_Union_Up_Cluster2_Genes <- subset(FR_Union_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

FR_Union_Up_Cluster3_Genes <- subset(FR_Union_Up, Group == "Cluster 3 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

FR_Union_Down_Cluster1_Genes <- subset(FR_Union_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

FR_Union_Down_Cluster2_Genes <- subset(FR_Union_Down, Group == "Cluster 2 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD56_Up_Cluster1_Genes <- subset(SD56_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD56_Up_Cluster2_Genes <- subset(SD56_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD56_Up_Cluster3_Genes <- subset(SD56_Up, Group == "Cluster 3 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD56_Down_Cluster1_Genes <- subset(SD56_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD56_Down_Cluster2_Genes <- subset(SD56_Down, Group == "Cluster 2 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD56_Down_Cluster3_Genes <- subset(SD56_Down, Group == "Cluster 3 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD3_SD56_RS2_Up_Cluster1_Genes <- subset(SD3_SD56_RS2_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD3_SD56_RS2_Up_Cluster2_Genes <- subset(SD3_SD56_RS2_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SR_Union_Up_Cluster1_Genes <- subset(SR_Union_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SR_Union_Up_Cluster2_Genes <- subset(SR_Union_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SR_Union_Up_Cluster3_Genes <- subset(SR_Union_Up, Group == "Cluster 3 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SR_Union_Down_Cluster1_Genes <- subset(SR_Union_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SR_Union_Down_Cluster2_Genes <- subset(SR_Union_Down, Group == "Cluster 2 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SR_Union_Down_Cluster3_Genes <- subset(SR_Union_Down, Group == "Cluster 3 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD56_RS2_Up_Cluster1_Genes <- subset(SD56_RS2_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD56_RS2_Up_Cluster2_Genes <- subset(SD56_RS2_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD56_RS2_Down_Cluster1_Genes <- subset(SD56_RS2_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD56_RS2_Down_Cluster2_Genes <- subset(SD56_RS2_Down, Group == "Cluster 2 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD56_RS2_Down_Cluster3_Genes <- subset(SD56_RS2_Down, Group == "Cluster 3 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD56_RS2_Down_Cluster4_Genes <- subset(SD56_RS2_Down, Group == "Cluster 4 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD3_SD56_RS2_RS6_Up_Cluster1_Genes <- subset(SD3_SD56_RS2_RS6_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD3_SD56_RS2_RS6_Up_Cluster2_Genes <- subset(SD3_SD56_RS2_RS6_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD3_SD56_RS2_RS6_Down_Cluster1_Genes <- subset(SD3_SD56_RS2_RS6_Down,
                                               Group == "Cluster 1 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

NR_Union_Up_Cluster1_Genes <- subset(NR_Union_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

NR_Union_Up_Cluster2_Genes <- subset(NR_Union_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

NR_Union_Down_Cluster1_Genes <- subset(NR_Union_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

SD56_RS2_RS6_Down_Cluster1_Genes <- subset(SD56_RS2_RS6_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

RS2_Up_Cluster1_Genes <- subset(RS2_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

RS2_Up_Cluster2_Genes <- subset(RS2_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()

RS2_Down_Cluster1_Genes <- subset(RS2_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(Gene_Names) %>% pull() %>% str_split(",") %>% unlist() %>%
  noquote() %>% trimws() %>% extract(. != "") %>% unique()


### Condense Clusters ----

## Steps: extract fold enrichment, transform values into numerics, and use psych package to calculate geometric mean of fold enrichment and p-values

SD3_SD56_Up_Cluster1_FGM <- subset(SD3_SD56_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD3_SD56_Up_Cluster1_PGM <- subset(SD3_SD56_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD3_SD56_Up_Cluster2_FGM <- subset(SD3_SD56_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD3_SD56_Up_Cluster2_PGM <- subset(SD3_SD56_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD3_SD56_Up_Cluster3_FGM <- subset(SD3_SD56_Up, Group == "Cluster 3 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD3_SD56_Up_Cluster3_PGM <- subset(SD3_SD56_Up, Group == "Cluster 3 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

FR_Union_Up_Cluster1_FGM <- subset(FR_Union_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

FR_Union_Up_Cluster1_PGM <- subset(FR_Union_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

FR_Union_Up_Cluster2_FGM <- subset(FR_Union_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

FR_Union_Up_Cluster2_PGM <- subset(FR_Union_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

FR_Union_Up_Cluster3_FGM <- subset(FR_Union_Up, Group == "Cluster 3 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

FR_Union_Up_Cluster3_PGM <- subset(FR_Union_Up, Group == "Cluster 3 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

FR_Union_Down_Cluster1_FGM <- subset(FR_Union_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

FR_Union_Down_Cluster1_PGM <- subset(FR_Union_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

FR_Union_Down_Cluster2_FGM <- subset(FR_Union_Down, Group == "Cluster 2 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

FR_Union_Down_Cluster2_PGM <- subset(FR_Union_Down, Group == "Cluster 2 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_Up_Cluster1_FGM <- subset(SD56_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_Up_Cluster1_PGM <- subset(SD56_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_Up_Cluster2_FGM <- subset(SD56_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_Up_Cluster2_PGM <- subset(SD56_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_Up_Cluster3_FGM <- subset(SD56_Up, Group == "Cluster 3 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_Up_Cluster3_PGM <- subset(SD56_Up, Group == "Cluster 3 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_Down_Cluster1_FGM <- subset(SD56_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_Down_Cluster1_PGM <- subset(SD56_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_Down_Cluster2_FGM <- subset(SD56_Down, Group == "Cluster 2 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_Down_Cluster2_PGM <- subset(SD56_Down, Group == "Cluster 2 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_Down_Cluster3_FGM <- subset(SD56_Down, Group == "Cluster 3 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_Down_Cluster3_PGM <- subset(SD56_Down, Group == "Cluster 3 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD3_SD56_RS2_Up_Cluster1_FGM <- subset(SD3_SD56_RS2_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD3_SD56_RS2_Up_Cluster1_PGM <- subset(SD3_SD56_RS2_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD3_SD56_RS2_Up_Cluster2_FGM <- subset(SD3_SD56_RS2_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD3_SD56_RS2_Up_Cluster2_PGM <- subset(SD3_SD56_RS2_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SR_Union_Up_Cluster1_FGM <- subset(SR_Union_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SR_Union_Up_Cluster1_PGM <- subset(SR_Union_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SR_Union_Up_Cluster2_FGM <- subset(SR_Union_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SR_Union_Up_Cluster2_PGM <- subset(SR_Union_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SR_Union_Up_Cluster3_FGM <- subset(SR_Union_Up, Group == "Cluster 3 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SR_Union_Up_Cluster3_PGM <- subset(SR_Union_Up, Group == "Cluster 3 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SR_Union_Down_Cluster1_FGM <- subset(SR_Union_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SR_Union_Down_Cluster1_PGM <- subset(SR_Union_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SR_Union_Down_Cluster2_FGM <- subset(SR_Union_Down, Group == "Cluster 2 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SR_Union_Down_Cluster2_PGM <- subset(SR_Union_Down, Group == "Cluster 2 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SR_Union_Down_Cluster3_FGM <- subset(SR_Union_Down, Group == "Cluster 3 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SR_Union_Down_Cluster3_PGM <- subset(SR_Union_Down, Group == "Cluster 3 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_RS2_Up_Cluster1_FGM <- subset(SD56_RS2_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_RS2_Up_Cluster1_PGM <- subset(SD56_RS2_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_RS2_Up_Cluster2_FGM <- subset(SD56_RS2_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_RS2_Up_Cluster2_PGM <- subset(SD56_RS2_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_RS2_Down_Cluster1_FGM <- subset(SD56_RS2_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_RS2_Down_Cluster1_PGM <- subset(SD56_RS2_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_RS2_Down_Cluster2_FGM <- subset(SD56_RS2_Down, Group == "Cluster 2 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_RS2_Down_Cluster2_PGM <- subset(SD56_RS2_Down, Group == "Cluster 2 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_RS2_Down_Cluster3_FGM <- subset(SD56_RS2_Down, Group == "Cluster 3 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_RS2_Down_Cluster3_PGM <- subset(SD56_RS2_Down, Group == "Cluster 3 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_RS2_Down_Cluster4_FGM <- subset(SD56_RS2_Down, Group == "Cluster 4 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_RS2_Down_Cluster4_PGM <- subset(SD56_RS2_Down, Group == "Cluster 4 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD3_SD56_RS2_RS6_Up_Cluster1_FGM <- subset(SD3_SD56_RS2_RS6_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD3_SD56_RS2_RS6_Up_Cluster1_PGM <- subset(SD3_SD56_RS2_RS6_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD3_SD56_RS2_RS6_Up_Cluster2_FGM <- subset(SD3_SD56_RS2_RS6_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD3_SD56_RS2_RS6_Up_Cluster2_PGM <- subset(SD3_SD56_RS2_RS6_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD3_SD56_RS2_RS6_Down_Cluster1_FGM <- subset(SD3_SD56_RS2_RS6_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD3_SD56_RS2_RS6_Down_Cluster1_PGM <- subset(SD3_SD56_RS2_RS6_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

NR_Union_Up_Cluster1_FGM <- subset(NR_Union_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

NR_Union_Up_Cluster1_PGM <- subset(NR_Union_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

NR_Union_Up_Cluster2_FGM <- subset(NR_Union_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

NR_Union_Up_Cluster2_PGM <- subset(NR_Union_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

NR_Union_Down_Cluster1_FGM <- subset(NR_Union_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

NR_Union_Down_Cluster1_PGM <- subset(NR_Union_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_RS2_RS6_Down_Cluster1_FGM <- subset(SD56_RS2_RS6_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

SD56_RS2_RS6_Down_Cluster1_PGM <- subset(SD56_RS2_RS6_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

RS2_Up_Cluster1_FGM <- subset(RS2_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

RS2_Up_Cluster1_PGM <- subset(RS2_Up, Group == "Cluster 1 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

RS2_Up_Cluster2_FGM <- subset(RS2_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

RS2_Up_Cluster2_PGM <- subset(RS2_Up, Group == "Cluster 2 Up") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()

RS2_Down_Cluster1_FGM <- subset(RS2_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(Fold_Enrichment) %>% pull() %>% as.numeric() %>% geometric.mean()

RS2_Down_Cluster1_PGM <- subset(RS2_Down, Group == "Cluster 1 Down") %>%
  dplyr::select(p_Value) %>% pull() %>% as.numeric() %>% geometric.mean()


### Assemble Data Frames of Clustered Terms ----

SD3_SD56_Up_Clustered_DF <- data.frame(
  Order = 1:3,
  Gene_Count = c(length(SD3_SD56_Up_Cluster1_Genes), length(SD3_SD56_Up_Cluster2_Genes),
                 length(SD3_SD56_Up_Cluster3_Genes)),
  Enrichment = c(SD3_SD56_Up_Cluster1_FGM, SD3_SD56_Up_Cluster2_FGM, SD3_SD56_Up_Cluster3_FGM),
  p_Value = c(SD3_SD56_Up_Cluster1_PGM, SD3_SD56_Up_Cluster2_PGM, SD3_SD56_Up_Cluster3_PGM),
  Terms = c("Transcription regulation (BP, MF)", "MAPK signaling pathway (MF, KEGG)",
            "Pathways in cancer (KEGG)"),
  row.names = c("Cluster1", "Cluster2", "Cluster3"))

FR_Union_Up_Clustered_DF <- data.frame(
  Order = 1:3,
  Gene_Count = c(length(FR_Union_Up_Cluster1_Genes), length(FR_Union_Up_Cluster2_Genes),
                 length(FR_Union_Up_Cluster3_Genes)),
  Enrichment = c(FR_Union_Up_Cluster1_FGM, FR_Union_Up_Cluster2_FGM, FR_Union_Up_Cluster3_FGM),
  p_Value = c(FR_Union_Up_Cluster1_PGM, FR_Union_Up_Cluster2_PGM, FR_Union_Up_Cluster3_PGM),
  Terms = c("Transcription regulation (BP, MF)", "Phosphorylation (MF)",
            "Pathways in cancer (KEGG)"),
  row.names = c("Cluster1", "Cluster2", "Cluster3"))

FR_Union_Down_Clustered_DF <- data.frame(
  Order = 1:2,
  Gene_Count = c(length(FR_Union_Down_Cluster1_Genes), length(FR_Union_Down_Cluster2_Genes)),
  Enrichment = c(FR_Union_Down_Cluster1_FGM, FR_Union_Down_Cluster2_FGM),
  p_Value = c(FR_Union_Down_Cluster1_PGM, FR_Union_Down_Cluster2_PGM),
  Terms = c("DNA repair (BP, KEGG)", "Metabolism (KEGG)"),
  row.names = c("Cluster1", "Cluster2"))

SD56_Up_Clustered_DF <- data.frame(
  Order = 1:3,
  Gene_Count = c(length(SD56_Up_Cluster1_Genes), length(SD56_Up_Cluster2_Genes),
                 length(SD56_Up_Cluster3_Genes)),
  Enrichment = c(SD56_Up_Cluster1_FGM, SD56_Up_Cluster2_FGM, SD56_Up_Cluster3_FGM),
  p_Value = c(SD56_Up_Cluster1_PGM, SD56_Up_Cluster2_PGM, SD56_Up_Cluster3_PGM),
  Terms = c("Chromatin regulation (BP, MF)", "Phosphorylation (MF)",
            "Insulin Regulation/FoxO signaling pathway (KEGG)"),
  row.names = c("Cluster1", "Cluster2", "Cluster3"))

SD56_Down_Clustered_DF <- data.frame(
  Order = 1:3,
  Gene_Count = c(length(SD56_Down_Cluster1_Genes), length(SD56_Down_Cluster2_Genes),
                 length(SD56_Down_Cluster3_Genes)),
  Enrichment = c(SD56_Down_Cluster1_FGM, SD56_Down_Cluster2_FGM, SD56_Down_Cluster3_FGM),
  p_Value = c(SD56_Down_Cluster1_PGM, SD56_Down_Cluster2_PGM, SD56_Down_Cluster3_PGM),
  Terms = c("Oxidative phosphorylation (KEGG)", "Glutathione metabolism (KEGG)",
            "Response to infection (KEGG)"),
  row.names = c("Cluster1", "Cluster2", "Cluster3"))

SD3_SD56_RS2_Up_Clustered_DF <- data.frame(
  Order = 1:2,
  Gene_Count = c(length(SD3_SD56_RS2_Up_Cluster1_Genes), length(SD3_SD56_RS2_Up_Cluster2_Genes)),
  Enrichment = c(SD3_SD56_RS2_Up_Cluster1_FGM, SD3_SD56_RS2_Up_Cluster2_FGM),
  p_Value = c(SD3_SD56_RS2_Up_Cluster1_PGM, SD3_SD56_RS2_Up_Cluster2_PGM),
  Terms = c("RNA-binding (BP, MF, KEGG)", "Stress response (BP, MF, KEGG)"),
  row.names = c("Cluster1", "Cluster2"))

SR_Union_Up_Clustered_DF <- data.frame(
  Order = 1:3,
  Gene_Count = c(length(SR_Union_Up_Cluster1_Genes), length(SR_Union_Up_Cluster2_Genes),
                 length(SR_Union_Up_Cluster3_Genes)),
  Enrichment = c(SR_Union_Up_Cluster1_FGM, SR_Union_Up_Cluster2_FGM, SR_Union_Up_Cluster3_FGM),
  p_Value = c(SR_Union_Up_Cluster1_PGM, SR_Union_Up_Cluster2_PGM, SR_Union_Up_Cluster3_PGM),
  Terms = c("Stress response (BP, MF, KEGG)", "RNA-binding (BP, MF, KEGG)",
            "Longevity regulating pathways (KEGG)"),
  row.names = c("Cluster1", "Cluster2", "Cluster3"))

SR_Union_Down_Clustered_DF <- data.frame(
  Order = 1:3,
  Gene_Count = c(length(SR_Union_Down_Cluster1_Genes), length(SR_Union_Down_Cluster2_Genes),
                 length(SR_Union_Down_Cluster3_Genes)),
  Enrichment = c(SR_Union_Down_Cluster1_FGM, SR_Union_Down_Cluster2_FGM,
                 SR_Union_Down_Cluster3_FGM),
  p_Value = c(SR_Union_Down_Cluster1_PGM, SR_Union_Down_Cluster2_PGM,
              SR_Union_Down_Cluster3_PGM),
  Terms = c("Cholesterol biosynthesis (BP, KEGG)", "Ion transport (BP, MF)",
            "Neuroendocrine signaling pathways (KEGG)"),
  row.names = c("Cluster1", "Cluster2", "Cluster3"))

SD56_RS2_Up_Clustered_DF <- data.frame(
  Order = 1:2,
  Gene_Count = c(length(SD56_RS2_Up_Cluster1_Genes), length(SD56_RS2_Up_Cluster2_Genes)),
  Enrichment = c(SD56_RS2_Up_Cluster1_FGM, SD56_RS2_Up_Cluster2_FGM),
  p_Value = c(SD56_RS2_Up_Cluster1_PGM, SD56_RS2_Up_Cluster2_PGM),
  Terms = c("mTOR signaling pathway (KEGG)", "Cell adhesion (KEGG)"),
  row.names = c("Cluster1", "Cluster2"))

SD56_RS2_Down_Clustered_DF <- data.frame(
  Order = 1:4,
  Gene_Count = c(length(SD56_RS2_Down_Cluster1_Genes), length(SD56_RS2_Down_Cluster2_Genes),
                 length(SD56_RS2_Down_Cluster3_Genes), length(SD56_RS2_Down_Cluster4_Genes)),
  Enrichment = c(SD56_RS2_Down_Cluster1_FGM, SD56_RS2_Down_Cluster2_FGM,
                 SD56_RS2_Down_Cluster3_FGM, SD56_RS2_Down_Cluster4_FGM),
  p_Value = c(SD56_RS2_Down_Cluster1_PGM, SD56_RS2_Down_Cluster2_PGM,
              SD56_RS2_Down_Cluster3_PGM, SD56_RS2_Down_Cluster4_PGM),
  Terms = c("Cholesterol biosynthesis (BP, KEGG)", "Ion transport (BP, MF)",
            "Cell signaling (KEGG)", "Glucose metabolism (KEGG)"),
  row.names = c("Cluster1", "Cluster2", "Cluster3", "Cluster4"))

SD3_SD56_RS2_RS6_Up_Clustered_DF <- data.frame(
  Order = 1:2,
  Gene_Count = c(length(SD3_SD56_RS2_RS6_Up_Cluster1_Genes), length(SD3_SD56_RS2_RS6_Up_Cluster2_Genes)),
  Enrichment = c(SD3_SD56_RS2_RS6_Up_Cluster1_FGM, SD3_SD56_RS2_RS6_Up_Cluster2_FGM),
  p_Value = c(SD3_SD56_RS2_RS6_Up_Cluster1_PGM, SD3_SD56_RS2_RS6_Up_Cluster2_PGM),
  Terms = c("Stress response (BP, MF, KEGG)", "Protein processing in endoplasmic reticulum (KEGG)"),
  row.names = c("Cluster1", "Cluster2"))

SD3_SD56_RS2_RS6_Down_Clustered_DF <- data.frame(
  Order = 1,
  Gene_Count = length(SD3_SD56_RS2_RS6_Down_Cluster1_Genes),
  Enrichment = SD3_SD56_RS2_RS6_Down_Cluster1_FGM,
  p_Value = SD3_SD56_RS2_RS6_Down_Cluster1_PGM,
  Terms = "Metabolic pathways (MF, KEGG)",
  row.names = c("Cluster1"))

NR_Union_Up_Clustered_DF <- data.frame(
  Order = 1:2,
  Gene_Count = c(length(NR_Union_Up_Cluster1_Genes), length(NR_Union_Up_Cluster2_Genes)),
  Enrichment = c(NR_Union_Up_Cluster1_FGM, NR_Union_Up_Cluster2_FGM),
  p_Value = c(NR_Union_Up_Cluster1_PGM, NR_Union_Up_Cluster2_PGM),
  Terms = c("Stress response (BP, MF, KEGG)",
            "MAPK and PI3K/AKT signaling pathways (KEGG)"),
  row.names = c("Cluster1", "Cluster2"))

NR_Union_Down_Clustered_DF <- data.frame(
  Order = 1,
  Gene_Count = length(NR_Union_Down_Cluster1_Genes),
  Enrichment = NR_Union_Down_Cluster1_FGM,
  p_Value = NR_Union_Down_Cluster1_PGM,
  Terms = "Cholesterol biosynthesis (BP)",
  row.names = c("Cluster1"))

SD56_RS2_RS6_Down_Clustered_DF <- data.frame(
  Order = 1,
  Gene_Count = length(SD56_RS2_RS6_Down_Cluster1_Genes),
  Enrichment = SD56_RS2_RS6_Down_Cluster1_FGM,
  p_Value = SD56_RS2_RS6_Down_Cluster1_PGM,
  Terms = "Cholesterol biosynthesis (BP)",
  row.names = c("Cluster1"))

RS2_Up_Clustered_DF <- data.frame(
  Order = 1:2,
  Gene_Count = c(length(RS2_Up_Cluster1_Genes), length(RS2_Up_Cluster2_Genes)),
  Enrichment = c(RS2_Up_Cluster1_FGM, RS2_Up_Cluster2_FGM),
  p_Value = c(RS2_Up_Cluster1_PGM, RS2_Up_Cluster2_PGM),
  Terms = c("RNA-binding (BP, MF, KEGG)", "Electron transport (affected in neurodegeneration) (BP, KEGG)"),
  row.names = c("Cluster1", "Cluster2"))

RS2_Down_Clustered_DF <- data.frame(
  Order = 1,
  Gene_Count = length(RS2_Down_Cluster1_Genes),
  Enrichment = RS2_Down_Cluster1_FGM,
  p_Value = RS2_Down_Cluster1_PGM,
  Terms = "Phosphorylation (MF)",
  row.names = "Cluster1")


### Obtain Key Word for Each Term ----

## Loop through unclustered terms and automatically determine KW category

# create a new column titled KW (short for keyword)
SD3_SD56_Up$KW <- NA
z <- SD3_SD56_Up$Category

for (i in seq_along(z)) {
  if (SD3_SD56_Up$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    SD3_SD56_Up$KW[i] <- "(BP)"}
  if (SD3_SD56_Up$Category[i] == "KEGG_PATHWAY") {
  SD3_SD56_Up$KW[i] <- "(KEGG)" }
  if (SD3_SD56_Up$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    SD3_SD56_Up$KW[i] <- "(MF)"}
}

SD3_SD56_Down$KW <- NA
z <- SD3_SD56_Down$Category

for (i in seq_along(z)) {
  if (SD3_SD56_Down$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    SD3_SD56_Down$KW[i] <- "(BP)"}
  if (SD3_SD56_Down$Category[i] == "KEGG_PATHWAY") {
    SD3_SD56_Down$KW[i] <- "(KEGG)" }
  if (SD3_SD56_Down$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    SD3_SD56_Down$KW[i] <- "(MF)"}
}

SD56_Up$KW <- NA
z <- SD56_Up$Category

for (i in seq_along(z)) {
  if (SD56_Up$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    SD56_Up$KW[i] <- "(BP)"}
  if (SD56_Up$Category[i] == "KEGG_PATHWAY") {
    SD56_Up$KW[i] <- "(KEGG)" }
  if (SD56_Up$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    SD56_Up$KW[i] <- "(MF)"}
}

SD56_Down$KW <- NA
z <- SD56_Down$Category

for (i in seq_along(z)) {
  if (SD56_Down$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    SD56_Down$KW[i] <- "(BP)"}
  if (SD56_Down$Category[i] == "KEGG_PATHWAY") {
    SD56_Down$KW[i] <- "(KEGG)" }
  if (SD56_Down$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    SD56_Down$KW[i] <- "(MF)"}
}

FR_Union_Up$KW <- NA
z <- FR_Union_Up$Category

for (i in seq_along(z)) {
  if (FR_Union_Up$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    FR_Union_Up$KW[i] <- "(BP)"}
  if (FR_Union_Up$Category[i] == "KEGG_PATHWAY") {
    FR_Union_Up$KW[i] <- "(KEGG)" }
  if (FR_Union_Up$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    FR_Union_Up$KW[i] <- "(MF)"}
}

FR_Union_Down$KW <- NA
z <- FR_Union_Down$Category

for (i in seq_along(z)) {
  if (FR_Union_Down$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    FR_Union_Down$KW[i] <- "(BP)"}
  if (FR_Union_Down$Category[i] == "KEGG_PATHWAY") {
    FR_Union_Down$KW[i] <- "(KEGG)" }
  if (FR_Union_Down$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    FR_Union_Down$KW[i] <- "(MF)"}
}

SD3_SD56_RS2_Up$KW <- NA
z <- SD3_SD56_RS2_Up$Category

for (i in seq_along(z)) {
  if (SD3_SD56_RS2_Up$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    SD3_SD56_RS2_Up$KW[i] <- "(BP)"}
  if (SD3_SD56_RS2_Up$Category[i] == "KEGG_PATHWAY") {
    SD3_SD56_RS2_Up$KW[i] <- "(KEGG)" }
  if (SD3_SD56_RS2_Up$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    SD3_SD56_RS2_Up$KW[i] <- "(MF)"}
}

SD3_SD56_RS2_Down$KW <- NA
z <- SD3_SD56_RS2_Down$Category

for (i in seq_along(z)) {
  if (SD3_SD56_RS2_Down$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    SD3_SD56_RS2_Down$KW[i] <- "(BP)"}
  if (SD3_SD56_RS2_Down$Category[i] == "KEGG_PATHWAY") {
    SD3_SD56_RS2_Down$KW[i] <- "(KEGG)" }
  if (SD3_SD56_RS2_Down$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    SD3_SD56_RS2_Down$KW[i] <- "(MF)"}
}

SD56_RS2_Up$KW <- NA
z <- SD56_RS2_Up$Category

for (i in seq_along(z)) {
  if (SD56_RS2_Up$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    SD56_RS2_Up$KW[i] <- "(BP)"}
  if (SD56_RS2_Up$Category[i] == "KEGG_PATHWAY") {
    SD56_RS2_Up$KW[i] <- "(KEGG)" }
  if (SD56_RS2_Up$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    SD56_RS2_Up$KW[i] <- "(MF)"}
}

SD56_RS2_Down$KW <- NA
z <- SD56_RS2_Down$Category

for (i in seq_along(z)) {
  if (SD56_RS2_Down$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    SD56_RS2_Down$KW[i] <- "(BP)"}
  if (SD56_RS2_Down$Category[i] == "KEGG_PATHWAY") {
    SD56_RS2_Down$KW[i] <- "(KEGG)" }
  if (SD56_RS2_Down$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    SD56_RS2_Down$KW[i] <- "(MF)"}
}

SR_Union_Up$KW <- NA
z <- SR_Union_Up$Category

for (i in seq_along(z)) {
  if (SR_Union_Up$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    SR_Union_Up$KW[i] <- "(BP)"}
  if (SR_Union_Up$Category[i] == "KEGG_PATHWAY") {
    SR_Union_Up$KW[i] <- "(KEGG)" }
  if (SR_Union_Up$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    SR_Union_Up$KW[i] <- "(MF)"}
}

SR_Union_Down$KW <- NA
z <- SR_Union_Down$Category

for (i in seq_along(z)) {
  if (SR_Union_Down$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    SR_Union_Down$KW[i] <- "(BP)"}
  if (SR_Union_Down$Category[i] == "KEGG_PATHWAY") {
    SR_Union_Down$KW[i] <- "(KEGG)" }
  if (SR_Union_Down$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    SR_Union_Down$KW[i] <- "(MF)"}
}

SD3_SD56_RS2_RS6_Up$KW <- NA
z <- SD3_SD56_RS2_RS6_Up$Category

for (i in seq_along(z)) {
  if (SD3_SD56_RS2_RS6_Up$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    SD3_SD56_RS2_RS6_Up$KW[i] <- "(BP)"}
  if (SD3_SD56_RS2_RS6_Up$Category[i] == "KEGG_PATHWAY") {
    SD3_SD56_RS2_RS6_Up$KW[i] <- "(KEGG)" }
  if (SD3_SD56_RS2_RS6_Up$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    SD3_SD56_RS2_RS6_Up$KW[i] <- "(MF)"}
}

SD3_SD56_RS2_RS6_Down$KW <- NA
z <- SD3_SD56_RS2_RS6_Down$Category

for (i in seq_along(z)) {
  if (SD3_SD56_RS2_RS6_Down$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    SD3_SD56_RS2_RS6_Down$KW[i] <- "(BP)"}
  if (SD3_SD56_RS2_RS6_Down$Category[i] == "KEGG_PATHWAY") {
    SD3_SD56_RS2_RS6_Down$KW[i] <- "(KEGG)" }
  if (SD3_SD56_RS2_RS6_Down$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    SD3_SD56_RS2_RS6_Down$KW[i] <- "(MF)"}
}

NR_Union_Up$KW <- NA
z <- NR_Union_Up$Category

for (i in seq_along(z)) {
  if (NR_Union_Up$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    NR_Union_Up$KW[i] <- "(BP)"}
  if (NR_Union_Up$Category[i] == "KEGG_PATHWAY") {
    NR_Union_Up$KW[i] <- "(KEGG)" }
  if (NR_Union_Up$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    NR_Union_Up$KW[i] <- "(MF)"}
}

NR_Union_Down$KW <- NA
z <- NR_Union_Down$Category

for (i in seq_along(z)) {
  if (NR_Union_Down$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    NR_Union_Down$KW[i] <- "(BP)"}
  if (NR_Union_Down$Category[i] == "KEGG_PATHWAY") {
    NR_Union_Down$KW[i] <- "(KEGG)" }
  if (NR_Union_Down$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    NR_Union_Down$KW[i] <- "(MF)"}
}

SD56_RS2_RS6_Up$KW <- NA
z <- SD56_RS2_RS6_Up$Category

for (i in seq_along(z)) {
  if (SD56_RS2_RS6_Up$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    SD56_RS2_RS6_Up$KW[i] <- "(BP)"}
  if (SD56_RS2_RS6_Up$Category[i] == "KEGG_PATHWAY") {
    SD56_RS2_RS6_Up$KW[i] <- "(KEGG)" }
  if (SD56_RS2_RS6_Up$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    SD56_RS2_RS6_Up$KW[i] <- "(MF)"}
}

SD56_RS2_RS6_Down$KW <- NA
z <- SD56_RS2_RS6_Down$Category

for (i in seq_along(z)) {
  if (SD56_RS2_RS6_Down$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    SD56_RS2_RS6_Down$KW[i] <- "(BP)"}
  if (SD56_RS2_RS6_Down$Category[i] == "KEGG_PATHWAY") {
    SD56_RS2_RS6_Down$KW[i] <- "(KEGG)" }
  if (SD56_RS2_RS6_Down$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    SD56_RS2_RS6_Down$KW[i] <- "(MF)"}
}

RS2_Up$KW <- NA
z <- RS2_Up$Category

for (i in seq_along(z)) {
  if (RS2_Up$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    RS2_Up$KW[i] <- "(BP)"}
  if (RS2_Up$Category[i] == "KEGG_PATHWAY") {
    RS2_Up$KW[i] <- "(KEGG)" }
  if (RS2_Up$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    RS2_Up$KW[i] <- "(MF)"}
}

RS2_Down$KW <- NA
z <- RS2_Down$Category

for (i in seq_along(z)) {
  if (RS2_Down$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    RS2_Down$KW[i] <- "(BP)"}
  if (RS2_Down$Category[i] == "KEGG_PATHWAY") {
    RS2_Down$KW[i] <- "(KEGG)" }
  if (RS2_Down$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    RS2_Down$KW[i] <- "(MF)"}
}

RS6_Up$KW <- NA
z <- RS6_Up$Category

for (i in seq_along(z)) {
  if (RS6_Up$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    RS6_Up$KW[i] <- "(BP)"}
  if (RS6_Up$Category[i] == "KEGG_PATHWAY") {
    RS6_Up$KW[i] <- "(KEGG)" }
  if (RS6_Up$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    RS6_Up$KW[i] <- "(MF)"}
}

RS6_Down$KW <- NA
z <- RS6_Down$Category

for (i in seq_along(z)) {
  if (RS6_Down$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    RS6_Down$KW[i] <- "(BP)"}
  if (RS6_Down$Category[i] == "KEGG_PATHWAY") {
    RS6_Down$KW[i] <- "(KEGG)" }
  if (RS6_Down$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    RS6_Down$KW[i] <- "(MF)"}
}


### Assemble Data Frames of Unclustered Terms ----

SD3_SD56_Up_Unclustered_DF <- subset(SD3_SD56_Up, Group == "Unclustered Up") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 4:8) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  # use regex to extract function from term name
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(SD3_SD56_Up_Unclustered_DF) <- SD3_SD56_Up_Unclustered_DF$Terms

SD3_SD56_Down_Unclustered_DF <- subset(SD3_SD56_Down, Group == "Unclustered Down") %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(SD3_SD56_Down_Unclustered_DF) <- SD3_SD56_Down_Unclustered_DF$Terms

FR_Union_Up_Unclustered_DF <- subset(FR_Union_Up, Group == "Unclustered Up") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 4:12) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(FR_Union_Up_Unclustered_DF) <- FR_Union_Up_Unclustered_DF$Terms

FR_Union_Down_Unclustered_DF <- subset(FR_Union_Down, Group == "Unclustered Down") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 3:20) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(FR_Union_Down_Unclustered_DF) <- FR_Union_Down_Unclustered_DF$Terms

SD56_Up_Unclustered_DF <- subset(SD56_Up, Group == "Unclustered Up") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 4:12) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(SD56_Up_Unclustered_DF) <- SD56_Up_Unclustered_DF$Terms

SD56_Down_Unclustered_DF <- subset(SD56_Down, Group == "Unclustered Down") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 4:25) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(SD56_Down_Unclustered_DF) <- SD56_Down_Unclustered_DF$Terms

SD3_SD56_RS2_Up_Unclustered_DF <- subset(SD3_SD56_RS2_Up, Group == "Unclustered Up") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 3:5) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(SD3_SD56_RS2_Up_Unclustered_DF) <- SD3_SD56_RS2_Up_Unclustered_DF$Terms

SD3_SD56_RS2_Down_Unclustered_DF <- subset(SD3_SD56_RS2_Down, Group == "Unclustered Down") %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(SD3_SD56_RS2_Down_Unclustered_DF) <- SD3_SD56_RS2_Down_Unclustered_DF$Terms

SR_Union_Up_Unclustered_DF <- subset(SR_Union_Up, Group == "Unclustered Up") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 4:8) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(SR_Union_Up_Unclustered_DF) <- SR_Union_Up_Unclustered_DF$Terms

SR_Union_Down_Unclustered_DF <- subset(SR_Union_Down, Group == "Unclustered Down") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 4:16) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(SR_Union_Down_Unclustered_DF) <- SR_Union_Down_Unclustered_DF$Terms

SD56_RS2_Up_Unclustered_DF <- subset(SD56_RS2_Up, Group == "Unclustered Up") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 3:7) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(SD56_RS2_Up_Unclustered_DF) <- SD56_RS2_Up_Unclustered_DF$Terms

SD56_RS2_Down_Unclustered_DF <- subset(SD56_RS2_Down, Group == "Unclustered Down") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 5:13) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(SD56_RS2_Down_Unclustered_DF) <- SD56_RS2_Down_Unclustered_DF$Terms

SD3_SD56_RS2_RS6_Up_Unclustered_DF <- subset(SD3_SD56_RS2_RS6_Up, Group == "Unclustered Up") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 3:5) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW)%>% as.data.frame()
rownames(SD3_SD56_RS2_RS6_Up_Unclustered_DF) <- SD3_SD56_RS2_RS6_Up_Unclustered_DF$Terms

SD3_SD56_RS2_RS6_Down_Unclustered_DF <- subset(SD3_SD56_RS2_RS6_Down,
                                               Group == "Unclustered Down") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 2:3) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(SD3_SD56_RS2_RS6_Down_Unclustered_DF) <- SD3_SD56_RS2_RS6_Down_Unclustered_DF$Terms

NR_Union_Up_Unclustered_DF <- subset(NR_Union_Up, Group == "Unclustered Up") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 3:6) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(NR_Union_Up_Unclustered_DF) <- NR_Union_Up_Unclustered_DF$Terms

NR_Union_Down_Unclustered_DF <- subset(NR_Union_Down, Group == "Unclustered Down") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 2:7) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(NR_Union_Down_Unclustered_DF) <- NR_Union_Down_Unclustered_DF$Terms

SD56_RS2_RS6_Up_Unclustered_DF <- subset(SD56_RS2_RS6_Up, Group == "Unclustered Up") %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(SD56_RS2_RS6_Up_Unclustered_DF) <- SD56_RS2_RS6_Up_Unclustered_DF$Terms

SD56_RS2_RS6_Down_Unclustered_DF <- subset(SD56_RS2_RS6_Down, Group == "Unclustered Down") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 2:5) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(SD56_RS2_RS6_Down_Unclustered_DF) <- SD56_RS2_RS6_Down_Unclustered_DF$Terms

RS2_Up_Unclustered_DF <- subset(RS2_Up, Group == "Unclustered Up") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 3:11) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(RS2_Up_Unclustered_DF) <- RS2_Up_Unclustered_DF$Terms

RS2_Down_Unclustered_DF <- subset(RS2_Down, Group == "Unclustered Down") %>%
  dplyr::select(Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  mutate(Order = 2:7) %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(RS2_Down_Unclustered_DF) <- RS2_Down_Unclustered_DF$Terms

RS6_Up_Unclustered_DF <- subset(RS6_Up, Group == "Unclustered Up") %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(RS6_Up_Unclustered_DF) <- RS6_Up_Unclustered_DF$Terms

RS6_Down_Unclustered_DF <- subset(RS6_Down, Group == "Unclustered Down") %>%
  dplyr::select(Order, Gene_Count, Fold_Enrichment, p_Value, Term, KW) %>%
  rename(Enrichment = Fold_Enrichment, Terms = Term) %>%
  mutate(Terms = sub(".*[~:]", "", Terms)) %>%
  mutate(Terms = paste(Terms, KW)) %>%
  dplyr::select(-KW) %>%
  as.data.frame()
rownames(RS6_Down_Unclustered_DF) <- RS6_Down_Unclustered_DF$Terms


### Theme Elements for Plotting ----

upregulated_text_theme <- theme(text = element_text(colour = "#000000", size=7), 
                                axis.text.y =  element_text(colour = "#000000", size=7),
                                legend.text = element_text(colour = "#000000", size=7),
                                axis.text.x = element_blank(),
                                axis.ticks.x = element_blank())

downregulated_text_theme <- theme(text = element_text(colour = "#000000", size=7), 
                                  axis.text.y =  element_text(colour = "#000000", size=7), 
                                  legend.text = element_text(colour = "#000000", size=7), 
                                  axis.text.x = element_text(colour = "#000000", size=7))


### SD3 + SD5-6 Bubble Plot ----

SD3_SD56_Up_Plot <- rbind(SD3_SD56_Up_Clustered_DF,
                          SD3_SD56_Up_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Cluster2", "Cluster3", "Unclustered"),
                     times = c(1, 1, 1, 5)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

SD3_SD56_Up_Plot_Condensed <- ggplot(
  data = SD3_SD56_Up_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                # size of bubbles based on the number of unique genes within term or cluster
                size = Gene_Count,
                # color gradient is based on significance (p-Value)
                color = p_Value)) +
  geom_point(alpha = 1) + # alpha is the opacity of the bubbles
  # limits sets the min and max Gene_Count to be included, breaks defines the size of bubbles
  # included on the legend, range sets the min and max size of bubbles after transformation
  scale_size_continuous(limits = c(1,200), breaks = c(50,100,150,200), range = c(1,4.345)) +
  # red color palette for upregulated terms
  shades::lightness(scale_colour_distiller(palette = "Greens",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,8)) + # set min and max for x-axis
  theme_light() + # set the background color 
  upregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "", y = NULL)

SD3_SD56_Down_Plot <- SD3_SD56_Down_Unclustered_DF %>%
  mutate(Group = rep("Unclustered", times = 4),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

SD3_SD56_Down_Plot_Condensed <- ggplot(
  data = SD3_SD56_Down_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,200), breaks = c(50,100,150,200), range = c(1,4.345)) + 
  # blue color palette for downregulated terms
  shades::lightness(scale_colour_distiller(palette = "Purples",
                                           limits = c(0,0.05)), scalefac(.7)) +
  coord_cartesian(xlim = c(1,8)) +
  theme_light() +
  downregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL)

# Save overlaid plot of up/down
# pdf("SD3_SD56_Bubble_Plot_Updated.pdf", width = 7.139, height = 2.733)
# grid.newpage()
# grid.draw(rbind(ggplotGrob(SD3_SD56_Up_Plot_Condensed),
#                 ggplotGrob(SD3_SD56_Down_Plot_Condensed),
#                 size = "first"))
# dev.off()


### SD5-6 Bubble Plot ----

SD56_Up_Plot <- rbind(SD56_Up_Clustered_DF,
                      SD56_Up_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Cluster2", "Cluster3", "Unclustered"),
                     times = c(1, 1, 1, 9)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

SD56_Up_Plot_Condensed <- ggplot(
  data = SD56_Up_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,200), breaks = c(50,100,150,200), range = c(1,4.345)) +
  shades::lightness(scale_colour_distiller(palette = "Greens",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,8)) +
  theme_light() +
  upregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "", y = NULL)

SD56_Down_Plot <- rbind(SD56_Down_Clustered_DF,
                        SD56_Down_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Cluster2", "Cluster3", "Unclustered"),
                     times = c(1, 1, 1, 22)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

SD56_Down_Plot_Condensed <- ggplot(
  data = SD56_Down_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,200), breaks = c(50,100,150,200), range = c(1,4.345)) +
  shades::lightness(scale_colour_distiller(palette = "Purples",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,8)) +
  theme_light() +
  downregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL)

# Save overlaid plot of up/down
# pdf("SD56_Bubble_Plot_Updated.pdf", width = 7.7435, height = 6.395)
# grid.newpage()
# grid.draw(rbind(ggplotGrob(SD56_Up_Plot_Condensed),
#                 ggplotGrob(SD56_Down_Plot_Condensed),
#                 size = "first"))
# dev.off()


### (SD3 + SD5-6) OR SD5-6 Bubble Plot ----

FR_Union_Up_Plot <- rbind(FR_Union_Up_Clustered_DF,
                          FR_Union_Up_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Cluster2", "Cluster3", "Unclustered"),
                     times = c(1, 1, 1, 9)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

FR_Union_Up_Plot_Condensed <- ggplot(
  data = FR_Union_Up_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,300), breaks = c(70,140,210,280), range = c(1,6)) +
  shades::lightness(scale_colour_distiller(palette = "Greens",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,6)) +
  theme_light() +
  upregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "", y = NULL)

FR_Union_Down_Plot <- rbind(FR_Union_Down_Clustered_DF,
                            FR_Union_Down_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Cluster2", "Unclustered"),
                     times = c(1, 1, 18)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

FR_Union_Down_Plot_Condensed <- ggplot(
  data = FR_Union_Down_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,300), breaks = c(70,140,210,280), range = c(1,6)) +
  shades::lightness(scale_colour_distiller(palette = "Purples",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,6)) +
  theme_light() +
  downregulated_text_theme+
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL)

# Save overlaid plot of up/down
# pdf("Figures/Fast_Recoverers_Bubble_Plot.pdf", width = 7.14127, height = 7.27)
# grid.newpage()
# grid.draw(rbind(ggplotGrob(FR_Union_Up_Plot_Condensed),
#                 ggplotGrob(FR_Union_Down_Plot_Condensed),
#                 size = "first"))
# dev.off()


### (SD3 + SD5-6) + RS2 Bubble Plot ----

SD3_SD56_RS2_Up_Plot <- rbind(SD3_SD56_RS2_Up_Clustered_DF,
                              SD3_SD56_RS2_Up_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Cluster2", "Unclustered"),
                     times = c(1, 1, 3)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

SD3_SD56_RS2_Up_Plot_Condensed <- ggplot(
  data = SD3_SD56_RS2_Up_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,200), breaks = c(35,70,105,140), range = c(1,7.7)) +
  shades::lightness(scale_colour_distiller(palette = "Greens",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  scale_x_continuous(breaks = seq(1, 13, by = 4)) +
  coord_cartesian(xlim = c(1,13)) +
  theme_light() +
  upregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "", y = NULL)

SD3_SD56_RS2_Down_Plot <- SD3_SD56_RS2_Down_Unclustered_DF %>%
  mutate(Group = rep("Unclustered", times = 4),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

SD3_SD56_RS2_Down_Plot_Condensed <- ggplot(
  data = SD3_SD56_RS2_Down_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,200), breaks = c(35,70,105,140), range = c(1,7.7)) +
  shades::lightness(scale_colour_distiller(palette = "Purples",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  scale_x_continuous(breaks = seq(1, 13, by = 4)) +
  coord_cartesian(xlim = c(1,13)) +
  theme_light() +
  downregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL)

# Save overlaid plot of up/down
# pdf("SD3_SD56_RS2_Bubble_Plot_Updated.pdf", width = 6.53691, height = 2.711)
# grid.newpage()
# grid.draw(rbind(ggplotGrob(SD3_SD56_RS2_Up_Plot_Condensed),
#                 ggplotGrob(SD3_SD56_RS2_Down_Plot_Condensed),
#                 size = "first"))
# dev.off()


### SD5-6 + RS2 Bubble Plot ----

SD56_RS2_Up_Plot <- rbind(SD56_RS2_Up_Clustered_DF,
                          SD56_RS2_Up_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Cluster2", "Unclustered"),
                     times = c(1, 1, 5)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

SD56_RS2_Up_Plot_Condensed <- ggplot(
  data = SD56_RS2_Up_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,200), breaks = c(35,70,105,140), range = c(1,7.7)) +
  shades::lightness(scale_colour_distiller(palette = "Greens",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  scale_x_continuous(breaks = seq(1, 13, by = 4)) +
  coord_cartesian(xlim = c(1,13)) +
  theme_light() +
  upregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "", y = NULL)

SD56_RS2_Down_Plot <- rbind(SD56_RS2_Down_Clustered_DF,
                            SD56_RS2_Down_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Cluster2", "Cluster3", "Cluster4",
                       "Unclustered"),
                     times = c(1, 1, 1, 1, 9)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

SD56_RS2_Down_Plot_Condensed <- ggplot(
  data = SD56_RS2_Down_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,200), breaks = c(35,70,105,140), range = c(1,7.7)) +
  shades::lightness(scale_colour_distiller(palette = "Purples",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  scale_x_continuous(breaks = seq(1, 13, by = 4)) +
  coord_cartesian(xlim = c(1,13)) +
  theme_light() +
  downregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL)

# Save overlaid plot of up/down
# pdf("SD56_RS2_Bubble_Plot_Updated.pdf", width = 6.70, height = 5.18)
# grid.newpage()
# grid.draw(rbind(ggplotGrob(SD56_RS2_Up_Plot_Condensed),
#                 ggplotGrob(SD56_RS2_Down_Plot_Condensed),
#                 size = "first"))
# dev.off()


### ((SD3 + SD5-6) + RS2) OR (SD5-6 + RS2) Bubble Plot ----

SR_Union_Up_Plot <- rbind(SR_Union_Up_Clustered_DF,
                          SR_Union_Up_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Cluster2", "Cluster3", "Unclustered"),
                     times = c(1, 1, 1, 5)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

SR_Union_Up_Plot_Condensed <- ggplot(
  data = SR_Union_Up_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,200), breaks = c(50,100,150,200), range = c(1,6)) +
  shades::lightness(scale_colour_distiller(palette = "Greens",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,6)) +
  theme_light() +
  upregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "", y = NULL)

SR_Union_Down_Plot <- rbind(SR_Union_Down_Clustered_DF,
                            SR_Union_Down_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Cluster2", "Cluster3", "Unclustered"),
                     times = c(1, 1, 1, 13)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

SR_Union_Down_Plot_Condensed <- ggplot(
  data = SR_Union_Down_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,200), breaks = c(50,100,150,200), range = c(1,6)) +
  shades::lightness(scale_colour_distiller(palette = "Purples",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,6)) +
  theme_light() +
  downregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL)

# Save overlaid plot of up/down
# pdf("Figures/Slow_Recoverers_Bubble_Plot.pdf", width = 7.0239, height = 5.915)
# grid.newpage()
# grid.draw(rbind(ggplotGrob(SR_Union_Up_Plot_Condensed),
#                 ggplotGrob(SR_Union_Down_Plot_Condensed),
#                 size = "first"))
# dev.off()


### (SD3 + SD56) + (RS2 + RS6) Bubble Plot ----

SD3_SD56_RS2_RS6_Up_Plot <- rbind(SD3_SD56_RS2_RS6_Up_Clustered_DF,
                                  SD3_SD56_RS2_RS6_Up_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Cluster2", "Unclustered"),
                     times = c(1, 1, 3)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

SD3_SD56_RS2_RS6_Up_Plot_Condensed <- ggplot(
  data = SD3_SD56_RS2_RS6_Up_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,28), breaks = c(7,14,21,28), range = c(1,6)) +
  shades::lightness(scale_colour_distiller(palette = "Greens",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,10)) +
  theme_light() +
  upregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "", y = NULL)

SD3_SD56_RS2_RS6_Down_Plot <- rbind(SD3_SD56_RS2_RS6_Down_Clustered_DF,
                                    SD3_SD56_RS2_RS6_Down_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Unclustered"), times = c(1, 2)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

SD3_SD56_RS2_RS6_Down_Plot_Condensed <- ggplot(
  data = SD3_SD56_RS2_RS6_Down_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,28), breaks = c(7,14,21,28), range = c(1,6)) +
  shades::lightness(scale_colour_distiller(palette = "Purples",
                                           limits = c(0,0.05)), scalefac(.8)) +
  coord_cartesian(xlim = c(1,10)) +
  theme_light() +
  downregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL)

# Save overlaid plot of up/down
# pdf("SD3_SD56_RS2_RS6_Bubble_Plot_Updated.pdf", width = 6.6985, height = 2.663)
# grid.newpage()
# grid.draw(rbind(ggplotGrob(SD3_SD56_RS2_RS6_Up_Plot_Condensed),
#                 ggplotGrob(SD3_SD56_RS2_RS6_Down_Plot_Condensed),
#                 size = "first"))
# dev.off()


### SD56 + (RS2 + RS6) Bubble Plot ----

SD56_RS2_RS6_Up_Plot <- SD56_RS2_RS6_Up_Unclustered_DF %>%
  mutate(Group = rep("Unclustered", times = 7),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

SD56_RS2_RS6_Up_Plot_Condensed <- ggplot(
  data = SD56_RS2_RS6_Up_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,28), breaks = c(7,14,21,28), range = c(1,6)) +
  shades::lightness(scale_colour_distiller(palette = "Greens",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,10)) +
  theme_light() +
  upregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "", y = NULL)

SD56_RS2_RS6_Down_Plot <- rbind(SD56_RS2_RS6_Down_Clustered_DF,
                                SD56_RS2_RS6_Down_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Unclustered"), times = c(1, 4)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

SD56_RS2_RS6_Down_Plot_Condensed <- ggplot(
  data = SD56_RS2_RS6_Down_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,28), breaks = c(7,14,21,28), range = c(1,6)) +
  shades::lightness(scale_colour_distiller(palette = "Purples",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,10)) +
  theme_light() +
  downregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL)

# Save overlaid plot of up/down
# pdf("SD56_RS2_RS6_Bubble_Plot_Updated.pdf", width = 7.39, height = 3.18)
# grid.newpage()
# grid.draw(rbind(ggplotGrob(SD56_RS2_RS6_Up_Plot_Condensed),
#                 ggplotGrob(SD56_RS2_RS6_Down_Plot_Condensed),
#                 size = "first"))
# dev.off()


### ((SD3 + SD56) + (RS2 + RS6)) OR (SD56 + (RS2 + RS6)) Bubble Plot ----

NR_Union_Up_Plot <- rbind(NR_Union_Up_Clustered_DF,
                          NR_Union_Up_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Cluster2", "Unclustered"),
                     times = c(1, 1, 4)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

NR_Union_Up_Plot_Condensed <- ggplot(
  data = NR_Union_Up_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,60), breaks = c(15,30,45,60), range = c(1,6.2)) +
  shades::lightness(scale_colour_distiller(palette = "Greens",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,11)) +
  theme_light() +
  upregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "", y = NULL)

NR_Union_Down_Plot <- rbind(NR_Union_Down_Clustered_DF,
                            NR_Union_Down_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Unclustered"),
                     times = c(1, 6)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

NR_Union_Down_Plot_Condensed <- ggplot(
  data = NR_Union_Down_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,60), breaks = c(15,30,45,60), range = c(1,6.2)) +
  shades::lightness(scale_colour_distiller(palette = "Purples",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,11)) +
  theme_light() +
  downregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL)

# Save overlaid plot of up/down
# pdf("Figures/Non_Recoverers_Bubble_Plot.pdf", width = 7.251, height = 3.557)
# grid.newpage()
# grid.draw(rbind(ggplotGrob(NR_Union_Up_Plot_Condensed),
#                 ggplotGrob(NR_Union_Down_Plot_Condensed),
#                 size = "first"))
# dev.off()


### RS2 Bubble Plot ----

RS2_Up_Plot <- rbind(RS2_Up_Clustered_DF,
                     RS2_Up_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Cluster2", "Unclustered"),
                     times = c(1, 1, 9)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

RS2_Up_Plot_Condensed <- ggplot(
  data = RS2_Up_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,60), breaks = c(15,30,45,60), range = c(1,6)) +
  shades::lightness(scale_colour_distiller(palette = "Greens",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,10)) +
  theme_light() +
  upregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "", y = NULL)

RS2_Down_Plot <- rbind(RS2_Down_Clustered_DF,
                       RS2_Down_Unclustered_DF) %>%
  mutate(Group = rep(c("Cluster1", "Unclustered"),
                     times = c(1, 6)),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

RS2_Down_Plot_Condensed <- ggplot(
  data = RS2_Down_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,60), breaks = c(15,30,45,60), range = c(1,6)) +
  shades::lightness(scale_colour_distiller(palette = "Purples",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,10)) +
  theme_light() +
  downregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL)

# Save overlaid plot of up/down
# pdf("RS2_Bubble_Plot_Updated.pdf", width = 6.6, height = 4.522)
# grid.newpage()
# grid.draw(rbind(ggplotGrob(RS2_Up_Plot_Condensed),
#                 ggplotGrob(RS2_Down_Plot_Condensed),
#                 size = "first"))
# dev.off()


### RS6 Bubble Plot ----

RS6_Up_Plot <- RS6_Up_Unclustered_DF %>%
  mutate(Group = rep("Unclustered", times = 15),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

RS6_Up_Plot_Condensed <- ggplot(
  data = RS6_Up_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,60), breaks = c(15,30,45,60), range = c(1,6)) +
  shades::lightness(scale_colour_distiller(palette = "Greens",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,10)) +
  theme_light() +
  upregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "", y = NULL)

RS6_Down_Plot <- RS6_Down_Unclustered_DF %>%
  mutate(Group = rep("Unclustered", times = 2),
         Order = as.integer(Order),
         Gene_Count = as.numeric(Gene_Count),
         Enrichment = as.numeric(Enrichment),
         p_Value = as.numeric(p_Value))

RS6_Down_Plot_Condensed <- ggplot(
  data = RS6_Down_Plot,
  mapping = aes(x = Enrichment,
                y = reorder(Terms, -Order),
                size = Gene_Count,
                color = p_Value)) +
  geom_point(alpha = 1) +
  scale_size_continuous(limits = c(1,60), breaks = c(15,30,45,60), range = c(1,6)) +
  shades::lightness(scale_colour_distiller(palette = "Purples",
                                           limits = c(0,0.05)), scalefac(0.7)) +
  coord_cartesian(xlim = c(1,10)) +
  theme_light() +
  downregulated_text_theme +
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL)

# Save overlaid plot of up/down
# pdf("RS6_Bubble_Plot_Updated.pdf", width = 8.0386, height = 3.9999)
# grid.newpage()
# grid.draw(rbind(ggplotGrob(RS6_Up_Plot_Condensed),
#                 ggplotGrob(RS6_Down_Plot_Condensed),
#                 size = "first"))
# dev.off()
