### UpSetPlot.R
## Author: Alexander Popescu (July 2023, last updated May 2024)

## Adapted from Elena Zuin's code for the scRNA-seq pipeline

## Use ComplexUpset (https://github.com/krassowski/complex-upset) to visualize intersections of differentially expressed genes (DEGs) across experimental conditions with UpSet plots, intersecting upregulated and downregulated genes separately to maximize biological relevance


### Packages ----

library(ggplot2) # v 3.5.1
library(UpSetR) # v 1.4.0
library(ComplexUpset) # v 1.3.3
library(patchwork) # v 1.2.0
library(readxl) # v 1.4.3
library(gtools) # v 3.9.5

# sink('sessionInfo_UpSetPlot.txt')
# sessionInfo()
# sink()


### Import DGE Lists ----

## DEGs meet the double threshold: FDR < 0.05, log2 CPM > 0

file <- "Supplemental_Table_S3_log2CPM_Cut.xlsx"

SD3_Up <- read_xlsx(file, sheet = 1)
deg.sd3.up <- SD3_Up$`ENSEMBL ID`  # 939

SD3_Down <- read_xlsx(file, sheet = 2)
deg.sd3.down <- SD3_Down$`ENSEMBL ID`  # 1373

SD56_Up <- read_xlsx(file, sheet = 3)
deg.sd56.up <- SD56_Up$`ENSEMBL ID`  # 3095

SD56_Down <- read_xlsx(file, sheet = 4)
deg.sd56.down <- SD56_Down$`ENSEMBL ID`  # 4390

RS2_Up <- read_xlsx(file, sheet = 5)
deg.rs2.up <- RS2_Up$`ENSEMBL ID`  # 1782

RS2_Down <- read_xlsx(file, sheet = 6)
deg.rs2.down <- RS2_Down$`ENSEMBL ID`  # 2123

RS6_Up <- read_xlsx(file, sheet = 7)
deg.rs6.up <- RS6_Up$`ENSEMBL ID`  # 1103

RS6_Down <- read_xlsx(file, sheet = 8)
deg.rs6.down <- RS6_Down$`ENSEMBL ID`  # 885

all_Up <- rbind(SD3_Up, SD56_Up, RS2_Up, RS6_Up) # 6919
all_Down <- rbind(SD3_Down, SD56_Down, RS2_Down, RS6_Down) # 8771

all_DEG <- rbind(all_Up, all_Down)
id <- unique(all_DEG$`ENSEMBL ID`)
duplicated_rows <- duplicated(all_DEG$`ENSEMBL ID`)

all_DEG <- all_DEG[!duplicated_rows |
                     (duplicated_rows & !all_DEG$`ENSEMBL ID` %in% id), ] # 9015

## Create two data frames representing intersection matrices

lt1 <- list("SD3" = deg.sd3.up, "SD5-6" = deg.sd56.up,
            "RS2" = deg.rs2.up, "RS6" = deg.rs6.up)
DEG_Up <- fromList(lt1)

lt2 <- list("SD3" = deg.sd3.down, "SD5-6" = deg.sd56.down,
            "RS2" = deg.rs2.down, "RS6" = deg.rs6.down)
DEG_Down <- fromList(lt2)


### Upregulated DEGs, Main Figure ----

plot_up_main <- upset(data = DEG_Up,
            intersect = c("RS6", "RS2", "SD5-6", "SD3"),
            sort_sets = FALSE,  # keep bottom-top order specified by intersect
            min_size = 100,  # set min size for intersection to be shown
            
            # remove SD56/RS6 and specify order of intersections
            sort_intersections = FALSE,
            intersections = list(
              "SD5-6", c("SD3", "SD5-6"),
              c("SD5-6", "RS2"), c("SD3", "SD5-6", "RS2"),
              c("SD5-6", "RS2", "RS6"), c("SD3", "SD5-6", "RS2", "RS6"),
              "RS2", "RS6"),
            name = "",
            keep_empty_groups = TRUE,
            
            # selectively modify aesthetics
            queries = list(upset_query(set = "SD3", fill = "#FF6666"),
                           upset_query(set = "SD5-6", fill = "#E61919"),
                           upset_query(set = "RS2", fill = "#6495ED"),
                           upset_query(set = "RS6", fill = "#336699")),
            
            base_annotations = list("Intersection Size" = (intersection_size(
              bar_number_threshold = 1,  # show all numbers on top of bars
              width = 0.5,   # reduce width of the bars
              text = list(size = 7, vjust = -0.7))   # enlarge numbers on top of the bars
              + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                   limits = c(0, 2000))  # set y-axis limits
              + theme(panel.grid.major = element_blank(),  # hide grid lines
                      panel.grid.minor = element_blank(),  # show axis lines
                      axis.line = element_line(colour = 'black'),
                      axis.title = element_text(size = 20)))),

            # customize stripes
            stripes = upset_stripes(geom = geom_segment(linewidth = 12),  # enlarge stripes
                                    colors = c("grey95", "white")),  # color stripes
            # customize intersection matrix
            matrix = intersection_matrix(geom = geom_point(shape = "circle filled",
                                                           size = 5,
                                                           stroke = 0)),
            
            # customize horizontal bar plot of set size
            set_sizes = (upset_set_size(geom = geom_bar(width = 0.4))
                         + geom_text(aes(label = after_stat(count)),
                                     hjust = -0.3, vjust = -1.5, stat = "count",
                                     color = "black", size = 6)
                         + theme(axis.line.x = element_line(colour = 'black'),
                                 axis.ticks.x = element_line(),
                                 axis.text = element_text(size = 16),
                                 axis.title = element_text(size = 20),
                                 plot.title = element_text(size = 20,
                                                           face = "bold"))
                         + ylab("set size"))
            + labs(title = "Upregulated Genes"),
            
            themes = upset_default_themes(text = element_text(size = 25, color = "black")),
            wrap = TRUE)

# Save plot
ggsave(plot_up_main, file = "Figure3_Up_Main.pdf",
       width = 30, height = 30, units = "cm")


### Upregulated DEGs, Supplemental Figure ----

plot_up_supp <- upset(data = DEG_Up,
            intersect = c("RS6", "RS2", "SD5-6", "SD3"),
            sort_sets = FALSE,
            min_size = 1,  # all non-empty intersections
            sort_intersections_by = c("degree", "cardinality"),
            sort_intersections = "ascending",
            name = "",
            keep_empty_groups = TRUE,
            queries = list(upset_query(set = "SD3", fill = "#FF6666"),
                           upset_query(set = "SD5-6", fill = "#E61919"),
                           upset_query(set = "RS2", fill = "#6495ED"),
                           upset_query(set = "RS6", fill = "#336699")),
            base_annotations = list("Intersection Size" = (intersection_size(
              bar_number_threshold = 1,
              width = 0.5,
              text = list(size = 7, vjust = -0.7))
              + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                   limits = c(0, 2000))
              + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = 'black'),
                      axis.title = element_text(size = 20)))),
            stripes = upset_stripes(geom = geom_segment(linewidth = 12),
                                    colors = c("grey95", "white")),
            matrix = intersection_matrix(geom = geom_point(shape = "circle filled", 
                                                           size = 5,
                                                           stroke = 0)),
            set_sizes = (upset_set_size(geom = geom_bar(width = 0.4))
                         + geom_text(aes(label = after_stat(count)),
                                     hjust = -0.3, vjust = -1.5, stat = "count",
                                     color = "black", size = 6)
                         + theme(axis.line.x = element_line(colour = 'black'),
                                 axis.ticks.x = element_line(),
                                 axis.text = element_text(size = 16),
                                 axis.title = element_text(size = 20),
                                 plot.title = element_text(size = 20,
                                                           face = "bold"))
                         + ylab("set size"))
            + labs(title = "Upregulated Genes"),
            themes = upset_default_themes(text = element_text(size = 25, color = "black")),
            wrap = TRUE)

# Save plot
# ggsave(plot_up_supp, file = "Figures/Figure3_Up_Supplement.pdf",
#        width = 30, height = 30, units = "cm")


### Downregulated DEGs, Main Figure ----

plot_down_main <- upset(DEG_Down,
            intersect = c("RS6", "RS2", "SD5-6", "SD3"),
            sort_sets = FALSE,
            min_size = 100,
            sort_intersections = FALSE,
            intersections = list(
              "SD5-6", c("SD3", "SD5-6"),
              c("SD5-6", "RS2"), c("SD3", "SD5-6", "RS2"),
              c("SD5-6", "RS2", "RS6"), c("SD3", "SD5-6", "RS2", "RS6"),
              "RS2", "RS6"),
            name = "",
            keep_empty_groups = TRUE,
            queries = list(upset_query(set = "SD3", fill = "#FF6666"),
                           upset_query(set = "SD5-6", fill = "#E61919"),
                           upset_query(set = "RS2", fill = "#6495ED"),
                           upset_query(set = "RS6", fill = "#336699")),
            base_annotations = list("Intersection Size" = (intersection_size(
              bar_number_threshold = 1,
              width = 0.5,
              text = list(size = 7, vjust = -0.7))
              + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                   limits = c(0, 2000))
              + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = 'black'),
                      axis.title = element_text(size = 20)))),
            stripes = upset_stripes(geom = geom_segment(linewidth = 12),
                                    colors = c("grey95", "white")),
            matrix = intersection_matrix(geom = geom_point(shape = "circle filled", 
                                                           size = 5,
                                                           stroke = 0)),
            set_sizes = (upset_set_size(geom = geom_bar(width = 0.4))
                         + geom_text(aes(label = after_stat(count)),
                                     hjust = -0.3, vjust = -1.5, stat = "count",
                                     color = "black", size = 6)
                         + theme(axis.line.x = element_line(colour = 'black'),
                                 axis.ticks.x = element_line(),
                                 axis.text = element_text(size = 16),
                                 axis.title = element_text(size = 20),
                                 plot.title = element_text(size = 20,
                                                           face = "bold"))
                         + ylab("set size"))
            + labs(title = "Downregulated Genes"),
            themes = upset_default_themes(text = element_text(size = 25, color = "black")),
            wrap = TRUE)

# Save plot
ggsave(plot_down_main, file = "Figure3_Down_Main.pdf",
       width = 30, height = 30, units = "cm")


### Downregulated DEGs, Supplemental Figure ----

plot_down_supp <- upset(DEG_Down,
            intersect = c("RS6", "RS2", "SD5-6", "SD3"),
            sort_sets = FALSE,
            min_size = 1,
            sort_intersections = FALSE,
            intersections = list(
              "SD3", "RS2", "RS6", "SD5-6", c("SD3", "RS2"), c("SD3", "RS6"),
              c("RS2", "RS6"), c("SD5-6", "RS6"), c("SD3", "SD5-6"),
              c("SD5-6", "RS2"), c("SD3", "SD5-6", "RS6"),
              c("SD5-6", "RS2", "RS6"), c("SD3", "SD5-6", "RS2"),
              c("SD3", "SD5-6", "RS2", "RS6"), c("SD3", "RS2", "RS6")),
            name = "",
            keep_empty_groups = TRUE,
            queries = list(upset_query(set = "SD3", fill = "#FF6666"),
                           upset_query(set = "SD5-6", fill = "#E61919"),
                           upset_query(set = "RS2", fill = "#6495ED"),
                           upset_query(set = "RS6", fill = "#336699")),
            base_annotations = list("Intersection Size" = (intersection_size(
              bar_number_threshold = 1,
              width = 0.5,
              text = list(size = 7, vjust = -0.7))
              + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                   limits = c(0, 2000))
              + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = 'black'),
                      axis.title = element_text(size = 20)))),
            stripes = upset_stripes(geom = geom_segment(linewidth = 12),
                                    colors = c("grey95", "white")),
            matrix = intersection_matrix(geom = geom_point(shape = "circle filled", 
                                                           size = 5,
                                                           stroke = 0)),
            set_sizes = (upset_set_size(geom = geom_bar(width = 0.4))
                         + geom_text(aes(label = after_stat(count)),
                                     hjust = -0.3, vjust = -1.5, stat = "count",
                                     color = "black", size = 6)
                         + theme(axis.line.x = element_line(colour = 'black'),
                                 axis.ticks.x = element_line(),
                                 axis.text = element_text(size = 16),
                                 axis.title = element_text(size = 20),
                                 plot.title = element_text(size = 20,
                                                           face = "bold"))
                         + ylab("set size"))
            + labs(title = "Downregulated Genes"),
            themes = upset_default_themes(text = element_text(size = 25, color = "black")),
            wrap = TRUE)

# Save plot
# ggsave(plot_down_supp, file = "Figures/Figure3_Down_Supplement.pdf",
#        width = 30, height = 30, units = "cm")


### Save Gene Lists for All Intersections ----

lst.u1 <- setdiff(deg.sd3.up, c(deg.sd56.up, deg.rs2.up, deg.rs6.up))  # 69
lst.u2 <- setdiff(deg.rs2.up, c(deg.sd3.up, deg.sd56.up, deg.rs6.up))  # 423
lst.u3 <- setdiff(deg.rs6.up, c(deg.sd3.up, deg.sd56.up, deg.rs2.up))  # 463
lst.u4 <- setdiff(deg.sd56.up, c(deg.sd3.up, deg.rs2.up, deg.rs6.up))  # 1149
lst.u5 <- setdiff(intersect(deg.sd3.up, deg.rs2.up), c(deg.sd56.up, deg.rs6.up))  # 3
lst.u6 <- setdiff(intersect(deg.sd3.up, deg.rs6.up), c(deg.sd56.up, deg.rs2.up))  # 15
lst.u7 <- setdiff(intersect(deg.rs2.up, deg.rs6.up), c(deg.sd3.up, deg.sd56.up))  # 52
lst.u8 <- setdiff(intersect(deg.sd56.up, deg.rs6.up), c(deg.sd3.up, deg.rs2.up))  # 182
lst.u9 <- setdiff(intersect(deg.sd3.up, deg.sd56.up), c(deg.rs2.up, deg.rs6.up))  # 407
lst.u10 <- setdiff(intersect(deg.sd56.up, deg.rs2.up), c(deg.sd3.up, deg.rs6.up))  # 692
lst.u11 <- setdiff(Reduce(intersect, list(deg.sd3.up, deg.sd56.up, deg.rs6.up)), deg.rs2.up)  # 53
lst.u12 <- setdiff(Reduce(intersect, list(deg.sd56.up, deg.rs2.up, deg.rs6.up)), deg.sd3.up)  # 220
lst.u13 <- setdiff(Reduce(intersect, list(deg.sd3.up, deg.sd56.up, deg.rs2.up)), deg.rs6.up)  # 274
lst.u14 <- Reduce(intersect, list(deg.sd3.up, deg.sd56.up, deg.rs2.up, deg.rs6.up))  # 118

lst.d1 <- setdiff(deg.sd3.down, c(deg.sd56.down, deg.rs2.down, deg.rs6.down))  # 159
lst.d2 <- setdiff(deg.rs2.down, c(deg.sd3.down, deg.sd56.down, deg.rs6.down))  # 356
lst.d3 <- setdiff(deg.rs6.down, c(deg.sd3.down, deg.sd56.down, deg.rs2.down))  # 335
lst.d4 <- setdiff(deg.sd56.down, c(deg.sd3.down, deg.rs2.down, deg.rs6.down))  # 1898
lst.d5 <- setdiff(intersect(deg.sd3.down, deg.rs2.down), c(deg.sd56.down, deg.rs6.down))  # 5
lst.d6 <- setdiff(intersect(deg.sd3.down, deg.rs6.down), c(deg.sd56.down, deg.rs2.down))  # 2
lst.d7 <- setdiff(intersect(deg.rs2.down, deg.rs6.down), c(deg.sd3.down, deg.sd56.down))  # 18
lst.d8 <- setdiff(intersect(deg.sd56.down, deg.rs6.down), c(deg.sd3.down, deg.rs2.down))  # 170
lst.d9 <- setdiff(intersect(deg.sd3.down, deg.sd56.down), c(deg.rs2.down, deg.rs6.down))  # 522
lst.d10 <- setdiff(intersect(deg.sd56.down, deg.rs2.down), c(deg.sd3.down, deg.rs6.down))  # 931
lst.d11 <- setdiff(Reduce(intersect, list(deg.sd3.down, deg.sd56.down, deg.rs6.down)), deg.rs2.down)  # 57
lst.d12 <- setdiff(Reduce(intersect, list(deg.sd56.down, deg.rs2.down, deg.rs6.down)), deg.sd3.down)  # 185
lst.d13 <- setdiff(Reduce(intersect, list(deg.sd3.down, deg.sd56.down, deg.rs2.down)), deg.rs6.down)  # 510
lst.d14 <- Reduce(intersect, list(deg.sd3.down, deg.sd56.down, deg.rs2.down, deg.rs6.down))  # 117
lst.d15 <- setdiff(Reduce(intersect, list(deg.sd3.down, deg.rs2.down, deg.rs6.down)), deg.sd56.down)  # 1

## Add back in annotations from original data frame

lst1.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u1,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst2.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u2,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst3.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u3,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst4.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u4,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst5.up.annot <- unique(
  all_Up[(all_Up$`ENSEMBL ID` %in% lst.u5),
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst6.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u6,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst7.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u7,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst8.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u8,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst9.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u9,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst10.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u10,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst11.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u11,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst12.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u12,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst13.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u13,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst14.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u14,
         c("ENSEMBL ID", "Gene description", "Gene name")])

lst1.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d1,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst2.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d2,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst3.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d3,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst4.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d4,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst5.down.annot <- unique(
  all_Down[(all_Down$`ENSEMBL ID` %in% lst.d5),
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst6.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d6,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst7.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d7,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst8.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d8,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst9.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d9,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst10.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d10,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst11.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d11,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst12.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d12,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst13.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d13,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst14.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d14,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst15.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d15,
           c("ENSEMBL ID", "Gene description", "Gene name")])

## Save intersections

# write.table(x = lst1.up.annot, file = "lst1-up.txt", sep = "\t")
# write.table(x = lst2.up.annot, file = "lst2-up.txt", sep = "\t")
# write.table(x = lst3.up.annot, file = "lst3-up.txt", sep = "\t")
# write.table(x = lst4.up.annot, file = "lst4-up.txt", sep = "\t")
# write.table(x = lst5.up.annot, file = "lst5-up.txt", sep = "\t")
# write.table(x = lst6.up.annot, file = "lst6-up.txt", sep = "\t")
# write.table(x = lst7.up.annot, file = "lst7-up.txt", sep = "\t")
# write.table(x = lst8.up.annot, file = "lst8-up.txt", sep = "\t")
# write.table(x = lst9.up.annot, file = "lst9-up.txt", sep = "\t")
# write.table(x = lst10.up.annot, file = "lst10-up.txt", sep = "\t")
# write.table(x = lst11.up.annot, file = "lst11-up.txt", sep = "\t")
# write.table(x = lst12.up.annot, file = "lst12-up.txt", sep = "\t")
# write.table(x = lst13.up.annot, file = "lst13-up.txt", sep = "\t")
# write.table(x = lst14.up.annot, file = "lst14-up.txt", sep = "\t")

# write.table(x = lst1.down.annot, file = "lst1-down.txt", sep = "\t")
# write.table(x = lst2.down.annot, file = "lst2-down.txt", sep = "\t")
# write.table(x = lst3.down.annot, file = "lst3-down.txt", sep = "\t")
# write.table(x = lst4.down.annot, file = "lst4-down.txt", sep = "\t")
# write.table(x = lst5.down.annot, file = "lst5-down.txt", sep = "\t")
# write.table(x = lst6.down.annot, file = "lst6-down.txt", sep = "\t")
# write.table(x = lst7.down.annot, file = "lst7-down.txt", sep = "\t")
# write.table(x = lst8.down.annot, file = "lst8-down.txt", sep = "\t")
# write.table(x = lst9.down.annot, file = "lst9-down.txt", sep = "\t")
# write.table(x = lst10.down.annot, file = "lst10-down.txt", sep = "\t")
# write.table(x = lst11.down.annot, file = "lst11-down.txt", sep = "\t")
# write.table(x = lst12.down.annot, file = "lst12-down.txt", sep = "\t")
# write.table(x = lst13.down.annot, file = "lst13-down.txt", sep = "\t")
# write.table(x = lst14.down.annot, file = "lst14-down.txt", sep = "\t")
# write.table(x = lst15.down.annot, file = "lst15-down.txt", sep = "\t")

