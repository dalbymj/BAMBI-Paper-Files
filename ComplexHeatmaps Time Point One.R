
# source('http://bioconductor.org/biocLite.R')
# biocLite("tidyverse")
# biocLite("stringi")
# biocLite("survival")
# biocLite("phyloseq")
# biocLite("cluster")
# biocLite("ComplexHeatmap")
# biocLite("dplyr")

# install.packages("ggplot2")
# install.packages("vegan")

library(stringi)
library(survival)
library(phyloseq)
library(cluster)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(vegan)
library(circlize)
library(tidyverse)


###################################################################
# Heatmap for time point 1
###################################################################
# First merge, select, and filter the data into the correct format.

# Set the working directory to import this file.
setwd("C:/Users/dalby/OneDrive/Documents/BAMBI/Test data")
# Import the recoded metadata.
meta_data <- read.csv("recoded_metadata.csv", header=TRUE, sep=",", strip.white = TRUE, na.strings=c("NR", "NA", "?")) 
# Remove the first colomn of numbers that are added when importing .csv files.
meta_data <- meta_data[c(2:ncol(meta_data))]

# Reanme the Treament factors Control and Probiotic as "Control" and "Bif/Lacto".
levels(meta_data$Treatment) <- c("Control", "Bif/Lacto")

# Set the working directory to import this file.
setwd("C:/Users/dalby/OneDrive/Documents/BAMBI/All Plate Files")
# Import the normalised sequence data.
seq_data <- read.csv("mergedplates_transposed_subsampled.csv", header=TRUE, sep=",", strip.white = TRUE, na.strings=c("NR", "NA", "?"))
# Remove the first two colomns from this file.
seq_data1 <- seq_data[c(3:ncol(seq_data))]
# Convert the data in the file into proportions based on the calculated row totals.
seq_data.prop <- seq_data1/rowSums(seq_data1)

# Selct the Lane_ID column from the sequence data file into its own data frame.
lane_id <- subset(seq_data, select=c(Lane_ID))
# Combine the Land_ID column with the data frame containing the sequence data as proportions.
lane_id_seq_data.prop <- cbind(lane_id, seq_data.prop)
# Merge the sequence and metadata together useing the shared Lane_ID colulmn.
prop_all_metadata_merged <- merge(meta_data, lane_id_seq_data.prop, by="Lane_ID")

# Select data
seqs_subsampled_all <- dplyr::select(filter(prop_all_metadata_merged,
                                            sample_age %in% c(1)),
                                     c(Lane_ID, Acinetobacter:ncol(prop_all_metadata_merged)))

##Selct the Lane_ID column from the sequence data file into its own data frame.
Lane_ID <- subset(seqs_subsampled_all, select=c(Lane_ID))
seqs_subsampled_all <- seqs_subsampled_all[c(2:ncol(seqs_subsampled_all))]


# Transpose the selected data so that the genus names become row names.
seqs_subsampled_all_t <- t(seqs_subsampled_all)

# Convert back to a dataframe.
seqs_subsampled_all_t <- as.data.frame(seqs_subsampled_all_t)

seqs_subsampled_all_t <- tibble::rownames_to_column(seqs_subsampled_all_t, var = "rowname")

rownames <- subset(seqs_subsampled_all_t, select=c(rowname))

seqs_subsampled_all_t <- seqs_subsampled_all_t[c(2:ncol(seqs_subsampled_all_t))]

# Calculate the sum of each row for each genus.
seqs_subsampled_all_t <- seqs_subsampled_all_t %>% mutate(rowsum = rowSums(.))

rowsums <- dplyr::select(seqs_subsampled_all_t, c(rowsum))

seqs_subsampled_totals <- cbind(rownames, rowsums)

seqs_subsampled_sorted <- seqs_subsampled_totals[order(seqs_subsampled_totals$rowsum, decreasing=TRUE),] 

##Set the first column as rownames.
rownames(seqs_subsampled_sorted) <- seqs_subsampled_sorted[,1]
seqs_subsampled_sorted$rowname = NULL

# Transpose the selected data so that the genus names become row names.
seqs_subsampled_sorted_t <- t(seqs_subsampled_sorted)
seqs_subsampled_sorted_t <- as.data.frame(seqs_subsampled_sorted_t)

top_15 <- dplyr::select(seqs_subsampled_sorted_t, c(1:10))

n1 <- names(top_15)
seqs_subsampled_15 <- seqs_subsampled_all[, which(names(seqs_subsampled_all) %in% n1)]

##Combine the Land_ID column with the data frame containing the sequence data as proportions.
seqs_subsampled_15 <- cbind(Lane_ID, seqs_subsampled_15)


prop_metadata_merged <- merge(meta_data, seqs_subsampled_15, by="Lane_ID")

# Select the data
seqs_subsampled <- dplyr::select(filter(prop_metadata_merged,
                                        sample_age %in% c(1)),
                                 c(sample_ID, 33:ncol(prop_metadata_merged)))

seqs_subsampled <- column_to_rownames(seqs_subsampled, var = "sample_ID")

seqs_subsampled <- seqs_subsampled[, 1:ncol(seqs_subsampled)]*100

meta_subsampled <- dplyr::select(filter(prop_metadata_merged,
                                        sample_age %in% c(1)),
                                 c(sample_ID, Treatment))
meta_subsampled <- column_to_rownames(meta_subsampled, var = "sample_ID")

##################################################################
# Use the merged, selected, and filtered data to make the heatmap.

# Use all sequence data to calculate the Bray Curtis matrix.
seqs_subsampled_all2 <- seqs_subsampled_all[c(2:ncol(seqs_subsampled_all))]
# Calculate the Bray-Curtis dissimilarity matrix on the sample rows.
seqs_subsampled_all_dist <- vegdist(seqs_subsampled_all2, method = "bray")
# Cluster the rows by hierarchical clustering.
row_clus <- hclust(seqs_subsampled_all_dist, "average")
# Calculate the Bray-Curtis dissimilarity matrix on the genus columns.
tp_one_dist_g <- vegdist(t(seqs_subsampled), method = "bray")
# Cluster the columns by hierarchical clustering.
col_clus <- hclust(tp_one_dist_g, "average")

# Annotation data frame
# Define the groups for the colour bar.
annot_df <- data.frame(Treatment = meta_subsampled$Treatment)
# Define colors for each groups.
col = list(Treatment = c("Bif/Lacto" = "lightskyblue1", "Control" = "red1"))

# Select the column of Bifidobacteria proportions from the sequence table.
Bifidobacteria <- data.frame(Bifidobacterium = seqs_subsampled$Bifidobacterium)

# Define the colour palette for the heatmap colours.
scalered <- colorRampPalette(colors=c("white",
                                      "orange",
                                      "darkorange",
                                      "red",
                                      "darkred",
                                      "black"))(100)

# Create the coloured side bar group annotation for the heatmap.
annotation1 <- rowAnnotation(df = annot_df, col = col,
                             annotation_width = unit(c(.25), "cm"),
                             annotation_legend_param = list(title = "Group",
                                                            title_gp = gpar(fontsize = 7),
                                                            labels_gp = gpar(fontsize = 7)))


# Create the Bifidobacteria barplot annotation for the heatmap.
annotation2 <- rowAnnotation("Bifidobacterium (%)" = anno_barplot(as.matrix(Bifidobacteria),
                                                   which = "row",
                                                   axis = TRUE,
                                                   bar_width = 0.05,
                                                   title_gp = gpar(fontsize = 7),
                                                   axis_gp = gpar(fontsize = 7)),
                     show_annotation_name = TRUE,
                     annotation_name_gp = gpar(fontsize = 7, fontface = "italic"),
                     annotation_name_rot = c(0),
                     annotation_name_offset = unit(10, "mm"),
                     annotation_width = unit(c(1.25), "cm"))

# Create the heatmap.
heatmap1 <- Heatmap(seqs_subsampled,
                    name = "Proportion",
                    col = scalered,
               cluster_rows = row_clus,
               cluster_columns = col_clus,
               show_row_names = FALSE,
               row_names_gp = gpar(fontsize = 2),
               column_names_gp = gpar(fontsize = 7, fontface = "italic"),
               column_title_gp = gpar(fontsize = 7),
               row_title = "Infant samples",
               row_title_gp = gpar(fontsize = 7),
               column_title = "Genus",
               column_title_side = "bottom",
               heatmap_legend_param = list(title = "Proportional\nabundance (%)",
                                           at = c(0,20,40,60,80,100),
                                           labels = c("0","20","40","60","80","100"),
                                           title_gp = gpar(fontsize = 7),
                                           labels_gp = gpar(fontsize = 7)))

# Combine the heatmap and the annotation together in the order in which they are to appear.
p <- heatmap1 + annotation1 + annotation2
p


setwd("C:/Users/dalby/OneDrive/Documents/BAMBI/ComplexHeatmaps")
pdf("Time point one 26.10.2018.pdf", width = 3.7, height = 6)
print(p)
dev.off()

