# install.packages("ComplexHeatmap")
library(ComplexHeatmap)
# install.packages("vegan")
library(vegan)
# install.packages("dplyr")
library(dplyr)
# install.packages("tidyverse")
library(tidyverse)

# Heatmap for time point 3 (samples from 30-49 days of age).

# First merge, select, and filter the data into the correct format.

# Set the working directory to import the metadata
setwd("C:/Users/dalby/OneDrive/Documents/BAMBI/Test data")
# Import the recoded metadata
meta_data <- read.csv("recoded_metadata 01.12.2018.csv", header=TRUE, sep=",", strip.white = TRUE, na.strings=c("NR", "NA", "?")) 
# Remove the first colomn of numbers that are added when importing .csv files.
meta_data <- meta_data[c(2:ncol(meta_data))]

# Reanme the Treament factors Control and Probiotic as "Control" and "Bif/Lacto".
levels(meta_data$Treatment) <- c("Control", "Bif/Lacto")

# Set the working directory to import the sequence data file.
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
prop_all_metadata_merged <- merge(meta_data,
                                  lane_id_seq_data.prop,
                                  by="Lane_ID")

# Select data
seqs_subsampled_all <- dplyr::select(filter(prop_all_metadata_merged,
                                            sample_age %in% c(3)), # Select only rows from Time point 2
                                     c(Lane_ID, # Selet the Lane_ID column
                                       Acinetobacter:ncol(prop_all_metadata_merged))) # Select all sequence columns

##Selct the Lane_ID column from the sequence data file into its own data frame
Lane_ID <- subset(seqs_subsampled_all, select=c(Lane_ID))
seqs_subsampled_all <- seqs_subsampled_all[c(2:ncol(seqs_subsampled_all))]

# Transpose the selected data so that the genus names become row names
seqs_subsampled_all_t <- t(seqs_subsampled_all)
# Convert the matrix back into a dataframe
seqs_subsampled_all_t <- as.data.frame(seqs_subsampled_all_t)
# Make the row names into a column
seqs_subsampled_all_t <- tibble::rownames_to_column(seqs_subsampled_all_t, var = "rowname")
# Subset the rownames colum into its own dataframe
rownames <- subset(seqs_subsampled_all_t, select=c(rowname))
# Subset the sequence data into its own dataframe
seqs_subsampled_all_t <- seqs_subsampled_all_t[c(2:ncol(seqs_subsampled_all_t))]

# Select the top 10 genus by abundance at time point 1

# Calculate the sum of each row for each genus
seqs_subsampled_all_t <- seqs_subsampled_all_t %>% mutate(rowsum = rowSums(.))
# Subset the rowsum colum into its own dataframe
rowsums <- dplyr::select(seqs_subsampled_all_t, c(rowsum))
# Bind together the rownames and rowsums columns
seqs_subsampled_totals <- cbind(rownames, rowsums)
# Sort the rowsums by decreasing size to use this to select the top ten genus
seqs_subsampled_sorted <- seqs_subsampled_totals[order(seqs_subsampled_totals$rowsum, decreasing=TRUE),] 

# Set the first rownames column as rownames
rownames(seqs_subsampled_sorted) <- seqs_subsampled_sorted[,1]
# Remove the rownames colum
seqs_subsampled_sorted$rowname = NULL

# Transpose the selected data so that the genus names become row names.
seqs_subsampled_sorted_t <- t(seqs_subsampled_sorted)
# Convert the transposed data back into a dataframe
seqs_subsampled_sorted_t <- as.data.frame(seqs_subsampled_sorted_t)

# Select the first 10 columns into a new dataframe (these are the top 10 genus by abundance)
top_10 <- dplyr::select(seqs_subsampled_sorted_t, c(1:10))
# Select just the names of the top 10 genus
names <- names(top_10)
# Use the top 10 genus names to select these genus from the main sequence data file
seqs_subsampled_10 <- seqs_subsampled_all[, which(names(seqs_subsampled_all) %in% names)]

# Combine the Lane_ID column with the data frame containing the top 10 genus sequence data
seqs_subsampled_10 <- cbind(Lane_ID, seqs_subsampled_10)

# Merge the Meta data and top 10 genus sequence sequence data by the shared Lane_ID column
prop_metadata_merged <- merge(meta_data, seqs_subsampled_10, by="Lane_ID")

# Select the sequence data to use
seqs_subsampled_10 <- dplyr::select(filter(prop_metadata_merged,
                                        sample_age %in% c(3)), # Select only rows from Time point 2
                                 c(sample_ID, 36:ncol(prop_metadata_merged))) # Select all sequence columns

# Make the sample_ID column into rownames
seqs_subsampled_10 <- column_to_rownames(seqs_subsampled_10, var = "sample_ID")
seqs_subsampled_10 <- seqs_subsampled_10[, 1:ncol(seqs_subsampled_10)]*100

# Select the metadata to use
meta_subsampled <- dplyr::select(filter(prop_metadata_merged,
                                        sample_age %in% c(3)), #Select only rows from Time point 2
                                 c(sample_ID, Treatment)) # Select the sample_ID and Treatment columns

# Make the sample_ID column into rownames
meta_subsampled <- column_to_rownames(meta_subsampled, var = "sample_ID")


# Use the merged, selected, and filtered data to calculate the heatmap

# Use all sequence data to calculate the Bray Curtis matrix
seqs_subsampled_all2 <- seqs_subsampled_all[c(2:ncol(seqs_subsampled_all))]
# Calculate the Bray-Curtis dissimilarity matrix on the sample rows
seqs_subsampled_all_dist <- vegdist(seqs_subsampled_all2, method = "bray")
# Cluster the rows by hierarchical clustering
row_clustering <- hclust(seqs_subsampled_all_dist, "average")
# Calculate the Bray-Curtis dissimilarity matrix on the genus columns
tp_one_dist_g <- vegdist(t(seqs_subsampled_10), method = "bray")
# Cluster the columns by hierarchical clustering
col_clustering <- hclust(tp_one_dist_g, "average")

# ComplexHeatmap Annotation

# Define the groups for the colour bar
annot_df <- data.frame(Treatment = meta_subsampled$Treatment)
# Define colors for each groups
col = list(Treatment = c("Bif/Lacto" = "lightskyblue1", "Control" = "red1"))

# Select the column of Bifidobacteria proportions from the sequence table
Bifidobacteria <- data.frame(Bifidobacterium = seqs_subsampled_10$Bifidobacterium)

# Define the colour palette for the heatmap colours
colour_palette <- colorRampPalette(colors=c("white",
                                            "orange",
                                            "darkorange",
                                            "red",
                                            "darkred",
                                            "black"))(100)

# Create the coloured side bar group annotation for the heatmap
sidebar_annotation <- rowAnnotation(df = annot_df, # Dataframe containing treatment groups
                                    col = col, # The list of treatment groups and their assigned colours
                                    annotation_width = unit(c(.25), "cm"), # Set the width of the side bar
                                    annotation_legend_param = list(title = "Group", # Sidebar legend title
                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size

# Create the Bifidobacteria barplot annotation for the heatmap
barplot_annotation <- rowAnnotation("Bifidobacterium (%)" = anno_barplot(as.matrix(Bifidobacteria), # Make the Bif dataframe into a matrix
                                    which = "row", # Use rows in the matrix to make the plot
                                    axis = TRUE, # Turns on axis lables
                                    bar_width = 0.05, # Width of each bar in barplot
                                    title_gp = gpar(fontsize = 7), # Axis title font size
                                    axis_gp = gpar(fontsize = 7)), # Axis labels font size
                                    show_annotation_name = TRUE, # Turns on axis name
                                    annotation_name_gp = gpar(fontsize = 7, # Axis title font size
                                                              fontface = "italic"), # Axis title in italics
                                    annotation_name_rot = c(0), # Rotation of axis title
                                    annotation_name_offset = unit(10, "mm"), # Axis title distance from barplot
                                    annotation_width = unit(c(1.25), "cm")) # Width of the barplot

# Create the heatmap
heatmap <- Heatmap(seqs_subsampled_10, # The dataframe containing the heatmap data
                   name = "Proportion", # Name is used as the title of the heatmap legend if shown
                   col = colour_palette, # The predefined colour palette
                   cluster_rows = row_clustering, # Cluster the rows using the predefined clustering
                   cluster_columns = col_clustering, # Cluster the columns using the predefined clustering
                   show_row_names = FALSE, # Show or hide the row names, TRUE to show rownames
                   row_names_gp = gpar(fontsize = 2), # Row name font size
                   column_names_gp = gpar(fontsize = 7, # Column name font size
                                          fontface = "italic"), # Column names in italics
                   column_title_gp = gpar(fontsize = 7), # Column title font size
                   row_title = "Infant samples", # Set row title
                   row_title_gp = gpar(fontsize = 7), # Set row title font size
                   column_title = "Genus", # Set column title
                   column_title_side = "bottom", # Set column title font size
                   heatmap_legend_param = list(title = "Proportional\nabundance (%)", # Set legend title
                                               at = c(0,20,40,60,80,100), # Set legend scale breaks
                                               labels = c("0","20","40","60","80","100"), # Set legend scale labels
                                               title_gp = gpar(fontsize = 7), # Set legend title font size
                                               labels_gp = gpar(fontsize = 7))) # Set legend label font size

# Combine the heatmap and the annotation together in the order in which they are to appear
p <- heatmap + sidebar_annotation + barplot_annotation
p

# Set the file to save the figure into
setwd("C:/Users/dalby/OneDrive/Documents/BAMBI/ComplexHeatmaps")
# Set the saved pdf file name and define the size of the plot
pdf("ComplexHeatmap time point three 30.11.2018.pdf", width = 3.7, height = 6)
# pring(p) saves the figure into a file
print(p)
dev.off()
