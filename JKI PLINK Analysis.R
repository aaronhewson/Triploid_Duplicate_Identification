# This script integrates data from the Julius Kuhn Institute with data from the Jim Dunckley Orchard and Plant & Food Research, to then find duplicates with PLINK.
# Prior to using this script, JKI data was formatted to match PLINK specifications, but still in need of transposing.

#Load packages
library(tibble)
library(igraph)

library(dplyr)
library(tidyr)

#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/Triploid Duplicate Identification/Inputs/JKI_Inputs")

#Load JKI data and transpose
JKI <- read.delim("JKI_Samples_Fixed.txt", header = FALSE, stringsAsFactors = FALSE)
JKI_t <- as.data.frame(t(JKI))

#Add empty columns for PLINK formatting
JKI_t <- add_column(JKI_t, Fa = 0, Mo = 0, Se = 0, Ph = 0, .before = "V2" )
JKI_t <- add_column(JKI_t, Fid = 0, .before = "V1" )

#Remove header row 
JKI_t = JKI_t[-1, ]

#Remove row and column names
colnames(JKI_t) <- NULL
rownames(JKI_t) <- NULL

#Load in PLINK .ped
ped <- read.csv("JD_PFR_PLINK.ped", header = FALSE,sep = "\t")

#Bind JKI data onto PLINK .ped
names(JKI_t) <- names(ped)
combined_ped <- bind_rows(ped, JKI_t)

#Load .map file
map <- read.csv("JD_PFR_PLINK.map", header = FALSE, sep ="\t")

#Save .map file and combined .ped file with same base
write.table(map, "JKI_PLINK.map", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(combined_ped, "JKI_PLINK.ped", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


##Running PLINK

#clear workspace
rm(list=ls())

#set working directory [must contain plink.exe and files for analysis]
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/Triploid Duplicate Identification/Inputs/JKI_Inputs")

#Run PLINK
system("plink --file JKI_PLINK --missing-genotype 0 --genome full ")

#Read genome file
genome <- read.table("plink.genome", header = TRUE, sep = "", stringsAsFactors = FALSE)
write.table(genome, "C:/Users/curly/Desktop/Apple Genotyping/Results/Triploid Duplicates/JKI_Duplicates/PLINK_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

##Grouping duplicates

#Filter for PI_HAT >0.96 (duplicate threshold)
genome <- genome[!(genome$PI_HAT < 0.96), ]
genome <- subset(genome, select = c("IID1","IID2"))

#Group duplicates with igraph
graph <- graph_from_data_frame(genome, directed = FALSE)
components <- components(graph)

#Sort groupings by number of duplicates
group_sizes <- table(components$membership)
sorted_group_ids <- order(group_sizes)
new_ids <- match(components$membership, sorted_group_ids)
V(graph)$group <- new_ids
grouped_samples <- split(names(components$membership), new_ids)

#Pad group with length less than max length with NA's
max_len <- max(sapply(grouped_samples, length))
padded_list <- lapply(grouped_samples, function(x) {c(x, rep(" ", max_len - length(x)))})

#Write groupings to a dataframe
dd <- as.data.frame(do.call(rbind, padded_list))

#Add a number for each group
dd <- cbind(Group = seq_len(nrow(dd)), dd)

# Add the number of duplicates in each grouping
sample_counts <- rowSums(dd[, -1] != " ")
dd <- add_column(dd, SampleCount = sample_counts, .after = "Group")

#Rename columns
colnames(dd) <- c("Group", "SampleCount", "ID1","ID2","ID3","ID4","ID5","ID6","ID7","ID8","ID9","ID10","ID11","ID12")

#Save .csv of duplicate groupings
write.csv(dd, "C:/Users/curly/Desktop/Apple Genotyping/Results/Triploid Duplicates/JKI_Duplicates/Grouped_Duplicates.csv", row.names = FALSE)

