#This script integrates 480K SNP data from Denance et al. 2019, with the Jim Dunckley and Plant & Food Research samples to find duplicates with PLINK.
#The "JD_PFR_All.ped" and "JD_PFR_All.map" from the Howard_2021 analysis are reused; these are filtered for diploids only and have SNP locations added.

#Load packages
library(tibble)
library(igraph)

library(dplyr)
library(tidyr)

#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/20K_480K PLINK Duplicate Identification/Inputs/Denance_2019")

#Use PLINK to extract overlapping SNPs from JD and PFR samples
system("plink --file JD_PFR_All --extract 50K_480K_extract.txt --make-bed --out 480K_JD_PFR")
system("plink --bfile 480K_JD_PFR --recode --tab --out 480K_JD_PFR")

#Use PLINK to extract overlapping SNPs from 480K Denance et al. 2019 samples, recode SNP positions, and output as .ped and .map
system("plink --bfile FruitBreedomics_apple_320K_SNP --extract 50K_480K_extract.txt --recode tab --out 480K_Samples")
system("plink --file 480K_Samples --update-chr 480K_chr.txt --update-cm 480K_cm.txt --update-map 480K_map.txt --make-bed --out 480K_Samples")
system("plink --bfile 480K_Samples --recode tab --out 480K_Samples")



##Combining Denance et al. 2019 480K data with JD_PFR samples

#clear workspace
rm(list=ls())

#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/20K_480K PLINK Duplicate Identification/Inputs/Denance_2019")

#Load JD_PFR data
JD_ped <- read.csv("480K_JD_PFR.ped", header = FALSE,sep = "\t")

#Load Denance 480K data, add prefix to names, set irrelevant tags to 0
DN_ped <- read.csv("480K_Samples.ped", header = FALSE,sep = "\t")
DN_ped$V2 <- paste0("DN_", DN_ped$V2) 
DN_ped$V1 <- 0
DN_ped$V3 <- 0
DN_ped$V4 <- 0

#Bind Denance 480K data onto PLINK .ped, set other variables to 0
names(DN_ped) <- names(JD_ped)
combined_ped <- bind_rows(JD_ped, DN_ped)

#Load .map file
map <- read.csv("480K_JD_PFR.map", header = FALSE, sep ="\t")

#Save .map file and combined .ped file with same base
write.table(map, "480K_PLINK.map", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(combined_ped, "480K_PLINK.ped", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


##Running PLINK Duplicate Analysis

#clear workspace
rm(list=ls())

#set working directory [must contain plink.exe and files for analysis]
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/20K_480K PLINK Duplicate Identification/Inputs/Denance_2019")

#Run PLINK
system("plink --file 480K_PLINK --missing-genotype 0 --genome full ")

#Read genome file
genome <- read.table("plink.genome", header = TRUE, sep = "", stringsAsFactors = FALSE)
write.table(genome, "C:/Users/curly/Desktop/Apple Genotyping/Results/20K_480K PLINK Duplicate Identification/480K Results/PLINK_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

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
colnames(dd) <- c("Group", "SampleCount", "ID1","ID2","ID3","ID4","ID5","ID6","ID7","ID8","ID9","ID10","ID11","ID12","ID13", "ID14")

#Save .csv of duplicate groupings
write.csv(dd, "C:/Users/curly/Desktop/Apple Genotyping/Results/20K_480K PLINK Duplicate Identification/480K Results/Grouped_Duplicates.csv", row.names = FALSE)


