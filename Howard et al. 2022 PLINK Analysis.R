#This script integrates 20K SNP data from Howard et al. 2022, with the Jim Dunckley and Plant & Food Research samples to find duplicates with PLINK.
#JD_PFR genotype files are re-exported from Axiom to have all SNPs, and the Howard 2022 data formatted as close to PLINK format as possible in Excel first, then saved as .txt.


#Load packages
library(tibble)
library(igraph)

library(dplyr)
library(tidyr)

#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/20K_480K PLINK Duplicate Identification/Inputs/Howard_2022")

# Initial JD_PFR .ped file curation ---------------------------------------
#Load .ped file and remove .CEL from sample filenames
ped <- read.csv("JD_PFR_Raw.ped", header = FALSE,sep = "\t")
ped[[1]] <- sub("\\.CEL$", "", ped[[1]])

#Remove triploids from .ped file
ids_to_remove <- read.delim("TriploidSampleNames.txt", header = FALSE, stringsAsFactors = FALSE)
ids_to_remove <- ids_to_remove$V1
ped <- ped[!(ped$V1 %in% ids_to_remove), ]

#Add empty columns for PLINK formatting
ped <- add_column(ped, Fa = 0, Mo = 0, Se = 0, Ph = 0, .before = "V2" )
ped <- add_column(ped, Fid = 0, .before = "V1" )

#Remove header row and save .ped file
ped = ped[-1, ]
write.table(ped, "JD_PFR_All.ped", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




# Initial JD_PFR .map file curation ---------------------------------------

#Load .map file and BLAST results
map <- read.csv("JD_PFR_Raw.map", header = FALSE, sep ="\t")
BLAST <- read.csv("BLAST results.tsv", header = FALSE, sep = "\t")

#Match Marker IDs
match_ids <- match(map$V2, BLAST$V2)
matched <- !is.na(match_ids)
map[matched, ] <- BLAST[match_ids[matched], ]

#Remove header row, and set any chromosome numbers larger than 17 to 0
map = map[-1, ]
map$V1 <- as.numeric(as.character(map$V1))
map[map$V1 > 17, c("V1", "V4")] <- 0

#Save .map file (with locations)
write.table(map, "JD_PFR_All.map", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




# Initial Howard .ped file curation ---------------------------------------
#Load .ped file
HW_ped <- read.csv("Howard_2022_Raw.ped", header = FALSE,sep = "\t")


#Add empty columns for PLINK formatting
HW_ped <- add_column(HW_ped, Fa = 0, Mo = 0, Se = 0, Ph = 0, .before = "V2" )
HW_ped <- add_column(HW_ped, Fid = 0, .before = "V1" )

#Remove header row and save .ped file
HW_ped = HW_ped[-1, ]
write.table(HW_ped, "Howard_All.ped", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




# Initial Howard .map file curation ---------------------------------------

#Load .map file and BLAST results
map <- read.csv("Howard_2022_Raw.map", header = FALSE, sep ="\t")
BLAST <- read.csv("BLAST results Howard.tsv", header = FALSE, sep = "\t")

#Match Marker IDs
match_ids <- match(map$V2, BLAST$V2)
matched <- !is.na(match_ids)
map[matched, ] <- BLAST[match_ids[matched], ]

#Remove header row, and set any chromosome numbers larger than 17 to 0
map = map[-1, ]
map$V1 <- as.numeric(as.character(map$V1))
map[map$V1 > 17, c("V1", "V4")] <- 0

#Save .map file (with locations)
write.table(map, "Howard_All.map", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# PLINK formatting --------------------------------------------------------


#Use PLINK to extract SNPs and recode alleles for JD_PFR samples
system("plink --file JD_PFR_All --extract Extract_SNP.txt --update-alleles Recode_50K.txt --make-bed --out JD_PFR_Ready")
system("plink --bfile JD_PFR_Ready --recode --tab --out JD_PFR_Ready")


#Use PLINK to rename SNPs for Howard 2022 data, and extract SNPS
system("plink --file Howard_All --update-name Howard_rename.txt --make-bed --out Howard_All")
system("plink --bfile Howard_All --extract Extract_SNP.txt  --recode tab --out Howard_Ready")



# Combine JD_PFR and Howard data ------------------------------------------

#Load in JD and Howard .ped files
JDped <- read.csv("JD_PFR_Ready.ped", header = FALSE,sep = "\t")
HWped <- read.csv("Howard_Ready.ped", header = FALSE,sep = "\t")

#Bind Howard 20K data onto PLINK .ped
names(HWped) <- names(JDped)
combined_ped <- bind_rows(JDped, HWped)

#Load .map file
map <- read.csv("JD_PFR_Ready.map", header = FALSE, sep ="\t")

#Save .map file and combined .ped file with same base
write.table(map, "Ready_PLINK.map", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(combined_ped, "Ready_PLINK.ped", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# Running PLINK duplicate analysis ----------------------------------------

#clear workspace
rm(list=ls())

#set working directory [must contain plink.exe and files for analysis]
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/20K_480K PLINK Duplicate Identification/Inputs/Howard_2022")

#Run PLINK
system("plink --file Ready_PLINK --missing-genotype 0 --genome full ")

#Read genome file
genome <- read.table("plink.genome", header = TRUE, sep = "", stringsAsFactors = FALSE)
write.table(genome, "C:/Users/curly/Desktop/Apple Genotyping/Results/20K_480K PLINK Duplicate Identification/Howard et al. 2022 20K Results/PLINK_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

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
colnames(dd) <- c("Group", "SampleCount", "ID1","ID2","ID3","ID4","ID5","ID6","ID7","ID8","ID9","ID10","ID11","ID12","ID13")

#Save .csv of duplicate groupings
write.csv(dd, "C:/Users/curly/Desktop/Apple Genotyping/Results/20K_480K PLINK Duplicate Identification/Howard et al. 2022 20K Results/Grouped_Duplicates.csv", row.names = FALSE)


