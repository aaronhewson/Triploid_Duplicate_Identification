# PLINK duplicate analysis for Jim Dunckley Orchard/Plant & Food Research samples genotyped in 2025 for Aaron Hewson's MSc project, alongside samples from Julius Kuhn Institute genotyped on 50K SNP array in 2025.
# Prior to using this script, JKI data provided as .xlsx was formatted in Excel to match PLINK style, and only SNPs selected as "PolyHighRes" in JD/PFR genotyping were kept.


# Load Packages -----------------------------------------------------------

#Load packages
library(tibble)
library(igraph)

library(dplyr)
library(tidyr)

# Set Working Directory ---------------------------------------------------

#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/Triploid Duplicate Identification/Inputs/JKI_2025")

# Format JD/PFR .ped File --------------------------------------------------------

#Load .ped file and remove .CEL from sample filenames
ped <- read.csv("JD_PFR_Genotyping.ped", header = FALSE,sep = "\t")
ped[[1]] <- sub("\\.CEL$", "", ped[[1]])

#Remove triploids from .ped file
ids_to_keep <- read.delim("TriploidSampleNames.txt", header = FALSE, stringsAsFactors = FALSE)
ids_to_keep <- ids_to_keep$V1
ped <- ped[(ped$V1 %in% ids_to_keep), ]

#Add empty columns for PLINK formatting
ped <- add_column(ped, Fa = 0, Mo = 0, Se = 0, Ph = 0, .before = "V2" )
ped <- add_column(ped, Fid = 0, .before = "V1" )

#Remove header row and save .ped file
ped = ped[-1, ]
write.table(ped, "JD_PFR_PLINK.ped", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Format JD/PFR .map File --------------------------------------------------------

#Load .map file and BLAST results
map <- read.csv("JD_PFR_Genotyping.map", header = FALSE, sep ="\t")
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
write.table(map, "JD_PFR_PLINK.map", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Load and Format JKI Samples ---------------------------------------------

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

# Combine JKI and JD/PFR Samples ------------------------------------------

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


# Run PLINK Duplicate Analysis --------------------------------------------

#clear workspace
rm(list=ls())

#set working directory [must contain plink.exe and files for analysis]
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/Triploid Duplicate Identification/Inputs/JKI_2025")

#Run PLINK
system("plink --file JKI_PLINK --missing-genotype 0 --genome full ")


# Load and Save PLINK .genome File --------------------------------------

#Read genome file
genome <- read.table("plink.genome", header = TRUE, sep = "", stringsAsFactors = FALSE)
write.table(genome, "C:/Users/curly/Desktop/Apple Genotyping/Results/Triploid Duplicates/JKI_Duplicates/PLINK_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Grouping Duplicate IDs --------------------------------------------------

#Filter for PI_HAT >0.90 (duplicate threshold)
genome <- genome[!(genome$PI_HAT < 0.90), ]
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
