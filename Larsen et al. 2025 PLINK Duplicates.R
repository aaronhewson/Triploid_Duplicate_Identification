# PLINK duplicate analysis for Jim Dunckley Orchard/Plant & Food Research samples genotyped in 2025 for Aaron Hewson's MSc project, alongside samples from Larsen et al. 2025 (https://doi.org/10.1007/s10722-024-02104-1).
# Prior to using this script, Larsen et al. 2025 data downloaded as a .vcf file was formatted in Excel to match PLINK style. Triploid genotype calls were replaced with biallelic genotype calls.
# Genotype files for JD/PFR are exported from Axiom Analysis Suite with all SNPs present.
# High-quality SNPS overlapping JD/PFR and Larsen et al. 2025 were selected using a list provided by Nick Howard.


# Load Packages -----------------------------------------------------------

#Load packages
library(tibble)
library(igraph)

library(dplyr)
library(tidyr)

# Set Working Directory ---------------------------------------------------
#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/Triploid Duplicate Identification/Inputs/Larsen_2025")


# Initial JD/PFR .ped file curation ---------------------------------------
#Load .ped file and remove .CEL from sample filenames
ped <- read.csv("JD_PFR_Raw.ped", header = FALSE,sep = "\t")
ped[[1]] <- sub("\\.CEL$", "", ped[[1]])

#Keep triploids only in .ped file
ids_to_keep <- read.delim("TriploidSampleNames.txt", header = FALSE, stringsAsFactors = FALSE)
ids_to_keep <- ids_to_keep$V1
ped <- ped[(ped$V1 %in% ids_to_keep), ]

#Add empty columns for PLINK formatting
ped <- add_column(ped, Fa = 0, Mo = 0, Se = 0, Ph = 0, .before = "V2" )
ped <- add_column(ped, Fid = 0, .before = "V1" )

#Remove header row and save .ped file
ped = ped[-1, ]
write.table(ped, "JD_PFR_All.ped", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Initial JD/PFR .map file curation ---------------------------------------

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

# Format Files and Extract SNPs with PLINK --------------------------------

#Convert Larsen et al. 2025 data from VCF to PLINK .ped and .map
system("plink --vcf Danish_20K.vcf --const-fid 0 --allow-extra-chr --recode tab --out Danish_20K")

#Use PLINK to extract overlapping SNPs from JD and PFR samples
system("plink --file JD_PFR_All --extract JD_PFR_Danish_ExtractList.txt --make-bed --out Danish_JD_PFR")
system("plink --bfile Danish_JD_PFR --recode --tab --out Danish_JD_PFR")

#Use PLINK to extract overlapping SNPs from 20K Larsen et al. 2025 samples, recode SNP positions, and output as .ped and .map
system("plink --file Danish_20K --extract Danish_ExtractList.txt --allow-extra-chr --update-chr Danish_chr.txt --update-cm Danish_cm.txt --update-map Danish_map.txt --make-bed --out Danish_20K")
system("plink --bfile Danish_20K  --update-name Danish_name.txt --recode tab --out Danish_Ready")


# Combine JD/PFR and Larsen et al. 2025 Samples ---------------------------

#clear workspace
rm(list=ls())

#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/Triploid Duplicate Identification/Inputs/Larsen_2025")

#Load JD_PFR data
JD_ped <- read.csv("Danish_JD_PFR.ped", header = FALSE,sep = "\t")

#Load Larsen 2025 Danish 20K data, add prefix to names, set irrelevant tags to 0
LD_ped <- read.csv("Danish_Ready.ped", header = FALSE,sep = "\t")
LD_ped$V2 <- paste0("LD_", LD_ped$V2) 
LD_ped$V1 <- 0
LD_ped$V3 <- 0
LD_ped$V4 <- 0

#Bind Larsen 2025 Danish 20K data onto PLINK .ped, set other variables to 0
names(LD_ped) <- names(JD_ped)
combined_ped <- bind_rows(JD_ped, LD_ped)

#Load .map file
map <- read.csv("Danish_JD_PFR.map", header = FALSE, sep ="\t")

#Save .map file and combined .ped file with same base
write.table(map, "Danish_PLINK.map", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(combined_ped, "Danish_PLINK.ped", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



# Run PLINK Duplicate Analysis --------------------------------------------

#clear workspace
rm(list=ls())

#set working directory [must contain plink.exe and files for analysis]
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/Triploid Duplicate Identification/Inputs/Larsen_2025")

#Run PLINK
system("plink --file Danish_PLINK --missing-genotype 0 --genome full ")

# Load and Save PLINK .genome File ----------------------------------------

#Read genome file
genome <- read.table("plink.genome", header = TRUE, sep = "", stringsAsFactors = FALSE)
write.table(genome, "C:/Users/curly/Desktop/Apple Genotyping/Results/Triploid Duplicates/Larsen_2025_Duplicates/PLINK_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

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
write.csv(dd, "C:/Users/curly/Desktop/Apple Genotyping/Results/Triploid Duplicates/Larsen_2025_Duplicates/Grouped_Duplicates.csv", row.names = FALSE)

