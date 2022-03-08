# Analysis microbial communities from total RNA data processed using phyloFlash.


# Load libraries
library(dplyr)
library(microbiome)
library(microViz)
library(viridis)
library(vegan)
library(here)


# Create scale function to normalize library sizes
scale_reads <- function(physeq, n) {
  physeq.scale <-
    transform_sample_counts(physeq, function(x) {
      (n * x/sum(x))
    })
  otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}


# Make a list of paths of all .csv files from the phyloFlash output (*NTUabundance.csv)
all_paths <-list.files(path = "path/to/phyloFlash/TotalRNA/NTUabundance/files", pattern = "*.csv",full.names = TRUE)

#import all files
all_content <-lapply(all_paths, read.table,sep = ",", encoding = "UTF-8")

#combine all files
NTUtable <- Reduce(function(x, y) merge(x, y, by = 'V1', all = TRUE), all_content)

#get list of filenames and subtract path and extentions to get sample names
all_filenames <- all_paths %>% basename() %>% as.list()
all_filenames <-sub('\\.phyloFlash.NTUabundance.csv$', '', all_filenames)

# add entry for first column and replace column names with sample names
all_filenames <- c(L1 = "",all_filenames)
colnames(NTUtable)<-all_filenames

#create tax-table
path_split <- strsplit(NTUtable[,1], ";")
silva <- read.table("tax_slv_ssu_138.1.txt", h = F, sep = "\t", stringsAsFactors = F) # import taxa map from version used to annotate with phyloFLash

silva_map <- data.frame(  # prepare taxa map in right format for parsing function
  path = gsub(";$", "", silva$V1),
  node = sapply(strsplit(silva$V1, ";"), function(x) x[length(x)]),
  rank = silva$V3,
  stringsAsFactors = T
)


# Create new function for SILVA taxonomy
SILVAtaxopath <- function(tax, SILVA){    # parsing function provided by Christiane Hassenrück @chassenr on github
  output <- matrix(NA, nrow = length(tax), ncol = length(levels(SILVA$rank)))
  colnames(output) <- levels(SILVA$rank)
  for (i in 1:length(tax)) {
    for (j in 1:length(levels(SILVA$rank))) {
      if (paste(tax[[i]][1:j], collapse = ";") %in% SILVA$path) {
        output[i, as.character(SILVA[SILVA$path == paste(tax[[i]][1:j], collapse = ";"), "rank"])] <- as.character(SILVA[SILVA$path == paste(tax[[i]][1:j], collapse = ";"), "node"])
      }
    }
  }
  return(output)
}

# Use function to organize taxonomic calls
TAXmat <- SILVAtaxopath(path_split,silva_map)
prefix<- "NTU"  #create rownames corresponding to NTUs
suffix<- seq(1:nrow(NTUtable))
NTU.names<- paste(prefix,suffix, sep = "")
row.names(TAXmat)<-NTU.names
TAXmat <- TAXmat[,c("domain","major_clade","kingdom","phylum","class","order","family","genus")]


# create phyloseq tax table
TAX<-tax_table(TAXmat) 

#create OTU matrix for phyloseq import
df<-NTUtable[,-1]
sample_prefix<-"IS19"
colnames(df)<-paste(sample_prefix,colnames(df),sep = "_")
colnames(df)[16:18]
colnames(df)[16:18]<-c("Field-BL","Lab-BL","NEG")
rownames(df)<-NTU.names
df[is.na(df)] <- 0
NTUmat<-as.matrix(df)
NTU<-otu_table(NTUmat, taxa_are_rows = TRUE)


# create sample data for phyloseq
SAMPLE_table<-read.csv(file = "path/to/SAMPLE_data.csv", header = TRUE, row.names = 1,sep = ",")
SAMPLE_DATA<-sample_data(SAMPLE_table)


# Check OTU table names vs mapping file (sample_data) names (if they don't match, samples won't be included in final phyloseq object).
all(colnames(NTU) %in% rownames(SAMPLE_DATA))


# Create phyloseq object using the created tables
TotalRNA<-phyloseq(NTU,TAX,SAMPLE_DATA)


# fix taxa table to replace unknowns etc. (using microViz package)
TotalRNA_fixed<-tax_fix(TotalRNA, min_length = 4)

TotalRNA_fixed<- TotalRNA_fixed %>%
  tax_fix(
    min_length = 4,
    unknowns = c("uncultured", "uncultured class", "uncultured order", "uncultured family", "uncultured phylum", "Unknown Family family","Unknown", "endosymbionts", "Incertae Sedis", "Incertae Sedis class", "Incertae Sedis order","Incertae Sedis family"),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )
TotalRNA_fixed<-subset_taxa(TotalRNA_fixed,genus != "Unknown Family family")

# Save phyloseq objects for posterity
saveRDS(TotalRNA, here("TotalRNA.rds"))
saveRDS(TotalRNA_fixed, here("TotalRNA_fixed.rds"))