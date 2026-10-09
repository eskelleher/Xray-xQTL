#load packages
library(tidyverse)
library(DESeq2)
library (ggplot2)
library(readr)
library(sva)
library(WGCNA)
library(flashClust)
library(curl)
library (dplyr)

#Load Kallisto transcript-level abundance files for all samples at 0, 1, and 4 hours.
# 0 hr
line96093_R1_0hr <- read.delim("96093-R1-0hr.tsv", header = TRUE, sep = "\t")
line96093_R2_0hr <- read.delim("96093-R2-0hr.tsv", header = TRUE, sep = "\t")
line96093_R3_0hr <- read.delim("96093-R3-0hr.tsv", header = TRUE, sep = "\t")
line96197_R1_0hr <- read.delim("96197-R1-0hr.tsv", header = TRUE, sep = "\t")
line96197_R2_0hr <- read.delim("96197-R2-0hr.tsv", header = TRUE, sep = "\t")
line96197_R3_0hr <- read.delim("96197-R3-0hr.tsv", header = TRUE, sep = "\t")
line96943_R1_0hr <- read.delim("96943-R1-0hr.tsv", header = TRUE, sep = "\t")
line96943_R2_0hr <- read.delim("96943-R2-0hr.tsv", header = TRUE, sep = "\t")
line96943_R3_0hr <- read.delim("96943-R3-0hr.tsv", header = TRUE, sep = "\t")
line96977_R1_0hr <- read.delim("96977-R1-0hr.tsv", header = TRUE, sep = "\t")
line96977_R2_0hr <- read.delim("96977-R2-0hr.tsv", header = TRUE, sep = "\t")
line96977_R3_0hr <- read.delim("96977-R3-0hr.tsv", header = TRUE, sep = "\t")

# 1 hr
line96093_R1_1hr <- read.delim("96093-R1-1hr.tsv", header = TRUE, sep = "\t")
line96093_R2_1hr <- read.delim("96093-R2-1hr.tsv", header = TRUE, sep = "\t")
line96093_R3_1hr <- read.delim("96093-R3-1hr.tsv", header = TRUE, sep = "\t")
line96197_R1_1hr <- read.delim("96197-R1-1hr.tsv", header = TRUE, sep = "\t")
line96197_R2_1hr <- read.delim("96197-R2-1hr.tsv", header = TRUE, sep = "\t")
line96197_R3_1hr <- read.delim("96197-R3-1hr.tsv", header = TRUE, sep = "\t")
line96943_R1_1hr <- read.delim("96943-R1-1hr.tsv", header = TRUE, sep = "\t")
line96943_R2_1hr <- read.delim("96943-R2-1hr.tsv", header = TRUE, sep = "\t")
line96943_R3_1hr <- read.delim("96943-R3-1hr.tsv", header = TRUE, sep = "\t")
line96977_R1_1hr <- read.delim("96977-R1-1hr.tsv", header = TRUE, sep = "\t")
line96977_R2_1hr <- read.delim("96977-R2-1hr.tsv", header = TRUE, sep = "\t")
line96977_R3_1hr <- read.delim("96977-R3-1hr.tsv", header = TRUE, sep = "\t")

# 4 hr
line96093_R1_4hr <- read.delim("96093-R1-4hr.tsv", header = TRUE, sep = "\t")
line96093_R2_4hr <- read.delim("96093-R2-4hr.tsv", header = TRUE, sep = "\t")
line96093_R3_4hr <- read.delim("96093-R3-4hr.tsv", header = TRUE, sep = "\t")
line96197_R1_4hr <- read.delim("96197-R1-4hr.tsv", header = TRUE, sep = "\t")
line96197_R2_4hr <- read.delim("96197-R2-4hr.tsv", header = TRUE, sep = "\t")
line96197_R3_4hr <- read.delim("96197-R3-4hr.tsv", header = TRUE, sep = "\t")
line96943_R1_4hr <- read.delim("96943-R1-4hr.tsv", header = TRUE, sep = "\t")
line96943_R2_4hr <- read.delim("96943-R2-4hr.tsv", header = TRUE, sep = "\t")
line96943_R3_4hr <- read.delim("96943-R3-4hr.tsv", header = TRUE, sep = "\t")
line96977_R1_4hr <- read.delim("96977-R1-4hr.tsv", header = TRUE, sep = "\t")
line96977_R2_4hr <- read.delim("96977-R2-4hr.tsv", header = TRUE, sep = "\t")
line96977_R3_4hr <- read.delim("96977-R3-4hr.tsv", header = TRUE, sep = "\t")


#Load the FlyBase transcript-to-gene annotation used to convert transcript IDs to gene IDs.
index <- read.csv("FlyBase-tx2gn-FULL.csv", header = TRUE, stringsAsFactors = FALSE)

names(index) <- c("gene_id", "target_id")
head(index)





# Define a function to summarize transcript-level Kallisto estimated counts into gene-level counts.
tr2gn <- function(kallisto_df, sample_name, tx2gene) {
    merged <- dplyr::left_join(kallisto_df, tx2gene, by = "target_id")
  
  # Now sum transcript counts per gene
  gene_counts <- merged %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(est_counts = sum(est_counts, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  # Rename column to sample name
  colnames(gene_counts)[2] <- sample_name
  
  return(gene_counts)
}


# Use tr2gn() to summarize transcript-level estimated counts into gene-level counts for each sample.

line96093_R1_0hr.counts <- tr2gn(line96093_R1_0hr, "line96093_R1_0hr", index)
line96093_R3_0hr.counts <- tr2gn(line96093_R3_0hr, "line96093_R3_0hr", index)
line96197_R2_0hr.counts <- tr2gn(line96197_R2_0hr, "line96197_R2_0hr", index)
line96943_R1_0hr.counts <- tr2gn(line96943_R1_0hr, "line96943_R1_0hr", index)
line96943_R3_0hr.counts <- tr2gn(line96943_R3_0hr, "line96943_R3_0hr", index)
line96977_R2_0hr.counts <- tr2gn(line96977_R2_0hr, "line96977_R2_0hr", index)


line96093_R1_1hr.counts <- tr2gn(line96093_R1_1hr, "line96093_R1_1hr", index)
line96093_R3_1hr.counts <- tr2gn(line96093_R3_1hr, "line96093_R3_1hr", index)
line96197_R2_1hr.counts <- tr2gn(line96197_R2_1hr, "line96197_R2_1hr", index)
line96943_R1_1hr.counts <- tr2gn(line96943_R1_1hr, "line96943_R1_1hr", index)
line96943_R3_1hr.counts <- tr2gn(line96943_R3_1hr, "line96943_R3_1hr", index)
line96977_R2_1hr.counts <- tr2gn(line96977_R2_1hr, "line96977_R2_1hr", index)

line96093_R1_4hr.counts <- tr2gn(line96093_R1_4hr, "line96093_R1_4hr", index)
line96093_R3_4hr.counts <- tr2gn(line96093_R3_4hr, "line96093_R3_4hr", index)
line96197_R2_4hr.counts <- tr2gn(line96197_R2_4hr, "line96197_R2_4hr", index)
line96943_R1_4hr.counts <- tr2gn(line96943_R1_4hr, "line96943_R1_4hr", index)
line96943_R3_4hr.counts <- tr2gn(line96943_R3_4hr, "line96943_R3_4hr", index)
line96977_R2_4hr.counts <- tr2gn(line96977_R2_4hr, "line96977_R2_4hr", index)

line96093_R2_0hr.counts <- tr2gn(line96093_R2_0hr, "line96093_R2_0hr", index)
line96197_R1_0hr.counts <- tr2gn(line96197_R1_0hr, "line96197_R1_0hr", index)
line96197_R3_0hr.counts <- tr2gn(line96197_R3_0hr, "line96197_R3_0hr", index)
line96943_R2_0hr.counts <- tr2gn(line96943_R2_0hr, "line96943_R2_0hr", index)
line96977_R1_0hr.counts <- tr2gn(line96977_R1_0hr, "line96977_R1_0hr", index)
line96977_R3_0hr.counts <- tr2gn(line96977_R3_0hr, "line96977_R3_0hr", index)

line96093_R2_1hr.counts <- tr2gn(line96093_R2_1hr, "line96093_R2_1hr", index)
line96197_R1_1hr.counts <- tr2gn(line96197_R1_1hr, "line96197_R1_1hr", index)
line96197_R3_1hr.counts <- tr2gn(line96197_R3_1hr, "line96197_R3_1hr", index)
line96943_R2_1hr.counts <- tr2gn(line96943_R2_1hr, "line96943_R2_1hr", index)
line96977_R1_1hr.counts <- tr2gn(line96977_R1_1hr, "line96977_R1_1hr", index)
line96977_R3_1hr.counts <- tr2gn(line96977_R3_1hr, "line96977_R3_1hr", index)


line96093_R2_4hr.counts <- tr2gn(line96093_R2_4hr, "line96093_R2_4hr", index)
line96197_R1_4hr.counts <- tr2gn(line96197_R1_4hr, "line96197_R1_4hr", index)
line96197_R3_4hr.counts <- tr2gn(line96197_R3_4hr, "line96197_R3_4hr", index)
line96943_R2_4hr.counts <- tr2gn(line96943_R2_4hr, "line96943_R2_4hr", index)
line96977_R1_4hr.counts <- tr2gn(line96977_R1_4hr, "line96977_R1_4hr", index)
line96977_R3_4hr.counts <- tr2gn(line96977_R3_4hr, "line96977_R3_4hr", index)

#Combine all gene-level count tables into a single count matrix containing all samples.
count_data <- full_join(line96093_R1_0hr.counts,line96093_R3_0hr.counts) 
count_data <- full_join(count_data,line96197_R2_0hr.counts) 
count_data <- full_join(count_data,line96943_R1_0hr.counts) 
count_data <- full_join(count_data,line96943_R3_0hr.counts)
count_data <- full_join(count_data,line96977_R2_0hr.counts)

count_data <- full_join(count_data,line96093_R1_1hr.counts)
count_data <- full_join(count_data,line96093_R3_1hr.counts)
count_data <- full_join(count_data,line96197_R2_1hr.counts)
count_data <- full_join(count_data,line96943_R1_1hr.counts)
count_data <- full_join(count_data,line96943_R3_1hr.counts)
count_data <- full_join(count_data,line96977_R2_1hr.counts)

count_data <- full_join(count_data,line96093_R1_4hr.counts)
count_data <- full_join(count_data,line96093_R3_4hr.counts)
count_data <- full_join(count_data,line96197_R2_4hr.counts)
count_data <- full_join(count_data,line96943_R1_4hr.counts)
count_data <- full_join(count_data,line96943_R3_4hr.counts)
count_data <- full_join(count_data,line96977_R2_4hr.counts)


count_data <- full_join(count_data,line96093_R2_0hr.counts)
count_data <- full_join(count_data,line96197_R1_0hr.counts)
count_data <- full_join(count_data,line96197_R3_0hr.counts)
count_data <- full_join(count_data,line96943_R2_0hr.counts)
count_data <- full_join(count_data,line96977_R1_0hr.counts)
count_data <- full_join(count_data,line96977_R3_0hr.counts)



count_data <- full_join(count_data,line96093_R2_1hr.counts)
count_data <- full_join(count_data,line96197_R1_1hr.counts)
count_data <- full_join(count_data,line96197_R3_1hr.counts)
count_data <- full_join(count_data,line96943_R2_1hr.counts)
count_data <- full_join(count_data,line96977_R1_1hr.counts)
count_data <- full_join(count_data,line96977_R3_1hr.counts)

count_data <- full_join(count_data,line96093_R2_4hr.counts)
count_data <- full_join(count_data,line96197_R1_4hr.counts)
count_data <- full_join(count_data,line96197_R3_4hr.counts)
count_data <- full_join(count_data,line96943_R2_4hr.counts)
count_data <- full_join(count_data,line96977_R1_4hr.counts)
count_data <- full_join(count_data,line96977_R3_4hr.counts)


# Replace missing gene/sample combinations with zero counts.
count_data[is.na(count_data)] <- 0

# Convert the gene-level count data frame to a count matrix.
countData <- as.matrix(count_data[, -which(names(count_data) == "gene_id")])

# Set row names to gene IDs.
rownames(countData) <- count_data$gene_id


#make colData object with information about every sample
dataframe <- read_csv("gene_expression_dataframe.csv")
names(dataframe)
colData <- column_to_rownames(dataframe, var = "...1")
colnames(countData) <- gsub("4h$", "4hr", colnames(countData))
rownames(colData) <- gsub("4h$", "4hr", rownames(colData))
colData$Replicate <- factor(colData$Replicate)
colData$Time <- factor(colData$Time)

colData <- colData[colnames(countData), ]
all(colnames(countData) == rownames(colData))
countData <- round(countData)
storage.mode(countData) <- "integer"


# Confirm dimensions of the gene-level count matrix
dim(countData)

# Create a DESeq2 object for count normalization and variance-stabilizing transformation prior to WGCNA.
d <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~1)

# Remove genes with fewer than 10 counts in at least three samples to reduce noise from very lowly expressed genes.
smallestGroupSize <- 3

keep <- rowSums(counts(d) >= 10) >= smallestGroupSize

d <- d[keep, ]

# Number of genes retained after low-count filtering.
sum(keep)

# Normalize for differences in sequencing depth and estimate gene-level dispersion before variance stabilization.
d <- estimateSizeFactors(d)

d <- estimateDispersions(d)

# Apply a variance-stabilizing transformation to reduce the dependence of variance on mean expression for correlation-based network analysis.
counttable <- data.frame(getVarianceStabilizedData(d))

dim(counttable)

edata <- as.matrix(counttable)



#Calculate the variance of each gene across all samples to identify the most variable genes for WGCNA.

gene_variance <- apply(
  edata,
  1,
  var
)

summary(gene_variance)


#Description: Rank genes from highest to lowest expression variance to select the most variable genes for network construction.
variance_rank <- order(
  gene_variance,
  decreasing = TRUE
)

high_var_genes_10000 <- names(gene_variance)[
  variance_rank[1:10000]
]

length(high_var_genes_10000)


#Restrict the expression matrix to the 10,000 most variable genes selected for WGCNA.

vsd_highvar_10000 <- edata[
  high_var_genes_10000,
]

dim(vsd_highvar_10000)


#Estimate the number of hidden sources of variation under different biological models to assess potential unwanted variation in the dataset.
edata_10000 <- vsd_highvar_10000

mod <- model.matrix(~ Replicate:Background + Background + Allele + Time, data = colData)
num.sv(edata_10000, mod,method = "leek")
#this suggests that there are no significant sources of unaccounted variation in the model


mod <- model.matrix(~ Background + Allele + Time,data = colData)
num.sv(edata_10000,mod,method = "leek")
#if we remove replicate we can see the batch effects


mod <- model.matrix(~ Replicate:Background + Background + Time,data = colData)
num.sv(edata_10000,mod,method = "leek")
#if we remove allele 


mod <- model.matrix(~ Replicate:Background + Background + Allele,data = colData)
num.sv(edata_10000,mod,method = "leek")
#if we remove time



# Use ComBat to adjust for replicate-by-background batch effects before WGCNA.
modcombat_10000 <- model.matrix(~1,data = colData)

combat_edata_10000 <- ComBat(dat = vsd_highvar_10000, batch = interaction(colData$Replicate,colData$Background),mod = modcombat_10000,par.prior = TRUE, prior.plots = FALSE)


#Transpose the expression matrix so samples are rows and genes are columns, as required by WGCNA.
datExpr_10000 <- as.data.frame(t(combat_edata_10000)
)

datExpr2_10000 <- as.matrix(datExpr_10000)

storage.mode(datExpr2_10000) <- "numeric"


datExpr2_10000.df <- as_tibble(datExpr2_10000,rownames = NA)

#Remove genes with constant expression across samples because they cannot contribute to a correlation-based co-expression network.
datExpr3_10000.df <- datExpr2_10000.df %>%dplyr::select(dplyr::where(~ dplyr::n_distinct(.) > 1)) %>%as.data.frame()

rownames(datExpr3_10000.df) <- rownames(datExpr2_10000.df)

datExpr3_10000 <- as.matrix(datExpr3_10000.df)

#Evaluate candidate soft-thresholding powers to identify a power that produces an approximately scale-free network.

sft_10000 <- pickSoftThreshold(
  datExpr3_10000,
  dataIsExpr = TRUE,
  corFnc = bicor,
  networkType = "signed"
)

sft_10000_df <- data.frame(sft_10000$fitIndices) %>%dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(
  sft_10000_df,
  aes(x = Power, y = model_fit, label = Power)
) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.80, col = "red") +
  ylim(c(min(sft_10000_df$model_fit), 1.05)) +
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R²") +
  ggtitle("Scale independence - top 10,000 genes") +
  theme_classic()



#Construct a signed weighted gene co-expression network and identify groups of co-expressed genes as modules.

temp_cor <- cor
cor <- WGCNA::cor

maxBlockSize <- 5000

netwk_10000 <- blockwiseModules(
  datExpr3_10000,
  power = 14,  # Soft-thresholding power selected based on the scale-free topology analysis.
  networkType = "signed", #Genes are grouped based on the direction and strength of their expression correlations.
  corType = "bicor", #Use biweight midcorrelation, which is more robust to potential outliers than Pearson correlation.
  TOMType = "signed", #Construct a signed topological overlap matrix.
  randomSeed = 1234,
  deepSplit = 2, #Controls how finely the gene dendrogram is split into initial modules.
  pamRespectsDendro = FALSE,
  minModuleSize = 30, #Modules must contain at least 30 genes.
  maxBlockSize = maxBlockSize,
  blockSizePenaltyPower = 5,
  nPreclusteringCenters = as.integer(
    min(
      ncol(datExpr3_10000) / 20,
      100 * ncol(datExpr3_10000) / maxBlockSize
    )
  ),
  reassignThreshold = 0.2, # Allows genes to be reassigned when their module membership provides stronger evidence for another module. 
  mergeCutHeight = 0.25, #Merge modules with highly similar expression profiles.
  saveTOMs = TRUE,
  saveTOMFileBase = "ER_10000",
  numericLabels = FALSE,
  verbose = 3
)

cor <- temp_cor



# Convert WGCNA numeric module assignments to color labels and create a table linking each gene to its module.
mergedColors_10000 <- data.frame(
  gene_id = names(netwk_10000$colors),
  colors = labels2colors(netwk_10000$colors)
)

table(mergedColors_10000$colors)

mergedColors_10000