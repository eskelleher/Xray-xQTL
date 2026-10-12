#load required packages for GO analysis 

library(clusterProfiler) # GO enrichment analysis
library(org.Dm.eg.db) # Drosophila melanogaster gene annotations
library(AnnotationDbi) # Gene ID mapping
library(ggplot2) # Plotting
library(enrichplot) # GO enrichment visualization



# Load a CSV file containing the FlyBase gene IDs for the  
# The input file must contain a column named "gene_id".

# Choose a CSV file according to the module of interest.
midnightblue <- read.csv("module_midnightblue.csv", stringsAsFactors = FALSE)

# Convert FlyBase gene IDs to Entrez Gene IDs, which are used as input for the GO enrichment analysis below.
entrez_midnightblue<- mapIds(org.Dm.eg.db,
                             keys = midnightblue$gene_id,
                             column = "ENTREZID",
                             keytype = "FLYBASE", # keytype = "FLYBASE" specifies the identifier type in the input gene list.
                             multiVals = "first") # multiVals = "first" retains the first mapped Entrez ID if multiple mappings are available for a FlyBase ID.

# Perform GO enrichment

ego_midnightblue<- enrichGO(gene = entrez_midnightblue,
                            OrgDb = org.Dm.eg.db,
                            keyType = "ENTREZID",
                            ont = "MF", # ont specifies the GO category:: MF= Molecular function, BP= Biological process, CC= Cellular component
                            pAdjustMethod = "BH", # pAdjustMethod = "BH" applies the Benjamini-Hochberg procedure to adjust p-values for multiple testing.
                            
                            qvalueCutoff = 0.05)

# Plot top 10 enriched Molecular function terms
# Note: The number of terms displayed can be changed using showCategory.
dotplot(ego_midnightblue, showCategory = 10, title = "Module midnightblue")

# We can also plot these as cnet plots, showing genes alongside the enriched GO terms
#Convert the enrichment results to use readable gene symbols instead of Entrez IDs where annotations are available. This makes gene labels easier to interpret in downstream plots.


ego_midnightblue <- setReadable(
  ego_midnightblue,
  OrgDb = org.Dm.eg.db,
  keyType = "ENTREZID"
)

# Generate a category network (cnet) plot for the top 3 enriched GO terms. The plot links enriched terms to the genes associated with each term, helping visualize overlap among gene sets.

cnetplot(
  ego_midnightblue,
  showCategory = 3, 
  node_label = "all" # node_label = "all" labels both GO term and gene nodes.
) 
