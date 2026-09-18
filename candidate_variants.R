setwd("~/Documents/manuscripts/XQTL/G3_submission/github/variant_identification")
library(tidyverse)

#full list of 3609 alternate alleles
genos <- read_delim("split.filtered.B6neB7.vcf",delim="\t",col_names=TRUE,skip=22)
names(genos)[4] <- "Ref"

## determine the length of the alternate allele, this is helpful for categorizing the type of variant
genos$reflength <- nchar(genos$Ref)
genos$altlength <- nchar(genos$ALT)

#determine the class of variant based on length
genos$variant_type <- "unclear"
genos$variant_type[genos$reflength == 1 & genos$altlength== 1] <- "SNP"
genos$variant_type[genos$variant_type != "SNP" & genos$reflength < 10 & abs(genos$reflength-genos$altlength)<= 100] <- "small INDEL"
genos$variant_type[genos$variant_type != "SNP" & genos$altlength < 10 & abs(genos$reflength-genos$altlength) <= 100] <- "small INDEL"
genos$variant_type[genos$reflength < 10 & abs(genos$reflength-genos$altlength)> 100] <- "large INDEL"
genos$variant_type[genos$altlength < 10 & abs(genos$reflength-genos$altlength)> 100] <- "large INDEL"
genos$variant_type[genos$reflength > 10*genos$altlength & genos$reflength> 100] <- "large INDEL"
genos$variant_type[genos$altlength > 10*genos$reflength & genos$altlength> 100] <- "large INDEL"
genos$variant_type[genos$reflength >100 & genos$altlength> 100] <- "complex"

write_csv(genos,"genos.annotated.vcf")

#manually inspect the "complex" and "unclear" variants, often these reflect positions where there is an indel and then nested variation in the longer allele. Also delete the alleles where neither B2 nor B6 is ALT. THe curated alleles are then run through VEP to relate them to gene function.
#import my cleaned clean version of the data
genos <- read_csv("genos.annotated.curated.csv")
##there are 2132 lines, i.e. variants in the file, which correspond to 2046 variable positions in dm6
pull(genos,ID) %>% unique() %>% length()
genos$variant_type <- factor(genos$variant_type, levels = c("SNP","small INDEL","large INDEL","small ambiguous","large ambiguous"))


# read in VEP output with variant effects
VEP <- read_csv("VEP.csv",col_names=TRUE) 
names(VEP)[1] <- "ID"

#merge VEP with molecular descriptions and genotypes. Note that the merged table is longer than any individual table because variant positions can have multiple effects (in VEP) and also multiple alleles (in genos)
merged <-left_join(VEP,genos)

#make a new column "functional" describing whether a variant may impact ANY genes (i.e. is it within a gene body or 5Kb of any encoded transcript). If it is not is labelled as intergenic by VEP
merged$functional <- "genic"
merged$functional[merged$Consequence == "intergenic_variant"] <- "intergenic"

#now make a new table where very variable position is described a single time as either "genic" or "intergenic", note this table has 2050 variable positions because a handfull of variants have alleles of different types

variant.tbl <- select(merged,ID,variant_type,functional) %>% distinct()

pull(variant.tbl,ID) %>% table() %>% sort(decreasing=TRUE)

#make the graph
ggplot(variant.tbl,aes(variant_type,fill=functional)) + geom_bar() +
  scale_fill_manual(values=c("#18dbbe","#2d3d3b")) + theme_classic() +
  xlab("") + theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ylab("candidate variant count\n")
ggsave("variant_effects.pdf",height=3,width=3)

#make a table describing whether different types of effects are observed in each gene 
gene.impacts <- filter(merged, Consequence != "intergenic_variant") %>% select(SYMBOL,Consequence)  %>% distinct() %>% mutate(dummy="yes") %>% pivot_wider(names_from = Consequence,values_from = dummy, values_fill = "no") 

#intron_variant column should read yes for the non-coding RNAs. We can drop the column intron_variant&non_coding_transcript_variant after we make this change
gene.impacts$intron_variant[gene.impacts$`intron_variant&non_coding_transcript_variant` =="yes"] <- "yes"
gene.impacts<- select(gene.impacts, -`intron_variant&non_coding_transcript_variant`)

#fix gene name for sph242
gene.impacts$SYMBOL[gene.impacts$SYMBOL == "CG40160"] <- "Sph242"

#add row for CR41623
names(gene.impacts)

gene.impacts <- rbind(gene.impacts, tibble(SYMBOL = "CR46123",intron_variant = "no", upstream_gene_variant = "no", 
                 downstream_gene_variant = "no", synonymous_variant = "no", non_coding_transcript_exon_variant = "no",
                 `5_prime_UTR_variant` = "no", `3_prime_UTR_variant` = "no", missense_variant = "no"))

#add information about differential expression
gene.impacts$DEG <- "not DEG"
gene.impacts$DEG[gene.impacts$SYMBOL == "Pzl"] <- "up in B7"
gene.impacts$DEG[gene.impacts$SYMBOL == "l(3)80Fg"] <- "up in B7"
gene.impacts$DEG[gene.impacts$SYMBOL == "spok"] <- "up in B7"
gene.impacts$DEG[gene.impacts$SYMBOL == "lovit"] <- "up in B7"
gene.impacts$DEG[gene.impacts$SYMBOL == "CR42650"] <- "up in B7"
gene.impacts$DEG[gene.impacts$SYMBOL == "Pzl"] <- "up in B7"
gene.impacts$DEG[gene.impacts$SYMBOL == "Myo81F"] <- "up in B6"
gene.impacts$DEG[gene.impacts$SYMBOL == "Dbp80"] <- "up in B6"
gene.impacts$DEG[gene.impacts$SYMBOL == "Sph242"] <- "up in B6"
gene.impacts$DEG[gene.impacts$SYMBOL == "CR46123"] <- "up in B6"

#now we have to pivot long again to make the graph
gene.impacts <- pivot_longer(gene.impacts, !SYMBOL, names_to ="impact",values_to = "fill")

#cleanup and reorder names
gene.impacts$impact <- gsub("3_prime","3'",gene.impacts$impact)
gene.impacts$impact <- gsub("5_prime","5'",gene.impacts$impact)
gene.impacts$impact <- gsub("_variant","",gene.impacts$impact)
gene.impacts$impact <- gsub("_gene","",gene.impacts$impact)
gene.impacts$impact <- gsub("_UTR","-UTR",gene.impacts$impact)
gene.impacts$impact <- gsub("non_coding_transcript_exon","ncRNA-transcript",gene.impacts$impact)
levels(as.factor(gene.impacts$impact))
gene.impacts$impact <- factor(gene.impacts$impact, levels=c("missense","synonymous","5'-UTR","3'-UTR","ncRNA-transcript","intron","upstream","downstream","DEG"))
gene.impacts$fill <- factor(gene.impacts$fill, levels=c("no","yes","not DEG","up in B6","up in B7"))

#make a variable to describe facets
gene.impacts$facet <- "variants"
gene.impacts$facet[gene.impacts$impact=="DEG"] <- "DEG"
gene.impacts$facet <- factor(gene.impacts$facet, levels=c("variants","DEG"))

#now we can make a tile graph relating the presence of different types of variants to genes within the QTL
ggplot(gene.impacts, aes(x=impact,y=SYMBOL,fill=fill))+  
  geom_tile() + theme_classic() +
  facet_wrap(~facet,scales="free_x",space="free_x")+
  theme(axis.text.x=element_text(angle=90,hjust=1),axis.text.y=element_text(face="italic")) +
  scale_fill_manual(values=c("#332f2c","#18dbbe","#2d3d3b","#338a7d","#Db2f18")) +
  xlab("") + ylab("")
ggsave("impacted_genes.pdf",height=4.5,width=4)


