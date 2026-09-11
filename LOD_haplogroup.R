library(viridis)
library(tidyverse)
library(ggplotify)	
library(gridExtra)

#import haplotypes for 200 Kb, 500 Kb and 1Mb windows
hap200 <- read_delim("allhaps.200Kb.txt",col_names=TRUE) 
hap500 <- read_delim("allhaps.500Kb.txt",col_names=TRUE) 
hap1000 <- read_delim("allhaps.1Mb.txt",col_names=TRUE) 

#make AB8 values into B8 so they are sorted last.
hap200$founder[hap200$founder == "AB8"] <- "B8"
hap500$founder[hap500$founder == "AB8"] <- "B8"
hap1000$founder[hap1000$founder == "AB8"] <- "B8"

#make chr into a factor
hap200$chr <- as.factor(hap200$chr)
hap500$chr <- as.factor(hap500$chr)
hap1000$chr <- as.factor(hap1000$chr)
###extract windows from both sides of the 3rd chromosome centromere, exclude windows were AF is not estimated
###we only need values from one replicate:treatment. We will use identical AF estimates to detect pairs and groups of haplotypes that cannot be differentiated from each other within a window.



#the function below will take a list of 8 haplotype frequencies in the same position, sort them, keep the unique values first founder name (in ascending order 1-8)
#the each uniquely identifiable haplotypes is then assigned a "haplogroup" name based o the first founder that has this haplotype
haplogroup_name <- function(df){
  df <- data.frame(df)
  left_join(df, arrange(df, freq,founder) %>% distinct(freq,.keep_all=TRUE) %>% rename(c("haplogroup"="founder"))) 
}

cent.200 <- rbind(filter(hap200, chr == "chr3R", pos < 7307440, !is.na(freq), pool== "A.C") %>% select(-pool),
                 filter(hap200, chr == "chr3L", pos > 22586280, !is.na(freq), pool== "A.C") %>% select(-pool))
## groups data in cent.3R by chromosome and position and assigns a founder name to each unique haplogroup
cent.200.haplogroups <- cent.200 %>% 
  group_by(chr,pos) %>% 
  nest() %>% 
  mutate(allele = map(data,haplogroup_name)) %>% 
  unnest() 

###now isolate all the haplotype frequency data in the window of interest
cent.200 <- rbind(filter(hap200, chr == "chr3R", pos < 7307440, !is.na(freq)),
                  filter(hap200, chr == "chr3L", pos > 22586280, !is.na(freq)))

### sum the haplotype frequencies of haplotypes that cannot be distinguished to generate the frequency of the haplogroup.
cent.200 <- left_join(cent.200,select(cent.200.haplogroups,-freq)) %>% group_by(chr,pos,pool,haplogroup) %>%
  summarise(freq_sum = sum(freq, na.rm = TRUE))


#repeat for haplotypes based on 500 Kb windows
cent.500 <- rbind(filter(hap500, chr == "chr3R", pos < 7307440, !is.na(freq), pool== "A.C") %>% select(-pool),
                  filter(hap500, chr == "chr3L", pos > 22586280, !is.na(freq), pool== "A.C") %>% select(-pool))

cent.500.haplogroups = cent.500 %>% 
  group_by(chr,pos) %>% 
  nest() %>% 
  mutate(allele = map(data,haplogroup_name)) %>% 
  unnest() 

cent.500 <- rbind(filter(hap500, chr == "chr3R", pos < 7307440, !is.na(freq)),
                  filter(hap500, chr == "chr3L", pos > 22586280, !is.na(freq)))

cent.500 <- left_join(cent.500,select(cent.500.haplogroups,-freq)) %>% group_by(chr,pos,pool,haplogroup) %>%
  summarise(freq_sum = sum(freq, na.rm = TRUE))


#repeat for haplotypes based on 1M windows
cent.1000 <- rbind(filter(hap1000, chr == "chr3R", pos < 7307440, !is.na(freq), pool== "A.C") %>% select(-pool),
                   filter(hap1000, chr == "chr3L", pos > 22586280, !is.na(freq), pool== "A.C") %>% select(-pool))

cent.1000.haplogroups = cent.1000 %>% 
  group_by(chr,pos) %>% 
  nest() %>% 
  mutate(allele = map(data,haplogroup_name)) %>% 
  unnest() 

cent.1000 <- rbind(filter(hap1000, chr == "chr3R", pos < 7307440, !is.na(freq)),
                  filter(hap1000, chr == "chr3L", pos > 22586280, !is.na(freq)))

cent.1000 <- left_join(cent.1000,select(cent.1000.haplogroups,-freq)) %>% group_by(chr,pos,pool,haplogroup) %>%
  summarise(freq_sum = sum(freq, na.rm = TRUE))

#now we are ready to reanalyze the data, with haplogroup as a predictor instead of haplotype the anova function, significance is determined from the interaction term between haplogroup and treatment
anova_full <- function(df){
  df = data.frame(df)
  df$Nrep = as.factor(df$rep)
  df$trt = as.factor(df$trt)
  df$haplogroup = as.factor(df$haplogroup)
  df<- df %>% 
    mutate(afreq=asin(sqrt(freq_sum))) %>%
    select(-c(Nrep,freq_sum))
  tt <- tapply(df$afreq,df$haplogroup,mean) #get the mean of afreq for each haplotype
  tt2 <- names(tt)[tt > 0.14] #keep haplogroups whose afreq greater than 0.14 (haplotypes > 0.02)
  df2 <-  df %>% filter(haplogroup %in% tt2) %>% droplevels() #remove haplotypes which are not in tt2
  out = anova(lm(afreq~haplogroup+haplogroup:trt+haplogroup:rep,data=df2))
  -pf(out[2,3]/out[4,3],out[2,1],out[4,1],lower.tail=FALSE,log.p=TRUE)/log(10)	
}

haplogroup.LODS.200 <- cent.200  %>%  
  separate(pool,c("rep", "trt")) %>% 
  group_by(chr,pos) %>% 
  nest() %>% 
  mutate(LOD_full = map(data,anova_full)) %>% 
  select(-data) %>%
  unnest(cols = c(LOD_full)) %>%
  mutate(Ichr=recode(chr,'chrX'=1,'chr2L'=2,'chr2R'=3,'chr3L'=4,'chr3R'=5)) %>%
  unite("ID", chr:pos, remove = FALSE)

haplogroup.LODS.200$window.size <- "200 Kb"


haplogroup.LODS.500 <- cent.500  %>%  
  separate(pool,c("rep", "trt")) %>% 
  group_by(chr,pos) %>% 
  nest() %>% 
  mutate(LOD_full = map(data,anova_full)) %>% 
  select(-data) %>%
  unnest(cols = c(LOD_full)) %>%
  mutate(Ichr=recode(chr,'chrX'=1,'chr2L'=2,'chr2R'=3,'chr3L'=4,'chr3R'=5)) %>%
  unite("ID", chr:pos, remove = FALSE)

haplogroup.LODS.500$window.size <- "500 Kb"


haplogroup.LODS.1000 <- cent.1000  %>%  
  separate(pool,c("rep", "trt")) %>% 
  group_by(chr,pos) %>% 
  nest() %>% 
  mutate(LOD_full = map(data,anova_full)) %>% 
  select(-data) %>%
  unnest(cols = c(LOD_full)) %>%
  mutate(Ichr=recode(chr,'chrX'=1,'chr2L'=2,'chr2R'=3,'chr3L'=4,'chr3R'=5)) %>%
  unite("ID", chr:pos, remove = FALSE)

haplogroup.LODS.1000$window.size <- "1000 Kb"

# NEXT we loess smooth lod scores for all 3 haplotype windows

#drop chromosomes not represented from the data
haplogroup.LODS.200$chr <- droplevels(haplogroup.LODS.200$chr)

#for 200 Kb windows
# SET "span" value for loess smoothing
# - using a chromosome arm-specific value
#   so the # of markers for the fit is always
#   approx the same (at ~21)
# - assumes number of rows per chromosome arm is given by the counts in this table:
table(haplogroup.LODS.200$chr)
chr_ids <- c("chr3L","chr3R")
chr_loess_spans <- matrix(c(0.09606987, 0.09361702),ncol=1)
rownames(chr_loess_spans) <- chr_ids
colnames(chr_loess_spans) <- "loess_span_val"
chr_loess_spans

###now we loess smooth the LOD scores
outLODS.200 <- NULL
for(chr in levels(haplogroup.LODS.200$chr)){
  temp = haplogroup.LODS.200[haplogroup.LODS.200$chr==chr,] # EXTRACT 1 chromosome arm of data
  
  # Next line loess smoothing 
  # - degree = 2 (degree of polynomial; can be 1 or 2)
  # - span = a proportion (If # of points is N and span=0.5, then for given
  #   position (x) loess will use the 0.5 * N closest datapoints to x for the fit)
  fit1 <- loess(temp$LOD_full ~ temp$pos, degree=2, span = chr_loess_spans[rownames(chr_loess_spans)==chr,], family="gaussian")
  temp$smooLOD = fit1$fitted
  # KEEP the newly-smoothed dataset
  outLODS.200 <- rbind(outLODS.200,temp)
}

outLODS.200$window.size <- "200 Kb"


#for 500 Kb windows
haplogroup.LODS.500$chr <- droplevels(haplogroup.LODS.500$chr)
table(haplogroup.LODS.500$chr)
chr_ids <- c("chr3L","chr3R")
chr_loess_spans <- matrix(c(0.0876494, 0.09243697),ncol=1)
rownames(chr_loess_spans) <- chr_ids
colnames(chr_loess_spans) <- "loess_span_val"
chr_loess_spans
outLODS.500 <- NULL
for(chr in levels(haplogroup.LODS.500$chr)){
  temp = haplogroup.LODS.500[haplogroup.LODS.500$chr==chr,] 
  fit1 <- loess(temp$LOD_full ~ temp$pos, degree=2, span = chr_loess_spans[rownames(chr_loess_spans)==chr,], family="gaussian")
  temp$smooLOD = fit1$fitted
  outLODS.500 <- rbind(outLODS.500,temp)
}
outLODS.500$window.size <- "500 Kb"

#for 1000 Kb windows
haplogroup.LODS.1000$chr <- droplevels(haplogroup.LODS.1000$chr)
table(haplogroup.LODS.1000$chr)
chr_ids <- c("chr3L","chr3R")
chr_loess_spans <- matrix(c(0.088, 0.088),ncol=1)
rownames(chr_loess_spans) <- chr_ids
colnames(chr_loess_spans) <- "loess_span_val"
chr_loess_spans
outLODS.1000 <- NULL
for(chr in levels(haplogroup.LODS.1000$chr)){
  temp = haplogroup.LODS.1000[haplogroup.LODS.1000$chr==chr,] 
  fit1 <- loess(temp$LOD_full ~ temp$pos, degree=2, span = chr_loess_spans[rownames(chr_loess_spans)==chr,], family="gaussian")
  temp$smooLOD = fit1$fitted
  outLODS.1000 <- rbind(outLODS.1000,temp)
}
outLODS.1000$window.size <- "1 Mb"


#join the data into a single object for graphing
cent.3 <- rbind(outLODS.200,outLODS.500,outLODS.1000)
#order window size factor levels so they are logical
cent.3$window.size <- factor(cent.3$window.size, levels = c("200 Kb", "500 Kb", "1 Mb"))

### the delta 3LOD window can be determined by manual inspection of the the outLODS objects. Here we draw a box that will represent the interval graphically. For 500 Kb windows is is 3L:25576280
rect.tbl <- tibble(x1 = c(25586280,2307440), x2=c(27866280,3467440), y1 = rep(0,2), y2 = rep(5,2),chr=c("chr3L","chr3R"))

#build plot
ggplot(cent.3, aes(pos,smooLOD,color=window.size)) +
  geom_rect(data=rect.tbl,aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2),fill="#866e3c", alpha=.2,inherit.aes = FALSE) +geom_line() + 
   theme_classic() + facet_grid(~chr,scales="free") +
   scale_color_manual(values=c("#Db2f18","#db751b","#86443c")) + 
  geom_hline(yintercept = 4, linetype = "dashed", colour = "#866e3c") +  
   theme_classic() + ylab("-log10(p)\n") + theme(legend.position="bottom")+
   xlab("\nPhysical Location (Mb)") + scale_x_continuous(breaks = c(25586280,3467440))

ggsave("haplogroup_mapping.pdf",width=3.5,height=3)

### in the second part of this script we evaluate how haplotype frequencies differ between experimental and control samples 

###calculate allele frequency changes (ratio of treatment to control) for each replicate
Fs.200 = cent.200 %>% filter(!is.na(freq_sum)) %>% 
  separate(pool,c("rep","trt")) %>%
  filter(chr != "chr4") %>%
  droplevels() %>%
  pivot_wider(names_from = trt, values_from = freq_sum) %>%
  mutate(freqdff = T/C) %>%
  select(chr,pos,haplogroup,freqdff,rep) 

##isolate frequency changes at the lod peak in 200 Kb mapping is chr3R 3327440
peak.200.Fs <- filter(Fs.200, chr == "chr3R", pos == 3327440)


###inspect the raw haplotype frequencies at this position to determine who should be dropped and which haplotypes are grouped
filter(hap200, chr == "chr3R", pos == 3327440)
#B1/B5 are basically absent freq << 0.02
#B2/B3/B4/B6 are a haplogroup, we can't tell them apart




###drop haplogroup B1 from the plot because this haplotype is rare/absent not included in the calculation of the LOD score

ggplot(filter(peak.200.Fs,haplogroup != "B1") , aes(haplogroup, freqdff, colour = haplogroup)) +
  theme_classic() +
  stat_summary(geom="bar",fun="mean",fill = "white",size=1.5) +
  geom_point(size=1.5) +
  scale_color_manual(values=palette2) +
  scale_x_discrete(labels=c("B2" = "B2,B3,B4,B6"))+ 
  geom_hline(yintercept=1,linetype=2) +
  theme(legend.position = "none") +
  xlab("\nfounder haplotype") + ylab("relative haplotype frequency (T/C)\n") 


###repeat for 500 Kb windows
Fs.500 = cent.500 %>% filter(!is.na(freq_sum)) %>% 
  separate(pool,c("rep","trt")) %>%
  filter(chr != "chr4") %>%
  droplevels() %>%
  pivot_wider(names_from = trt, values_from = freq_sum) %>%
  mutate(freqdff = T/C) %>%
  select(chr,pos,haplogroup,freqdff,rep) 

##lod peak in 200 Kb mapping is chr3R 3317440
peak.500.Fs <- filter(Fs.500, chr == "chr3R", pos == 3317440)

###inspect the raw haplotype frequencies at this position to determine who should be dropped and which haplotypes are grouped
filter(hap500, chr == "chr3R", pos == 3317440)
#B1/B5 are basically absent freq << 0.02
#B3/B4 are a haplogroup, we can't tell them apart
#B2/B6 are a haplogroup, we can't tell them apart.


ggplot(filter(peak.500.Fs,haplogroup != "B1", haplogroup != "B5"), aes(haplogroup, freqdff, colour = haplogroup)) +
  theme_classic() +
  stat_summary(geom="bar",fun="mean",fill = "white",size=1.5) +
  geom_point(size=1.5) +
  scale_color_manual(values=c("#338a7d","#866e3c","#Db2f18","#db715b")) +theme(legend.position = "none") +
  scale_x_discrete(labels=c("B2" = "B2,B6","B3"="B3,B4"))+ 
  geom_hline(yintercept=1,linetype=2) +
  xlab("\nfounder haplotype") + ylab("relative haplotype frequency (T/C)\n") 
ggsave("deltafreq.pdf", height =3, width=3)


#reanalyze the data to see if some haplotypes are more important than others 
#drop B1/B5 because not present
peak.500 <- filter(cent.500, chr == "chr3R", pos == 3317440,haplogroup != "B1",haplogroup != "B5")
peak.500$afreq <- asin(sqrt(peak.500$freq_sum))
peak.500 <- peak.500 %>% separate(pool,c("rep","trt"))
peak.500.lm <- lm(afreq~haplogroup+haplogroup:trt+haplogroup:rep,peak.500)
-pf(anova(peak.500.lm)[2,3]/anova(peak.500.lm)[4,3],anova(peak.500.lm)[2,1],anova(peak.500.lm)[4,1],lower.tail=FALSE,log.p=TRUE)/log(10)


###drop B2/B6 LOD score drops 
peak.500.lm <- lm(afreq~haplogroup+haplogroup:trt+haplogroup:rep,filter(peak.500,haplogroup != "B2"))
anova(peak.500.lm)
-pf(anova(peak.500.lm)[2,3]/anova(peak.500.lm)[4,3],anova(peak.500.lm)[2,1],anova(peak.500.lm)[4,1],lower.tail=FALSE,log.p=TRUE)/log(10)


###drop B3/B4 LOD score increases
peak.500.lm <- lm(afreq~haplogroup+haplogroup:trt+haplogroup:rep,filter(peak.500,haplogroup != "B3"))
anova(peak.500.lm)
-pf(anova(peak.500.lm)[2,3]/anova(peak.500.lm)[4,3],anova(peak.500.lm)[2,1],anova(peak.500.lm)[4,1],lower.tail=FALSE,log.p=TRUE)/log(10)

###drop B7 LOD score drops 
peak.500.lm <- lm(afreq~haplogroup+haplogroup:trt+haplogroup:rep,filter(peak.500,haplogroup != "B7"))
anova(peak.500.lm)
-pf(anova(peak.500.lm)[2,3]/anova(peak.500.lm)[4,3],anova(peak.500.lm)[2,1],anova(peak.500.lm)[4,1],lower.tail=FALSE,log.p=TRUE)/log(10)

###drop B8 LOD score drops, still significant
peak.500.lm <- lm(afreq~haplogroup+haplogroup:trt+haplogroup:rep,filter(peak.500,haplogroup != "B8"))
anova(peak.500.lm)
-pf(anova(peak.500.lm)[2,3]/anova(peak.500.lm)[4,3],anova(peak.500.lm)[2,1],anova(peak.500.lm)[4,1],lower.tail=FALSE,log.p=TRUE)/log(10)


# let's look at the haplogroup frequencies in the control samples
filter(peak.500, trt == "C",chr=="chr3R",pos==3317440)

ggplot(filter(peak.500, trt == "C",chr=="chr3R"), aes(x=rep,y=freq_sum,fill=haplogroup)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c("#338a7d","#866e3c","#Db2f18","#db715b")) +theme(legend.position = "none") +
  theme_classic() + xlab("\nreplicate") +
  ylab("haplogroup frequency\n")
ggsave("haplogroup_freq.pdf", height = 3, width = 3 )



