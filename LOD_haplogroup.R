library(viridis)
library(tidyverse)
library(ggplotify)	
library(gridExtra)

#import haplotypes for 200 Kb and 500 Kb windows
hap200 <- read_delim("allhaps.200Kb.txt",col_names=TRUE) 
hap500 <- read_delim("allhapsXRAY.500Kb.txt",col_names=TRUE) 
hap1000 <- read_delim("allhapsXray500K.txt",col_names=TRUE) 

#make AB8 values into B8 so they are sorted last.
hap200$founder[hap200$founder == "AB8"] <- "B8"
hap500$founder[hap500$founder == "AB8"] <- "B8"
hap1000$founder[hap1000$founder == "AB8"] <- "B8"
###extract windows from the right side of the 3rd chromosome centromere, exclude windows were AF is not estimated
###we only need values from one replicate:treatment. We will use identical AF estimates to detect pairs and groups of haplotypes that cannot be differentiated from each other within a window.
cent.3R.200 <- filter(hap200, chr == "chr3R", pos < 5000000, !is.na(freq), pool== "A.C") %>% select(-pool)




#the function below will take a list of 8 haplotype frequencies in the same position, sort them, keep the unique values first founder name (in ascending order 1-8)
#the each uniquely identifiable haplotypes is then assigned a "haplogroup" name based o the first founder that has this haplotype
haplogroup_name <- function(df){
  df <- data.frame(df)
  left_join(df, arrange(df, freq,founder) %>% distinct(freq,.keep_all=TRUE) %>% rename(c("haplogroup"="founder"))) 
}


## groups data in cent.3R by chromosome and position and assigns a founder name to each unique haplogroup
cent.3R.200 <- cent.3R.200 %>% 
  group_by(chr,pos) %>% 
  nest() %>% 
  mutate(allele = map(data,haplogroup_name)) %>% 
  unnest() 
###describe the window size so that the data sites can be compared 
cent.3R.200$window.size <- "200Kb"

#repeat for haplotypes based on 500 Kb windows
cent.3R.500 <- filter(hap500, chr == "chr3R", pos < 5000000, !is.na(freq), pool== "A.C") %>% select(-pool)
cent.3R.500 = cent.3R.500 %>% 
  group_by(chr,pos) %>% 
  nest() %>% 
  mutate(allele = map(data,haplogroup_name)) %>% 
  unnest() 

cent.3R.500$window.size <- "500Kb"

#repeat for haplotypes based on 1M windows
cent.3R.1000 <- filter(hap1000, chr == "chr3R", pos < 5000000, !is.na(freq), pool== "A.C") %>% select(-pool)
cent.3R.1000 = cent.3R.1000 %>% 
  group_by(chr,pos) %>% 
  nest() %>% 
  mutate(allele = map(data,haplogroup_name)) %>% 
  unnest() 

cent.3R.1000$window.size <- "1 Mb"

#join the data into a single object for graphing
cent.3R <- rbind(cent.3R.200,cent.3R.500,cent.3R.1000)

#set the palette for the graph:
palette2 <- c("#18dbbe","#db751b","#86443c","#332f2c","#2d3d3b","#338a7d","#Db2f18","#866e3c")

#make a graph to compare the window sizes:

rectanges <- xmin


ggplot(filter(cent.3R,founder!="B1",founder!="B5"), aes(pos,founder,fill=haplogroup)) +
  geom_tile() +geom_segment(aes(x=3127440,xend=3347440, y=2,yend=2))+
  geom_segment(aes(x=2677440,xend=3437440, y=1,yend=1)) + 
  theme_classic() + facet_grid(window.size~.) +
  scale_fill_manual(values=palette2) + theme_classic() + xlab("\n3R position") + ylab("founder haplotype") 


#now we need to make new haplotype table that sums frequencies across the haplogroup at a given position

cent.haplogroups.200 <- left_join(filter(hap200, chr == "chr3R",!is.na(freq), pos< 5000000),select(cent.3R.200,-freq)) %>% group_by(chr,pos,pool,haplogroup) %>%
  summarise(freq_sum = sum(freq, na.rm = TRUE))


#the anova function, significance is determined from the interaction term between haplogroup and treatment
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

haplogroup.LODS.200 <- cent.haplogroups.200  %>%  
  separate(pool,c("rep", "trt")) %>% 
  group_by(chr,pos) %>% 
  nest() %>% 
  mutate(LOD_full = map(data,anova_full)) %>% 
  select(-data) %>%
  unnest(cols = c(LOD_full)) %>%
  mutate(Ichr=recode(chr,'chrX'=1,'chr2L'=2,'chr2R'=3,'chr3L'=4,'chr3R'=5)) %>%
  unite("ID", chr:pos, remove = FALSE)

#repeat for the larger window.size
cent.haplogroups.500 <- left_join(filter(hap500, chr == "chr3R",!is.na(freq), pos< 5000000),select(cent.3R.500,-freq)) %>% group_by(chr,pos,pool,haplogroup) %>%
  summarise(freq_sum = sum(freq, na.rm = TRUE))

haplogroup.LODS.500 <- cent.haplogroups.500  %>%  
  separate(pool,c("rep", "trt")) %>% 
  group_by(chr,pos) %>% 
  nest() %>% 
  mutate(LOD_full = map(data,anova_full)) %>% 
  select(-data) %>%
  unnest(cols = c(LOD_full)) %>%
  mutate(Ichr=recode(chr,'chrX'=1,'chr2L'=2,'chr2R'=3,'chr3L'=4,'chr3R'=5)) %>%
  unite("ID", chr:pos, remove = FALSE)

filter(haplogroup.LODS.200,LOD_full>4)
filter(haplogroup.LODS.500,LOD_full>4)

