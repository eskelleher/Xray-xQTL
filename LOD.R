library(tidyverse)


# import the estimated haplotype frequencies. The data set below uses 200 Kb windows and 20 Kb step sizes, we also provide raw data 500 Kb and 1 Mb windows. 
hap <- read_delim("allhaps.200Kb.txt",col_names=TRUE) 

#functions to fit linear models for the full experiment and extract the appropriate P-value
## in this model, haplotypes with frequency below 0.02 are dropped from the analysis. Rare haplotypes will not respond much to a single generation of selection anyways, so this gives more robust mapping results.
anova_full <- function(df){
  df = data.frame(df)
  df$Nrep = as.factor(df$rep)
  df$trt = as.factor(df$trt)
  df$founder = as.factor(df$founder)
  df<- df %>% 
    mutate(afreq=asin(sqrt(freq))) %>%
    select(-c(Nrep,freq))
  tt <- tapply(df$afreq,df$founder,mean) #get the mean of afreq for each haplotype
  tt2 <- names(tt)[tt > 0.14] #keep haplotypes whose afreq greater than 0.14 (haplotypes > 0.02)
  df2 <-  df %>% filter(founder %in% tt2) %>% droplevels() #remove haplotypes which are not in tt2
  out = anova(lm(afreq~founder+founder:trt+founder:rep,data=df2))
  -pf(out[2,3]/out[4,3],out[2,1],out[4,1],lower.tail=FALSE,log.p=TRUE)/log(10)	
}

LODS = hap %>% filter(!is.na(freq)) %>% 
  separate(pool,c("rep", "trt")) %>% 
  filter(chr != "chr4") %>% 
  droplevels() %>%
  group_by(chr,pos) %>% 
  nest() %>% 
  mutate(LOD_full = map(data,anova_full)) %>% 
  select(-data) %>%
  unnest(cols = c(LOD_full)) %>%
  mutate(Ichr=recode(chr,'chrX'=1,'chr2L'=2,'chr2R'=3,'chr3L'=4,'chr3R'=5)) %>%
  unite("ID", chr:pos, remove = FALSE)


###view significant positions and manually inspect for 3 LOD peaks
filter(LODS, LOD_full>4)



##  function to return a df with cumulative bp positions
linearize_genome = function(df){
  df %>% 
    # Compute chromosome size
    group_by(Ichr) %>% 
    summarise(chr_len=max(pos)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(df, ., by=c("Ichr"="Ichr")) %>%
    # Add a cumulative position of each SNP
    arrange(Ichr, pos) %>%
    mutate( BPcum=pos+tot)
}


plot.LODs <- linearize_genome(LODS)

# change the value of Ichr to 7 for LOD scores above threshold so that they appear as a different color.
plot.LODs$Ichr[plot.LODs$LOD_full>4] <- 7

## function to make plot	
make.Manhattan = function(df,mychar,mylab,threshold,ylimit){
  mychar=sym(mychar)
  myaxis = df %>%
    group_by(Ichr) %>% 
    summarize(center=( max(BPcum) + min(BPcum) ) / 2, chrlab = chr[1] )
  
  ggplot(df, aes_string(x="BPcum", y=mychar)) +
    ylab("-log10(p)\n") +
    xlab("\nPhysical Location (Mb)") +
    ggtitle(mylab) + 
    theme(plot.title = element_text(vjust = - 10, hjust=0.025, size=10)) +
     geom_point( aes(color=as.factor(Ichr),size = as.factor(Ichr), alpha = as.factor(Ichr))) +
    scale_color_manual(values = c("grey30", "grey70", "grey30", "grey70", "grey30","#Db2f18")) +
    scale_size_manual(values = c(rep(0.5,5),1))+
    scale_alpha_manual(values = c(rep(0.3,5),1))+
    # threshold
    geom_hline(yintercept = threshold, linetype = "dashed", colour = "#866e3c") +  
    # custom X axis:
    scale_x_continuous(label = myaxis$chrlab, breaks= myaxis$center ) +
    scale_y_continuous(expand = c(0, 0), limits=c(0,ylimit) ) +     # remove space between plot area and x axis
    # Custom the theme:
    theme_classic() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
    theme(axis.text=element_text(size=9),axis.title=element_text(size=10)) + 
    theme(legend.position = "none") 
}

make.Manhattan(plot.LODs,"LOD_full","200 Kb Window Size (Standard)",4,12)
ggsave("LOD.plot.200.pdf",width=3.5,height=2.75)

#filter out only replicate A to identify haplotypes that can't be resolved, since this is a function of the genetic variation (or lack thereof) present in the window. 
#In these cases two founder haplotypes will have the AF and hence the same AF ratio



one.rep <- filter(hap, pool =="A.C")

###count the number of unique alleles at each position. 
hapcount <- tibble(
  "pos" = tapply(one.rep$freq, one.rep$pos, unique) %>% lapply(length) %>% t() %>% colnames(), 
  "unique" = tapply(one.rep$freq, one.rep$pos, unique) %>% lapply(length) %>% unlist()
) 

hapcount$pos <- as.numeric(hapcount$pos)

# append hapcounts to LOD scores
plot.LODs <- left_join(plot.LODs,hapcount)

###make a plot showing the number of unique haplotypes at each position

make.Manhattan = function(df,mychar,mylab,threshold,ylimit){
  mychar=sym(mychar)
  myaxis = df %>%
    group_by(Ichr) %>% 
    summarize(center=( max(BPcum) + min(BPcum) ) / 2, chrlab = chr[1] )
  
  ggplot(df, aes_string(x="BPcum", y=mychar)) +
    ylab("-log10(p)\n") +
    xlab("\nPhysical Location (Mb)") +
    ggtitle(mylab) + 
    theme(plot.title = element_text(vjust = - 10, hjust=0.025, size=10)) +
     geom_point(size=.3, aes(color=unique, alpha = as.factor(Ichr))) +
    scale_color_gradient(low= "#Db2f18", high = "#332f2c") +
     scale_alpha_manual(values = c(rep(0.3,5),1))+
        # threshold
    geom_hline(yintercept = threshold, linetype = "dashed", colour = "#866e3c") +  
    # custom X axis:
    scale_x_continuous(label = myaxis$chrlab, breaks= myaxis$center ) +
    scale_y_continuous(expand = c(0, 0), limits=c(0,ylimit) ) +     # remove space between plot area and x axis
    # Custom the theme:
    theme_classic() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
    theme(axis.text=element_text(size=9),axis.title=element_text(size=10)) +
    theme(legend.position="none")
  
  
}

make.Manhattan(plot.LODs,"LOD_full","200 Kb Window Size (Standard)",4,12)
ggsave("LOD.plot.hap.pdf",width=3.5,height=2.75)
