runscan <- function(pool, OutName, folder, SNPtable, foundernames){

	#  pool = pool to get frequencies for
	#  OutName = prefix for output file
	#  folder = output folder for results
	#  SNPtable = file of SNP frequencies in pool and founders.  The structure of the file matters
	#     the 1st column is a chromosome name, and the 2nd a position, there can be additional columns
	#     eventually you want paired columns of the 1st is the frequency in some sample and the 2nd the coverage
	#     I have a python script that makes these files from "vcf" files, vcf are a very standard format
	#     this file should be preprocessed with GrannySNP.R
	#  foundernames = file with the names of the founder haplotypes (and a short name) 

	poolroot=pool
	yy = read.table(file=foundernames,header=FALSE)
	founderroots = as.character(yy[,1])
	foundershortnames = as.character(yy[,2])
	filename=paste(folder,"/",OutName,"_hap_freq.txt",sep='')
	
	if (!file.exists(filename)){
		xx = read.table(file=SNPtable)
		# the idea here is that freq columns can be preceeded with "freq_" or "f"
		Npool <-  match(paste("freq_",pool,sep=''),names(xx))
		founders <- match(paste("freq_",founderroots,sep=''),names(xx))
		Npool_alt <-  match(paste("f",pool,sep=''),names(xx))
		founders_alt <- match(paste("f",founderroots,sep=''),names(xx))
		if(sum(!is.na(Npool_alt))>sum(!is.na(Npool))){Npool=Npool_alt}
		if(sum(!is.na(founders_alt))>sum(!is.na(founders))){founders=founders_alt}

		Nf = length(founders)
		stepSize = 20000
		WINSIZE = 200000
		names = paste(foundershortnames, collapse=";")
		lppp = 0
		for(chr in c("chrX","chr2L","chr2R","chr3L","chr3R")){
			Maxpos = max(xx$POS[xx$CHROM == chr],na.rm=TRUE)
			Minpos = min(xx$POS[xx$CHROM == chr],na.rm=TRUE)
			ppp = seq(Minpos+WINSIZE,Maxpos-WINSIZE,stepSize)
			lppp = lppp + length(ppp)
			}
		# the data frame is going more tidy-like
		lppp = lppp * Nf
		ddd = data.frame(chr=rep("",lppp),pos=rep(0,lppp),pool=rep("",lppp),founder=rep("",lppp),freq=rep(NA,lppp),stringsAsFactors=FALSE)
		i=1
		for(chr in c("chrX","chr2L","chr2R","chr3L","chr3R")){
			Maxpos = max(xx$POS[xx$CHROM == chr],na.rm=TRUE)
			Minpos = min(xx$POS[xx$CHROM == chr],na.rm=TRUE)
			ppp = seq(Minpos+WINSIZE,Maxpos-WINSIZE,stepSize)
			for(pos in ppp){
				#pos=9700
				#pos=50700
				# median non-recombined block size is 14kb
				#  so step size should be much smaller than that (500bp = 1/30th)
				#  I would guess optimum window size is +/- 7kb (= 14 total) ish.  
				#  but markers are limiting...
				MynumberSNP<- NA
				Groups = rep(NA,Nf)
				groups = NA
				nGroups = NA
				Freqs = rep(NA,Nf)
				GroupAvFreqs = rep(NA,Nf)

				if(sum(xx$CHROM == chr & xx$POS > pos - WINSIZE & xx$POS < pos + WINSIZE)>100){
					predictors <- subset(xx, CHROM == chr & POS > pos - WINSIZE & POS < pos + WINSIZE, select=founders)
					Y <- subset(xx, CHROM == chr & POS > pos - WINSIZE & POS < pos + WINSIZE, select=Npool)
					MynumberSNP <- nrow(predictors)
					# pick a sigma, near the ends of chromosomes sigma has to be bigger than in the middle
					predictNotMissing = apply(predictors,1,function(x) sum(is.na(x))==0) & (!is.na(Y))
					predictors = predictors[predictNotMissing,]
					Y = Y[predictNotMissing,]

					Groups <- cutree(hclust(dist(t(predictors))),h=5)
					nGroups = nlevels(as.factor(Groups))
					groups = paste(as.numeric(Groups),collapse=";")
					d = ncol(predictors)
					A = predictors
					E = t(matrix(rep(1,d)))
					F = 1
					G = diag(rep(1,d))
					H = matrix(rep(0.0003,d))
					if(nrow(A)>=100){
						out = lsei(A=A,B=Y,E=E,F=F,G=G,H=H,verbose=TRUE)
						if(!out$IsError){
							Freqs <- out$X
							meanbygroup <- tapply(Freqs,as.numeric(Groups),mean)
							GroupAvFreqs <- meanbygroup[as.numeric(Groups)]      # freqs replaced by group means
							} # model could be fitted
						} # nrow >= 100
				} # >100 raw SNPs present
				
				temp = cbind(rep(chr,Nf), rep(pos,Nf), rep(pool,Nf), as.character(foundershortnames), round(as.numeric(GroupAvFreqs),4))
				ddd[i:(i+Nf-1),] = temp
				i=i+Nf
					
				}  # over positions
			}  # chr
		write.table(ddd,filename,row.names=FALSE,sep="\t",quote=FALSE)
		} # file does not exist
	} # function


