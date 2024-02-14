### Script to identify a specific read matching dwell-time and current intensity from the 
### nanocompore database with the counterparts from nanopolish eventalign.

## All the conditions of interest es:
samples <- c("WT1","WT2","WT3")

## Loading of nanopolish eventalign tsv files
eventalignReports <- lapply(samples,function(i)
{
	reportTmp <- read.table(paste0("path/to/nanopolish/results/",i,"/nanopolishcomp/out_eventalign_collapse.tsv")
	,col.names=c("ref_pos","ref_kmer","num_events","dwell_time","NNNNN_dwell_time","mismatch_dwell_time","mean","median","num_signals")
	,fill=TRUE,header=FALSE,comment.char="")

	# From the eventalign format to a matrix of sites for TARDBP
	reportTmp <- reportTmp[-(which(reportTmp[,1]=="ref_pos")),]
	readsIdxTmp <- which(reportTmp[,2]=="ENSMUST00000084125::4:148612381-148627019")+1
	readsIdxTmp <- c(readsIdxTmp,(nrow(reportTmp)+1))
	readsTmp <- reportTmp[which(reportTmp[,2]=="ENSMUST00000084125::4:148612381-148627019"),1]

	reportListTmp <- lapply(seq_along(readsIdxTmp)[-1],function(k)
	{
		reportTmp[readsIdxTmp[[(k-1)]]:(readsIdxTmp[[(k)]]-2),]
	})
	
	readsTmp <- rep(readsTmp,sapply(reportListTmp,nrow))
	reportListTmp <- do.call("rbind", reportListTmp)
	cbind(readsTmp,reportListTmp)
})
names(eventalignReports) <- samples

## Probabilities files loading - from singleMolecule.py
sites <- list.files(".",pattern="singleDwellTimes")
sites <- gsub("singleDwellTimes_","",sites)

sitesData <- lapply(sites,function(site)
{
	# Sites data formatting
	siteData <- cbind("sample"=read.table(paste0("singleMoleculeCond_",site),header=FALSE)
					 ,"class"=read.table(paste0("singleMoleculeClass_",site),header=FALSE)
					 ,"p"=read.table(paste0("singleMoleculeProba_",site),header=FALSE)
					 ,"dwell"=read.table(paste0("singleDwellTimes_",site),header=FALSE)
					 ,"current"=read.table(paste0("singleIntensityTimes_",site),header=FALSE))
	siteData[,1] <- gsub("_","",siteData[,1])
	
	# Site to read match
	siteReads <- sapply(1:nrow(siteData),function(i)
	{
		j <- siteData[i,]
		reportTmp <- eventalignReports[j[[1]]][[1]]
		# Identical current intensity
		reportTmp <- reportTmp[which(as.numeric(reportTmp$median) == as.numeric(j[[5]])),]
		# Identical position
		reportTmp <- reportTmp[which(reportTmp$ref_pos == strsplit(site,"_")[[1]][[1]]),]
		# Closer dwell time
		reportTmp <- reportTmp[which.min(abs(as.numeric(reportTmp$dwell_time) - as.numeric(j[[4]]))),]
		reportTmp$readsTmp
	})
	
	siteData <- cbind(siteReads,siteData)
	colnames(siteData) <- c("read_id","sample","class","proba")
	siteData
})
names(sitesData) <- gsub(".csv","",sites)
sitesData <- sitesData[sort(names(sitesData))]
saveRDS(sitesData,"TARDBP_sitesData.rds")

## Reads spanning all the sites
commonReads <- intersect(sitesData[[1]][,1],sitesData[[2]][,1])
for(i in 3:length(sitesData))
{
	commonReads <- intersect(commonReads,sitesData[[i]][,1])
}

## Selection of the m6A cluster
sitesP <- sapply(sitesData,function(i)
{
	rownames(i) <- i[,1]
	i <- i[commonReads,]
	classWTtmp <- table(i[grepl("WT",i[,2]),3])/length(grep("WT",i[,2]))
	classKDtmp <- table(i[grepl("KD",i[,2]),3])/length(grep("KD",i[,2]))

	classIdxTmp <- which(classKDtmp<classWTtmp)

	pTmp <- as.numeric(sapply(strsplit(i[,"proba"],","),"[[",classIdxTmp))
})
rownames(sitesP) <- commonReads

saveRDS(sitesP,"TARDBP_sitesProbabilities.rds")