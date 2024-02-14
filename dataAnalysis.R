####################################################################################
### Code for Nanocompore post-processing analyses, and gene-expression analyses .###
####################################################################################

### Libraries
  tryCatch(library(tidyverse),error=function(e){print("Missing tidyverse")}) # Basic libraries for graphic outputs and data input
  tryCatch(library(motifStack),error=function(e){print("Missing motifStack")}) # For graphic representation of the motif
  tryCatch(library(Guitar),error=function(e){print("Missing Guitar")}) # For guitar plot
  tryCatch(library(cowplot),error=function(e){print("Missing cowplot")}) # For guitar plot
  tryCatch(library(lemon),error=function(e){print("Missing lemon")}) # For guitar plot
  tryCatch(library(clusterProfiler),error=function(e){print("Missing clusterProfiler")}) # For enrichments
  tryCatch(library(biomaRt),error=function(e){print("Missing biomaRt")}) # For enrichments
  tryCatch(library(org.Mm.eg.db),error=function(e){print("Missing org.Mm.eg.db")}) # For enrichments
  tryCatch(library(TxDb.Mmusculus.UCSC.mm10.knownGene),error=function(e){print("Missing TxDb.Mmusculus.UCSC.mm10.knownGene")}) # For guitar plot
  tryCatch(library(compEpiTools),error=function(e){print("Missing compEpiTools")}) # For enrichments
  tryCatch(library(GenomicAlignments),error=function(e){print("Missing GenomicAlignments")}) # For genes counts
  tryCatch(library(pheatmap),error=function(e){print("Missing pheatmap")}) # For replicates comparison
  tryCatch(library(DESeq2),error=function(e){print("Missing DESeq2")}) # For counts estimation
  tryCatch(library(BSgenome.Mmusculus.UCSC.mm10),error=function(e){print("Missing BSgenome.Mmusculus.UCSC.mm10")}) # For motif analysis
  tryCatch(library(readxl),error=function(e){print("Missing readxl")}) # For motif analysis
  tryCatch(library(msigdbr),error=function(e){print("Missing msigdbr")}) # For enrichments

### Data
  ## Nanocompore
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

  nanocomporeOutput <- "/path/to/out_nanocompore_results.tsv"
  mouseAnnotation <- read_tsv("/path/to/Mus_musculus.GRCm38.96.info.txt", col_names=c("GeneId", "GeneName", "TxId", "TxName", "TxBiotype"))

  ## Association of genes and transcripts from the annotation to the k-mers reported by nanocompore
  nanocomporeResults <- read_tsv(nanocomporeOutput, col_types="icicccddddddddcicdd") %>%
                        mutate(TxId=gsub("::.+", "", ref_id)) %>%
                        left_join(mouseAnnotation)

  drach <- c("AAACA", "AAACC", "AAACT", "AGACA", "AGACC", "AGACT",
             "GAACA", "GAACC", "GAACT", "GGACA", "GGACC", "GGACT", 
             "TAACA", "TAACC", "TAACT", "TGACA", "TGACC", "TGACT")

  ## Gene expression
  # We load in R each bam file obtained through the standard Nanocompore pipeline
  # (alignment to the transcriptome thorugh minimap2 genome mode), and we estimate the
  # number of reads per transcript.
  #
  # bam <- readGAlignments("path/to/minimap.filt.sort.bam")
  # counts <- table(grep("ENSMUST",unlist(strsplit(as.character(seqnames(bam)),"::")),value=TRUE))

  # saveRDS(list("WT1"=countsWT1,"WT2"=countsWT2,"WT3"=countsWT3,"KD1"=countsKD1,"KD2"=countsKD2,"KD3"=countsKD3),"counts.rds")
  counts <- readRDS("counts.rds")

### DESeq2
  commonTranscripts <- intersect(intersect(intersect(names(counts[["WT1"]]),names(counts[["WT2"]])),names(counts[["WT3"]]))
                                ,intersect(intersect(names(counts[["KD1"]]),names(counts[["KD2"]])),names(counts[["KD3"]])))

  counts <- sapply(counts,function(i)i[commonTranscripts])

  ## DESeq2 - standard pipeline from the vignette to quantify normalized counts and genes differentially expressed
  cts <- as.matrix(counts)
  coldata <- data.frame("condition"=c("untreated","untreated","untreated","treated","treated","treated"))
  rownames(coldata) <- c("WT1","WT2","WT3","KD1","KD2","KD3")
  coldata$condition <- factor(coldata$condition)

  all(rownames(coldata) == colnames(counts))

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = ~ condition)

  dds$condition <- relevel(dds$condition, ref = "untreated")

  dds <- DESeq(dds)
  
  DESeq2_normalized_counts <- counts(dds, normalized=TRUE)

  DESeq2_results <- results(dds, alpha=0.05)

  ## Clustering between replicated based on normalized counts pearson correlation
  pdf("replicatesClustering_DESeq2.pdf", width=4, height=4)
  plot(hclust(as.dist(1-cor(DESeq2_normalized_counts))),xlab="", ylab="1 - Pearson Correlation")
  dev.off()

  ## Normalized counts pearson correlation matrix
  pdf("replicatesCorrelation_DESeq2.pdf", width=4, height=4)
  pheatmap(cor(DESeq2_normalized_counts), display_numbers=TRUE)
  dev.off()

  DESeq2_normalized_counts <- cbind(DESeq2_normalized_counts
								   ,"WT"=apply(DESeq2_normalized_counts[,1:3],1,sum)
								   ,"KD"=apply(DESeq2_normalized_counts[,4:6],1,sum))

  ## TARDBP expression analysis - LESS THAN 6% OF THE TRANSCRIPTS ARE MORE EXPRESSED THAN TARDBP IN WT
  table(WT=DESeq2_normalized_counts[,"WT"]>DESeq2_normalized_counts["ENSMUST00000084125","WT"]
  	   ,KD=DESeq2_normalized_counts[,"KD"]>DESeq2_normalized_counts["ENSMUST00000084125","KD"])

  round(table(WT=DESeq2_normalized_counts[,"WT"]>DESeq2_normalized_counts["ENSMUST00000084125","WT"]
             ,KD=DESeq2_normalized_counts[,"KD"]>DESeq2_normalized_counts["ENSMUST00000084125","KD"])/nrow(DESeq2_normalized_counts),2)

  ## TARDBP targets expression analysis
  # Ensembl genes
  ensembl <- useMart("ensembl",host = "https://dec2021.archive.ensembl.org/")
  ensemblH = useDataset("hsapiens_gene_ensembl",mart=ensembl)

  ensemblGenesH <- getBM(attributes = c("ensembl_gene_id","start_position","end_position","strand","chromosome_name"), mart = ensemblH)
  ensemblGenesH[ensemblGenesH[,"strand"]==1,"strand"] <- "+"
  ensemblGenesH[ensemblGenesH[,"strand"]==(-1),"strand"] <- "-"

  colnames(ensemblGenesH) <- c("id","start","end","strand","chr")
  
  ensemblGenesH <- makeGRangesFromDataFrame(ensemblGenesH,keep.extra.columns=TRUE)
  seqlevelsStyle(ensemblGenesH) <- "UCSC"

  ensemblGenesHwidth <- width(ensemblGenesH)
  names(ensemblGenesHwidth) <- ensemblGenesH$id

  # iCLIP sites from https://doi.org/10.1016/j.cell.2021.07.018
  TARDBPsites <- read.table("tardbp-293flp-2-20181129-ju_trimmed_reads_single.bed")
  colnames(TARDBPsites) <- c("chr","start","end","","score","strand")

  TARDBPsites <- makeGRangesFromDataFrame(TARDBPsites,keep.extra.columns=TRUE )

  TARDBPgenes <- findOverlaps(TARDBPsites,ensemblGenesH)
  TARDBPgenes <- TARDBPgenes[isUnique(TARDBPgenes@from),]
  
  TARDBPsites <- TARDBPsites[TARDBPgenes@from,]
  TARDBPcounts <- TARDBPsites$score
  names(TARDBPcounts) <- ensemblGenesH$id[TARDBPgenes@to]

  # Top ranking genes for iCLIP signal
  TARDBPcounts <- split(TARDBPcounts,names(TARDBPcounts))

  TARDBPnormalizedRatio <- sort(sapply(TARDBPcounts,length)/ensemblGenesHwidth[names(TARDBPcounts)],decreasing=TRUE)

  plot(x=seq_along(TARDBPnormalizedRatio),y=log10(TARDBPnormalizedRatio),pch=".",ylab="Log10(iCLIP sites/length)",xlab="index")
  points(which(names(TARDBPnormalizedRatio)=="ENSG00000120948"),log10(TARDBPnormalizedRatio["ENSG00000120948"]),col=2,pch=16)

  abline(v=500*1:5,col="lightgrey",lty=2)

  # Top 500 genes as TARDBP targets
  TARDBPtargetsH <- names(TARDBPnormalizedRatio)[1:1000]

  # From human to mouse
  ensemblM = useDataset("mmusculus_gene_ensembl",mart=ensembl)
  TARDBPtargetsH2M <- getLDS(attributes = c("ensembl_gene_id"),
							filters = "ensembl_gene_id",
							values = TARDBPtargetsH,
							mart = ensemblH,
							attributesL = c("ensembl_gene_id"),
							martL = ensemblM)
  TARDBPtargetsM <- unique(TARDBPtargetsH2M[,2])

  ## TARDBP targets expression analysis
  ensemblGenes2TxM <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id"), mart = ensemblM)
  tmp <- ensemblGenes2TxM$ensembl_transcript_id
  ensemblGenes2TxM <- ensemblGenes2TxM[,1]
  names(ensemblGenes2TxM) <- tmp

  DESeq2_normalized_counts_genes_WT <- DESeq2_normalized_counts[,"WT"]
  tmp <- ensemblGenes2TxM[names(DESeq2_normalized_counts_genes_WT)]
  names(DESeq2_normalized_counts_genes_WT) <- unname(tmp)
  DESeq2_normalized_counts_genes_WT <- DESeq2_normalized_counts_genes_WT[!is.na(names(DESeq2_normalized_counts_genes_WT))]
  DESeq2_normalized_counts_genes_WT <- split(DESeq2_normalized_counts_genes_WT,names(DESeq2_normalized_counts_genes_WT))
  DESeq2_normalized_counts_genes_WT <- sapply(DESeq2_normalized_counts_genes_WT,sum)

  DESeq2_normalized_counts_genes_KD <- DESeq2_normalized_counts[,"KD"]
  tmp <- ensemblGenes2TxM[names(DESeq2_normalized_counts_genes_KD)]
  names(DESeq2_normalized_counts_genes_KD) <- unname(tmp)
  DESeq2_normalized_counts_genes_KD <- DESeq2_normalized_counts_genes_KD[!is.na(names(DESeq2_normalized_counts_genes_KD))]
  DESeq2_normalized_counts_genes_KD <- split(DESeq2_normalized_counts_genes_KD,names(DESeq2_normalized_counts_genes_KD))
  DESeq2_normalized_counts_genes_KD <- sapply(DESeq2_normalized_counts_genes_KD,sum)

  pdf("TARDBPtargetsExpression.pdf",width=5,height=4)
  boxplot(list("TARBDP targets"=DESeq2_normalized_counts_genes_WT[names(DESeq2_normalized_counts_genes_WT)%in%TARDBPtargetsM]
              ,"Others"=DESeq2_normalized_counts_genes_WT[!(names(DESeq2_normalized_counts_genes_WT)%in%TARDBPtargetsM)])
		 ,outline=FALSE,ylab="Normalized gene counts",varwidth=TRUE)

  text(2,2500,paste0("Wilcoxon p-value: "
  					,signif(wilcox.test(DESeq2_normalized_counts_genes_WT[names(DESeq2_normalized_counts_genes_WT)%in%TARDBPtargetsM]
  									   ,DESeq2_normalized_counts_genes_WT[!(names(DESeq2_normalized_counts_genes_WT)%in%TARDBPtargetsM)])
  					$p.value,1)))
  dev.off()

  ## DEGs identification
  DESeq2_padj <- DESeq2_results$padj
  DESeq2_FC <- DESeq2_results$log2FoldChange
  names(DESeq2_padj) <-  names(DESeq2_FC) <- rownames(DESeq2_results)
  DESeq2_padj <- DESeq2_padj[is.finite(DESeq2_padj)]
  DESeq2_FC <- DESeq2_FC[names(DESeq2_padj)]

  DEGs <- names(DESeq2_padj)[DESeq2_padj<0.05]
  DEGsUP <- names(which(DESeq2_FC[DEGs]>0))
  DEGsDOWN <- names(which(DESeq2_FC[DEGs]<0))

  ## From DEG transcripts to DEG genes
  AllGenesDESeq2 <- unique(getBM(attributes = c("ensembl_transcript_id","ensembl_gene_id")
                               , filters = "ensembl_transcript_id"
                               , values = names(DESeq2_padj)
                               , mart = ensemblM))
  rownames(AllGenesDESeq2) <- AllGenesDESeq2[,"ensembl_transcript_id"]

  DESeq2_padj_genes <- DESeq2_padj
  names(DESeq2_padj_genes) <- AllGenesDESeq2[names(DESeq2_padj_genes),2]
  DESeq2_padj_genes <- DESeq2_padj_genes[!is.na(DESeq2_padj_genes)]
  DESeq2_padj_genes <- sapply(split(DESeq2_padj_genes,names(DESeq2_padj_genes)),min)

  AllGenesDESeq2 <- unique(AllGenesDESeq2[,2])

  DEGsGenes <- unique(getBM(attributes = c("ensembl_gene_id"), filters = "ensembl_transcript_id", values = DEGs, mart = ensemblM))[[1]]
  DEGsUPGenes <- unique(getBM(attributes = c("ensembl_gene_id"), filters = "ensembl_transcript_id", values = DEGsUP, mart = ensemblM))[[1]]
  DEGsDOWNGenes <- unique(getBM(attributes = c("ensembl_gene_id"), filters = "ensembl_transcript_id", values = DEGsDOWN, mart = ensemblM))[[1]]

  DEGsGenes <- DEGsGenes[!(DEGsGenes%in%DEGsUPGenes&DEGsGenes%in%DEGsDOWNGenes)]
  DEGsUPGenes <- DEGsUPGenes[DEGsUPGenes%in%DEGsGenes]
  DEGsDOWNGenes <- DEGsDOWNGenes[DEGsDOWNGenes%in%DEGsGenes]

  table(TARDBPtargetsM%in%AllGenesDESeq2)

  ## TARDBP targets in DESeq2 database
  TARDBPtargetsDESeq2 <- intersect(TARDBPtargetsM,AllGenesDESeq2)

  ## DEG TARDBP targets
  table(TARDBPtargetsDESeq2%in%DEGsGenes)
  table(TARDBPtargetsDESeq2%in%DEGsUPGenes)
  table(TARDBPtargetsDESeq2%in%DEGsDOWNGenes)

  ## FCs
  DESeq2_FC_UP <- DESeq2_FC[DEGsUP]
  DESeq2_FC_DOWN <- DESeq2_FC[DEGsDOWN]

  DESeq2_FC_UP_genes <- unique(getBM(attributes = c("ensembl_transcript_id","ensembl_gene_id")
                                   , filters = "ensembl_transcript_id"
                                   , values = names(DESeq2_FC_UP)
                                   , mart = ensemblM))
  rownames(DESeq2_FC_UP_genes) <- DESeq2_FC_UP_genes[,"ensembl_transcript_id"]
  DESeq2_FC_UP <- DESeq2_FC_UP[rownames(DESeq2_FC_UP_genes)]
  names(DESeq2_FC_UP) <- DESeq2_FC_UP_genes[,2]
  DESeq2_FC_UP <- sapply(split(DESeq2_FC_UP,names(DESeq2_FC_UP)),mean)
  DESeq2_FC_UP_TARDBP <- DESeq2_FC_UP[TARDBPtargetsDESeq2[TARDBPtargetsDESeq2%in%DEGsUPGenes]]

  DESeq2_FC_DOWN_genes <- unique(getBM(attributes = c("ensembl_transcript_id","ensembl_gene_id")
                                   , filters = "ensembl_transcript_id"
                                   , values = names(DESeq2_FC_DOWN)
                                   , mart = ensemblM))
  rownames(DESeq2_FC_DOWN_genes) <- DESeq2_FC_DOWN_genes[,"ensembl_transcript_id"]
  DESeq2_FC_DOWN <- DESeq2_FC_DOWN[rownames(DESeq2_FC_DOWN_genes)]
  names(DESeq2_FC_DOWN) <- DESeq2_FC_DOWN_genes[,2]
  DESeq2_FC_DOWN <- sapply(split(DESeq2_FC_DOWN,names(DESeq2_FC_DOWN)),mean)
  DESeq2_FC_DOWN_TARDBP <- DESeq2_FC_DOWN[TARDBPtargetsDESeq2[TARDBPtargetsDESeq2%in%DEGsDOWNGenes]]

  # Enrichments
  fisher.test(table(AllGenesDESeq2%in%DEGsGenes,AllGenesDESeq2%in%TARDBPtargetsDESeq2))$p.value
  fisher.test(table(AllGenesDESeq2%in%DEGsUPGenes,AllGenesDESeq2%in%TARDBPtargetsDESeq2))$p.value
  fisher.test(table(AllGenesDESeq2%in%DEGsDOWNGenes,AllGenesDESeq2%in%TARDBPtargetsDESeq2))$p.value

  ## Nanocompore profiled genes
  nanocomporeResultsPvalue <- unlist(nanocomporeResults[,9])
  names(nanocomporeResultsPvalue) <- as.character(unlist(nanocomporeResults[,"GeneId"]))
  nanocomporeResultsPvalue <- nanocomporeResultsPvalue[is.finite(nanocomporeResultsPvalue)]
  nanocomporeResultsPvalue <- split(nanocomporeResultsPvalue,names(nanocomporeResultsPvalue))

  nanocomporeResultsNSites <- sapply(nanocomporeResultsPvalue,function(i)sum(i<0.05,na.rm=TRUE))

  nanocomporeResultsPvalue <- sapply(nanocomporeResultsPvalue,min,na.rm=TRUE)

  AllGenesDESeq2Nanocompore <- intersect(AllGenesDESeq2, names(nanocomporeResultsPvalue))
  DEGsGenesNanocompore <- intersect(DEGsGenes, names(nanocomporeResultsPvalue))
  DEGsUPGenesNanocompore <- intersect(DEGsUPGenes, names(nanocomporeResultsPvalue))
  DEGsDOWNGenesNanocompore <- intersect(DEGsDOWNGenes, names(nanocomporeResultsPvalue))

  nanocomporeResultsPvalue <- nanocomporeResultsPvalue[AllGenesDESeq2Nanocompore]

  ## Annotations
  DESeq2_FC_UP_TARDBP <- sort(DESeq2_FC_UP_TARDBP)
  DESeq2_FC_DOWN_TARDBP <- sort(DESeq2_FC_DOWN_TARDBP,decreasing=TRUE)

  DESeq2_PVALUE_UP_TARDBP <- nanocomporeResultsPvalue[names(DESeq2_FC_UP_TARDBP)]
  names(DESeq2_PVALUE_UP_TARDBP) <- names(DESeq2_FC_UP_TARDBP)

  DESeq2_PVALUE_DOWN_TARDBP <- nanocomporeResultsPvalue[names(DESeq2_FC_DOWN_TARDBP)]
  names(DESeq2_PVALUE_DOWN_TARDBP) <- names(DESeq2_FC_DOWN_TARDBP)
  
  DESeq2_N_UP_TARDBP <- nanocomporeResultsNSites[names(DESeq2_FC_UP_TARDBP)]
  names(DESeq2_N_UP_TARDBP) <- names(DESeq2_FC_UP_TARDBP)

  DESeq2_N_DOWN_TARDBP <- nanocomporeResultsNSites[names(DESeq2_FC_DOWN_TARDBP)]
  names(DESeq2_N_DOWN_TARDBP) <- names(DESeq2_FC_DOWN_TARDBP)

  DESeq2_WTCOUNTS_UP_TARDBP <- DESeq2_normalized_counts_genes_WT[names(DESeq2_FC_UP_TARDBP)]
  DESeq2_WTCOUNTS_DOWN_TARDBP <- DESeq2_normalized_counts_genes_WT[names(DESeq2_FC_DOWN_TARDBP)]
  
  DESeq2_KDCOUNTS_UP_TARDBP <- DESeq2_normalized_counts_genes_KD[names(DESeq2_FC_UP_TARDBP)]
  DESeq2_KDCOUNTS_DOWN_TARDBP <- DESeq2_normalized_counts_genes_KD[names(DESeq2_FC_DOWN_TARDBP)]

  ensembl2symbolUP <- unique(getBM(attributes = c("ensembl_gene_id","external_gene_name")
                                   , filters = "ensembl_gene_id"
                                   , values = names(DESeq2_FC_UP_TARDBP)
                                   , mart = ensemblM)) 
  rownames(ensembl2symbolUP) <- ensembl2symbolUP[,1]
  ensembl2symbolUP <- ensembl2symbolUP[names(DESeq2_FC_UP_TARDBP),]

  names(DESeq2_FC_UP_TARDBP) <- names(DESeq2_PVALUE_UP_TARDBP) <- names(DESeq2_N_UP_TARDBP) <- 
  names(DESeq2_WTCOUNTS_UP_TARDBP) <- names(DESeq2_KDCOUNTS_UP_TARDBP) <- ensembl2symbolUP[,2]

  ensembl2symbolDOWN <- unique(getBM(attributes = c("ensembl_gene_id","external_gene_name")
                                   , filters = "ensembl_gene_id"
                                   , values = names(DESeq2_FC_DOWN_TARDBP)
                                   , mart = ensemblM)) 
  rownames(ensembl2symbolDOWN) <- ensembl2symbolDOWN[,1]
  ensembl2symbolDOWN <- ensembl2symbolDOWN[names(DESeq2_FC_DOWN_TARDBP),]

  names(DESeq2_FC_DOWN_TARDBP) <- names(DESeq2_PVALUE_DOWN_TARDBP) <- names(DESeq2_N_DOWN_TARDBP) <- 
  names(DESeq2_WTCOUNTS_DOWN_TARDBP) <- names(DESeq2_KDCOUNTS_DOWN_TARDBP) <- ensembl2symbolDOWN[,2]

  annotationColors = list(m6A_min_pValue=c("(0.05, Inf]"="lightgrey","(1e-05,0.05]"="sienna2","(1e-10,1e-05]"="sienna3","(-Inf,1e-10]"="sienna4")
                         ,m6A_Sites=c("(-Inf,0]"="lightgrey","(0,1]"="skyblue2","(1,2]"="skyblue3","(2, Inf]"="skyblue4"))

  annotationUP <- data.frame(normCounts=round(log10((DESeq2_WTCOUNTS_UP_TARDBP+DESeq2_KDCOUNTS_UP_TARDBP)/2))
                            ,m6A_min_pValue=cut(DESeq2_PVALUE_UP_TARDBP,c(-Inf,1e-10,1e-5,0.05,Inf))
                            ,m6A_Sites=cut(DESeq2_N_UP_TARDBP,c(-Inf,0,1,2,Inf)))

  annotationDOWN <- data.frame(normCounts=round(log10((DESeq2_WTCOUNTS_DOWN_TARDBP+DESeq2_KDCOUNTS_DOWN_TARDBP)/2))
                            ,m6A_min_pValue=cut(DESeq2_PVALUE_DOWN_TARDBP,c(-Inf,1e-10,1e-5,0.05,Inf))
                            ,m6A_Sites=cut(DESeq2_N_DOWN_TARDBP,c(-Inf,0,1,2,Inf)))

  pheatmap(data.frame(c(DESeq2_FC_DOWN_TARDBP,DESeq2_FC_UP_TARDBP))
          ,cluster_rows=FALSE
          ,cluster_cols=FALSE
          ,annotation_row=rbind(annotationDOWN,annotationUP)
          ,breaks=seq(-1.5,1.5,length.out=25)
          ,col=colorRampPalette(c("blue","white","red"))(26)
          ,annotation_colors=annotationColors
          ,show_colnames=FALSE
          ,filename="TARDBP_DEGS.pdf"
          ,width=4
          ,height=16)

  # m6A DEGs overlaps
  fisher.test(table(AllGenesDESeq2Nanocompore%in%DEGsGenesNanocompore
  				   ,AllGenesDESeq2Nanocompore%in%names(which(nanocomporeResultsPvalue<0.05))))
  fisher.test(table(AllGenesDESeq2Nanocompore%in%DEGsUPGenesNanocompore
				   ,AllGenesDESeq2Nanocompore%in%names(which(nanocomporeResultsPvalue<0.05))))
  fisher.test(table(AllGenesDESeq2Nanocompore%in%DEGsDOWNGenesNanocompore
  				   ,AllGenesDESeq2Nanocompore%in%names(which(nanocomporeResultsPvalue<0.05))))

  TARDBPtargetsDESeq2Nanocompore <- intersect(TARDBPtargetsDESeq2,AllGenesDESeq2Nanocompore)

  # m6A TARDBP targets overlaps
  fisher.test(table(AllGenesDESeq2Nanocompore%in%names(which(nanocomporeResultsPvalue<0.05))
				   ,AllGenesDESeq2Nanocompore%in%TARDBPtargetsDESeq2Nanocompore))

  # m6A DEGs overlaps for TARDBP targets
  fisher.test(table(TARDBPtargetsDESeq2Nanocompore%in%names(which(nanocomporeResultsPvalue<0.05))
				   ,TARDBPtargetsDESeq2Nanocompore%in%DEGsGenesNanocompore))
  fisher.test(table(TARDBPtargetsDESeq2Nanocompore%in%names(which(nanocomporeResultsPvalue<0.05))
				   ,TARDBPtargetsDESeq2Nanocompore%in%DEGsUPGenesNanocompore))
  fisher.test(table(TARDBPtargetsDESeq2Nanocompore%in%names(which(nanocomporeResultsPvalue<0.05))
				   ,TARDBPtargetsDESeq2Nanocompore%in%DEGsDOWNGenesNanocompore))

  # Number of m6A sites for TARDBP targets
  nanocomporeResultsPvalue <- unlist(nanocomporeResults[,9])
  names(nanocomporeResultsPvalue) <- as.character(unlist(nanocomporeResults[,"GeneId"]))
  nanocomporeResultsPvalue <- nanocomporeResultsPvalue[is.finite(nanocomporeResultsPvalue)]
  nanocomporeResultsPvalue <- split(nanocomporeResultsPvalue,names(nanocomporeResultsPvalue))

  boxplot(list("TARDBP_targets"=sapply(nanocomporeResultsPvalue[TARDBPtargetsDESeq2Nanocompore],function(i)sum(i<0.05))
              ,"Other_genes"=sapply(nanocomporeResultsPvalue[setdiff(names(nanocomporeResultsPvalue),TARDBPtargetsDESeq2Nanocompore)],function(i)sum(i<0.05))
         ),ylab="# m6A sites",main="",outline=FALSE)

  wilcox.test(sapply(nanocomporeResultsPvalue[TARDBPtargetsDESeq2Nanocompore],function(i)sum(i<0.05))
             ,sapply(nanocomporeResultsPvalue[setdiff(names(nanocomporeResultsPvalue),TARDBPtargetsDESeq2Nanocompore)],function(i)sum(i<0.05)))

### General statistics
  ## All k-mers
  dim(nanocomporeResults)
  # [1] 4866607      24

  ## p-values per gene
  pValues <- as.numeric(nanocomporeResults$GMM_logit_pvalue)
  names(pValues) <- nanocomporeResults$GeneId

  ## K-mers with finite p-value
  table(is.finite(pValues))
  #  FALSE    TRUE 
  # 655849 4210758 

  ## K-mers with significant p-value
  table(is.finite(pValues)&pValues<=0.05)
  #   FALSE    TRUE 
  # 4861121    5486 

  ## Profiled genes
  length(as.numeric(table(as.character(names(pValues)))))
  # 4037

  ## m6A+ genes
  length(as.numeric(table(as.character(names(pValues[pValues<=0.05])))))
  # 1777

  ## Significant sites per m6A+ gene
  summary(as.numeric(table(as.character(names(pValues[pValues<=0.05])))))
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 1.00    1.00    2.00    3.08    4.00   49.00 

  ## Enrichment in TARDBP targets
  fisher.test(table("m6A+"=unique(names(pValues))%in%names(table(as.character(names(pValues[pValues<=0.05])))),"TARDBP"=unique(names(pValues))%in%TARDBPtargetsM))

  ## Stoichiometry for TARDBP
  x <- (unlist(nanocomporeResults[,"GMM_logit_pvalue"])<=0.05)
  y <- grepl("ENSMUSG00000041459",unlist(nanocomporeResults[,"GeneId"]))
  TARDBPcounts <- unlist(nanocomporeResults[which(x&y),"cluster_counts"])

  TARDBPsto <- unname(sapply(TARDBPcounts,function(i)
  {
    split1 <- strsplit(i,"_")[[1]]
    split2 <- split1[grep("[:]",split1)]
    split3 <- sapply(split2,function(j)strsplit(j,"[:]")[[1]][[2]])
    split4 <- sapply(split3,function(j)strsplit(j,"/")[[1]])

    out <- type.convert(rbind(split4,apply(split4,2,function(z)sum(as.numeric(z)))),"numeric")
    out <- cbind(apply(out[,1:3],1,sum),apply(out[,4:6],1,sum))
    rownames(out) = c("1","2","total")
    colnames(out) = c("WT","KD")

    # Selection of m6A cluster    
    WT <- out[1:2,1]/out[3,1]
    KD <- out[1:2,2]/out[3,2]

    out[which((WT-KD)>0),]/out[3,]
  }))

  TARDBPstoSingleSample <- unname(sapply(TARDBPcounts,function(i)
  {
    split1 <- strsplit(i,"_")[[1]]
    split2 <- split1[grep("[:]",split1)]
    split3 <- sapply(split2,function(j)strsplit(j,"[:]")[[1]][[2]])
    split4 <- sapply(split3,function(j)strsplit(j,"/")[[1]])

    out <- type.convert(rbind(split4,apply(split4,2,function(z)sum(as.numeric(z)))),"numeric")
    outSplit <- out


    out <- cbind(apply(out[,1:3],1,sum),apply(out[,4:6],1,sum))
    rownames(out) = c("1","2","total")
    colnames(out) = c("WT","KD")

    # Selection of m6A cluster    
    WT <- out[1:2,1]/out[3,1]
    KD <- out[1:2,2]/out[3,2]

    outSplit[which((WT-KD)>0),]/outSplit[3,]
  }))

  TARDBPstom6ASites <- TARDBPsto[,c(1,3,4,5,7)]
  colnames(TARDBPstom6ASites) <- paste0("site",1:5)
  
  TARDBPstoSingleSample <- TARDBPstoSingleSample[,c(1,3,4,5,7)]
  colnames(TARDBPstoSingleSample) <- paste0("site",1:5)

  pdf("stoichiometryTARDBP5-mers.pdf", width=6, height=8)
  barplot(TARDBPsto,ylim=c(0,1),beside=T,col=c("red","blue"),ylab="m6A stoichiometry",main="TARDBP",xlab="",las=2,angle=30)
  legend("topright",col=c("red","blue"),pch=15,legend=c("WT","Mettl3-KO"))
  dev.off()

  pdf("stoichiometryTARDBPSites.pdf", width=6, height=4)
  barplot(TARDBPstom6ASites,ylim=c(0,1),beside=T,col=c("red","blue"),ylab="m6A stoichiometry",main="TARDBP",xlab="",las=2)
  legend("topright",col=c("red","blue"),pch=15,legend=c("WT","Mettl3-KO"))
  dev.off()

  pdf("stoichiometryTARDBPSitesWT.pdf", width=3, height=4)
  barplot(TARDBPstom6ASites[1,],ylim=c(0,1),beside=T,col="red",ylab="m6A stoichiometry",main="TARDBP",xlab="",las=2)
  dev.off()

  ## Stoichiometry - on probability threshold 0.5
  pdf("stoichiometryTARDBPSitesWTReplicates.pdf", width=3, height=4)
  ytmp <- TARDBPstom6ASites[1,]
  xtmp <- barplot(ytmp,las=2,ylim=c(0,1),col="red",ylab="TARDBP m6A stoichiometry")

  matplot(xtmp,t(TARDBPstoSingleSample[1:3,]),pch=16,col=1,add=TRUE)
  legend("top",pch=c(16,15),col=c(1,"red"),legend=c("Single sample", "All samples"),bty="n")
  dev.off()


  
  # Number of sites - LESS THAN 8% OF THE GENES HAVE MORE m6A SITES THAN TARDBP
  GenepValues <- unlist(nanocomporeResults[,9])
  names(GenepValues) <- as.character(unlist(nanocomporeResults[,21]))
  GenepValues[!is.finite(GenepValues)] <- 1

  Genem6ASites <- sort(table(as.character(names(GenepValues[GenepValues<=0.05]))),decreasing=TRUE)

  m6Aquantile <- round((table(Genem6ASites>Genem6ASites["ENSMUSG00000041459"])/length(Genem6ASites))["TRUE"],2)

  # Significance - LESS THAN 7% OF THE GENES HAVE A LOWER MEDIAN SIGNIFICANCE THAN TARDBP
  nanocomporeSignificantResults <- nanocomporeResults[which(nanocomporeResults[,9]<=0.05),]
  GenemedianPvalue <- sort(sapply(split(nanocomporeSignificantResults[,9],nanocomporeSignificantResults[,21]),function(i)median(unlist(i))),decreasing=FALSE)
  medianPvaluequantile <- round((table(GenemedianPvalue<GenemedianPvalue["ENSMUSG00000041459"])/length(GenemedianPvalue))["TRUE"],2)

  # Joint analysis - ONLY 8 GENES HAVE SIMILAR MEDIAN SIGNIFICANCE AND NUMBER OF m6A SITES
  x <- log10(GenemedianPvalue)
  y <- as.numeric(Genem6ASites[names(GenemedianPvalue)])
  names(y) <- names(GenemedianPvalue)

  x1 <- x[x>=(-25)&y<=15]
  y1 <- y[names(x1)]

  xth <- 0.95*x1["ENSMUSG00000041459"]
  yth <- 0.95*y1["ENSMUSG00000041459"]

  pdf("m6AStatistics.pdf",width=8,height=5)
  plot(x1,y1,pch=20,xlab="Log10 median p-value",ylab="Number of m6A 5-mers")

  abline(v=xth,col="lightgrey",lty=2,lwd=2)
  abline(h=yth,col="lightgrey",lty=2,lwd=2)

  x2 <- x[x<xth&y>yth]
  y2 <- y[names(x2)]

  points(x2,y2,pch=16,col=2)

  topGenesSymbols <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                                  filter = "ensembl_gene_id",
                                  values = names(x2),
                                  mart = ensemblM)

  rownames(topGenesSymbols) <- topGenesSymbols[,1]
  topGenesSymbols <- topGenesSymbols[names(x2),]
  names(x2) <- names(y2) <- topGenesSymbols[,2]

  x2 <- x2[sort(names(x2),decreasing=TRUE)]
  y2 <- y2[names(x2)]

  text(x2+1*(3**(seq_along(x2)%%2)-2),y2+0.4,names(x2))

  x3 <- x[x<(-25)|y>15]
  y3 <- y[names(x3)]

  x3[x3<(-25)] <- (-25)
  y3[y3>(15)] <- 15

  points(x3,y3,pch=17,col=3)

  text(-24,15,paste0("N=",length(x2)))
  text(-2,15,paste0("N=",length(x[x>=xth&y>yth])))
  text(-24,5.5,paste0("N=",length(x[x<xth&y<=yth])))
  text(-2,5.5,paste0("N=",length(x[x>=xth&y<=yth])))
  dev.off()

  ## TARDBP sites probabilities heatmap
  # Loading files from singleMolecule.R
  TARDBP_sites <- readRDS("TARDBP_sitesData.rds")
  TARDBP_proba <- readRDS("TARDBP_sitesProbabilities.rds")

  WT_reads <- unique(unlist(sapply(TARDBP_sites,function(i)
  {
  	i[grep("WT",i[,2]),1]
  })))
  KD_reads <- unique(unlist(sapply(TARDBP_sites,function(i)
  {
  	i[grep("KD",i[,2]),1]
  })))

  sitesNames <- c("148617470"="site1","148617386"="site2","148617261"="site3","148617241"="site4","148617069"="site5")

  names(TARDBP_sites) <- sapply(strsplit(names(TARDBP_sites),"_"),"[[",2)
  names(TARDBP_sites) <- sitesNames[names(TARDBP_sites)]

  colnames(TARDBP_proba) <- sapply(strsplit(colnames(TARDBP_proba),"_"),"[[",2)
  colnames(TARDBP_proba) <- sitesNames[colnames(TARDBP_proba)]
  
  TARDBP_proba <- t(TARDBP_proba)

  pheatmap(TARDBP_proba[,colnames(TARDBP_proba)%in%WT_reads],color=colorRampPalette(c("blue","white","red"))(length(seq(0,1,by=0.01)))
  		  ,breaks=seq(0,1,by=0.01),cluster_rows=FALSE,show_colnames=FALSE,treeheight_col=0
  		  ,main="Modification probability",filename="WT_proba.pdf",width=6,height=3)

  pheatmap(TARDBP_proba[,colnames(TARDBP_proba)%in%KD_reads],color=colorRampPalette(c("blue","white","red"))(length(seq(0,1,by=0.01)))
  		  ,breaks=seq(0,1,by=0.01),cluster_rows=FALSE,show_colnames=FALSE,treeheight_col=0
  		  ,main="Modification probability",filename="KD_proba.pdf",width=6,height=3)

  ## Residual signal
  round(sum(TARDBP_proba[,colnames(TARDBP_proba)%in%WT_reads]>0.5)/length(TARDBP_proba[,colnames(TARDBP_proba)%in%WT_reads]),2) # 0.53
  round(sum(TARDBP_proba[,colnames(TARDBP_proba)%in%KD_reads]>0.5)/length(TARDBP_proba[,colnames(TARDBP_proba)%in%KD_reads]),2) # 0.24

  round(apply(TARDBP_proba[,colnames(TARDBP_proba)%in%WT_reads]>0.5,1,sum)/ncol(TARDBP_proba[,colnames(TARDBP_proba)%in%WT_reads]>0.5),2)
  # site1 site2 site3 site4 site5 
  #  0.61  0.70  0.38  0.44  0.54 

  ## Heatmap per sample
  WT1_reads <- unique(unlist(sapply(TARDBP_sites,function(i){i[grep("WT1",i[,2]),1]})))
  WT2_reads <- unique(unlist(sapply(TARDBP_sites,function(i){i[grep("WT2",i[,2]),1]})))
  WT3_reads <- unique(unlist(sapply(TARDBP_sites,function(i){i[grep("WT3",i[,2]),1]})))

  KD1_reads <- unique(unlist(sapply(TARDBP_sites,function(i){i[grep("KD1",i[,2]),1]})))
  KD2_reads <- unique(unlist(sapply(TARDBP_sites,function(i){i[grep("KD2",i[,2]),1]})))
  KD3_reads <- unique(unlist(sapply(TARDBP_sites,function(i){i[grep("KD3",i[,2]),1]})))

  WT1_reads_ord <- sort(apply(TARDBP_proba[,colnames(TARDBP_proba)%in%WT1_reads],2,median),index.return = TRUE)[[2]]
  WT2_reads_ord <- sort(apply(TARDBP_proba[,colnames(TARDBP_proba)%in%WT2_reads],2,median),index.return = TRUE)[[2]]
  WT3_reads_ord <- sort(apply(TARDBP_proba[,colnames(TARDBP_proba)%in%WT3_reads],2,median),index.return = TRUE)[[2]]

  KD1_reads_ord <- sort(apply(TARDBP_proba[,colnames(TARDBP_proba)%in%KD1_reads],2,median),index.return = TRUE)[[2]]
  KD2_reads_ord <- sort(apply(TARDBP_proba[,colnames(TARDBP_proba)%in%KD2_reads],2,median),index.return = TRUE)[[2]]
  KD3_reads_ord <- sort(apply(TARDBP_proba[,colnames(TARDBP_proba)%in%KD3_reads],2,median),index.return = TRUE)[[2]]

  pheatmap(TARDBP_proba[,colnames(TARDBP_proba)%in%WT1_reads][,WT1_reads_ord],color=colorRampPalette(c("blue","white","red"))(length(seq(0,1,by=0.01)))
        ,breaks=seq(0,1,by=0.01),cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=FALSE,treeheight_col=0
        ,main="Modification probability",filename="WT1_proba.pdf",width=2,height=3)

  pheatmap(TARDBP_proba[,colnames(TARDBP_proba)%in%KD1_reads][,KD1_reads_ord],color=colorRampPalette(c("blue","white","red"))(length(seq(0,1,by=0.01)))
        ,breaks=seq(0,1,by=0.01),cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=FALSE,treeheight_col=0
        ,main="Modification probability",filename="KD1_proba.pdf",width=2,height=3)

  pheatmap(TARDBP_proba[,colnames(TARDBP_proba)%in%WT2_reads][,WT2_reads_ord],color=colorRampPalette(c("blue","white","red"))(length(seq(0,1,by=0.01)))
        ,breaks=seq(0,1,by=0.01),cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=FALSE,treeheight_col=0
        ,main="Modification probability",filename="WT2_proba.pdf",width=2,height=3)

  pheatmap(TARDBP_proba[,colnames(TARDBP_proba)%in%KD2_reads][,KD2_reads_ord],color=colorRampPalette(c("blue","white","red"))(length(seq(0,1,by=0.01)))
        ,breaks=seq(0,1,by=0.01),cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=FALSE,treeheight_col=0
        ,main="Modification probability",filename="KD2_proba.pdf",width=2,height=3)

  pheatmap(TARDBP_proba[,colnames(TARDBP_proba)%in%WT3_reads][,WT3_reads_ord],color=colorRampPalette(c("blue","white","red"))(length(seq(0,1,by=0.01)))
        ,breaks=seq(0,1,by=0.01),cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=FALSE,treeheight_col=0
        ,main="Modification probability",filename="WT3_proba.pdf",width=2,height=3)

  pheatmap(TARDBP_proba[,colnames(TARDBP_proba)%in%KD3_reads][,KD3_reads_ord],color=colorRampPalette(c("blue","white","red"))(length(seq(0,1,by=0.01)))
        ,breaks=seq(0,1,by=0.01),cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=FALSE,treeheight_col=0
        ,main="Modification probability",filename="KD3_proba.pdf",width=2,height=3)

  ## Residual signal
  round(sum(TARDBP_proba[,colnames(TARDBP_proba)%in%WT1_reads]>0.5)/length(TARDBP_proba[,colnames(TARDBP_proba)%in%WT1_reads]),2) # 0.58
  round(sum(TARDBP_proba[,colnames(TARDBP_proba)%in%KD1_reads]>0.5)/length(TARDBP_proba[,colnames(TARDBP_proba)%in%KD1_reads]),2) # 0.23

  round(sum(TARDBP_proba[,colnames(TARDBP_proba)%in%WT2_reads]>0.5)/length(TARDBP_proba[,colnames(TARDBP_proba)%in%WT2_reads]),2) # 0.51
  round(sum(TARDBP_proba[,colnames(TARDBP_proba)%in%KD2_reads]>0.5)/length(TARDBP_proba[,colnames(TARDBP_proba)%in%KD2_reads]),2) # 0.24

  round(sum(TARDBP_proba[,colnames(TARDBP_proba)%in%WT3_reads]>0.5)/length(TARDBP_proba[,colnames(TARDBP_proba)%in%WT3_reads]),2) # 0.51
  round(sum(TARDBP_proba[,colnames(TARDBP_proba)%in%KD3_reads]>0.5)/length(TARDBP_proba[,colnames(TARDBP_proba)%in%KD3_reads]),2) # 0.24

  ## TARDBP sites co-occurrence
  # Classification matrix
  x <- (TARDBP_proba[,colnames(TARDBP_proba)%in%WT_reads])>0.75
  xpos <- c(1.5,4.5,7.5,10.5,13.5,16.5,19.5,22.5,25.5,28.5)

  # Second order interactions
  cooccurrenceStatistics2 <- apply(combn(rownames(x),m=2),2,function(comb)
  {
  	i <- comb[[1]]
  	j <- comb[[2]]
    
    y1 <- x[i,]
    y2 <- x[j,]
    # Experimental co-occurrences
    cooccurTmp <- table(y1&y2)["TRUE"]

    # Null model with 1k random configurations
    nullModel <- sapply(1:1000,function(k)
    {
    	set.seed(k)
    	y1R <- sample(y1,length(y1))
    	y2R <- sample(y2,length(y2))
    	table(y1R&y2R)["TRUE"]
    })
    
    c(cooccurTmp,nullModel)
  })

  rownames(cooccurrenceStatistics2) <- NULL
  colnames(cooccurrenceStatistics2) <- apply(combn(rownames(x),m=2),2,paste0,collapse="_")

  # Distribution momenta
  meanTmp2 <- apply(cooccurrenceStatistics2[-1,],2,mean)
  sdTmp2 <- apply(cooccurrenceStatistics2[-1,],2,sd)

  ylim2 <- (max(meanTmp2,cooccurrenceStatistics2[1,]))
  ylim2 <- ylim2*1.1

  pdf("cooccurrenceBarplots2.pdf",width=6,height=8)
  par(mfrow=c(2,1))
  barplot(rbind(meanTmp2,cooccurrenceStatistics2[1,]),las=2,beside=TRUE,col=c("red","blue"),main="Co-occurrences of\ndegree 2",ylab="Co-occurrences",ylim=c(0,ylim2))
  
  for(i in seq_along(meanTmp2))
  {
	arrows(x0=xpos[[i]],y0=meanTmp2[[i]],x1=xpos[[i]],y1=(meanTmp2[[i]]+sdTmp2[[i]]),angle=90,length=0.05)
	arrows(x0=xpos[[i]],y0=meanTmp2[[i]],x1=xpos[[i]],y1=(meanTmp2[[i]]-sdTmp2[[i]]),angle=90,length=0.05)

	text((xpos[[i]]+0.5),max(cooccurrenceStatistics2[1,i],(meanTmp2[[i]]+sdTmp2[[i]]))+(0.05*ylim2)
		,format(round(table(cooccurrenceStatistics2[-1,i]>cooccurrenceStatistics2[1,i])["TRUE"]/length(cooccurrenceStatistics2[-1,i]),2), nsmall = 2))
  }

  # legend("topright",col=c("red","blue"),pch=15,legend=c("Null model","Experimental"))
  dev.off()

  # Third order interactions
  cooccurrenceStatistics3 <- apply(combn(rownames(x),m=3),2,function(comb)
  {
  	i <- comb[[1]]
  	j <- comb[[2]]
  	l <- comb[[3]]
    
    y1 <- x[i,]
    y2 <- x[j,]
    y3 <- x[l,]
    # Experimental co-occurrences
    cooccurTmp <- table(y1&y2&y3)["TRUE"]

    # Null model with 1k random configurations
    nullModel <- sapply(1:1000,function(k)
    {
    	set.seed(k)
    	y1R <- sample(y1,length(y1))
    	y2R <- sample(y2,length(y2))
    	y3R <- sample(y3,length(y3))
    	table(y1R&y2R&y3R)["TRUE"]
    })
    
    c(cooccurTmp,nullModel)
  })

  rownames(cooccurrenceStatistics3) <- NULL
  colnames(cooccurrenceStatistics3) <- apply(combn(rownames(x),m=3),2,paste0,collapse="_")

  # Distribution momenta
  meanTmp3 <- apply(cooccurrenceStatistics3[-1,],2,mean)
  sdTmp3 <- apply(cooccurrenceStatistics3[-1,],2,sd)

  ylim3 <- (max(meanTmp3,cooccurrenceStatistics3[1,]))
  ylim3 <- ylim3*1.1

  pdf("cooccurrenceBarplots3.pdf",width=6,height=8)
  par(mfrow=c(2,1))
  barplot(rbind(meanTmp3,cooccurrenceStatistics3[1,]),las=2,beside=TRUE,col=c("red","blue"),main="Co-occurrences of\ndegree 3",ylab="Co-occurrences",ylim=c(0,(max(meanTmp3,cooccurrenceStatistics3[1,])+10)))
  
  for(i in seq_along(meanTmp3))
  {
	arrows(x0=xpos[[i]],y0=meanTmp3[[i]],x1=xpos[[i]],y1=(meanTmp3[[i]]+sdTmp3[[i]]),angle=90,length=0.05)
	arrows(x0=xpos[[i]],y0=meanTmp3[[i]],x1=xpos[[i]],y1=(meanTmp3[[i]]-sdTmp3[[i]]),angle=90,length=0.05)
	text((xpos[[i]]+0.5),max(cooccurrenceStatistics3[1,i],(meanTmp3[[i]]+sdTmp3[[i]]))+(0.05*ylim3)
		,format(round(table(cooccurrenceStatistics3[-1,i]>cooccurrenceStatistics3[1,i])["TRUE"]/length(cooccurrenceStatistics3[-1,i]),2), nsmall = 2))

  }

  # legend("topright",col=c("red","blue"),pch=15,legend=c("Null model","Experimental"))
  dev.off()

  # Fourth order interactions
  cooccurrenceStatistics4 <- apply(combn(rownames(x),m=4),2,function(comb)
  {
  	i <- comb[[1]]
  	j <- comb[[2]]
  	l <- comb[[3]]
  	p <- comb[[4]]
    
    y1 <- x[i,]
    y2 <- x[j,]
    y3 <- x[l,]
    y4 <- x[p,]
    # Experimental co-occurrences
    cooccurTmp <- table(y1&y2&y3&y4)["TRUE"]

    # Null model with 1k random configurations
    nullModel <- sapply(1:1000,function(k)
    {
    	set.seed(k)
    	y1R <- sample(y1,length(y1))
    	y2R <- sample(y2,length(y2))
    	y3R <- sample(y3,length(y3))
    	y4R <- sample(y4,length(y4))
    	table(y1R&y2R&y3R&y4R)["TRUE"]
    })
    
    c(cooccurTmp,nullModel)
  })

  rownames(cooccurrenceStatistics4) <- NULL
  colnames(cooccurrenceStatistics4) <- apply(combn(rownames(x),m=4),2,paste0,collapse="_")

  cooccurrenceStatistics4[!is.finite(cooccurrenceStatistics4)] <- 0

  # Distribution momenta
  meanTmp4 <- apply(cooccurrenceStatistics4[-1,],2,mean)
  sdTmp4 <- apply(cooccurrenceStatistics4[-1,],2,sd)

  ylim4 <- (max(meanTmp4,cooccurrenceStatistics4[1,]))
  ylim4 <- ylim4*1.1

  pdf("cooccurrenceBarplots4.pdf",width=6,height=8)
  par(mfrow=c(2,1))
  barplot(rbind(c(meanTmp4,rep(NaN,5)),c(cooccurrenceStatistics4[1,],rep(NaN,5))),las=2,beside=TRUE,col=c("red","blue"),main="Co-occurrences of\ndegree 4",ylab="Co-occurrences",ylim=c(0,(max(meanTmp4,cooccurrenceStatistics4[1,])+10)))
  
  for(i in seq_along(meanTmp4))
  {
	arrows(x0=xpos[[i]],y0=meanTmp4[[i]],x1=xpos[[i]],y1=(meanTmp4[[i]]+sdTmp4[[i]]),angle=90,length=0.05)
	arrows(x0=xpos[[i]],y0=meanTmp4[[i]],x1=xpos[[i]],y1=(meanTmp4[[i]]-sdTmp4[[i]]),angle=90,length=0.05)
	text((xpos[[i]]+0.5),max(cooccurrenceStatistics4[1,i],(meanTmp4[[i]]+sdTmp4[[i]]))+(0.05*ylim4)
		,format(round(table(cooccurrenceStatistics4[-1,i]>cooccurrenceStatistics4[1,i])["TRUE"]/length(cooccurrenceStatistics4[-1,i]),2), nsmall = 2))
  }

  # legend("topright",col=c("red","blue"),pch=15,legend=c("Null model","Experimental"))
  dev.off()

  # Fifth order interactions
  cooccurrenceStatistics5 <- apply(combn(rownames(x),m=5),2,function(comb)
  {
  	i <- comb[[1]]
  	j <- comb[[2]]
  	l <- comb[[3]]
  	p <- comb[[4]]
  	o <- comb[[5]]
    
    y1 <- x[i,]
    y2 <- x[j,]
    y3 <- x[l,]
    y4 <- x[p,]
    y5 <- x[o,]
    # Experimental co-occurrences
    cooccurTmp <- table(y1&y2&y3&y4&y5)["TRUE"]

    # Null model with 1k random configurations
    nullModel <- sapply(1:1000,function(k)
    {
    	set.seed(k)
    	y1R <- sample(y1,length(y1))
    	y2R <- sample(y2,length(y2))
    	y3R <- sample(y3,length(y3))
    	y4R <- sample(y4,length(y4))
    	y5R <- sample(y5,length(y5))
    	table(y1R&y2R&y3R&y4R&y5R)["TRUE"]
    })
    
    c(cooccurTmp,nullModel)
  })

  rownames(cooccurrenceStatistics5) <- NULL
  colnames(cooccurrenceStatistics5) <- apply(combn(rownames(x),m=5),2,paste0,collapse="_")

  cooccurrenceStatistics5[!is.finite(cooccurrenceStatistics5)] <- 0

  meanTmp5 <- mean(cooccurrenceStatistics5[-1,])
  sdTmp5 <- sd(cooccurrenceStatistics5[-1,])

  ylim5 <- (max(meanTmp5,cooccurrenceStatistics5[1,]))
  ylim5 <- ylim5*1.1

  pdf("cooccurrenceBarplots5.pdf",width=6,height=8)
  par(mfrow=c(2,1))
  barplot(rbind(c(meanTmp5,rep(NaN,9)),c(cooccurrenceStatistics5[1,],rep(NaN,9))),las=2,beside=TRUE,col=c("red","blue"),main="Co-occurrences of\ndegree 5",ylab="Co-occurrences",ylim=c(0,ylim5))
  
  for(i in seq_along(meanTmp5))
  {
	arrows(x0=xpos[[i]],y0=meanTmp5[[i]],x1=xpos[[i]],y1=(meanTmp5[[i]]+sdTmp5[[i]]),angle=90,length=0.05)
	arrows(x0=xpos[[i]],y0=meanTmp5[[i]],x1=xpos[[i]],y1=(meanTmp5[[i]]-sdTmp5[[i]]),angle=90,length=0.05)
	text((xpos[[i]]+0.5),max(cooccurrenceStatistics5[1,i],(meanTmp5[[i]]+sdTmp5[[i]]))+(0.05*ylim5)
		,format(round(table(cooccurrenceStatistics5[-1,i]>cooccurrenceStatistics5[1,i])["TRUE"]/length(cooccurrenceStatistics5[-1,i]),2), nsmall = 2))

  }

  # legend("topright",col=c("red","blue"),pch=15,legend=c("Null model","Experimental"))
  dev.off()

  ## Co-occurrences barplot
  input <- c(sort(cooccurrenceStatistics2[1,],decreasing=TRUE)
			,sort(cooccurrenceStatistics3[1,],decreasing=TRUE)
			,sort(cooccurrenceStatistics4[1,],decreasing=TRUE)
			,sort(cooccurrenceStatistics5[1,],decreasing=TRUE))
  names(input) <- gsub("_","&",names(input))

  probas <- c(sapply(seq_along(cooccurrenceStatistics2[1,]),function(i)table(cooccurrenceStatistics2[-1,i]>cooccurrenceStatistics2[1,i])["TRUE"]/length(cooccurrenceStatistics2[-1,i]))
  			, sapply(seq_along(cooccurrenceStatistics3[1,]),function(i)table(cooccurrenceStatistics3[-1,i]>cooccurrenceStatistics3[1,i])["TRUE"]/length(cooccurrenceStatistics3[-1,i]))
  			, sapply(seq_along(cooccurrenceStatistics4[1,]),function(i)table(cooccurrenceStatistics4[-1,i]>cooccurrenceStatistics4[1,i])["TRUE"]/length(cooccurrenceStatistics4[-1,i]))
  			, sapply(seq_along(cooccurrenceStatistics5[1,]),function(i)table(cooccurrenceStatistics5[-1,i]>cooccurrenceStatistics5[1,i])["TRUE"]/length(cooccurrenceStatistics5[-1,i])))

  names(probas) <- c(colnames(cooccurrenceStatistics2),colnames(cooccurrenceStatistics3),colnames(cooccurrenceStatistics4),colnames(cooccurrenceStatistics5))
  names(probas) <- gsub("_","&",names(probas))
  probas <- probas[names(input)]

  colorFunction <- function(p)
  {
  	pL <- p
  	colsTmp <- colorRampPalette(c("darkred","red","beige","white"))(21)
  	colsTmp[1+(round(pL,2)*20)]
  }

  pdf("cooccurrenceUpsetPlot.pdf",width=5,height=4)
  upset(fromExpression(input),keep.order=TRUE,sets=sort(colnames(fromExpression(input))))
  foe <- barplot(unname(input),xlab="",ylab="co-occurrences",col=colorFunction(probas),ylim=c(0,160))
  # text(foe[,1],input+25,round(probas,3),srt=90)
  barplot(as.matrix(rep(0.1,11),ncol=1),col=colorFunction(seq(0,1,by=0.1)),ylab="Co-occurence significance")
  dev.off()

### General plots
  ## Only sites k-mers
  nanocomporeResultsSignificant <- nanocomporeResults[is.finite(pValues)&pValues<=0.05,]

  ## Guitar plots - time demanding
  # Significant sites
  ncmp_results_sel <- nanocomporeResultsSignificant %>%
    mutate(start=genomicPos, end=genomicPos+5, name=ref_id, score=0) %>%
    dplyr::select(chr, start, end, name, score, strand)
  sites <- ncmp_results_sel %>% with(., GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), strand=strand))
  newStyle <- mapSeqlevels(seqlevels(sites), "UCSC")
  sites <- renameSeqlevels(sites, newStyle)
  diffSites <- list(sites)
  names(diffSites) <- c("Nanocompore peaks")

  pdf("nanocomporeProfileSignificantCI_GMM_logit_pvalue.pdf")
  GuitarPlot(stGRangeLists = diffSites, txTxdb=txdb, pltTxType=c("mrna"), headOrtail=FALSE, enableCI=TRUE)
  dev.off()

  # Random sites
  set.seed(1)
  ncmp_results_sel <- nanocomporeResults[sample(1:nrow(nanocomporeResults),1e5),] %>%
    mutate(start=genomicPos, end=genomicPos+5, name=ref_id, score=0) %>%
    dplyr::select(chr, start, end, name, score, strand)
  sites <- ncmp_results_sel %>% with(., GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), strand=strand))
  newStyle <- mapSeqlevels(seqlevels(sites), "UCSC")
  sites <- renameSeqlevels(sites, newStyle)
  diffSites <- list(sites)
  names(diffSites) <- c("Nanocompore peaks")

  pdf("nanocomporeProfileNotSignificantCI_GMM_logit_pvalue.pdf")
  GuitarPlot(stGRangeLists = diffSites, txTxdb=txdb, pltTxType=c("mrna"), headOrtail=FALSE, enableCI=TRUE)
  dev.off()

  ## A containing k-mers
  pdf("AContainingKmers_GMM_logit_pvalue.pdf",width=4,height=4)
  barplot(c(length(grep("A",nanocomporeResults$ref_kmer))/length(nanocomporeResults$ref_kmer)
           ,length(grep("A",nanocomporeResultsSignificant$ref_kmer))/length(nanocomporeResultsSignificant$ref_kmer))
         ,col=c("grey","black"),ylab="Frequency",xlab="",main="A containing k-mers",ylim=c(0,2))
  legend("top",col=c("grey","black","red"),lwd=c(NA,NA,2),pch=c(15,15,NA),legend=c("All","Significant","Null model"),bg="white",box.col="white",lty=c(NA,NA,2))
  abline(h=((4**5-3**5)/4**5),lwd=2,lty=2,col=2)
  dev.off()

  ## DRACH k-mers
  pdf("DrachKmers_GMM_logit_pvalue.pdf",width=4,height=4)
  barplot(c(tryCatch(table(nanocomporeResults$ref_kmer %in% drach)[["TRUE"]],error=function(e)0)/length(nanocomporeResults$ref_kmer)
           ,tryCatch(table(nanocomporeResultsSignificant$ref_kmer %in% drach)[["TRUE"]],error=function(e)0)/length(nanocomporeResultsSignificant$ref_kmer))
         ,col=c("grey","black"),ylab="Frequency",xlab="",main="DRACH k-mers",ylim=c(0,1))
  legend("top",col=c("grey","black","red"),lwd=c(NA,NA,2),pch=c(15,15,NA),legend=c("All","Significant","Null model"),bg="white",box.col="white",lty=c(NA,NA,2))
  abline(h=(length(drach)/4**5),lwd=2,lty=2,col=2)
  dev.off()

  ## Motif search with the MEME suite online (see the supplemental folder memeResults.zip for additional info)
  # Sequences for all the significant k-mers
  significantRanges <- makeGRangesFromDataFrame(data.frame("chr"=nanocomporeResultsSignificant$chr
                                                          ,"start"=nanocomporeResultsSignificant$genomicPos
                                                          ,"end"=nanocomporeResultsSignificant$genomicPos+5
                                                          ,"strand"=nanocomporeResultsSignificant$strand))
  significantRanges <- resize(significantRanges, width = 11, fix='center')

  newStyle <- mapSeqlevels(seqlevels(significantRanges), "UCSC")
  significantRanges <- renameSeqlevels(significantRanges, newStyle)
  names(significantRanges) <- 1:length(significantRanges)

  # Sequences for 5000 random k-mers
  allRanges <- makeGRangesFromDataFrame(data.frame("chr"=nanocomporeResults$chr
                                                  ,"start"=nanocomporeResults$genomicPos
                                                  ,"end"=nanocomporeResults$genomicPos+5
                                                  ,"strand"=nanocomporeResults$strand))
  allRanges <- resize(allRanges, width = 11, fix='center')

  newStyle <- mapSeqlevels(seqlevels(allRanges), "UCSC")
  allRanges <- renameSeqlevels(allRanges, newStyle)

  set.seed(2)
  backgroundRanges <- allRanges[sample(1:length(allRanges),5000),]
  names(backgroundRanges) <- 1:length(backgroundRanges)

  significantSequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10,significantRanges)
  backgroundSequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10,backgroundRanges)

  writeXStringSet(significantSequences,"significantSequences.fa")
  writeXStringSet(backgroundSequences,"backgroundSequences.fa")

### Gene ontology
  ## Ensembl 2 Entrez
  ensembl2entrez2symbol <- getBM(attributes = c("ensembl_transcript_id","ensembl_gene_id","entrezgene_id","external_gene_name"),
                                  filter = "ensembl_transcript_id",
                                  values = names(dds),
                                  mart = ensemblM)

  ## Gene IDs conversion, from ENS(G/T) to Entrez
  convertListToEntrez <- function(inputList,ensembl2entrez2symbol,key="g")
  {
    if(key=="g")
    {
      listTmp <- unlist(lapply(names(inputList),function(g)
      {
        entrezTmp <- ensembl2entrez2symbol[ensembl2entrez2symbol[,"ensembl_gene_id"]==g,"entrezgene_id"]
        pTmp <- rep(inputList[g],length(entrezTmp))
        names(pTmp) <- entrezTmp
        return(pTmp)
      }))    
    }else{
      listTmp <- unlist(lapply(names(inputList),function(g)
      {
        entrezTmp <- ensembl2entrez2symbol[ensembl2entrez2symbol[,"ensembl_transcript_id"]==g,"entrezgene_id"]
        pTmp <- rep(inputList[g],length(entrezTmp))
        names(pTmp) <- entrezTmp
        return(pTmp)
      }))
    }

    listTmp <- sort(listTmp,decreasing=FALSE)
    listTmp <- listTmp[is.finite(listTmp)]
    listTmp <- listTmp[!is.na(names(listTmp))]
    listTmp <- listTmp[unique(names(listTmp))]
    listTmp
  }

  pUS <- convertListToEntrez(pUS,ensembl2entrez2symbol)

  ## Expressed genes enrichments
  # BP
  BP_GO <- enrichGO(gene = names(pUS)
                    ,OrgDb = org.Mm.eg.db
                    ,ont = "BP"
                    ,readable = FALSE)

  BP_GO_simp <- simplifyGOterms(BP_GO[BP_GO$qvalue<=1e-40,"ID"],maxOverlap=0.15,ontology="BP",go2allEGs=org.Mm.egGO2ALLEGS)
  BP_GO[BP_GO[,"ID"]%in%BP_GO_simp,c("Description","qvalue")]

  #                                      Description       qvalue
  # GO:0042254                   ribosome biogenesis 4.137162e-50
  # GO:0000398        mRNA splicing, via spliceosome 9.669332e-49
  # GO:0022618    ribonucleoprotein complex assembly 9.620967e-47
  # GO:0007005            mitochondrion organization 4.720284e-44
  # GO:0010498 proteasomal protein catabolic process 1.349416e-41

  BP_GO[BP_GO$ID%in%c("GO:0007409","GO:0008088","GO:0010976","GO:0045666"),c("Description","qvalue")]
  #                                                     Description       qvalue
  # GO:0007409                                         axonogenesis 6.549498e-14
  # GO:0010976 positive regulation of neuron projection development 1.753411e-12
  # GO:0008088                              axo-dendritic transport 2.410731e-12
  # GO:0045666        positive regulation of neuron differentiation 3.251317e-11

  # MF
  MF_GO <- enrichGO(gene = names(pUS)
                    ,OrgDb = org.Mm.eg.db
                    ,ont = "MF"
                    ,readable = FALSE)

  MF_GO_simp <- simplifyGOterms(MF_GO[MF_GO$qvalue<=1e-40,"ID"],maxOverlap=0.15,ontology="MF",go2allEGs=org.Mm.egGO2ALLEGS)
  MF_GO[MF_GO[,"ID"]%in%MF_GO_simp,c("Description","qvalue")]

  #                                   Description       qvalue
  # GO:0003735 structural constituent of ribosome 1.288523e-68

  # CC
  CC_GO <- enrichGO(gene = names(pUS)
                    ,OrgDb = org.Mm.eg.db
                    ,ont = "CC"
                    ,readable = FALSE)

  CC_GO_simp <- simplifyGOterms(CC_GO[CC_GO$qvalue<=1e-40,"ID"],maxOverlap=0.15,ontology="CC",go2allEGs=org.Mm.egGO2ALLEGS)
  CC_GO[CC_GO[,"ID"]%in%CC_GO_simp,c("Description","qvalue")]

  #                                             Description       qvalue
  # GO:0098800 inner mitochondrial membrane protein complex 9.340954e-53
  # GO:0015934                      large ribosomal subunit 1.670651e-45
  # GO:0005761                       mitochondrial ribosome 1.670651e-45
  # GO:0005681                         spliceosomal complex 9.753158e-45
  # GO:0043209                                myelin sheath 6.936292e-44

  KEGG <- enrichKEGG(gene = names(pUS)
                    ,organism = "mmu"
                    ,keyType = "kegg"
                    ,universe = keys(org.Mm.eg.db)
                    ,pAdjustMethod = "BH")

  KEGG_simp <- KEGG[KEGG$qvalue<=1e-40,c("Description","qvalue")]

  enrichmentResults <- rbind(MF_GO[MF_GO[,"ID"]%in%MF_GO_simp,c("Description","qvalue")]
                            ,KEGG_simp
                            ,CC_GO[CC_GO[,"ID"]%in%CC_GO_simp,c("Description","qvalue")]
                            ,BP_GO[BP_GO[,"ID"]%in%BP_GO_simp,c("Description","qvalue")])

  pdf("enrichmentResults.pdf",width=11,height=8)
  par(mfrow=c(1,3),cex.lab=1.5)
  plot.new()
  x <- (-log10(enrichmentResults[,"qvalue"]))
  names(x) <- enrichmentResults[,"Description"]
  barplot(x,horiz=TRUE,main="",xlab="-Log10(Adjusted p-value)",las=2,cex.axis=1.5,cex.names=1.5,col=c(rep(3,length(MF_GO_simp)),rep(5,nrow(KEGG_simp)),rep(2,length(CC_GO_simp)),rep(4,length(BP_GO_simp))))
  plot.new()
  legend("center",legend=c("Biological Processes","Cellular Components","KEGG pathways","Molecular Functions"),col=c(3,5,2,4),pch=15,cex=1.5)
  dev.off()

### Output files reported in SupplementalTable1
## Sheet1 - m6A+_5-mers
write.csv(as.matrix(nanocomporeResults[is.finite(pValues)&pValues<=0.05,c("GeneId","TxId","chr","genomicPos","strand","ref_kmer","GMM_logit_pvalue")]),"m6A+_5-mers.csv")

## Sheet2 - DESeq2_results
write.csv(DESeq2_results,"DESeq2_results.csv")

## Sheet3 - DESeq2_normalizedCounts
write.csv(DESeq2_normalized_counts,"DESeq2_normalizedCounts.csv")

sessionInfo()
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-conda_cos6-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)

# Matrix products: default
# BLAS/LAPACK: /X/Y/Z/conda2/envs/mp_r3.6.1/lib/libopenblasp-r0.3.7.so

# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

# attached base packages:
#  [1] stats4    parallel  grid      stats     graphics  grDevices utils    
#  [8] datasets  methods   base     

# other attached packages:
#  [1] eulerr_6.1.1                             
#  [2] msigdbr_7.4.1                            
#  [3] readxl_1.3.1                             
#  [4] BSgenome.Mmusculus.UCSC.mm10_1.4.0       
#  [5] BSgenome_1.54.0                          
#  [6] DESeq2_1.26.0                            
#  [7] pheatmap_1.0.12                          
#  [8] GenomicAlignments_1.22.0                 
#  [9] Rsamtools_2.2.0                          
# [10] SummarizedExperiment_1.16.0              
# [11] DelayedArray_0.12.0                      
# [12] BiocParallel_1.20.0                      
# [13] matrixStats_0.55.0                       
# [14] ggrepel_0.8.1                            
# [15] compEpiTools_1.20.0                      
# [16] topGO_2.38.1                             
# [17] SparseM_1.77                             
# [18] GO.db_3.10.0                             
# [19] graph_1.64.0                             
# [20] TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0
# [21] org.Mm.eg.db_3.10.0                      
# [22] biomaRt_2.42.0                           
# [23] clusterProfiler_3.14.3                   
# [24] lemon_0.4.5                              
# [25] cowplot_1.1.1                            
# [26] Guitar_2.2.0                             
# [27] knitr_1.25                               
# [28] magrittr_2.0.1                           
# [29] rtracklayer_1.46.0                       
# [30] GenomicFeatures_1.38.0                   
# [31] AnnotationDbi_1.48.0                     
# [32] Biobase_2.46.0                           
# [33] motifStack_1.30.0                        
# [34] Biostrings_2.54.0                        
# [35] XVector_0.26.0                           
# [36] ade4_1.7-16                              
# [37] MotIV_1.42.0                             
# [38] GenomicRanges_1.38.0                     
# [39] GenomeInfoDb_1.22.0                      
# [40] IRanges_2.20.0                           
# [41] S4Vectors_0.24.0                         
# [42] BiocGenerics_0.32.0                      
# [43] grImport2_0.2-0                          
# [44] forcats_0.4.0                            
# [45] stringr_1.4.0                            
# [46] dplyr_1.0.2                              
# [47] purrr_0.3.3                              
# [48] readr_1.3.1                              
# [49] tidyr_1.0.0                              
# [50] tibble_2.1.3                             
# [51] ggplot2_3.3.3                            
# [52] tidyverse_1.2.1                          

# loaded via a namespace (and not attached):
#   [1] tidyselect_1.1.0         RSQLite_2.1.2            htmlwidgets_1.5.1       
#   [4] munsell_0.5.0            withr_2.1.2              colorspace_2.0-0        
#   [7] GOSemSim_2.12.1          rstudioapi_0.10          DOSE_3.12.0             
#  [10] urltools_1.7.3           GenomeInfoDbData_1.2.2   polyclip_1.10-0         
#  [13] bit64_0.9-7              farver_2.0.3             vctrs_0.3.6             
#  [16] generics_0.1.0           xfun_0.10                biovizBase_1.34.1       
#  [19] BiocFileCache_1.10.0     R6_2.4.0                 graphlayouts_0.7.1      
#  [22] locfit_1.5-9.1           AnnotationFilter_1.10.0  bitops_1.0-6            
#  [25] fgsea_1.12.0             gridGraphics_0.5-1       assertthat_0.2.1        
#  [28] scales_1.0.0             ggraph_2.0.4             nnet_7.3-12             
#  [31] enrichplot_1.6.1         gtable_0.3.0             tidygraph_1.2.0         
#  [34] ensembldb_2.10.0         seqLogo_1.52.0           rlang_0.4.10            
#  [37] genefilter_1.68.0        splines_3.6.1            lazyeval_0.2.2          
#  [40] acepack_1.4.1            dichromat_2.0-0          broom_0.5.2             
#  [43] europepmc_0.4            checkmate_1.9.4          BiocManager_1.30.10     
#  [46] reshape2_1.4.3           modelr_0.1.5             backports_1.2.1         
#  [49] qvalue_2.16.0            Hmisc_4.2-0              tools_3.6.1             
#  [52] ggplotify_0.0.5          gplots_3.0.1.1           RColorBrewer_1.1-2      
#  [55] ggridges_0.5.2           Rcpp_1.0.5               plyr_1.8.4              
#  [58] base64enc_0.1-3          progress_1.2.2           zlibbioc_1.32.0         
#  [61] RCurl_1.95-4.12          prettyunits_1.0.2        rpart_4.1-15            
#  [64] openssl_1.4.1            viridis_0.5.1            haven_2.1.1             
#  [67] cluster_2.1.0            data.table_1.13.6        DO.db_2.9               
#  [70] triebeard_0.3.0          ProtGenerics_1.18.0      xtable_1.8-4            
#  [73] hms_0.5.3                XML_3.98-1.20            jpeg_0.1-8.1            
#  [76] gridExtra_2.3            methylPipe_1.20.0        compiler_3.6.1          
#  [79] KernSmooth_2.23-18       crayon_1.3.4             htmltools_0.4.0         
#  [82] Formula_1.2-4            geneplotter_1.64.0       lubridate_1.7.4         
#  [85] DBI_1.1.0                tweenr_1.0.1             dbplyr_1.4.2            
#  [88] MASS_7.3-53              rappdirs_0.3.1           babelgene_21.4          
#  [91] rGADEM_2.34.1            Matrix_1.2-17            cli_1.1.0               
#  [94] marray_1.64.0            gdata_2.18.0             Gviz_1.14.2             
#  [97] igraph_1.2.4.1           pkgconfig_2.0.3          rvcheck_0.1.8           
# [100] foreign_0.8-72           xml2_1.2.2               annotate_1.64.0         
# [103] rvest_0.3.4              VariantAnnotation_1.32.0 digest_0.6.27           
# [106] cellranger_1.1.0         fastmatch_1.1-0          htmlTable_1.13.2        
# [109] curl_4.3                 gtools_3.8.2             lifecycle_0.2.0         
# [112] nlme_3.1-141             jsonlite_1.7.2           limma_3.42.0            
# [115] viridisLite_0.3.0        askpass_1.1              pillar_1.4.2            
# [118] lattice_0.20-41          httr_1.4.1               survival_2.44-1.1       
# [121] glue_1.4.2               png_0.1-7                bit_4.0.4               
# [124] ggforce_0.3.2            stringi_1.4.3            blob_1.2.1              
# [127] latticeExtra_0.6-28      caTools_1.18.0           memoise_1.1.0  