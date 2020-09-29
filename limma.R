## This is for limma
library(SplicingCompass); library(DEXSeq); library(limma)
library(readxl);library(tidyverse)
setwd("/blue/mateescu/lihe.liu/AltSplicing/DEXSeqPyScripts/dexseq_count")
sample_index =read_excel("/blue/mateescu/lihe.liu/Methylation_WGCNA/Samples_RNA-Seq.xlsx")
control_index = dplyr::filter(sample_index,TRT == "a") %>% dplyr::select('Tube ID') %>% unlist(use.names = F)
treatment_index = dplyr::filter(sample_index,TRT == "b") %>%  dplyr::select('Tube ID') %>% unlist(use.names = F)
conditions<-factor(c(rep("Control",9),rep("Meth",10)),levels=c("Control","Meth"))

expInf<-new("ExperimentInfo")
expInf<-setDescription(expInf,"NormalVsTumor")
expInf<-setGroupInfo(expInf, groupName1="Normal", sampleNumsGroup1=1:4,
                     groupName2="Tumor", sampleNumsGroup2=5:8)

covBedCountFilesNormal<-c(
  "/path_to_quantification_covBed_sim_scenario/SRR057649.covBed.counts",
  "/path_to_quantification_covBed_sim_scenario/SRR057650.covBed.counts",
  "/path_to_quantification_covBed_sim_scenario/SRR057651.covBed.counts",
  "/path_to_quantification_covBed_sim_scenario/SRR057652.covBed.counts")

covBedCountFilesTumor<-
  c("/path_to_quantification_covBed_sim_scenario/SRR057631.covBed.counts",
    "/path_to_quantification_covBed_sim_scenario/SRR057643.covBed.counts",
    "/path_to_quantification_covBed_sim_scenario/SRR057645.covBed.counts",
    "/path_to_quantification_covBed_sim_scenario/SRR057648.covBed.counts")

expInf<-setCovBedCountFiles(expInf, c(covBedCountFilesTumor,
                                      covBedCountFilesNormal))

junctionBedFilesNormal<-
  c("/path_to_alignment_genome_sim_scenario/SRR057649/junctions.bed",
    "/path_to_alignment_genome_sim_scenario/SRR057650/junctions.bed",
    "/path_to_alignment_genome_sim_scenario/SRR057651/junctions.bed",
    "/path_to_alignment_genome_sim_scenario/SRR057652/junctions.bed")
junctionBedFilesTumor<-
  c("/path_to_alignment_genome_sim_scenario/SRR057631/junctions.bed",
    "/path_to_alignment_genome_sim_scenario/SRR057643/junctions.bed",
    "/path_to_alignment_genome_sim_scenario/SRR057645/junctions.bed",
    "/path_to_alignment_genome_sim_scenario/SRR057648/junctions.bed")

expInf<-setJunctionBedFiles(expInf, c(junctionBedFilesTumor,junctionBedFilesNormal))
  expInf<-setReferenceAnnotation(expInf, "/path_to_annotation_files/flattened.splcmp.gtf")
  
  referenceAnnotationFormat<-list(IDFieldName="geneSymbol",idValSep=" ")
  expInf<-setReferenceAnnotationFormat(expInf,referenceAnnotationFormat)
  checkExperimentInfo(expInf)
  ## Constructing an object of class CountTable
  mycountTable<-new("CountTable")
  mycountTable<-constructCountTable(mycountTable,nCores=20, printDotPerGene=TRUE)
  sc<-new("SplicingCompass")
  sc<-constructSplicingCompass(sc, mycountTable, minOverallJunctionReadSupport=5,
                               nCores=20)
  # obataining significant DE genes
  sc<-initSigGenesFromResults(sc, adjusted=TRUE, threshold=0.05)
  sigGenes<-getSignificantGeneSymbols(sc)
  # obtaining a data frame with tested genes and correspoonding p-values
  resTab<-getResultTable(sc)
  resTab<-resTab[resTab$gene_id %in% iso_info$gene_id,]
  # tested gene ID
  genesTested<-getAllTestedGenes(sc)
  iso_info_SC<-iso_info[iso_info$gene_id %in% genesTested,]
  sigGenes<-getSignificantGeneSymbols(sc)
  sigGenes<-sigGenes[sigGenes %in% iso_info_SC$gene_id]
  save(resTab, file="SCresultTableNOv.RData")
  save(sc, file="SCobjectNOv.RData")
  save(mycountTable, file="SCCountTableNov.RData")
  
  #DEXSeq
  countfiles<-paste("path_to_quantification_DEXSeq_sim_scenario/htseq_sim_",
                    paste("SRR0576", c(49:52, 31, 43, 45, 48),".htseq.counts", sep="")
                    # building design matrix
                    sample_data<-data.frame(condition= conditions, levels=c("Normal","Tumor")))
  design<-~sample+exon+condition:exon
  row.names(sample_data)<-c(paste("C1R", 1:8, sep=""), paste("C2R", 1:8, sep=""))
  # build dexseq count matrix from HTseq output
  count_matrix<-DEXSeqDataSetFromHTSeq(countfiles, sample_data,design,
                                       flattenedfile="/path_to_annotation_files/flattened.dexseq.gtf")
  count_matrix<-estimateSizeFactors(count_matrix)
  count_matrix<-estimateDispersions(count_matrix, maxit=500,
                                    BPPARAM=MulticoreParam(workers=18))
  fullModel<- ~sample + exon + condition:exon
  reducedModel<- ~sample + exon
  count_matrix<-testForDEU(count_matrix, fullModel=fullModel,
                           reducedModel=reducedModel, BPPARAM=MulticoreParam(workers=20))
  count_matrix<-estimateExonFoldChanges(count_matrix, fitExpToVar="condition",
                                        BPPARAM=MulticoreParam(workers=20), denominator="Normal")
  myresults<-DEXSeqResults( count_matrix )
  perGeneQ<-perGeneQValue(myresults)
  myresultsDF<-as.data.frame(myresults)
  myresultsDF<-myresultsDF[!is.na(myresultsDF$padj) ,]
  myresultsDF$qvalGene<-do.call(c, lapply(1:nrow(myresultsDF), function(i){
    return(perGeneQ[names(perGeneQ) == myresultsDF$groupID[i]])
  }))
  myresultsDF<-myresultsDF[myresultsDF[,"groupID"]%in%unique(iso_info$gene_id), ]
  iso_info_dexseq<-iso_info[iso_info$gene_id %in%
                              unique(myresultsDF[,"groupID"]), ]
  DEXGenes<-unique(myresultsDF[myresultsDF$qvalGene < 0.05,"groupID"])
  save(count_matrix, file="DEXSeqCountMatrixSim.RData")
  
  
  #LimmaDS
  cm<-counts(count_matrix)[,1:8]
  cm2<-as.data.frame(rowRanges(count_matrix))
  geneInfo<-cm2[,c(1:7)]
  y.all <- DGEList(counts=cm, genes=geneInfo)
  isexpr <- rowSums(cpm(y.all) > 1) >=3
  y <- y.all[isexpr,,keep.lib.sizes=FALSE]
  save(y , file="expressCMNover.RData", compress="xz")
  y <- calcNormFactors(y)
  design <- model.matrix(~ conditions)
  v <- voom(y,design,plot=FALSE)
  fit <- lmFit(v, design)
  fit.de <- eBayes(fit, robust=TRUE)
  limmaResults<-data.frame(gene=fit.de$genes,
                           baseMean=exp(fit.de$coefficients[,1]), logFC=fit.de$coefficients[,2],
                           pval=fit.de$p.value[,2])
  limmaResults$padj<-p.adjust(limmaResults$pval, method="BH")
  ex <- diffSplice(fit[,"conditionTumor"], geneid = "groupID", exonid = "start")
  DSRes<-topSplice(ex, test="simes", n=length(ex))
  iso_info_limma<-iso_info[iso_info$gene_id %in% DSRes[,"groupID"],]
  DSRes<-DSRes[DSRes$groupID %in% unique(iso_info$gene_id ),]
  DELimma<-DSRes[DSRes$FDR < 0.05,"groupID"]
  save(ex, file="limmaDSNOver.RData", compress="xz")
  save(fit.de, file="fitdeLimmaNOver.RData", compress="xz")
  save(DSRes, file="LimmaDSRes.RData", compress="xz")