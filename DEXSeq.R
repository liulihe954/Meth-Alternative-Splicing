### prepare package and index 
suppressPackageStartupMessages(library("DEXSeq"))
library(tidyverse);library(readxl)
sample_index =read_excel("/blue/mateescu/lihe.liu/Methylation_WGCNA/Samples_RNA-Seq.xlsx")
control_index = dplyr::filter(sample_index,TRT == "a") %>% dplyr::select('Tube ID') %>% unlist(use.names = F)
treatment_index = dplyr::filter(sample_index,TRT == "b") %>%  dplyr::select('Tube ID') %>% unlist(use.names = F)
control_index_final = control_index[which(!(control_index %in% c('6268')))]
treatment_index_final = c(treatment_index[which(!(treatment_index %in% c('6228')))],'6228-1')
#conditions<-factor(c(rep("Control",9),rep("Meth",10)),levels=c("Control","Meth"))
# examples in: https://genomicsclass.github.io/book/pages/rnaseq_exon_usage.html

### read in count files
inDir = '/blue/mateescu/lihe.liu/AltSplicing/DEXSeqPyScripts/dexseq_count/'
countfiles = paste(inDir,paste("M",c(control_index_final,treatment_index_final),".count.txt",sep = ''),
                   sep='')
basename(countfiles)
flattenedFile = paste0('/blue/mateescu/lihe.liu/AltSplicing/ARS-UCD1.2/annotation/',
                       'flattened.clean.dexseq.gtf',
                       sep = '')
basename(flattenedFile)

# sample details
sampleTable = data.frame(
  row.names = c( paste0("M",c(control_index_final)),
                 paste0("M",c(treatment_index_final))),
  condition = c(rep('Control',length(control_index_final)),
                rep('Meth',length(treatment_index_final))))
  #,libType = rep('paired-end',length(c(control_index_final,treatment_index_final))))

# create test object
dxd_raw = DEXSeqDataSetFromHTSeq(
  countfiles,
  sampleData = sampleTable,
  design = ~ sample + exon + condition:exon,
  flattenedfile = flattenedFile
)

### filter-out low expr gene
all_gene = read.csv(countfiles[1],header = F,sep = '\t') %>%
  rename(Exon = V1, Count = V2) %>% 
  mutate(Exon = sub("\\:[0-9]*", "", Exon)) %>% 
  group_by(Exon) %>% 
  mutate(SampleSum = sum(Count)) %>% 
  slice_head() %>% 
  dplyr::select(-Count) %>% 
  mutate(Keep = ifelse(SampleSum >= 1, 1, 0)) %>% 
  dplyr::select(-SampleSum) 
head(all_gene)

#tmp = c()
for (i in c(2: length(basename(countfiles)))){
  all_gene_tmp = read.csv(countfiles[i],header = F,sep = '\t') %>% 
    rename(Exon = V1, Count = V2) %>% 
    mutate(Exon = sub("\\:[0-9]*", "", Exon)) %>% 
    group_by(Exon) %>% 
    mutate(SampleSum = sum(Count)) %>% 
    slice_head() %>% 
    dplyr::select(-Count) %>% 
    mutate(Keep = ifelse(SampleSum >= 1, 1, 0)) %>% 
    dplyr::select(-SampleSum) 
  #tmp[i] = all_gene_tmp[6,2]
  all_gene = all_gene_tmp %>% 
    left_join(all_gene, by=c('Exon' = 'Exon'))
}
#unlist(tmp,use.names = F)
#cm[1,]

all_gene_out = all_gene[-c(1:5),]
head(all_gene_out)

gene2keep = 
  all_gene[(which(rowSums(all_gene_out[,2:ncol(all_gene_out)])> 9)),1] %>% 
  unlist(use.names = F)

# keep only genes with high expr
dxd_filter = dxd_raw[geneIDs(dxd_raw) %in% gene2keep,]
colData(dxd_filter)

### Normalisation and Dispersion estimation
dxd_out1 = estimateSizeFactors(dxd_filter)
dxd_out2 = estimateDispersions(dxd_out1)
#plotDispEsts(dxd)

### models
fullModel<- ~ sample + exon + condition:exon 
reducedModel<- ~ sample + exon

###  test
dxd_DEU_DEXSeq = testForDEU(dxd_out2,
                 fullModel = fullModel,
                 reducedModel = reducedModel,
                 BPPARAM = MulticoreParam(workers=20))

dxd_DEU_DEXSeq_final = estimateExonFoldChanges(dxd_DEU_DEXSeq,fitExpToVar = "condition") # glm: count ∼ condition + exon + condition:exon.
dxd_DEU_DEXSeq_results = DEXSeqResults(dxd_DEU_DEXSeq_final)

#save(dxd,dxr1,file = 'DEXSeq_out.rda')
#save(dxd,dxr1,file = 'DEXSeq_out2.rda')
save(dxd_DEU_DEXSeq_final,dxd_DEU_DEXSeq_results,file = 'DEXSeq_out3.rda')
# head(dxr1,100)

load('DEXSeq_out3.rda')
#table ( dxd_DEU_DEXSeq_results$padj < 0.1 )
DEXSeq_final = data.frame(dxd_DEU_DEXSeq_results)
DEXSeq_out = DEXSeq_final[DEXSeq_final$pvalue < 0.01,]
dim(DEXSeq_out)
head(DEXSeq_out,100)
sig_gene = unique(DEXSeq_out$groupID)

length(sig_gene)

############################################################################################
# http://bioinf.wehi.edu.au/RNAseqCaseStudy/
library(limma)
library(edgeR)
cm<-counts(dxd_filter)[,1:19]
#head(cm)
cm2<-as.data.frame(rowRanges(dxd_filter))
#head(cm2)
geneInfo<-cm2#[,c(1:7)]
#dim(all_info)
#head(all_info)
#head(geneInfo)
y.all <- DGEList(counts=cm, genes=geneInfo)

isexpr <- rowSums(cpm(y.all) > 1) >= 9 # >= 9

y <- y.all[isexpr,keep.lib.sizes = FALSE]

save(y , file="expressCMNover_limma.RData", compress="xz")

y <- calcNormFactors(y)

conditions<-factor(c(rep("Control",9),
                     rep("Meth",10)),
                   levels=c("Control","Meth"))

design <- model.matrix( ~ conditions)

v <- voom(y,design,plot=FALSE)

fit <- lmFit(v, design)

fit.de <- eBayes(fit, robust=TRUE)

limmaResults<-data.frame(gene=fit.de$genes,
                         baseMean=exp(fit.de$coefficients[,1]), logFC=fit.de$coefficients[,2],
                         pval=fit.de$p.value[,2])
limmaResults$padj<-p.adjust(limmaResults$pval, method="BH")

ex <- diffSplice(fit[,"conditionsMeth"], geneid = "groupID", exonid = "start")
head(ex)
#dim(ex[[1]])

#F-tests for each gene.  ‘"t"’ gives t-tests for each exon. ‘
# "simes"’ gives genewise p-values derived from the t-tests after Simes adjustment for each gene.
#test_out = topSplice(ex,test="t")
DSRes<-topSplice(ex, test="simes", n = dim(ex)[1])
#DSRes_fdr<-topSplice(ex, test="simes", n = dim(ex)[1],FDR = 0.1) 

table(DSRes$P.Value <= 0.05)
table(DSRes$FDR <= 0.05)

# iso_info_limma<-iso_info[iso_info$gene_id %in% DSRes[,"groupID"],]
# DSRes<-DSRes[DSRes$groupID %in% unique(iso_info$gene_id ),]
# DELimma<-DSRes[DSRes$FDR < 0.05,"groupID"]

save(ex, file="limmaDSNOver.RData", compress="xz")
save(fit.de, file="fitdeLimmaNOver.RData", compress="xz")
save(DSRes, file="LimmaDSRes.RData", compress="xz")

load('limmaDSNOver.RData')
load('fitdeLimmaNOver.RData')
load('LimmaDSRes.RData')


