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
  mutate(Exon = sub("\\:(.*)", "", Exon)) %>% 
  group_by(Exon) %>% 
  mutate(SampleSum = sum(Count)) %>% 
  slice_head() %>% 
  dplyr::select(-Count) %>% 
  mutate(Keep = ifelse(SampleSum >= 1, 1, 0)) %>% 
  dplyr::select(-SampleSum) 

all_gene_test = read.csv(countfiles[1],header = F,sep = '\t') %>%
  rename(Exon = V1, Count = V2) %>% 
  mutate(Exon = sub("\\:(.*)", "", Exon)) %>% 
  group_by(Exon) %>% 
  mutate(SampleSum = sum(Count)) %>% 
  slice_head() %>% 
  dplyr::select(-Count) %>% 
  mutate(Keep = ifelse(SampleSum >= 1, 1, 0)) %>% 
  dplyr::select(-SampleSum) 
head(all_gene_test)

#tmp = c()
for (i in c(2: length(basename(countfiles)))){
  all_gene_tmp = read.csv(countfiles[i],header = F,sep = '\t') %>% 
    rename(Exon = V1, Count = V2) %>% 
    mutate(Exon = sub("\\:(.*)", "", Exon)) %>% 
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
head(all_gene)
# rm ambiguous / empty 
all_gene_out = all_gene[-c(1:5),]
head(all_gene_out)

# keep genes expressed in > 9 replicates
gene2keep = 
  all_gene_out[(which(rowSums(all_gene_out[,2:ncol(all_gene_out)])>9)),1] %>% 
  unlist(use.names = F)
length(gene2keep) # 12056: previously

# keep only genes with high expr
dxd_filter = dxd_raw[geneIDs(dxd_raw) %in% gene2keep,]
colData(dxd_filter)
save(dxd_filter,file = 'DEXSeq_dxd_filter_4count.rda')
###########################################
##############   DEXSeq   ################
############################################

### Normalisation and Dispersion estimation
dxd_out1 = estimateSizeFactors(dxd_filter)
dxd_out2 = estimateDispersions(dxd_out1)
save(dxd_out1,dxd_out2,file = 'DEXSeq_process.rda')

#plotDispEsts(dxd)

### models
fullModel<- ~ sample + exon + condition:exon 
reducedModel<- ~ sample + exon

###  test
dxd_DEU_DEXSeq = testForDEU(
  dxd_out2,
  fullModel = fullModel,
  reducedModel = reducedModel,
  BPPARAM = MulticoreParam(workers=20))

dxd_DEU_DEXSeq_final = estimateExonFoldChanges(dxd_DEU_DEXSeq,fitExpToVar = "condition") # glm: count ∼ condition + exon + condition:exon.
dxd_DEU_DEXSeq_results = DEXSeqResults(dxd_DEU_DEXSeq_final)

#save(dxd,dxr1,file = 'DEXSeq_out.rda') #save(dxd,dxr1,file = 'DEXSeq_out2.rda')
save(dxd_DEU_DEXSeq, dxd_DEU_DEXSeq_final,dxd_DEU_DEXSeq_results,file = 'DEXSeq_out_all.rda')
# head(dxr1,100)
