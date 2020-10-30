#prepare g2t
library(biomaRt);library(tidyverse)
library(readxl)
database = useMart('ensembl')
genome = useDataset("btaurus_gene_ensembl", mart = database)
#listAttributes(mart = genome)
#searchAttributes(mart = genome, pattern = "transcript")
attributes = c("ensembl_gene_id","ensembl_transcript_id")
gene = getBM(attributes,mart = genome)
write.csv(gene,file = 'g2t.txt',row.names=FALSE)
g2t = read.csv('g2t.txt') %>% rename(gene_id = 'ensembl_gene_id',
                                     transcript_id = 'ensembl_transcript_id')
head(g2t)

#
setwd("/blue/mateescu/lihe.liu/AltSplicing/DEXSeqPyScripts/dexseq_count")
sample_index =read_excel("/blue/mateescu/lihe.liu/Methylation_WGCNA/Samples_RNA-Seq.xlsx")
control_index = dplyr::filter(sample_index,TRT == "a") %>% dplyr::select('Tube ID') %>% unlist(use.names = F)
treatment_index = dplyr::filter(sample_index,TRT == "b") %>%  dplyr::select('Tube ID') %>% unlist(use.names = F)
#conditions<-factor(c(rep("Control",9),rep("Meth",10)),levels=c("Control","Meth"))
control_index_final = control_index[which(!(control_index %in% c('6268')))]
treatment_index_final = c(treatment_index[which(!(treatment_index %in% c('6228')))],'6228-1')
#

# isoform counts
iso_files<-paste("M",c(control_index_final,treatment_index_final), ".count.txt", sep="")
iso_res<-lapply(iso_files, function(iso_f){
  tmp_reads = read.delim(iso_f,header = F)
  names(tmp_reads) = c('gene_exon','counts')
  return(tmp_reads)
})

names(iso_res)<-sapply(iso_files,function(iso_f){
  return(strsplit(iso_f, split="[.]")[[1]][1])
})

#class(iso_res[[1]])
#head(iso_res[[1]])

iso_cm<-as.data.frame(do.call(cbind, lapply(1:length(iso_res), function(x){
  iso_sample <-iso_res[[x]]
  #return(names(iso_sample))
  return(iso_sample[,1])
})))

#names(iso_cm) = names(iso_res)

length(rownames(iso_cm))
       
rownames(iso_cm)<-iso_res[[1]]$transcript_id

colnames(iso_cm)<-c(paste("C1R", c(1:10), sep=""), paste("C2R", c(1:10), sep=""))

iso_tpm<-as.data.frame(do.call(cbind, lapply(1:length(iso_res), function(x){
  iso_sample<-iso_res[[x]]
  return(iso_sample$TPM)
})))

rownames(iso_tpm)<-iso_res[[1]]$transcript_id
colnames(iso_tpm)<-c(paste("C1R", c(1:10), sep=""), paste("C2R", c(1:10), sep=""))
dim(iso_cm)
# [1] 169481     20




