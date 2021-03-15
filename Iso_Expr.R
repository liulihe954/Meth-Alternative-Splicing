### pre
# setwd('../Iso_Expression/')
# iso_files = list.files()[-1]

#setwd('../counts/iso_counts/')
iso_files = list.files()[]
iso_res<-lapply(iso_files, function(iso_f){
  return(read.delim(iso_f))
})

names(iso_res)<-sapply(iso_files,function(iso_f){
  return(strsplit(iso_f, split="[.]")[[1]][1])
})

iso_cm<-as.data.frame(do.call(cbind, lapply(1:length(iso_res), function(x){
  iso_sample<-iso_res[[x]]
  return(iso_sample$expected_count)
})))

colnames(iso_cm) = names(iso_res)
rownames(iso_cm)<-iso_res[[1]]$transcript_id
length(unique(rownames(iso_cm))) # 43512 transcripts
# index
control_index_final = c('M6125','M6126','M6133','M6157','M6226','M6254','M6472','M6482','M6520')
treatment_index_final = c("M6152","M6153","M6197","M6412","M6415","M6419","M6438","M6453","M6485","M6228-1")
# order to 9-1o pattern
iso_cm_order = cbind(iso_cm[,(colnames(iso_cm)%in%control_index_final)],
                     iso_cm[,!(colnames(iso_cm)%in%control_index_final)])

head(iso_cm_order)

### DIE analysis
conditions<-factor(c(rep("Control",9),
                     rep("Meth",10)),
                   levels=c("Control","Meth"))

# Consider those isoforms expressed in at least one condition
#iso_cm_filter = iso_cm_order[(rowSums(iso_cm_order[,1:9]) >= 0 | rowSums(iso_cm_order[,10:19]) >= 0),]
to_keep = c()
for (i in seq_along(rownames(iso_cm_order))){
  print(i)
  tmp = iso_cm_order[i,]
  if (sum(tmp >= 1) >= 9){
    to_keep = append(to_keep,i)
  }
} #  19701

iso_cm_filter = iso_cm_order[to_keep,]

library(DESeq2)
## DESeq2
conditions<-factor(c(rep("Control",9),
                     rep("Meth",10)),
                   levels=c("Control","Meth"))

sample_data<-data.frame(condition=conditions)
rownames(sample_data) = names(iso_cm_filter)

dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(iso_cm_filter)),
                              colData = sample_data, 
                              design = formula(~condition,sample_data))
dds

# y<-estimateSizeFactors(y)
# y<-estimateDispersions(y)
# et_y_DES<-nbinomWaldTest(y)
out = DESeq(dds)
#out
#test = rowData(out)
#head(test)
#test= counts(dds) # original counts

DESeqres<-results(out,contrast = c("condition","Meth","Control"))

# Example from the manual
# dds <- makeExampleDESeqDataSet(m=4)
# dds <- DESeq(dds)
# res <- results(dds, contrast=c("condition","B","A"))
head(DESeqres)

#isoDESeq<-rownames(DESeqres[DESeqres$padj < 0.05 & !is.na(DESeqres$padj),])

save(DESeqres,out,file="DESeqRes1217.RData")

setwd('../Iso_Expression/')
load('DESeqRes1217.RData')

head(DESeqres)
table(DESeqres$pvalue <= 0.009) # 175
# 
# library(tidyverse)
# out = DESeqres %>%
#   data.frame() %>%
#   dplyr::filter(pvalue <= 0.009) %>%
#   arrange(pvalue)
# 
# library(openxlsx)
# 
# write.xlsx(out,row.names = T,
#            file = '/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/Manuscript/Supplementary_Files/Differential_Isoform_Expresion_Sig_out.xlsx')



# incoperate DIE and meth


# local: to generate sup file output for manuscript
load('/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/iso_with_count_in_exon.rda')
library(tidyverse)
head(out_format)
out_format_addE = out_format %>% 
  rowwise() %>% 
  mutate(exon_number = paste0('E00',exon_number))

data.frame(head(out_format_addE))
library(openxlsx)
write.xlsx(out_format_addE,row.names = F,
           file = '/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/Manuscript/Supplementary_Files/DIE_Meth.xlsx')


