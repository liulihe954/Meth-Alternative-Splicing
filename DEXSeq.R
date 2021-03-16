## Differential Exon Usage

# prepare package and index 
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
head(all_gene)
all_gene = read.csv(countfiles[1],header = F,sep = '\t') %>%
  dplyr::rename(Exon = V1, Count = V2) %>% 
  mutate(Exon = sub("\\:(.*)", "", Exon)) %>% 
  group_by(Exon) %>% 
  mutate(SampleSum = sum(Count)) %>% 
  slice_head() %>% 
  #dplyr::select(-Count) %>% 
  #mutate(Keep = ifelse(SampleSum >= 1, 1, 0)) %>% 
  dplyr::select(-SampleSum) 

#tmp = c()
for (i in c(2: length(basename(countfiles)))){
  all_gene_tmp = read.csv(countfiles[i],header = F,sep = '\t') %>% 
    dplyr::rename(Exon = V1, Count = V2) %>% 
    mutate(Exon = sub("\\:(.*)", "", Exon)) %>% 
    group_by(Exon) %>% 
    mutate(SampleSum = sum(Count)) %>% 
    slice_head() %>% 
    #dplyr::select(-Count) %>% 
    #mutate(Keep = ifelse(SampleSum >= 1, 1, 0)) %>% 
    dplyr::select(-SampleSum) 
  #tmp[i] = all_gene_tmp[6,2]
  all_gene = all_gene_tmp %>% 
    left_join(all_gene, by=c('Exon' = 'Exon'))
}
#unlist(tmp,use.names = F)
#cm[1,]

# rm ambiguous / empty 
all_gene_out = all_gene[-c(1:5),]
head(all_gene_out)
dim(all_gene_out)

all_gene_out = all_gene_out %>% 
  as_tibble() %>% 
  rowwise() %>% 
  mutate(rowsum = sum(get(colnames(all_gene_out)[2:20])))
table(all_gene_out$keep == 0)

# keep genes expressed in > 9 replicates
gene2keep = 
  all_gene_out[(which(rowSums(all_gene_out[,2:ncol(all_gene_out)])>9)),1] %>% 
  unlist(use.names = F)
length(gene2keep) # 12056: previously

#save(gene2keep,file = 'gene2keep.rda')
# keep only genes with high expr
dxd_filter = dxd_raw[geneIDs(dxd_raw) %in% gene2keep,]
#dim(dxd_raw)
#dim(dxd_filter)

colData(dxd_filter)

setwd('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R')
save(dxd_filter,file = 'DEXSeq_dxd_filter_4count.rda')

# test if genes kept overlap correctly
#load('DEXSeq_out_all.rda')
#DEXSeq_final_test = data.frame(dxd_DEU_DEXSeq_results)
#geneskept = unique(DEXSeq_final_test$groupID)

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

## PLOT1 - DEU  - Volcano



####################################################################################
############ Assemble: after count c proportion (multiple sbatch in hpg) ############
####################################################################################
suppressPackageStartupMessages(library("tidyverse"))
setwd('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/out/2b_assemble')
integrate_by_gene_filenames = list.files(path = ".",pattern=".rda")
#head(integrate_by_gene_filenames)
integrate_by_gene_out = data.frame()
for (i in seq_along(integrate_by_gene_filenames)){
  load(integrate_by_gene_filenames[i])
  print(i)
  #print(unique(DEXSeq_final_add_intr_out_final$groupID))
  integrate_by_gene_out = integrate_by_gene_out %>% 
     bind_rows(DEXSeq_final_add_intr_out_final)
  rm(DEXSeq_final_add_intr_out_final)
}
#length(unique(integrate_by_gene_out$groupID))
#
setwd('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/out/2b_assemble/assembled')
save(integrate_by_gene_out,file = 'integrate_by_gene_withCproportion_pval0.01.rda')
head(integrate_by_gene_out)

sum(integrate_by_gene_out$count_all) # 1311047
sum(integrate_by_gene_out$count_sig_p.01) # 51796

####################################################################################
############                local: Gather information                   ############
####################################################################################
# local: download output with three columns of propotions
load(
  '/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/DEXSeq_out_all.rda'
)
load(
  '/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/integrate_by_gene_withCproportion_pval0.01.rda'
)

#rm(integrate_by_gene_out)
#load('/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/integrate_by_gene_withCproportion.rda')

length(unique(integrate_by_gene_out$groupID)) # 12056 genes expressed
names(integrate_by_gene_out)
head(integrate_by_gene_out)
dim(integrate_by_gene_out)

integrate_by_gene_out_new = integrate_by_gene_out
integrate_by_gene_out_new_section = integrate_by_gene_out_new[,c(7:10)]

head(integrate_by_gene_out_new)
dim(integrate_by_gene_out_new)

# get pvalues from previous 
load('/Users/liulihe95/Desktop/Isoform-Expression/Iso_Expression/DIE_DEU_Universal_Info.rda')
view(head(Universal_exon_intron_transcript_final))
#Universal_exon_intron_transcript_final_test = Universal_exon_intron_transcript_final[1:30,]
Universal_exon_intron_transcript_final_new = Universal_exon_intron_transcript_final %>% 
  mutate(pvalue = ifelse(pvalue == 1,NA,pvalue)) %>% 
  dplyr::select(c(1:7))

Universal_DEU_info_sup_file = cbind(Universal_exon_intron_transcript_final_new,integrate_by_gene_out_new_section)
names(Universal_DEU_info_sup_file)

sum(Universal_DEU_info_sup_file$count_sig_p.01)
# load('../AltSplicing-R/DIE_selected_genes.rda')
# DIE_genes = all_gene_index
# DEU_genes = unique(Universal_DEU_info_sup_file$groupID)
# Overlap: F:940 T:11116 
head(DIE_genes[DIE_genes %in% DEU_genes])

sum(Universal_DEU_info_sup_file$count_sig_prom)

# Universal_DEU_info_sup_file_test = Universal_DEU_info_sup_file %>% 
#   mutate(EI_indic = substr(featureID,1,1)) %>% 
#   dplyr::filter(EI_indic == 'E')
# sum(Universal_DEU_info_sup_file_test$count_all)
# sum(Universal_DEU_info_sup_file_test$count_sig_p.01)

library(openxlsx)
write.xlsx(Universal_DEU_info_sup_file,row.names = F,
           file = '/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/Manuscript/Supplementary_Files/Universal_DEU_info.xlsx')



#### map summary

# missing
test_miss_chr = unique(as.character(Diff_C_all$chr))[!(unique(as.character(Diff_C_all$chr)) %in% unique(as.character(DEXSeq_final$genomicData.seqnames)))]
test_miss_chr 
# "chr21" "chr22" "chr23" "chr25" "chr26" 
#  WITH Condition
#  all_gene_out[(which(rowSums(all_gene_out[,2:ncol(all_gene_out)])>9)),1] 



load('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/out/DEXSeq_final.rda')
load('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/out/2b_assemble/assembled/integrate_by_gene_withCproportion.rda')
###
library(tidyverse)
aggrate_by_gene = integrate_by_gene_out %>% 
  #dplyr::slice(1:20) %>% 
  group_by(groupID) %>% 
  mutate(all_exon_prop = sum(count_sig[startsWith(featureID,'E')])/sum(count_all[startsWith(featureID,'E')])) %>% 
  mutate(all_intron_prop = sum(count_sig[startsWith(featureID,'I')])/sum(count_all[startsWith(featureID,'I')]))

aggrate_by_gene_to_join = aggrate_by_gene %>% 
  group_by(groupID) %>% 
  dplyr::select(all_exon_prop,all_intron_prop) %>% 
  slice(1)

####################################################################################
############  Merge1: add pvalue indicator columns (done seprately)       ###########
####################################################################################
load(
  '/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/aggrate_by_gene_withp.rda'
)
head(DEXSeq_final_add_intr_out_final)

####################################################################################
############ Merge2: add aggrate_by_gene_withp and integrate_by_gene_out    ############
####################################################################################
head(integrate_by_gene_out) # with all propotions
head(aggrate_by_gene_withp) # with all pvalues
names(aggrate_by_gene_withp)
library(tidyverse)

# start from all pvalues
aggrate_by_gene_withp_merge = aggrate_by_gene_withp %>% 
  dplyr::select(c(1:7)) %>% 
  mutate(Id_Univ = paste0(groupID,featureID))
head(aggrate_by_gene_withp_merge)
# left join id all propotions
integrate_by_gene_out_merge = integrate_by_gene_out %>% 
  as_tibble() %>% 
  dplyr::select(c(1:2,7:13)) %>% 
  mutate(Id_Univ = paste0(groupID,featureID)) %>% 
  dplyr::select(-groupID,-featureID) %>% 
  relocate(Id_Univ,.before=count_all)
head(integrate_by_gene_out_merge)

# left join all transcripts
names(DEXSeq_final)
DEXSeq_final_transcripts = DEXSeq_final %>% 
  as_tibble() %>% 
  dplyr::select(c(1,2,35)) %>% 
  mutate(Id_Univ = paste0(groupID,featureID)) %>% 
  dplyr::select(-groupID,-featureID) %>% 
  relocate(Id_Univ,.before=transcripts)
View(head(DEXSeq_final_transcripts))

###
Universal_exon_intron_transcript = aggrate_by_gene_withp_merge %>% 
  left_join(integrate_by_gene_out_merge, by =c('Id_Univ' = 'Id_Univ')) %>% 
  left_join(DEXSeq_final_transcripts, by =c('Id_Univ' = 'Id_Univ'))

head(Universal_exon_intron_transcript)

# loop to mark all significant ones
setwd('../Iso_Expression/')
load('DESeqRes1217.RData')
table(DESeqres$pvalue <= 0.009) # 175
transcript_thres = 0.009
sig_transcript_list = DESeqres %>% 
  as_tibble() %>% 
  mutate(transcript_id = rownames(DESeqres)) %>% 
  dplyr::filter(pvalue <= transcript_thres) %>% 
  pull(transcript_id)
length(sig_transcript_list)


sig_transcript = numeric(nrow(Universal_exon_intron_transcript))
for (i in seq_len(nrow(Universal_exon_intron_transcript))) {
  print(i)
  tmp_transcript = unlist(Universal_exon_intron_transcript[i,16],use.names = F)
  if (is.null(tmp_transcript)){sig_transcript[i] = NA} else {sig_transcript[i] = sum(tmp_transcript %in% sig_transcript_list)} 
}
Universal_exon_intron_transcript_final = Universal_exon_intron_transcript %>% 
  mutate(sig_transcript = sig_transcript) %>% 
  relocate(sig_transcript,.before=transcripts) %>% 
  dplyr::select(-Id_Univ) %>% 
  mutate(prop_sig.2 = count_sig.2/count_all,
         prop_sig.15 = count_sig.15/count_all,
         prop_sig.1 = count_sig.1/count_all)
View(head(Universal_exon_intron_transcript_final,20))

### output1 pvalue = 0.0009
Universal_exon_intron_transcript_final_sig = Universal_exon_intron_transcript_final %>% 
  dplyr::filter(pvalue <= 0.0009)

### output2 pvalue = 0.009
#
DE_Transcripts = data.frame(DESeqres) %>% 
  mutate(transcript_id = rownames(.)) %>% 
  relocate(transcript_id,.before = baseMean)
DE_Transcripts_sig = DE_Transcripts %>% 
  dplyr::filter(pvalue <= 0.009)

dim(DE_Transcripts_sig)

#
write.table(Universal_exon_intron_transcript_final,sep = ",",row.names = F,file = 'Universal_exon_intron_transcript_final.csv')
write.table(Universal_exon_intron_transcript_final_sig,sep = ",",row.names = F,file = 'Universal_exon_intron_transcript_final_sig.csv')
write.table(DE_Transcripts,row.names = F,sep = ",",file = 'DE_Transcripts.csv')
write.table(DE_Transcripts_sig,row.names = F,sep = ",",file = 'DE_Transcripts_sig.csv')

save(Universal_exon_intron_transcript_final,
     Universal_exon_intron_transcript_final_sig,
     DE_Transcripts,
     DE_Transcripts_sig,file = 'DIE_DEU_Universal_Info.rda')

sum(Universal_exon_intron_transcript_final$count_all)
sum(Universal_exon_intron_transcript_final$count_sig.2)

