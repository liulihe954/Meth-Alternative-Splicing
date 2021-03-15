# create gene index for loop - col1: geneid col2:chr

# test = read.table('GeneChr_Index.csv')
# test_1 = test[c(1:6000),]
# test_2 = test[c(6001:nrow(test)),]
# write.table(test_1,"GeneChr_Index_6000head.csv", col.names = F,row.names = F, quote = F)
# write.table(test_2,"GeneChr_Index_tail.csv", col.names = F,row.names = F, quote = F)

# read in results (to save as new rda; run just once)
# setwd('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R')
# load('DEXSeq_out_all.rda') # new
# load('Meth_Diff_CpG_Universe.rda')
# DEXSeq_final = data.frame(dxd_DEU_DEXSeq_results)
# DEXSeq_final$padj<-p.adjust(DEXSeq_final$pvalue, method="BH")

# load('gene2keep.rda') # check gene kept remain the same, run just once
# TE = read.csv('meth_prop/GeneChr_Index.csv',header = F,sep = '')
# head(TE)
# TEST_genes =as.character(unique(TE[,1]))
# table(TEST_genes %in% gene2keep)
# head(gene2keep)

# out put gene index (run only once)
# setwd('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop')
# GeneChr_Index = DEXSeq_final %>%
#   as_tibble() %>%
#   dplyr::select(-c(16:34)) %>%
#   group_by(groupID) %>%
#   dplyr::select(groupID,genomicData.seqnames) %>%
#   dplyr::distinct()
# write.table(GeneChr_Index,"GeneChr_Index.csv", col.names = F,row.names = F, quote = F)

args <- commandArgs(trailingOnly = TRUE)
#args = c('ENSBTAG00000000015','chr5')

# print(args[1])
# print(args[2])

suppressPackageStartupMessages(library("tidyverse"))

setwd('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/out')
# save(Diff_C_all,file = 'Diff_C_all.rda')
# save(DEXSeq_final,file = 'DEXSeq_final.rda')

load('Diff_C_all.rda')
load('DEXSeq_final.rda')

# thresholds
pvalue.thres = 0.01
#meth.diff.thres = 20
# test = Diff_C_all %>%
#   dplyr::filter(pvalue <= pvalue.thres) %>%
#   arrange(start)
# dim(test)


#qval.thres.2 = 0.2;qval.thres.15 = 0.15;qval.thres.1 = 0.1

# massge original dataset : add intro info
DEXSeq_final_add_intr = DEXSeq_final %>%
  as_tibble() %>%
  dplyr::select(c(1:2,11:14)) %>%
  dplyr::filter(groupID == args[1],genomicData.seqnames == args[2]) %>% 
  arrange(genomicData.start)

# get subset of diff C
# threshold

# all
Diff_C_all_subset = Diff_C_all %>% 
  dplyr::filter(chr == args[2]) %>% 
  arrange(start)

# pvalue <= 0.01
Diff_C_all_subset_sig_p.01 = Diff_C_all %>% 
  dplyr::filter(chr == args[2], 
                pvalue <= pvalue.thres) %>% 
  arrange(start)

# # .2
# Diff_C_all_subset_sig.15 = Diff_C_all %>% 
#   dplyr::filter(chr == args[2], 
#                 qvalue <= qval.thres.2, 
#                 abs(meth.diff) >= meth.diff.thres) %>% 
#   arrange(start)
# # .15
# Diff_C_all_subset_sig.15 = Diff_C_all %>% 
#   dplyr::filter(chr == args[2], 
#                 qvalue <= qval.thres.15, 
#                 abs(meth.diff) >= meth.diff.thres) %>% 
#   arrange(start)
# # .1
# Diff_C_all_subset_sig.1 = Diff_C_all %>% 
#   dplyr::filter(chr == args[2], 
#                 qvalue <= qval.thres.1, 
#                 abs(meth.diff) >= meth.diff.thres) %>% 
#   arrange(start)


# loop
if (nrow(DEXSeq_final_add_intr) <= 1){
  DEXSeq_final_add_intr_out_final = DEXSeq_final_add_intr %>% 
    mutate(count_all = 0,
           count_sig_p.01 = 0,
           count_all_prom = 0,
           count_sig_prom = 0
           # count_sig.2= 0,
           # count_sig.15= 0,
           # count_sig.1= 0,
           # 
           # prop_sig.2 = 0,
           # prop_sig.15 = 0,
           # prop_sig.1 = 0
           )
  save(DEXSeq_final_add_intr_out_final,file = paste0('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/out/2b_assemble/count_c_',args[1],'_p.01.rda'))
} else {
  tmp_row = DEXSeq_final_add_intr[1,]
  i = 1;n = 1
  Intro_out = data.frame()
  while (i < nrow(DEXSeq_final_add_intr)){
    last_exon_end = as.numeric(DEXSeq_final_add_intr[i,5])
    next_exon_start = as.numeric(DEXSeq_final_add_intr[(i+1),4])
    if ((next_exon_start-1) == last_exon_end){i = i + 1;next} 
    if (n > 3) {tmp_row[1,2] == paste0('I',paste0(rep("0",0), collapse = ""),n)}
    tmp_row[1,2] = paste0('I',paste0(rep("0",3-length(n)), collapse = ""),n)
    tmp_row[1,4] = last_exon_end + 1
    tmp_row[1,5] = next_exon_start-1
    tmp_row[1,6] = (next_exon_start-1) - (last_exon_end + 1) + 1 
    Intro_out = Intro_out %>% bind_rows(tmp_row)
    i = i + 1
    n = n + 1
  }
  
  DEXSeq_final_add_intr_out = DEXSeq_final_add_intr %>% 
    bind_rows(Intro_out) %>% 
    arrange(genomicData.start)
  
  # go over each row of exon/intron infomation
  
  # all sections
  count_all = numeric(nrow(DEXSeq_final_add_intr_out))
  count_sig_p.01 = numeric(nrow(DEXSeq_final_add_intr_out))
  
  # promoter
  count_all_prom = numeric(nrow(DEXSeq_final_add_intr_out))
  count_sig_prom = numeric(nrow(DEXSeq_final_add_intr_out))
  
  #count_sig.2 = numeric(nrow(DEXSeq_final_add_intr_out))
  #
  #count_sig.15 = numeric(nrow(DEXSeq_final_add_intr_out))
  #
  #count_sig.1 = numeric(nrow(DEXSeq_final_add_intr_out))
  #
  for (i in seq_along(count_all)){
    
    start = as.numeric(DEXSeq_final_add_intr_out[i,4])
    end = as.numeric(DEXSeq_final_add_intr_out[i,5])
    
    tmp_count_all = sum(Diff_C_all_subset$start %in% c(start:end))
    tmp_count_all.p.01 = sum(Diff_C_all_subset_sig_p.01$start %in% c(start:end))
    
    # tmp_count_all.2 = sum(Diff_C_all_subset_sig.2$start %in% c(start:end))
    # tmp_count_all.15 = sum(Diff_C_all_subset_sig.15$start %in% c(start:end))
    # tmp_count_all.1 = sum(Diff_C_all_subset_sig.1$start %in% c(start:end))
    
    # all sections
    count_all[i] = tmp_count_all
    count_sig_p.01[i] = tmp_count_all.p.01
    
    
    # promoters
    prom_end = as.numeric(DEXSeq_final_add_intr_out[1,4])
    prom_start = prom_end - 5000
    
    count_all_prom[i] = sum(Diff_C_all_subset$start %in% c(prom_start:prom_end))
    count_sig_prom[i] = sum(Diff_C_all_subset_sig_p.01$start %in% c(prom_start:prom_end))
    
    
    # count_sig.2[i] = tmp_count_all.2
    # count_sig.15[i] = tmp_count_all.15
    # count_sig.1[i] = tmp_count_all.1
  }
  
  #"/" <- function(x,y) ifelse(y==0,0,base:::"/"(x,y))
  
  DEXSeq_final_add_intr_out_final = DEXSeq_final_add_intr_out %>% 
    mutate(count_all = count_all, # count all in sections
           count_sig_p.01 = count_sig_p.01, # count sig in sections
           
           count_all_prom = count_all_prom, # count all in promoter
           count_sig_prom = count_sig_prom # count sig in promoter
           # count_sig.2 = count_sig.2,
           # count_sig.15 = count_sig.15,
           # count_sig.1 = count_sig.1,
           # 
           # prop_sig.2 = count_sig.2/count_all,
           # prop_sig.15 = count_sig.15/count_all,
           # prop_sig.1 = count_sig.1/count_all
           ) 
  save(DEXSeq_final_add_intr_out_final,file = paste0('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/out/2b_assemble/count_c_',args[1],'_p.01.rda'))
}
data.frame(DEXSeq_final_add_intr_out_final)
