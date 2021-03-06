# Prepare iso expr gene list - for proportion count

# format annotation files
setwd('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop')
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(tidyverse))
gtf_path = '/blue/mateescu/lihe.liu/AltSplicing/ARS-UCD1.2/annotation/Bos_taurus.ARS-UCD1.2.101.clean.v1.gtf'

gtf = rtracklayer::import(gtf_path)

gtf_df = as.data.frame(gtf) %>% as_tibble() %>% dplyr::select(c(1,2,3,10,12,14,19))

data.frame(head(gtf_df,20))

# gtf_2_test = '/blue/mateescu/lihe.liu/AltSplicing/ARS-UCD1.2/annotation/flattened.clean.dexseq.gtf'
# gtf_2_test = rtracklayer::import(gtf_2_test)
# gtf_2_test_df = as.data.frame(gtf_2_test) %>% as_tibble() %>% 
#   dplyr::filter(type != 'aggregate_gene')

#data.frame(head(gtf_df,10))
#save(gtf_df,file = 'bta_gtf_df.rda')

# get iso expr info
load('/blue/mateescu/lihe.liu/AltSplicing/rsem/iso_expr/DESeqRes1217.RData')

# gtf_df_reduce = gtf_df %>% 
#   dplyr::select(gene_id, transcript_id)
# 
# DIE_map_gene = DESeqres %>% 
#   as.data.frame() %>% 
#   mutate(TranscriptID = rownames(.)) %>% 
#   dplyr::left_join(gtf_df_reduce, by = c('TranscriptID' = 'transcript_id')) %>% 
#   relocate(TranscriptID,.before = baseMean) %>% 
#   relocate(gene_id,.before = TranscriptID)
# save(DIE_map_gene,file = 'DIE_map_gene.rda')


sig_transcript = DESeqres %>% 
  as.data.frame() %>% 
  mutate(TranscriptID = rownames(.)) %>% 
  dplyr::filter(pvalue <= 0.009) %>% 
  relocate(TranscriptID,.before = baseMean) %>% 
  pull(TranscriptID)
head(DESeqres)
length(sig_transcript)

# extract all gene id input in DIE
all_transcripts = unique(rownames(DESeqres))

all_gene_iso_expr = gtf_df %>% 
  dplyr::filter(transcript_id %in% all_transcripts) %>% 
  distinct(gene_id,seqnames) %>% 
  relocate(gene_id,.before = seqnames) %>% 
  arrange(gene_id)

dim(all_gene_iso_expr)
head(all_gene_iso_expr)
#write.table(sig_gene_iso_expr,"GeneChr_Index_iso_expr.csv", col.names = F,row.names = F, quote = F)

# extract gene id by transcript id 
sig_gene_iso_expr = gtf_df %>% 
  dplyr::filter(transcript_id %in% sig_transcript) %>% 
  distinct(gene_id,seqnames) %>% 
  relocate(gene_id,.before = seqnames) %>% 
  arrange(gene_id) 
dim(sig_gene_iso_expr)


#write.table(sig_gene_iso_expr,"GeneChr_Index_iso_expr.csv", col.names = F,row.names = F, quote = F)


# extract gene id by transcript id - full
target_gene_index = sig_gene_iso_expr %>% pull(gene_id)
all_gene_index = all_gene_iso_expr %>% pull(gene_id)
# #head(target_gene_index)
# save(target_gene_index,all_gene_index,file = 'DIE_selected_genes.rda')

sig_gene_iso_expr_full = gtf_df %>% 
  dplyr::filter(gene_id %in% target_gene_index) %>% 
  group_by(gene_id) #%>% mutate(count_isoform = length(unique(transcript_id)))

dim(sig_gene_iso_expr_full)
data.frame(head(sig_gene_iso_expr_full))

#length(unique(sig_gene_iso_expr_full$transcript_id)) # matched 175
length(unique(sig_gene_iso_expr_full$gene_id)) # matched 175

all_gene_iso_expr_full = gtf_df %>% 
  dplyr::filter(gene_id %in% all_gene_index) %>% 
  group_by(gene_id) 
#%>% mutate(count_isoform = length(unique(transcript_id)))
dim(all_gene_iso_expr_full)

data.frame(head(all_gene_iso_expr_full))
# # check if gene has multiple transcripts
# checkiso = sig_gene_iso_expr_full %>% 
#   group_by(gene_id) %>% 
#   slice(1)
# table(checkiso$count_isoform)
# 
# selected_transcript_full = sig_gene_iso_expr_full %>% 
#   dplyr::filter(count_isoform != 1)
# dim(selected_transcript_full)

###
load('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/out/Diff_C_all.rda')
load('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/out/DEXSeq_final.rda')

# thresholds
pvalue.thres = 0.01
#count_all = numeric(nrow(all_gene_iso_expr_full))
#count_sig = numeric(nrow(all_gene_iso_expr_full))

# for supplementary files - output

library(tidyverse)
#names(all_gene_iso_expr_full)
#table(all_gene_iso_expr_full$seqnames)
chr_idx = c(paste0('chr',seq(1:29)),'chrX')

######                                ######
###### create sup material DIE&Meth  ######
######                                ######
# 
# out = data.frame()
# for (i in seq_along(chr_idx)){
#   #i = 3
#   chr_tmp = chr_idx[i]
#   print(paste0('checking ',chr_idx[i]))
#   # all Cs to check
#   Diff_C_all_subset = Diff_C_all %>%
#     dplyr::filter(chr == chr_tmp)
#   
#   # Sig Cs to check: pvalue <= 0.01
#   Diff_C_all_subset_sig = Diff_C_all %>%
#     dplyr::filter(chr == chr_tmp ,
#                   pvalue <= pvalue.thres) 
#   #
#   out_tmp = all_gene_iso_expr_full %>%
#     data.frame() %>%
#     dplyr::filter(seqnames == chr_tmp) # %>% arrange(start)
# 
#   count_all_tmp = numeric(nrow(out_tmp))
#   count_sig_tmp = numeric(nrow(out_tmp))
#   
#   count_all_prom = numeric(nrow(out_tmp))
#   count_sig_prom = numeric(nrow(out_tmp))
#   
#   #summary(out_tmp$start)
#   #sum(Diff_C_all_subset %in% c(310346:120848565))
#   #c(310346,120848565) # all exons
#   
#   for (j in seq_len(nrow(out_tmp))){
#     if(j %% 100 == 0){print(j)}
#     #j = 3
#     start_tmp = unlist(out_tmp[j,2],use.names = F)
#     end_tmp = unlist(out_tmp[j,3],use.names = F)
#     
#     count_all_tmp[j] = sum(Diff_C_all_subset$start %in% c(start_tmp:end_tmp))
#     count_sig_tmp[j] = sum(Diff_C_all_subset_sig$start %in% c(start_tmp:end_tmp))
#     
#     count_all_prom[j] = sum(Diff_C_all_subset$start %in% c((start_tmp - 3000):start_tmp))
#     count_sig_prom[j] = sum(Diff_C_all_subset_sig$start %in% c((start_tmp - 3000):start_tmp))
#     
#   }
#   print(summary(count_sig_tmp))
#   out_tmp['count_all'] = count_all_tmp
#   out_tmp['count_sig'] = count_sig_tmp
#   
#   out_tmp['count_all_prom_x'] = count_all_prom
#   out_tmp['count_sig_prom_x'] = count_sig_prom
#   
#   out = rbind(out,out_tmp)
# }
# 
# out_tmp = out %>% 
#   group_by(gene_id) %>% 
#   arrange(start) %>% 
#   mutate(count_all_prom = count_all_prom_x[1]) %>% 
#   mutate(count_sig_prom = count_sig_prom_x[1])
# 
# sum(out_tmp$count_all_prom)
# sum(out_tmp$count_sig_prom)
# 
# 
# pval_container = DESeqres %>%
#   data.frame() %>%
#   mutate(transcript_id = rownames(.)) %>%
#   dplyr::select(transcript_id,pvalue)
# 
# out_format = out_tmp %>%
#   relocate(transcript_id,.before = seqnames) %>%
#   relocate(gene_id,.before = seqnames) %>%
#   relocate(exon_number,.before = seqnames) %>%
#   dplyr::select(-exon_id) %>%
#   group_by(transcript_id) %>%
#   arrange(transcript_id) %>%
#   left_join(pval_container, by = c('transcript_id'='transcript_id')) %>%
#   relocate(pvalue,.before = gene_id)
# 
# data.frame(head(out_format,10))
# 
# save(out_format,file = 'iso_with_count_in_exon.rda')


########################################################################################################################


# all Cs to check
Diff_C_all_subset = Diff_C_all %>% 
  dplyr::filter(chr == chr_tmp) %>% 
  arrange(start)
# Sig Cs to check: pvalue <= 0.01
Diff_C_all_subset_sig = Diff_C_all %>% 
  dplyr::filter(chr == chr_tmp , 
                pvalue <= pvalue.thres) %>% 
  arrange(start)

for (i in seq_len(nrow(selected_transcript_full))){
  chr_tmp_new = as.character(unlist(selected_transcript_full[i,1]))
  if (chr_tmp_new != chr_tmp){
    chr_tmp = as.character(unlist(selected_transcript_full[i,1]))
    # all Cs to check
    Diff_C_all_subset = Diff_C_all %>% 
      dplyr::filter(chr == chr_tmp) %>% 
      arrange(start)
    # Sig Cs to check: pvalue <= 0.01
    Diff_C_all_subset_sig = Diff_C_all %>% 
      dplyr::filter(chr == chr_tmp, 
                    pvalue <= pvalue.thres) %>% 
      arrange(start)
  }
  
  # count
  start = as.numeric(unlist(selected_transcript_full[i,2]))
  end = as.numeric(unlist(selected_transcript_full[i,3]))
  
  count_all[i] = sum(Diff_C_all_subset$start %in% c(start:end))
  count_sig[i] = sum(Diff_C_all_subset_sig$start %in% c(start:end))
  print(i)
}

"/" <- function(x,y) ifelse(y==0,0,base:::"/"(x,y))
selected_transcript_with_prop = selected_transcript_full %>% 
  ungroup() %>% 
  mutate(countall = count_all,
         countsig = count_sig) %>% 
  mutate(proportion = count_sig/count_all)

table(selected_transcript_with_prop$proportion)



######
proportions_sig_transcript = c()
proportions_other_transcript = c()

gene_list_multi_transcript = selected_transcript_full %>% 
  pull(gene_id) %>% 
  unique()
length(gene_list_multi_transcript)

for (i in seq_along(gene_list_multi_transcript)){
  
  temp = selected_transcript_with_prop %>% 
    dplyr::filter(gene_id == gene_list_multi_transcript[i]) %>% 
    mutate(ifsig = ifelse(transcript_id %in% sig_transcript, 'Sig','No')) %>% 
    ungroup() %>% 
    group_by(ifsig) %>% 
    mutate(prop = sum(countsig)/sum(countall))
  #data.frame(temp)
  #
  subset_sig = temp %>% 
    dplyr::filter(ifsig == 'Sig')
  
  chr_tmp = as.character(unlist(subset_sig[1,1]))
  # all Cs to check
  Diff_C_all_subset = Diff_C_all %>% 
    dplyr::filter(chr == chr_tmp) %>% 
    arrange(start)
  # Sig Cs to check: pvalue <= 0.01
  Diff_C_all_subset_sig = Diff_C_all %>% 
    dplyr::filter(chr == chr_tmp, 
                  pvalue <= pvalue.thres) %>% 
    arrange(start)
  
  countall_sig_tmp = countsig_sig_tmp = 0
  for (m in seq_len(nrow(subset_sig))){
    #data.frame(subset_sig)
    #m = 5
    start = as.numeric(unlist(subset_sig[m,2]))
    end = as.numeric(unlist(subset_sig[m,3]))
    
    countall_sig_tmp = countall_sig_tmp + sum(Diff_C_all_subset$start %in% c(start:end))
    countsig_sig_tmp = countsig_sig_tmp + sum(Diff_C_all_subset_sig$start %in% c(start:end))
  }
  
  proportions_sig_transcript = append(proportions_sig_transcript,
                                      countsig_sig_tmp/countall_sig_tmp)
  
  #
  subset_nonsig = temp %>% 
    dplyr::filter(ifsig == 'No') %>% 
    dplyr::distinct(exon_id,.keep_all = TRUE)
  
  countall_sig_tmp2 = countsig_sig_tmp2 = 0
  
  for (n in seq_len(nrow(subset_nonsig))){
    #
    start = as.numeric(unlist(subset_nonsig[n,2]))
    end = as.numeric(unlist(subset_nonsig[n,3]))
    
    countall_sig_tmp2 = countall_sig_tmp2 + sum(Diff_C_all_subset$start %in% c(start:end))
    countsig_sig_tmp2 = countsig_sig_tmp2 + sum(Diff_C_all_subset_sig$start %in% c(start:end))
  }
  proportions_other_transcript = append(proportions_other_transcript,
                                        countsig_sig_tmp2/countall_sig_tmp2)
  ###
  print(i)
}

ks.test(proportions_sig_transcript,
        proportions_other_transcript)


##### compare promoters

# extract all gene id input in DIE
all_transcripts = unique(rownames(DESeqres))

all_gene_iso_expr = gtf_df %>% 
  dplyr::filter(transcript_id %in% all_transcripts) %>% 
  distinct(gene_id,seqnames) %>% 
  relocate(gene_id,.before = seqnames) %>% 
  arrange(gene_id)
dim(all_gene_iso_expr)

#write.table(sig_gene_iso_expr,"GeneChr_Index_iso_expr.csv", col.names = F,row.names = F, quote = F)

# extract gene id by transcript id 
sig_gene_iso_expr = gtf_df %>% 
  dplyr::filter(transcript_id %in% sig_transcript) %>% 
  distinct(gene_id,seqnames) %>% 
  relocate(gene_id,.before = seqnames) %>% 
  arrange(gene_id) 
dim(sig_gene_iso_expr)
#write.table(sig_gene_iso_expr,"GeneChr_Index_iso_expr.csv", col.names = F,row.names = F, quote = F)
head(sig_gene_iso_expr)

# main
pvalue.thres=0.01
gene_list_to_sample = all_gene_iso_expr$gene_id[!(all_gene_iso_expr$gene_id %in% sig_gene_iso_expr$gene_id)]

#####
ratio_sig = c()
sig_gene_list = sig_gene_iso_expr$gene_id
"/" <- function(x,y) ifelse(y==0,0,base:::"/"(x,y))
for (i in seq_along(sig_gene_list)){
  tmp1 = gtf_df %>% 
    dplyr::filter(gene_id == sig_gene_list[i]) %>% 
    arrange(start)
  chr1 = as.character(unlist(tmp1[1,1]))
  prom_end1 = as.numeric(unlist(tmp1[1,2]))
  prom_start1 = prom_end1 - 3000
  # C info
  Diff_C_all_tmp1 = Diff_C_all %>% 
    dplyr::filter(chr == chr1) %>% 
    arrange(start)
  
  #
  Diff_C_sig_tmp1 = Diff_C_all %>% 
    dplyr::filter(chr == chr1, 
                  pvalue <= pvalue.thres) %>% 
    arrange(start)
  dim(Diff_C_sig_tmp1)
  # ratio 1
  ratio1 = sum(Diff_C_all_tmp1$start %in% c(prom_start1:prom_end1)) / sum(Diff_C_sig_tmp1$start %in% c(prom_start1:prom_end1))
  ratio_sig = append(ratio_sig,ratio1)
  print(i)
}



all_ks_pval = c()
"/" <- function(x,y) ifelse(y==0,0,base:::"/"(x,y))
for (n in c(1:100)){
  gene_iso_expr_sample = sample(gene_list_to_sample, nrow(sig_gene_iso_expr), replace = FALSE)
  ratio_bg = c()
  for (i in seq_along(gene_iso_expr_sample)){
    tmp1 = gtf_df %>% 
      dplyr::filter(gene_id == gene_iso_expr_sample[i]) %>% 
      arrange(start)
    chr1 = as.character(unlist(tmp1[1,1]))
    prom_end1 = as.numeric(unlist(tmp1[1,2]))
    prom_start1 = prom_end1 - 3000
    # C info
    Diff_C_all_tmp1 = Diff_C_all %>% 
      dplyr::filter(chr == chr1) %>% 
      arrange(start)
    
    #
    Diff_C_sig_tmp1 = Diff_C_all %>% 
      dplyr::filter(chr == chr1, 
                    pvalue <= pvalue.thres) %>% 
      arrange(start)
    
    # ratio 1
    ratio1 = sum(Diff_C_all_tmp1$start %in% c(prom_start1:prom_end1)) / sum(Diff_C_sig_tmp1$start %in% c(prom_start1:prom_end1))
    ratio_bg = append(ratio_bg,ratio1)
  }
  all_ks_pval[n] = ks.test(ratio_sig,ratio_bg)$p.value
  print(n)
}

summary(all_ks_pval)

save(all_ks_pval,file = 'all_ks_pval.rda')

load('all_ks_pval.rda')


png('ks_statistics.png')
x <- all_ks_pval
h<-hist(x, breaks=10, col="red", xlab='P-value',main="KS statistics")
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)
dev.off()

