# Prepare iso expr gene list - for proportion count

# format annotation files
setwd('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop')
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(tidyverse))
gtf_path = '/blue/mateescu/lihe.liu/AltSplicing/ARS-UCD1.2/annotation/Bos_taurus.ARS-UCD1.2.101.clean.v1.gtf'
gtf = rtracklayer::import(gtf_path)
gtf_df = as.data.frame(gtf) %>% as_tibble() %>% 
  dplyr::select(c(1,2,3,10,12,14,19))

head(gtf_df)
names(gtf_df)

save(gtf_df,file = 'bta_gtf_df.rda')

# get iso expr info
load('/blue/mateescu/lihe.liu/AltSplicing/rsem/iso_expr/DESeqRes1217.RData')
sig_transcript = DESeqres %>% 
  as.data.frame() %>% 
  mutate(TranscriptID = rownames(.)) %>% 
  dplyr::filter(pvalue <= 0.009) %>% 
  relocate(TranscriptID,.before = baseMean) %>% 
  pull(TranscriptID)

# extract gene id by transcript id 
sig_gene_iso_expr = gtf_df %>% 
  dplyr::filter(transcript_id %in% sig_transcript) %>% 
  distinct(gene_id,seqnames) %>% 
  relocate(gene_id,.before = seqnames) %>% 
  arrange(gene_id)
head(sig_gene_iso_expr)
#write.table(sig_gene_iso_expr,"GeneChr_Index_iso_expr.csv", col.names = F,row.names = F, quote = F)

# 