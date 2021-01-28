# DIE
# format annotation files
setwd('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop')
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(tidyverse))
gtf_path = '/blue/mateescu/lihe.liu/AltSplicing/ARS-UCD1.2/annotation/Bos_taurus.ARS-UCD1.2.101.clean.v1.gtf'
gtf = rtracklayer::import(gtf_path)
gtf_df = as.data.frame(gtf) %>% as_tibble() %>% 
  dplyr::select(c(1,2,3,10,12,14,19))

# get iso expr info
load('/blue/mateescu/lihe.liu/AltSplicing/rsem/iso_expr/DESeqRes1217.RData')
sig_transcript = DESeqres %>% 
  as.data.frame() %>% 
  mutate(TranscriptID = rownames(.)) %>% 
  dplyr::filter(pvalue <= 0.009) %>% 
  relocate(TranscriptID,.before = baseMean) %>% 
  pull(TranscriptID)
#head(DESeqres)

# extract all gene id input in DIE
all_transcripts = unique(rownames(DESeqres))

all_gene_iso_expr = gtf_df %>% 
  dplyr::filter(transcript_id %in% all_transcripts) %>% 
  distinct(gene_id,seqnames) %>% 
  relocate(gene_id,.before = seqnames) %>% 
  arrange(gene_id) %>% 
  pull(gene_id) %>% unique()
length(all_gene_iso_expr)

# extract gene id by transcript id 
sig_gene_iso_expr = gtf_df %>% 
  dplyr::filter(transcript_id %in% sig_transcript) %>% 
  distinct(gene_id,seqnames) %>% 
  relocate(gene_id,.before = seqnames) %>% 
  arrange(gene_id) %>% 
  pull(gene_id) %>% unique()
length(sig_gene_iso_expr)



#### DEU
load('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/out/DEXSeq_final.rda')
all_gene_deu = DEXSeq_final %>% 
  dplyr::pull(groupID) %>% 
  unique()
length(all_gene_deu)

sig_gene_deu = DEXSeq_final %>% 
  dplyr::filter(pvalue <= 0.009) %>% 
  dplyr::pull(groupID) %>% 
  unique()
length(sig_gene_deu)


# output gene list
#setwd('')
save(all_gene_iso_expr,sig_gene_iso_expr,all_gene_deu,sig_gene_deu,
     file = '/blue/mateescu/lihe.liu/AltSplicing/Enrichment_DIE_DEU/Enrich_Gene_List.rda')



# donwload from hpg
load('/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/Enrichment_DIE_DEU/Enrich_Gene_List.rda')
library(EnrichKit)
# input format

#all_gene_iso_expr
#sig_gene_iso_expr

#all_gene_deu
#sig_gene_deu

## DIE
# convert and orgnize
GeneInfo_DIE = convertNformatID(GeneSetNames=c("DIE"),
                            SigGene_list = list(sig_gene_iso_expr),
                            TotalGene_list = list(sig_gene_iso_expr),
                            IDtype = "ens") # Need to choose from c('ens','entrez','symbol')

# Resulting an integreted gene identifier object
setwd("/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/Enrichment_DIE_DEU/DIE")
HyperGEnrich(GeneSet = GeneInfo_DIE,
             Database = 'kegg', #c("go","kegg","interpro","mesh","msig","reactome")
             minOverlap = 4, # minimum overlap of pathway genes and total genes
             pvalue_thres = 0.05, # pvalue of fisher's exact test
             adj_pvalue_thres = 1, # adjusted pvalues based on multiple testing correction
             padj_method = "fdr", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
             NewDB = F)
HyperGEnrich(GeneSet = GeneInfo_DIE,
             Database = 'go', #c("go","kegg","interpro","mesh","msig","reactome")
             minOverlap = 4, # minimum overlap of pathway genes and total genes
             pvalue_thres = 0.05, # pvalue of fisher's exact test
             adj_pvalue_thres = 1, # adjusted pvalues based on multiple testing correction
             padj_method = "fdr", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
             NewDB = F)
HyperGEnrich(GeneSet = GeneInfo_DIE,
             Database = 'interpro', #c("go","kegg","interpro","mesh","msig","reactome")
             minOverlap = 4, # minimum overlap of pathway genes and total genes
             pvalue_thres = 0.05, # pvalue of fisher's exact test
             adj_pvalue_thres = 1, # adjusted pvalues based on multiple testing correction
             padj_method = "fdr", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
             NewDB = F)
HyperGEnrich(GeneSet = GeneInfo_DIE,
             Database = 'mesh', #c("go","kegg","interpro","mesh","msig","reactome")
             minOverlap = 4, # minimum overlap of pathway genes and total genes
             pvalue_thres = 0.05, # pvalue of fisher's exact test
             adj_pvalue_thres = 1, # adjusted pvalues based on multiple testing correction
             padj_method = "fdr", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
             NewDB = F)
HyperGEnrich(GeneSet = GeneInfo_DIE,
             Database = 'msig', #c("go","kegg","interpro","mesh","msig","reactome")
             minOverlap = 4, # minimum overlap of pathway genes and total genes
             pvalue_thres = 0.05, # pvalue of fisher's exact test
             adj_pvalue_thres = 1, # adjusted pvalues based on multiple testing correction
             padj_method = "fdr", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
             NewDB = F)
HyperGEnrich(GeneSet = GeneInfo_DIE,
             Database = 'reactome', #c("go","kegg","interpro","mesh","msig","reactome")
             minOverlap = 4, # minimum overlap of pathway genes and total genes
             pvalue_thres = 0.05, # pvalue of fisher's exact test
             adj_pvalue_thres = 1, # adjusted pvalues based on multiple testing correction
             padj_method = "fdr", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
             NewDB = F)



## DEU
# convert and orgnize
GeneInfo_DEU = convertNformatID(GeneSetNames=c("DIE"),
                                SigGene_list = list(sig_gene_deu),
                                TotalGene_list = list(all_gene_deu),
                                IDtype = "ens") # Need to choose from c('ens','entrez','symbol')


# Resulting an integreted gene identifier object
setwd("/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/Enrichment_DIE_DEU/DEU")
HyperGEnrich(GeneSet = GeneInfo_DEU,
             Database = 'kegg', #c("go","kegg","interpro","mesh","msig","reactome")
             minOverlap = 4, # minimum overlap of pathway genes and total genes
             pvalue_thres = 0.05, # pvalue of fisher's exact test
             adj_pvalue_thres = 1, # adjusted pvalues based on multiple testing correction
             padj_method = "fdr", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
             NewDB = F)
HyperGEnrich(GeneSet = GeneInfo_DEU,
             Database = 'go', #c("go","kegg","interpro","mesh","msig","reactome")
             minOverlap = 4, # minimum overlap of pathway genes and total genes
             pvalue_thres = 0.05, # pvalue of fisher's exact test
             adj_pvalue_thres = 1, # adjusted pvalues based on multiple testing correction
             padj_method = "fdr", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
             NewDB = F)
HyperGEnrich(GeneSet = GeneInfo_DEU,
             Database = 'interpro', #c("go","kegg","interpro","mesh","msig","reactome")
             minOverlap = 4, # minimum overlap of pathway genes and total genes
             pvalue_thres = 0.05, # pvalue of fisher's exact test
             adj_pvalue_thres = 1, # adjusted pvalues based on multiple testing correction
             padj_method = "fdr", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
             NewDB = F)
HyperGEnrich(GeneSet = GeneInfo_DEU,
             Database = 'mesh', #c("go","kegg","interpro","mesh","msig","reactome")
             minOverlap = 4, # minimum overlap of pathway genes and total genes
             pvalue_thres = 0.05, # pvalue of fisher's exact test
             adj_pvalue_thres = 1, # adjusted pvalues based on multiple testing correction
             padj_method = "fdr", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
             NewDB = F)
HyperGEnrich(GeneSet = GeneInfo_DEU,
             Database = 'msig', #c("go","kegg","interpro","mesh","msig","reactome")
             minOverlap = 4, # minimum overlap of pathway genes and total genes
             pvalue_thres = 0.05, # pvalue of fisher's exact test
             adj_pvalue_thres = 1, # adjusted pvalues based on multiple testing correction
             padj_method = "fdr", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
             NewDB = F)
HyperGEnrich(GeneSet = GeneInfo_DEU,
             Database = 'reactome', #c("go","kegg","interpro","mesh","msig","reactome")
             minOverlap = 4, # minimum overlap of pathway genes and total genes
             pvalue_thres = 0.05, # pvalue of fisher's exact test
             adj_pvalue_thres = 1, # adjusted pvalues based on multiple testing correction
             padj_method = "fdr", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
             NewDB = F)

# out put
DEU_out_list = list.files("/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/Enrichment_DIE_DEU/DEU")
require(openxlsx)
list_of_datasets =list()
for (i in seq_along(DEU_out_list)){
  name_tmp = strsplit(DEU_out_list[i],"-")[[1]][1]
  print(name_tmp)
  load(DEU_out_list[i])
  list_of_datasets[name_tmp] = results
  rm(results)
}
write.xlsx(list_of_datasets, file = "/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/Enrichment_DIE_DEU/DEU/Enrich_DEU.xlsx")

#
DIE_out_list = list.files("/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/Enrichment_DIE_DEU/DIE")
require(openxlsx)
list_of_datasets =list()
for (i in seq_along(DIE_out_list)){
  name_tmp = strsplit(DIE_out_list[i],"-")[[1]][1]
  print(name_tmp)
  load(DIE_out_list[i])
  list_of_datasets[name_tmp] = results
  rm(results)
}
write.xlsx(list_of_datasets, file = "/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/Enrichment_DIE_DEU/DIE/Enrich_DIE.xlsx")


