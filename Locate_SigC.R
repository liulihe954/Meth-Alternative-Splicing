######                                            ######
###### create sup material Sig C and everything  ######
######  
######
library(tidyverse)
load('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/out/Diff_C_all.rda')
load('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/out/DEXSeq_final.rda')

pvalue.thres = 0.01
chr_idx_ALL = c(paste0('chr',seq(1:29)),'chrX')

Diff_C_all_subset_sig = Diff_C_all %>%
  dplyr::filter(pvalue <= pvalue.thres) 
dim(Diff_C_all_subset_sig)
data.frame(head(Diff_C_all_subset_sig))


gtf_2_test = '/blue/mateescu/lihe.liu/AltSplicing/ARS-UCD1.2/annotation/flattened.clean.dexseq.gtf'
gtf_2_test = rtracklayer::import(gtf_2_test)
gtf_2_test_df = as.data.frame(gtf_2_test) %>% as_tibble() %>%
  dplyr::filter(type != 'aggregate_gene')


# test = gtf_2_test_df %>% 
#   dplyr::filter(gene_id == 'ENSBTAG00000006644') %>% 
#   dplyr::select(-strand,-source,-type,-score,-phase) %>% 
#   arrange(start) %>% 
#   mutate(exonic_part_number = paste0('E',exonic_part_number)) %>% 
#   rename(part_number = exonic_part_number)
# data.frame(test)

#exon_list = test
complete_gene = function(exon_list,prom_len = 3000){
#impute_intron_list = function(exon_list,prom_len = 3000){
  # GET prom row
  prom_start = as.numeric(exon_list[1,2]) - prom_len
  prom_end = as.numeric(exon_list[1,2]) - 1
  row_prom = exon_list[1,]
  row_prom['start'] = prom_start
  row_prom['end'] = prom_end
  row_prom['width'] = prom_len
  row_prom['transcripts'] = 'Upstream'
  row_prom['part_number'] = 'Upstream'
  
  # work on exon list
  if (length(exon_list)==1){
    out = dplyr::bind_rows(row_prom,exon_list)
  } else {
    # to do: get intron_list
    Intro_out = data.frame()
    n = 0
    for (i in seq_len(nrow(exon_list)-1)){
      row_intron_tmp = exon_list[i,]
      intron_start_tmp = as.numeric(exon_list[i,3]) + 1
      intron_end_tmp = as.numeric(exon_list[i+1,2]) - 1
      if (intron_start_tmp == (intron_end_tmp + 1)){
        next
      } else {
        n = n + 1
        row_intron_tmp['start'] = intron_start_tmp
        row_intron_tmp['end'] = intron_end_tmp
        row_intron_tmp['width'] = intron_end_tmp - intron_start_tmp + 1
        row_intron_tmp['transcripts'] = NA
        row_intron_tmp['part_number'] = paste0('I',paste0(rep("0",3-length(n)), collapse = ""),n)
        Intro_out = bind_rows(Intro_out,row_intron_tmp)
      }
    }
  # format output
  }
  out = dplyr::bind_rows(row_prom,exon_list,Intro_out) %>% 
    arrange(start)
  gene_head = as.numeric(out[1,2]) - 1
  gene_tail = as.numeric(out[nrow(out),3]) + 1
  output = list(out,gene_head,gene_tail)
  return(output)
}

#as_tibble(complete_gene(test)[[1]]) # works!

chr_idx_ALL

#chr_idx = chr_idx_ALL[1]

#as.numeric(test_out[3])
gene_info_full_all_chr = list()
for (chr in chr_idx_ALL){
  #
  print(paste0('checking - ',chr))
  df_per_chr = gtf_2_test_df %>% 
    dplyr::filter(seqnames == chr) %>% 
    dplyr::select(-strand,-source,-type,-score,-phase) %>% 
    mutate(exonic_part_number = paste0('E',exonic_part_number)) %>% 
    rename(part_number = exonic_part_number) %>% 
    arrange(start)

  #
  gene_list_per_chr = unique(df_per_chr$gene_id)

  print(paste0('total gene here - ',length(gene_list_per_chr)))
  #
  #
  intergenic_point_list = c(0)
  gene_info_full = data.frame()
  #
  for (gene in gene_list_per_chr){
    #print(paste0('checking gene - ',gene))
    gene_info_tmp = df_per_chr %>% 
      dplyr::filter(gene_id == gene)
    gene_info_full = bind_rows(gene_info_full,as_tibble(complete_gene(gene_info_tmp)[[1]]))
    intergenic_point_list = append(intergenic_point_list,c(complete_gene(gene_info_tmp)[[2]],complete_gene(gene_info_tmp)[[3]]))
  }
  intergenic_point_list = append(intergenic_point_list,intergenic_point_list[length(intergenic_point_list)])

  #data.frame(head(gene_info_full,20))
  #
  intergentic_rows = data.frame()
  for (point in seq_len(length(intergenic_point_list)/2)){
    start_tmp = intergenic_point_list[(2 * point)-1]
    end_tmp = intergenic_point_list[(2 * point)]
    intergentic_row = data.frame(
      seqnames = chr,
      start = start_tmp,
      end = end_tmp,
      width = end_tmp - start_tmp + 1,
      gene_id = 'IGR',
      transcripts = 'IGR',
      part_number = 'IGR'
    )
    intergentic_rows = bind_rows(intergentic_rows,intergentic_row)
  }
  gene_info_full = dplyr::bind_rows(gene_info_full,intergentic_rows) %>% 
    arrange(start)
  gene_info_full_all_chr[[chr]] = gene_info_full
  print(paste0('finished - ',chr))
}

length(gene_info_full_all_chr)
getwd()
save(chr_idx_ALL,gene_info_full_all_chr,file = '/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/matchC_with_exon/gene_info_full_all_chr.rda')
################################################################################
library(tidyverse)
load('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/matchC_with_exon/gene_info_full_all_chr.rda')
load('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/out/Diff_C_all.rda')

chr_idx_ALL_tmp = chr_idx_ALL[1]#
chr_idx_ALL_tmp = chr_idx_ALL[2]#
chr_idx_ALL_tmp = chr_idx_ALL[3] #
chr_idx_ALL_tmp = chr_idx_ALL[4] #
chr_idx_ALL_tmp = chr_idx_ALL[5] #
chr_idx_ALL_tmp = chr_idx_ALL[6]#
chr_idx_ALL_tmp = chr_idx_ALL[7]#
chr_idx_ALL_tmp = chr_idx_ALL[8] #
chr_idx_ALL_tmp = chr_idx_ALL[9] #
chr_idx_ALL_tmp = chr_idx_ALL[10] #

chr_idx_ALL_tmp = chr_idx_ALL[11]#
chr_idx_ALL_tmp = chr_idx_ALL[12] #
chr_idx_ALL_tmp = chr_idx_ALL[13] #
chr_idx_ALL_tmp = chr_idx_ALL[14] 
chr_idx_ALL_tmp = chr_idx_ALL[15] 
chr_idx_ALL_tmp = chr_idx_ALL[16] 
chr_idx_ALL_tmp = chr_idx_ALL[17]
chr_idx_ALL_tmp = chr_idx_ALL[18]
chr_idx_ALL_tmp = chr_idx_ALL[19] 
chr_idx_ALL_tmp = chr_idx_ALL[20] 

chr_idx_ALL_tmp = chr_idx_ALL[1]
chr_idx_ALL_tmp = chr_idx_ALL[2]
chr_idx_ALL_tmp = chr_idx_ALL[3] 
chr_idx_ALL_tmp = chr_idx_ALL[4] 
chr_idx_ALL_tmp = chr_idx_ALL[5] 
chr_idx_ALL_tmp = chr_idx_ALL[6]
chr_idx_ALL_tmp = chr_idx_ALL[7]
chr_idx_ALL_tmp = chr_idx_ALL[8] 
chr_idx_ALL_tmp = chr_idx_ALL[9] 
chr_idx_ALL_tmp = chr_idx_ALL[10] 

pvalue.thres = 0.01
Diff_C_all_subset_sig_for_sup = Diff_C_all %>%
  dplyr::filter(pvalue <= pvalue.thres) 
#dim(Diff_C_all_subset_sig_for_sup)
#head(Diff_C_all_subset_sig)

Diff_C_sig_with_region_list = list()
Diff_C_sig_with_region = data.frame()


for (chr_tmp in chr_idx_ALL_tmp){
  
  print(paste0('checking - ',chr_tmp))
  
  Diff_C_sig_with_region_tmp = Diff_C_all_subset_sig_for_sup %>% 
    dplyr::filter(chr == chr_tmp)
  
  loc_info_exon = gene_info_full_all_chr[[chr_tmp]] %>% 
    as_tibble() %>% 
    dplyr::filter(part_number != 'Upstream')
  
  loc_info_upstream = gene_info_full_all_chr[[chr_tmp]] %>% 
    as_tibble() %>% 
    dplyr::filter(part_number == 'Upstream')
  
  for (i in seq_len(nrow(Diff_C_sig_with_region_tmp))){
    print(paste0('checking - ',i))
    sig_c_loc = as.numeric(Diff_C_sig_with_region_tmp[i,2])
    find_index_exonlist = findInterval(sig_c_loc,loc_info_exon$start)
    
    if (sig_c_loc %in% c(as.numeric(loc_info_exon[find_index_exonlist,2]):as.numeric(loc_info_exon[find_index_exonlist,3]))){
      
      tmp_row = cbind(Diff_C_sig_with_region_tmp[i,],loc_info_exon[find_index_exonlist,5:7])
      
      Diff_C_sig_with_region = bind_rows(Diff_C_sig_with_region,tmp_row)
      
    } else {

      find_index_upstreamlist = findInterval(sig_c_loc,loc_info_upstream$start)
      
      if (sig_c_loc %in% c(as.numeric(loc_info_upstream[find_index_upstreamlist,2]):as.numeric(loc_info_exon[find_index_upstreamlist,3]))){
        
        tmp_row = cbind(Diff_C_sig_with_region_tmp[i,],loc_info_exon[find_index_upstreamlist,5:7])
        
        Diff_C_sig_with_region = bind_rows(Diff_C_sig_with_region,tmp_row)
      } else {
        
        tmp_addtion = data.frame(gene_id=NA, transcripts=NA, part_number = NA)
        
        tmp_row = cbind(Diff_C_sig_with_region_tmp[i,],tmp_addtion)
        
        Diff_C_sig_with_region = bind_rows(Diff_C_sig_with_region,tmp_row)
      }
    }
    Diff_C_sig_with_region_list[[chr_tmp]] = Diff_C_sig_with_region  
  }   
}

length(Diff_C_sig_with_region_list)
save(Diff_C_sig_with_region_list, file = paste0('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/matchC_with_exon/Diff_C_sig_with_region_list,',chr_idx_ALL_tmp,',.rda'))

# 
# test = Diff_C_sig_with_region_list[[chr]]
# data.frame(test[test$start == 2186932,])
# table(test$part_number)
# length(unique(test$start))

