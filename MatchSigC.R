# index = data.frame(Index = t(matrix(1:30,1)))
# write.table(index,
#             row.names = F,col.names = F,
#             file = '/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/matchC_with_exon/Chr_Idx.csv')
# 

args = commandArgs(trailingOnly = TRUE)

library(tidyverse)
load('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/matchC_with_exon/gene_info_full_all_chr.rda')
load('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/out/Diff_C_all.rda')

chr_idx_ALL_tmp = chr_idx_ALL[args[1]] 

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
save(Diff_C_sig_with_region_list, file = paste0('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/matchC_with_exon/Diff_C_sig_with_region_list_',chr_idx_ALL_tmp,'.rda'))

# 
# test = Diff_C_sig_with_region_list[[chr]]
# data.frame(test[test$start == 2186932,])
# table(test$part_number)
# length(unique(test$start))

