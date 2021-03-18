# index = data.frame(Index = t(matrix(1:30,1)))
# write.table(index,
#             row.names = F,col.names = F,
#             file = '/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/matchC_with_exon/Chr_Idx.csv')
# 
#args = commandArgs(trailingOnly = TRUE)
#print(args)
library(tidyverse)
load('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/matchC_with_exon/gene_info_full_all_chr.rda')
load('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/out/Diff_C_all.rda')
m = 30

#chr_idx_ALL_tmp = chr_idx_ALL[as.numeric(args[1])]
#print(chr_idx_ALL_tmp)

chr_idx_ALL_tmp = chr_idx_ALL[m]

pvalue.thres = 0.01

Diff_C_all_subset_sig_for_sup = Diff_C_all %>%
  dplyr::filter(pvalue <= pvalue.thres) 
#dim(Diff_C_all_subset_sig_for_sup)
#head(Diff_C_all_subset_sig)

Diff_C_sig_with_region_list = list()
Diff_C_sig_with_region = data.frame()

# Diff_C_sig_with_region_tmp = Diff_C_all_subset_sig_for_sup %>% 
#   dplyr::filter(chr == chr_idx_ALL_tmp)
# dim(Diff_C_sig_with_region_tmp)
for (chr_tmp in chr_idx_ALL_tmp){
  
  print(paste0('checking - ',chr_tmp))
  
  Diff_C_sig_with_region_tmp = Diff_C_all_subset_sig_for_sup %>% 
    dplyr::filter(chr == chr_tmp) #%>% 
    #slice_head(n = 9082) # chr 7
    #slice_tail(n = 280) # chr 30
  #dim(Diff_C_sig_with_region_tmp)
  
  loc_info_exon = gene_info_full_all_chr[[chr_tmp]] %>% 
    as_tibble() %>% 
    dplyr::filter(part_number != 'Upstream')
  
  loc_info_upstream = gene_info_full_all_chr[[chr_tmp]] %>% 
    as_tibble() %>% 
    dplyr::filter(part_number == 'Upstream')
  
  for (i in seq_len(nrow(Diff_C_sig_with_region_tmp))){
    #i = 38
    print(paste0('checking - ',i))
    
    sig_c_loc = as.numeric(Diff_C_sig_with_region_tmp[i,2])
    find_index_exonlist = findInterval(sig_c_loc,loc_info_exon$start)
    
    if (sig_c_loc %in% c(as.numeric(loc_info_exon[find_index_exonlist,2]):as.numeric(loc_info_exon[find_index_exonlist,3]))){
      
      tmp_row = cbind(Diff_C_sig_with_region_tmp[i,],loc_info_exon[find_index_exonlist,5:7])
      
      Diff_C_sig_with_region = bind_rows(Diff_C_sig_with_region,tmp_row)
      #print(tmp_row) #
      
    } else {
      
      find_index_upstreamlist = findInterval(sig_c_loc,loc_info_upstream$start)
      
      if (sig_c_loc %in% c(as.numeric(loc_info_upstream[find_index_upstreamlist,2]):as.numeric(loc_info_upstream[find_index_upstreamlist,3]))){
        
        tmp_row = cbind(Diff_C_sig_with_region_tmp[i,],loc_info_upstream[find_index_upstreamlist,5:7])
        
        Diff_C_sig_with_region = bind_rows(Diff_C_sig_with_region,tmp_row)
        
        print(tmp_row) #
      } else {
        
        tmp_addtion = data.frame(gene_id= '-', transcripts= '-', part_number = 'IGR')
        
        tmp_row = cbind(Diff_C_sig_with_region_tmp[i,],tmp_addtion)
        
        Diff_C_sig_with_region = bind_rows(Diff_C_sig_with_region,tmp_row)
        
        #print(tmp_row) #
      }
    }
    Diff_C_sig_with_region_list[[chr_tmp]] = Diff_C_sig_with_region  
  }   
}

length(Diff_C_sig_with_region_list)
name = paste0('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/matchC_with_exon/man2/Diff_C_sig_with_region_list_',chr_idx_ALL_tmp,'.rda')
save(Diff_C_sig_with_region_list, file = name)

# 
# #### Integrate everything ####

prefix = '/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/matchC_with_exon/man2/Diff_C_sig_with_region_list_chr'
suffix = '.rda'

sup_material_8_MatchSigC2Region = data.frame()
total = 0
library(dplyr)
for (index in c(1:29,'X')){
  if (exists('Diff_C_sig_with_region_list')){rm(Diff_C_sig_with_region_list)}
  # load
  load(paste0(prefix,index,suffix))
  # transfer
  tmp_file = Diff_C_sig_with_region_list[[1]]
  print(nrow(tmp_file))
  
  total = total + nrow(tmp_file)

  #print(table(tmp_file$part_number))

  # compile
  sup_material_8_MatchSigC2Region = dplyr::bind_rows(sup_material_8_MatchSigC2Region,tmp_file)
}


names(sup_material_8_MatchSigC2Region)

table(sup_material_8_MatchSigC2Region$part_number)

sup_material_8_MatchSigC2Region_out = sup_material_8_MatchSigC2Region %>%
  rowwise() %>%
  mutate(Transcript = ifelse( part_number %in% c('Upstream','IGR'),'-',transcripts))%>% 
  mutate(Gene = ifelse(part_number %in% c('IGR'),'-',gene_id)) %>% 
  mutate(Transcript = ifelse(is.na(Transcript),'-',Transcript))%>%
  mutate(Gene = ifelse(is.na(Gene),'-',Gene))#%>% 
  #mutate(part_number = ifelse(is.na(part_number),'IGR',part_number))

# 
# as.numeric(table(sup_material_8_MatchSigC2Region$part_number)['IGR'])
# 
# head(sup_material_8_MatchSigC2Region)
# 
# table(is.na(sup_material_8_MatchSigC2Region_out$gene_id))
# table(is.na(sup_material_8_MatchSigC2Region_out$transcripts))
# table(is.na(sup_material_8_MatchSigC2Region_out$part_number))
# dim(sup_material_8_MatchSigC2Region_out)

save(sup_material_8_MatchSigC2Region_out,file = '/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/matchC_with_exon/man/sup_material_8_MatchSigC2Region.rda')
library(openxlsx)
write.xlsx(sup_material_8_MatchSigC2Region_out,row.names = F,
           file = '/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/matchC_with_exon/man/sup_material_8_MatchSigC2Region.xlsx')
