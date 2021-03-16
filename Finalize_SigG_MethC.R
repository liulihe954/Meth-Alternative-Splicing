# library(openxlsx)
# write.xlsx(Universal_DEU_info_sup_file,row.names = F,
#            file = '/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/Manuscript/Supplementary_Files/Universal_DEU_info.xlsx')


# library(openxlsx)
# write.xlsx(out_format_addE,row.names = F,
#            file = '/Users/liulihe95/Desktop/Isoform-Expression/AltSplicing-R/Manuscript/Supplementary_Files/DIE_Meth.xlsx')


# make sure every DIE has meth record
dim(out_format_addE)

DIE_massage = out_format_addE %>% 
  group_by(gene_id) %>% 
  mutate(CountAll_sum = sum(count_all)) %>% 
  mutate(include = ifelse(CountAll_sum==0,0,1)) %>% 
  dplyr::select(gene_id,include) %>% 
  distinct(gene_id,.keep_all = T)

Gene_include_DIE = DIE_massage %>% 
  dplyr::filter(include == 1) %>% 
  pull(gene_id)
length(Gene_include_DIE)

Gene_exclude_DIE = DIE_massage %>% 
  dplyr::filter(include != 1) %>% 
  pull(gene_id)
length(Gene_exclude_DIE)

gene_list = out_format_addE %>% 
  group_by(gene_id) %>% 
  dplyr::filter(pvalue <= 0.009) %>% 
  dplyr::select(gene_id) %>% 
  unique() %>% unlist(use.names = F)

##
#names(out_format_addE)
DIE_massage_prom = out_format_addE %>% 
  group_by(gene_id) %>% 
  mutate(CountAll_sum_prom = sum(count_all_prom)) %>% 
  mutate(include = ifelse(CountAll_sum_prom == 0,0,1)) %>% 
  dplyr::select(gene_id,include) %>% 
  distinct(gene_id,.keep_all = T)

Gene_include_DIE_prom = DIE_massage_prom %>% 
  dplyr::filter(include == 1) %>% 
  pull(gene_id)
length(Gene_include_DIE_prom)

Gene_exclude_DIE_prom = DIE_massage_prom %>% 
  dplyr::filter(include != 1) %>% 
  pull(gene_id)
length(Gene_exclude_DIE_prom)

gene_list = out_format_addE %>% 
  group_by(gene_id) %>% 
  dplyr::filter(pvalue <= 0.009) %>% 
  dplyr::select(gene_id) %>% 
  unique() %>% unlist(use.names = F)

table(gene_list %in% Gene_include_DIE_prom)

test_include_gene_die = gene_list[gene_list %in% Gene_include_DIE_prom]
length(test_include_gene_die)


names(out_format_addE)

test_prom = out_format_addE %>% 
  group_by(gene_id) %>% 
  dplyr::filter(count_all_prom != 0) %>% 
  mutate(prom_prop = count_sig_prom/count_all_prom) %>% 
  dplyr::select(gene_id,prom_prop)


test_prom_include = test_prom %>% 
  dplyr::filter(gene_id %in% test_include_gene_die) %>% 
  distinct(gene_id,.keep_all = T) %>% 
  pull(prom_prop)

length(test_prom_include)

test_prom_exclude = test_prom %>% 
  dplyr::filter(!(gene_id %in% test_include_gene_die)) %>% 
  distinct(gene_id,.keep_all = T) %>% 
  pull(prom_prop) 
length(test_prom_exclude)

allp = c()
for (i in seq_len(1000)){
  g1 = test_prom_include
  g2 = test_prom_exclude[sample(1:length(test_prom_exclude),length(test_prom_include))]
  allp[i] = ks.test(g1,g2)$p.value
}
hist(allp)







# make sure every DEU has meth record

length(unique(out_format_addE$transcript_id))
head(Universal_DEU_info_sup_file)

DEU_massage = Universal_DEU_info_sup_file %>% 
  group_by(groupID) %>% 
  mutate(CountAll_sum = sum(count_all)) %>% 
  mutate(include = ifelse(CountAll_sum==0,0,1)) %>% 
  dplyr::select(groupID,include) %>% 
  distinct(groupID,.keep_all = T)

Gene_include_DEU = DEU_massage %>% 
  dplyr::filter(include == 1) %>% 
  pull(groupID)
length(Gene_include_DEU)

Gene_exclude_DEU = DEU_massage %>% 
  dplyr::filter(include != 1) %>% 
  pull(groupID)
length(Gene_exclude_DEU)


gene_list_DEU = Universal_exon_intron_transcript_final_new %>% 
  dplyr::filter(pvalue <= pthres)%>% 
  dplyr::select(groupID) %>% 
  unlist(use.names = F) %>% 
  unique()
table(gene_list_DEU %in% Gene_include_DEU)
## compare distribution of promoter area


#############
#names(Universal_DEU_info_sup_file)
DEU_massage_prom = Universal_DEU_info_sup_file %>% 
  group_by(groupID) %>% 
  mutate(CountAll_sum_prom = sum(count_all_prom)) %>% 
  mutate(include = ifelse(CountAll_sum_prom==0,0,1)) %>% 
  dplyr::select(groupID,include) %>% 
  distinct(groupID,.keep_all = T)

Gene_include_DEU_prom = DEU_massage_prom %>% 
  dplyr::filter(include == 1) %>% 
  pull(groupID)
length(Gene_include_DEU_prom)

Gene_exclude_DEU_prom = DEU_massage_prom %>% 
  dplyr::filter(include != 1) %>% 
  pull(groupID)
length(Gene_exclude_DEU_prom)


gene_list_DEU = Universal_exon_intron_transcript_final_new %>% 
  dplyr::filter(pvalue <= pthres)%>% 
  dplyr::select(groupID) %>% 
  unlist(use.names = F) %>% 
  unique()
length(gene_list_DEU)

table(gene_list_DEU %in% Gene_include_DEU_prom)
test_include_gene = gene_list_DEU[gene_list_DEU %in% Gene_include_DEU_prom]



test_prom = Universal_DEU_info_sup_file %>% 
  dplyr::filter(count_all_prom != 0) %>% 
  mutate(prom_prop = count_sig_prom/count_all_prom) %>% 
  dplyr::select(groupID,prom_prop)


test_prom_include = test_prom %>% 
  dplyr::filter(groupID %in% test_include_gene) %>% 
  distinct(groupID,.keep_all = T) %>% 
  pull(prom_prop)

test_prom_exclude = test_prom %>% 
  dplyr::filter(!(groupID %in% test_include_gene)) %>% 
  distinct(groupID,.keep_all = T) %>% 
  pull(prom_prop) 
  

allp = c()
for (i in seq_len(1000)){
  g1 = test_prom_include
  g2 = test_prom_exclude[sample(1:length(test_prom_exclude),length(test_prom_include))]
  allp[i] = ks.test(g1,g2)$p.value
}
hist(allp)


