library(tidyverse)
# local plotting relationship
load(
  'DEXSeq_final.rda'
)
load(
  'integrate_by_gene_withCproportion_pval0.01.rda'
)
dim(integrate_by_gene_out)

Universal_exon_intron_transcript_final_new = 
  cbind(Universal_exon_intron_transcript_final[,c(1:3)],
        integrate_by_gene_out[,c(4:10)]) %>% 
  as_tibble() %>% 
  mutate(prop_p.01 = count_sig_p.01/count_all)
head(Universal_exon_intron_transcript_final_new)

#
pthres = 9e-03
#
# gene_list = Universal_exon_intron_transcript_final %>% 
#   #group_by(groupID) %>% 
#   #mutate(sig_count = sum(pvalue <= pthres)) %>% 
#   #dplyr::filter(sig_count == 5)%>% 
#   dplyr::select(groupID) %>% unlist(use.names = F) %>% 
#   unique()
# length(gene_list)
# 1 158 - 2 19 - 3 2 - 4 1

test = Universal_exon_intron_transcript_final %>% 
  dplyr::filter(pvalue <= pthres) # 1745 exons
dim(test)

gene_list = Universal_exon_intron_transcript_final_new %>% 
  dplyr::filter(pvalue <= pthres)%>% 
  dplyr::select(groupID) %>% 
  unlist(use.names = F)
length(table((gene_list)))

#
most_adjacent = function(node,input_vector,select = 2){
  # node = 3
  # input_vector = c(1,2,5,7,8,9)
  if (is.na(node)){return(integer(0))}
  if (length(input_vector) == 1){
    return(input_vector)
  } else {
    tmp = (input_vector - node)
    if (sum(tmp > 0) == length(tmp)){
      return(input_vector[1])
    } else if (sum(tmp < 0) == length(tmp)){
      return(input_vector[length(input_vector)])
    } else {
      tmpout1 = which(tmp<0)
      out1 = tmpout1[length(tmpout1)]
      out2 = which(tmp>0)[1]
      return(c(input_vector[c(out1,out2)]))
    }
  }
}

# compare distr for exons: 
exon_prop_sig_overall = c()
exon_prop_sig_ave = c()
exon_prop_nonsig_overall = c()
exon_prop_nonsig_ave = c()

# compare distr for introns: 
intron_prop_surrd_overall = c()
intron_prop_surrd_ave = c()
intron_prop_other_overall = c()
intron_prop_other_ave = c()

gene_list_uniq = unique(gene_list)
###
for (i in seq_along(gene_list_uniq)){
  #i = 63
  tmp_file = Universal_exon_intron_transcript_final_new %>% 
    dplyr::filter(groupID == gene_list_uniq[i]) %>% 
    drop_na(pvalue)
  print(i)
  # differentiate intro and extron + rm na
  raw_index = c(1:nrow(tmp_file))
  Intro_dex = raw_index[str_detect(tmp_file$featureID,'I')]
  Intro_dex = Intro_dex[!is.na(Intro_dex)]
  Exon_dex = raw_index[str_detect(tmp_file$featureID,'E')]
  Exon_dex = Exon_dex[!is.na(Exon_dex)]
  
  # pull out temp columns
  col_total = tmp_file %>% pull(count_all)
  col_sig = tmp_file %>% pull(count_sig_p.01)
  col_prop = tmp_file %>% pull(prop_p.01)

  # EXONs
  # put threshold
  target = which(tmp_file$pvalue <= pthres)
  off_target = Exon_dex[!(Exon_dex %in% tareget)]
  
  # calculate props for "other"
  tmp_sig_prop_overall = sum(col_sig[target]) / sum(col_total[target])
  
  tmp_sig_prop_ave = mean(col_prop[target]/col_prop[target])
  
  tmp_other_overall_prop = sum(col_sig[off_target]) / sum(col_total[off_target])
  
  tmp_other_prop_ave = mean(col_prop[off_target]/col_prop[off_target])
  
  # append
  exon_prop_sig_overall = append(exon_prop_sig_overall,tmp_sig_prop_overall)
  
  exon_prop_sig_ave = append(exon_prop_sig_ave,tmp_sig_prop_ave)
  
  exon_prop_nonsig_overall = append(exon_prop_nonsig_overall,tmp_other_overall_prop)
  
  exon_prop_nonsig_ave = append(exon_prop_nonsig_ave,tmp_other_prop_ave)
  
  # INTRONs
  select_intron_index_tmp = c()
  for (m in seq_along(target)){
    tmp = most_adjacent(tareget[m],Intro_dex)
    select_intron_index_tmp = append(select_intron_index_tmp,tmp)
  }
  select_intron_index = unique(select_intron_index_tmp)
  drop_intron_index = Intro_dex[!(Intro_dex %in% select_intron_index)]
  #
  intron_prop_surrd_overall = append(intron_prop_surrd_overall,
                                     sum(col_sig[select_intron_index])/sum(col_total[select_intron_index]))
  
  intron_prop_surrd_ave = append(intron_prop_surrd_ave,
                                 mean(col_prop[select_intron_index]/col_prop[select_intron_index]))
  
  #
  intron_prop_other_overall = append(intron_prop_other_overall,
                                     sum(col_sig[drop_intron_index])/sum(col_total[drop_intron_index]))
  
  intron_prop_other_ave = append(intron_prop_other_ave,
                                 mean(col_prop[drop_intron_index]/col_prop[drop_intron_index]))
}

###
data <- data.frame(
  value=c(exon_prop_sig_overall,
          exon_prop_nonsig_overall),
  name =c(rep('sig exon - pval 0.01',length(exon_prop_sig_overall)),
          rep('non sig exon - pval 0.0',length(exon_prop_nonsig_overall))))

p1 =
  ggplot(data, aes(y=name, x=value, fill=name)) +
  geom_violin()+theme(legend.position="bottom")+
  coord_flip() 
  
p1

ks.test(exon_prop_sig_ave,
        exon_prop_nonsig_ave)

#####
data2 <- data.frame(
  value=c(intron_prop_surrd_overall,
          intron_prop_other_overall),
  name =c(rep('intron_prop_surrd_pval 0.01',length(intron_prop_surrd_overall)),
          rep('intron_prop_other_pval 0.01',length(intron_prop_other_overall))))

p2 = 
  ggplot(data2, aes(y=name, x=value, fill=name)) +
  geom_violin()+theme(legend.position="bottom")+
  coord_flip()


ks.test(intron_prop_surrd_overall,
        intron_prop_other_overall)




