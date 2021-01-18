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

# here we have prop information
head(aggrate_by_gene)

# here we have significance information
out = DEXSeq_final
out_reduce = out %>% 
  as_tibble() %>%
  group_by(groupID) %>% 
  dplyr::select(-c(3:5,7:35))
head(out_reduce)


load('all_pval1.rda')
load('all_pval2.rda')
load('all_pval3.rda')
load('all_pval4.rda')
load('all_pval5.rda')
all_pval = c(all_pval1,all_pval2,all_pval3,all_pval4,all_pval5)

aggrate_by_gene_withp = aggrate_by_gene %>% 
  data.frame() %>% 
  mutate(pvalue = all_pval) %>% 
  as_tibble() %>% 
  relocate(pvalue,.after = featureID)

save(aggrate_by_gene_withp,file = 'aggrate_by_gene_withp.rda')



# 1
all_pval1 = c()
all_geneid = (aggrate_by_gene %>% pull(groupID))[1:50000]
all_featureid = (aggrate_by_gene %>% pull(featureID))[1:50000]

#
for (i in seq_along(all_geneid)){
  if (i %% 100 == 0){print(i)}
  if (startsWith(all_featureid[i],'I')){
    all_pval1[i] = 1
  } else {
    all_pval1[i] = out_reduce %>% 
      dplyr::filter(groupID == all_geneid[i] & featureID == all_featureid[i]) %>% 
      pull(pvalue)
  }
}
save(all_pval1,file = '/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/combine_tmp/all_pval1.rda')



# 2
all_pval2 = c()
all_geneid = (aggrate_by_gene %>% pull(groupID))[50001:100000]
all_featureid = (aggrate_by_gene %>% pull(featureID))[50001:100000]

#
for (i in seq_along(all_geneid)){
  if (i %% 100 == 0){print(i)}
  if (startsWith(all_featureid[i],'I')){
    all_pval2[i] = 1
  } else {
    all_pval2[i] = out_reduce %>% 
      dplyr::filter(groupID == all_geneid[i] & featureID == all_featureid[i]) %>% 
      pull(pvalue)
  }
}
save(all_pval2,file = '/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/combine_tmp/all_pval2.rda')


# 3
all_pval3 = c()
all_geneid = (aggrate_by_gene %>% pull(groupID))[100001:150000]
all_featureid = (aggrate_by_gene %>% pull(featureID))[100001:150000]

#
for (i in seq_along(all_geneid)){
  if (i %% 100 == 0){print(i)}
  if (startsWith(all_featureid[i],'I')){
    all_pval3[i] = 1
  } else {
    all_pval3[i] = out_reduce %>% 
      dplyr::filter(groupID == all_geneid[i] & featureID == all_featureid[i]) %>% 
      pull(pvalue)
  }
}
save(all_pval3,file = '/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/combine_tmp/all_pval3.rda')


# 4
all_pval4 = c()
all_geneid = (aggrate_by_gene %>% pull(groupID))[150001:200000]
all_featureid = (aggrate_by_gene %>% pull(featureID))[150001:200000]

#
for (i in seq_along(all_geneid)){
  if (i %% 100 == 0){print(i)}
  if (startsWith(all_featureid[i],'I')){
    all_pval4[i] = 1
  } else {
    all_pval4[i] = out_reduce %>% 
      dplyr::filter(groupID == all_geneid[i] & featureID == all_featureid[i]) %>% 
      pull(pvalue)
  }
}
save(all_pval4,file = '/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/combine_tmp/all_pval4.rda')


# 5
all_pval5 = c()
all_geneid = (aggrate_by_gene %>% pull(groupID))[200001:271710]
all_featureid = (aggrate_by_gene %>% pull(featureID))[200001:271710]

#
for (i in seq_along(all_geneid)){
  if (i %% 100 == 0){print(i)}
  if (startsWith(all_featureid[i],'I')){
    all_pval5[i] = 1
  } else {
    all_pval5[i] = out_reduce %>% 
      dplyr::filter(groupID == all_geneid[i] & featureID == all_featureid[i]) %>% 
      pull(pvalue)
  }
}
save(all_pval5,file = '/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R/meth_prop/combine_tmp/all_pval5.rda')



