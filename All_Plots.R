####################################
################ DEU ###############
####################################

suppressPackageStartupMessages(library("DEXSeq"))
library(tidyverse);library(readxl)
sample_index =read_excel("/blue/mateescu/lihe.liu/Methylation_WGCNA/Samples_RNA-Seq.xlsx")
control_index = dplyr::filter(sample_index,TRT == "a") %>% dplyr::select('Tube ID') %>% unlist(use.names = F)
treatment_index = dplyr::filter(sample_index,TRT == "b") %>%  dplyr::select('Tube ID') %>% unlist(use.names = F)
control_index_final = control_index[which(!(control_index %in% c('6268')))]
treatment_index_final = c(treatment_index[which(!(treatment_index %in% c('6228')))],'6228-1')
#conditions<-factor(c(rep("Control",9),rep("Meth",10)),levels=c("Control","Meth"))
# examples in: https://genomicsclass.github.io/book/pages/rnaseq_exon_usage.html

### read in count files
inDir = '/blue/mateescu/lihe.liu/AltSplicing/DEXSeqPyScripts/dexseq_count/'
countfiles = paste(inDir,paste("M",c(control_index_final,treatment_index_final),".count.txt",sep = ''),
                   sep='')
basename(countfiles)
flattenedFile = paste0('/blue/mateescu/lihe.liu/AltSplicing/ARS-UCD1.2/annotation/',
                       'flattened.clean.dexseq.gtf',
                       sep = '')
basename(flattenedFile)

# sample details
sampleTable = data.frame(
  row.names = c( paste0("M",c(control_index_final)),
                 paste0("M",c(treatment_index_final))),
  condition = c(rep('Control',length(control_index_final)),
                rep('Meth',length(treatment_index_final))))
#,libType = rep('paired-end',length(c(control_index_final,treatment_index_final))))

# create test object
dxd_raw = DEXSeqDataSetFromHTSeq(
  countfiles,
  sampleData = sampleTable,
  design = ~ sample + exon + condition:exon,
  flattenedfile = flattenedFile
)

### filter-out low expr gene
all_gene = read.csv(countfiles[1],header = F,sep = '\t') %>%
  dplyr::rename(Exon = V1, Count = V2) %>% 
  mutate(Exon = sub("\\:(.*)", "", Exon)) %>% 
  group_by(Exon) %>% 
  mutate(SampleSum = sum(Count)) %>% 
  slice_head() %>% 
  dplyr::select(-Count) %>% 
  mutate(Keep = ifelse(SampleSum >= 1, 1, 0)) %>% 
  dplyr::select(-SampleSum) 


#tmp = c()
for (i in c(2: length(basename(countfiles)))){
  all_gene_tmp = read.csv(countfiles[i],header = F,sep = '\t') %>% 
    dplyr::rename(Exon = V1, Count = V2) %>% 
    mutate(Exon = sub("\\:(.*)", "", Exon)) %>% 
    group_by(Exon) %>% 
    mutate(SampleSum = sum(Count)) %>% 
    slice_head() %>% 
    dplyr::select(-Count) %>% 
    mutate(Keep = ifelse(SampleSum >= 1, 1, 0)) %>% 
    dplyr::select(-SampleSum) 
  #tmp[i] = all_gene_tmp[6,2]
  all_gene = all_gene_tmp %>% 
    left_join(all_gene, by=c('Exon' = 'Exon'))
}
#unlist(tmp,use.names = F)
#cm[1,]
head(all_gene)
# rm ambiguous / empty 
all_gene_out = all_gene[-c(1:5),]
head(all_gene_out)

# keep genes expressed in > 9 replicates
gene2keep = 
  all_gene_out[(which(rowSums(all_gene_out[,2:ncol(all_gene_out)])>9)),1] %>% 
  unlist(use.names = F)
length(gene2keep) # 12056: previously
#save(gene2keep,file = 'gene2keep.rda')
# keep only genes with high expr
dxd_filter = dxd_raw[geneIDs(dxd_raw) %in% gene2keep,]
#dim(dxd_raw)
#dim(dxd_filter)

colData(dxd_filter)

setwd('/blue/mateescu/lihe.liu/AltSplicing/AltSplicing-R')
save(dxd_filter,file = 'DEXSeq_dxd_filter_4count.rda')

# test if genes kept overlap correctly
#load('DEXSeq_out_all.rda')
#DEXSeq_final_test = data.frame(dxd_DEU_DEXSeq_results)
#geneskept = unique(DEXSeq_final_test$groupID)

###########################################
##############   DEXSeq   ################
############################################

### Normalisation and Dispersion estimation
dxd_out1 = estimateSizeFactors(dxd_filter)
dxd_out2 = estimateDispersions(dxd_out1)
save(dxd_out1,dxd_out2,file = 'DEXSeq_process.rda')

#plotDispEsts(dxd)

### models
fullModel<- ~ sample + exon + condition:exon 
reducedModel<- ~ sample + exon

###  test
dxd_DEU_DEXSeq = testForDEU(
  dxd_out2,
  fullModel = fullModel,
  reducedModel = reducedModel,
  BPPARAM = MulticoreParam(workers=20))

dxd_DEU_DEXSeq_final = estimateExonFoldChanges(dxd_DEU_DEXSeq,fitExpToVar = "condition") # glm: count ∼ condition + exon + condition:exon.

dxd_DEU_DEXSeq_results = DEXSeqResults(dxd_DEU_DEXSeq_final)

#save(dxd,dxr1,file = 'DEXSeq_out.rda') #save(dxd,dxr1,file = 'DEXSeq_out2.rda')
save(dxd_DEU_DEXSeq, dxd_DEU_DEXSeq_final,dxd_DEU_DEXSeq_results,file = 'DEXSeq_out_all.rda')

## PLOT1 - DEU  - Volcano
load('DEXSeq_out_all.rda')
names(DEXSeq_final)

# get gene name information and join list 
library(biomaRt)
biomart="ensembl";dataset="btaurus_gene_ensembl";attributes = c("ensembl_gene_id","external_gene_name")
database = useMart(biomart)
genome = useDataset(dataset, mart = database)
gene = getBM(attributes,mart = genome)
library(tidyverse)
df_DEU_Voc_plt = DEXSeq_final %>% 
  dplyr::select(c(1,2,6,7,10)) %>% 
  dplyr::left_join(gene, by = c('groupID' = 'ensembl_gene_id')) %>% 
  relocate(external_gene_name,.after = groupID) %>% 
  group_by(groupID) %>% 
  drop_na() %>% 
  mutate(count_sig_exon = as.numeric(table(pvalue <= 0.009)[2])) %>% 
  replace_na(list(count_sig_exon=0))

table(df_DEU_Voc_plt$count_sig_exon)
#
library(gplots)
choose_color = col2hex(c('grey','blue','brown','orange','pink','wheat','tan','red'))
match_color = data.frame(count = c(0:7),
                         color = choose_color)
#
df_DEU_Voc_plt_color = df_DEU_Voc_plt %>% 
  dplyr::left_join(match_color,by = c('count_sig_exon'='count')) %>% 
  #mutate(plot_label = ifelse(pvalue <= 0.009,external_gene_name,NA)) %>% 
  mutate(plot_label = ifelse(pvalue<=0.0001,
                             ifelse(external_gene_name == "",groupID,external_gene_name)
                             ,'test')) %>% 
  mutate(plot_color = ifelse(pvalue<=0.01,color,col2hex('grey'))) %>% 
  mutate(plot_color = factor(plot_color,levels = choose_color))
#
names(df_DEU_Voc_plt)
library(ggrepel)
number_ticks <- function(n) {function(limits) pretty(limits, n)}
Voc_DEU = 
ggplot(data = df_DEU_Voc_plt_color,
       aes(x = log2fold_Meth_Control,
           y = -log10(pvalue),
           #color= plot_color,
           label = plot_label)
       ) + 
  geom_point(size = 2,aes(color = plot_color)) + 
  #scale_color_identity(guide = "legend")+
  scale_color_manual(labels = c(0:7),
                     values = choose_color,
                     guide = T)+
  xlim(-30.1,30)+
  labs(y="-log10 (PValue)", x = "log2 (Fold Change)")+
  #scale_x_continuous(breaks = number_ticks(10))+
  #scale_fill_manual(values = c('grey','blue','brown','orange','pink','wheat','tan','red'))+
  geom_text_repel(size =3.5,#fontface = 'bold',
                  min.segment.length = unit(0, 'lines'),
                  segment.alpha = 0.7,
                  box.padding = 1, max.overlaps = Inf) +
  theme(legend.position = "bottom",
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=12,face = "bold"),
        legend.key=element_blank(),
        axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))+
  guides(color = guide_legend("Number of Significant Exons",nrow = 1))
Voc_DEU

tiff("Plots/Fig1b-Voc_DEU.tiff", width = 14, height = 14, units = 'in', res = 500)
print(Voc_DEU)
dev.off()


####################################
################ DIE ###############
####################################
### pre
setwd('../Iso_Expression/')
iso_files = list.files()[-1]

setwd('../counts/iso_counts/')
iso_files = list.files()
iso_res<-lapply(iso_files, function(iso_f){
  return(read.delim(iso_f))
})

names(iso_res)<-sapply(iso_files,function(iso_f){
  return(strsplit(iso_f, split="[.]")[[1]][1])
})

iso_cm<-as.data.frame(do.call(cbind, lapply(1:length(iso_res), function(x){
  iso_sample<-iso_res[[x]]
  return(iso_sample$expected_count)
})))

colnames(iso_cm) = names(iso_res)
rownames(iso_cm)<-iso_res[[1]]$transcript_id
length(unique(rownames(iso_cm))) # 43512 transcripts
# index
control_index_final = c('M6125','M6126','M6133','M6157','M6226','M6254','M6472','M6482','M6520')
treatment_index_final = c("M6152","M6153","M6197","M6412","M6415","M6419","M6438","M6453","M6485","M6228-1")
# order to 9-1o pattern
iso_cm_order = cbind(iso_cm[,(colnames(iso_cm)%in%control_index_final)],
                     iso_cm[,!(colnames(iso_cm)%in%control_index_final)])

### DIE analysis
conditions<-factor(c(rep("Control",9),
                     rep("Meth",10)),
                   levels=c("Control","Meth"))
# Consider those isoforms expressed in at least one condition
#iso_cm_filter = iso_cm_order[(rowSums(iso_cm_order[,1:9]) >= 0 | rowSums(iso_cm_order[,10:19]) >= 0),]
to_keep = c()
for (i in seq_along(rownames(iso_cm_order))){
  print(i)
  tmp = iso_cm_order[i,]
  if (sum(tmp >= 1) >= 9){
    to_keep = append(to_keep,i)
  }
} #  19701

iso_cm_filter = iso_cm_order[to_keep,]

library(DESeq2)
## DESeq2
conditions<-factor(c(rep("Control",9),
                     rep("Meth",10)),
                   levels=c("Control","Meth"))

sample_data<-data.frame(condition=conditions)
rownames(sample_data) = names(iso_cm_filter)

dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(iso_cm_filter)),
                              colData = sample_data, 
                              design = formula(~condition,sample_data))
#dds

# y<-estimateSizeFactors(y)
# y<-estimateDispersions(y)
# et_y_DES<-nbinomWaldTest(y)
out = DESeq(dds)
#out
#test = rowData(out)
#head(test)
#test= counts(dds) # original counts

DESeqres<-results(out,contrast = c("condition","Meth","Control"))

# Example from the manual
# dds <- makeExampleDESeqDataSet(m=4)
# dds <- DESeq(dds)
# res <- results(dds, contrast=c("condition","B","A"))
head(DESeqres)

#isoDESeq<-rownames(DESeqres[DESeqres$padj < 0.05 & !is.na(DESeqres$padj),])

save(DESeqres,out,file="DESeqRes1217.RData")

# local 
load('DIE_map_gene.rda')
table(DESeqres$pvalue <= 0.009) # 175

library(tidyverse)
df_DIE_Voc_plt = DIE_map_gene %>% 
  distinct(TranscriptID,.keep_all = T) %>% 
  dplyr::select(c(1,2,4,7)) %>% 
  dplyr::left_join(gene, by = c('gene_id' = 'ensembl_gene_id')) %>% 
  relocate(external_gene_name,.after = gene_id) %>% 
  mutate(plot_color = ifelse(pvalue<=0.01,col2hex('red'),col2hex('grey'))) %>% 
  mutate(plot_color = factor(plot_color,levels = c(col2hex('red'),col2hex('grey')))) %>% 
  mutate(plot_label = ifelse((pvalue <= 10^-2.5 | (pvalue <= 0.01 & abs(log2FoldChange)>= 0.8 )), ifelse(external_gene_name == "",gene_id,external_gene_name),NA))
  #mutate(plot_label2 = ifelse((pvalue <= 0.009 & is.na(plot_label)),gene_id,plot_label)) %>% 

head(df_DIE_Voc_plt)


#(df_DIE_Voc_plt[df_DIE_Voc_plt$log2FoldChange == max(df_DIE_Voc_plt$log2FoldChange),7])

library(ggrepel)
Voc_DIE = 
  ggplot(data = df_DIE_Voc_plt,
         aes(x = log2FoldChange,
             y = -log10(pvalue),
             #color= plot_color,
             label = plot_label)) + 
  geom_point(size = 2,aes(color = plot_color)) + 
  xlim(-5.6, 5.6) +
  scale_color_manual(values = c(col2hex('red'),col2hex('grey')),
                     guide = F)+
  theme(legend.position = "none")+
  labs(y="-log10 (PValue)", x = "log2 (Fold Change)")+
  geom_text_repel(size = 3.5,#fontface = 'bold',
                  min.segment.length = unit(0, 'lines'),
                  segment.alpha = 0.5,
                  box.padding =  .8, max.overlaps = 20)+
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))

#Voc_DIE

tiff("Plots/Fig1a-Voc_DIE.tiff", width = 14, height = 14, units = 'in', res = 500)
print(Voc_DIE)
dev.off()


#######################################
################ Proportion ###########
#####################################
# diff Meth Cs proportions are lower in Sig differentially used exons;
# Same trendin introns surrounding a differentially used exon

## skiped some lines

###
data = data.frame(
  value=c(exon_prop_sig_overall,
          exon_prop_nonsig_overall),
  name =c(rep('Significant',length(exon_prop_sig_overall)),
          rep('Non-significant',length(exon_prop_nonsig_overall))),
  type =c(rep('Exons',length(c(exon_prop_sig_overall,exon_prop_nonsig_overall)))))

ks.test(exon_prop_sig_ave,
        exon_prop_nonsig_ave)

#####
data2 <- data.frame(
  value=c(intron_prop_surrd_overall,
          intron_prop_other_overall),
  name =c(rep('Contiguous',length(intron_prop_surrd_overall)),
          rep('Non-contiguous',length(intron_prop_other_overall))),
  type =c(rep('Introns',length(c(intron_prop_surrd_overall,intron_prop_other_overall)))))


ks.test(intron_prop_surrd_overall,
        intron_prop_other_overall)


df_meth_prop = data %>% 
  bind_rows(data2) %>% 
  mutate(Significance = ifelse(name %in% c('Significant','Contiguous'),'DU','Non-DU')) %>% 
  mutate_at('Significance',as.factor) %>% 
  mutate_at('name',as.factor) %>% 
  mutate_at('type',as.factor) %>% 
  dplyr::filter(value<= 0.6)

df_meth_prop$name = factor(df_meth_prop$name,
                           levels = c('Significant','Non-significant','Contiguous','Non-contiguous'))

Violin_Prop = 
  ggplot(df_meth_prop,
       aes(y=value, x=name, fill=Significance)) +
  geom_violin() + theme(legend.position="bottom")+
  facet_wrap(~type,scales = "free")+
  ylab('Methylation Level')+
  xlab('')+
  guides(fill=guide_legend(title=""))+
  scale_fill_manual(values = c(col2hex('red'),col2hex('blue')))+
  theme(legend.position='none',
        legend.text=element_text(size=10),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.background = element_rect(color="black",linetype="solid"),
        legend.key.size = unit(.8,"line"),
        axis.text.x= element_text(color = 'black',face="bold", size = 12),
        axis.text.y= element_text(color = 'black',face="bold", size = 8))
Violin_Prop

tiff("Plots/Fig4-Violin-Prop.tiff", width = 12, height = 8, units = 'in', res = 500)
print(Violin_Prop)
dev.off()

#######################################
################ Enrichment ###########
######################################
library(openxlsx)
library(tidyverse)
plot_DEU = read.xlsx('./Enrichment_DIE_DEU/Enrichment_Summary_2plot_DEU.xlsx', sheet = 1) %>% 
  dplyr::filter(Include == 1) %>% 
  mutate(Tag = paste0(ID,'[',totalG,']')) %>% 
  relocate(Tag, .after = ID) %>% 
  mutate(Group = 'DEU')

plot_DIE = read.xlsx('./Enrichment_DIE_DEU/Enrichment_Summary_2plot_DIE.xlsx', sheet = 1) %>% 
  dplyr::filter(Include == 1) %>% 
  mutate(Tag = paste0(ID,'[',totalG,']')) %>% 
  relocate(Tag, .after = ID) %>% 
  mutate(Group = 'DIE')
names(plot_DIE)

plot_all = plot_DEU %>% 
  bind_rows(plot_DIE) %>% 
  mutate(Category = as.factor(Category)) %>% 
  mutate(log10pvalue = -log10(pvalue)) %>% 
  mutate(hitsPerc = - hitsPerc) %>% 
  group_by(Category)

View(plot_all[,c(2,9)])


plot_all$Tag =
  factor(plot_all$Tag,
         levels = rev(c(
           "GO:0005524[920]","GO:0016887[129]",#atp
           "GO:0030155[21]","GO:0045785[28]","GO:0051693[11]","bta04810[148]","IPR001589[18]",#Adherens junctions and Actin Cytoskeleton
           "GO:0030955[7]","GO:0036376[6]","M25065[8]",#Ion Binding and Transportation
           "R-BTA-8957322[25]","GO:0033539[6]","GO:0033211[5]",#Muscle Metabolism
           "M5909[192]","M18391[38]","GO:0003009[17]", #Muscle Function and Development
           "bta04512[57]","M27218[41]", #ECM-receptor Interaction
           "D019175[33]", #DNA METH
           "M8084[50]",# RNA SPLICING
           "GO:0046034[22]","GO:0042623[14]","M784[22]", #ATP
           "GO:0008270[538]","IPR036875[20]",#Ion Binding and Transportation
           "M26014[15]","M13561[17]","M15075[103]", #Postsynaptic Membrane
           "M2044[117]","D017398[1409]","M16965[19]",#RNA Splicing
           "GO:0030177[26]","bta04520[68]","GO:0032956[58]"#Adherens junctions and Actin Cytoskeleton
         )))


plot_all$Category =
  factor(plot_all$Category,
         levels = (c(
           "ATP",#atp
           "Adherens junctions and Actin Cytoskeleton",#
           "Ion Binding and Transportation",#
           "Muscle Metabolism Function and Development",#
           "ECM-receptor Interaction", #
           "DNA Methylation", #
           "RNA Splicing",# 
           #"Muscle Function and Development", #
           "Postsynaptic Membrane"
         )))


ggEnrich = 
  ggplot(plot_all) + coord_flip()+
  geom_bar(mapping = aes(x = Tag, y = hitsPerc,fill = Category),#
           stat = "identity", width = 0.2) +
  geom_point(mapping = aes(x = Tag, y = 10 * log10pvalue),size = 1) +
  #scale_y_continuous() 
  scale_y_continuous(name = expression(bold("Percentage of Significant Genes")),
                     sec.axis = sec_axis(~. * 0.1, breaks = c(0,5,10), 
                                         name = expression(bold("-log10(P-value)")))
                     ,limits = c(0,100), breaks = seq(0,100,by=5)
                     ,expand = c(0,1)) +
  guides(fill=guide_legend(title="",
                           keywidth = 1, nrow = 2, 
                           keyheight = 1))+
  theme(legend.position = "bottom",
        legend.text=element_text(size=10,face="bold",color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=8,face="bold",color = "black"),
        axis.text.x = element_text(size=10,face="bold",color = "black")) +
  #facet_wrap( ~ Group, nrow = 2 ,scales = "free",)
  facet_grid(Group~.,scales = "free_y",switch = 'y')


tiff("Plots/Fig2-Enrichment.tiff", width = 16, height = 8, units = 'in', res = 500)
#tiff("Fig3-Enrichment-Bar.tiff", width = 10, height = 6, units = 'in', res = 500)
print(ggEnrich)
dev.off()




