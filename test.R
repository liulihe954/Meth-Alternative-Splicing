# example
BiocManager::install("pasilla")
library(pasilla)
inDir = system.file("extdata", package="pasilla")
countFiles = list.files(inDir, pattern="fb.txt$", full.names=TRUE)

basename(countFiles)
flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
basename(flattenedFile)
sampleTable = data.frame(
  row.names = c( "treated1", "treated2", "treated3",
                 "untreated1", "untreated2", "untreated3", "untreated4" ),
  condition = c("knockdown", "knockdown", "knockdown",
                "control", "control", "control", "control" ),
  libType = c( "single-end", "paired-end", "paired-end",
               "single-end", "single-end", "paired-end", "paired-end" ) )

suppressPackageStartupMessages( library( "DEXSeq" ) )
dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

#isexpr <- rowSums(cpm(y.all) > 1) >= 9 # >= 9

genesForSubset = read.table(
  file.path(inDir, "geneIDsinsubset.txt"),
  stringsAsFactors=FALSE)[[1]]

dxd = dxd[geneIDs( dxd ) %in% genesForSubset,]
colData(dxd)


##########################################
BiocManager::install("Rsubread")
library(Rsubread)
library(limma)
library(edgeR)
targets <- readTargets()
celltype <- factor(targets$CellType)
design <- model.matrix(~celltype)


buildindex(basename="chr1",reference="hg19_chr1.fa")

align(index="chr1",readfile1=targets$InputFile,input_format="gzFASTQ",output_format="BAM",
      output_file=targets$OutputFile,unique=TRUE,indels=5)

fc <- featureCounts(files=targets$OutputFile,annot.inbuilt="hg19")
samplecounts = fc$counts
samplegenes = fc$annotation[,c("GeneID","Length")]
sampleout <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])




