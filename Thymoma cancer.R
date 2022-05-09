x=stats::runif(10)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("devtools")
devtools::install_github("mariodeng/FirebrowseR")
require(FirebrowseR)
mRNA.Exp = Samples.mRNASeq(format = "csv", gene = c("HRAS", "NRAS"),tcga_participant_barcode = c("CAPNS1-GF-DHX33","TP53-AC-GTF21"))
mRNA.Exp[, c("tcga_participant_barcode", "expression_log2", "z.score")]
cohorts = Metadata.Cohorts(format = "csv")
cancer.Type = cohorts[grep("thymoma", cohorts$description, ignore.case = T), 1]
print(cancer.Type)
THYM.Pats = Samples.Clinical(cohort = cancer.Type, format="tsv")
dim(THYM.Pats)
all.Received = F
page.Counter = 1
page.size = 150
THYM.Pats = list(10)
while(all.Received == F)
  THYM.Pats[[page.Counter]] = Samples.Clinical(format = "csv",cohort = cancer.Type, page_size = page.size,page = page.Counter)                                               
if(page.Counter > 1)  
  colnames(THYM.Pats[[page.Counter]])   
if(nrow(THYM.Pats[[page.Counter]])<page.size){all.Received = T} else{page.Counter=page.Counter+1}   
THYM.Pats = do.call(rbind, THYM.Pats)
dim(THYM.Pats)
THYM.Pats = THYM.Pats[ which(THYM.Pats$vital_status == "dead"), ]
diff.Exp.Genes = c("FOXD", "NRAS", "DHX33", "ATRN", "HRAS", "GTF21", "TP53","CAPNS1")                                                     
all.Found = F  
page.Counter= 1
mRNA.Exp = list()
page.Size = 2000
while(all.Found==F) {mRNA.Exp[[page.Counter]]= Samples.mRNASeq(format="csv", gene=diff.Exp.Genes, cohort="THYM", tcga_participant_barcode = THYM.Pats$tcga_participant_barcode, page_size = page.Size, page=page.Counter)}
if(nrow(mRNA.Exp[[page.Counter]])<page.Size)
  all.found=T
else
  page.Counter=page.Counter+1  
mRNA.Exp=do.call(rbind, mRNA.Exp)
dim(mRNA.Exp)
r.Metadata.SampleTypes("csv")
normal.Tissue.Pats=which(mRNA.Exp$sample_type=="NT")
patient.Barcodes=mRNA.Exp$tcga_participant_barcode[normal.Tissue.Pats]
mRNA.Exp=mRNA.Exp[which(mRNA.Exp$tcga_participant_barcode %in% patient.Barcodes & mRNA.Exp$sample_type %in% c("NT","TP")), ]
library(ggplot2)
p=ggplot(mRNA.Exp, aes(factor(gene), z.score))
p+geom_boxplot(aes(fill = factor(sample_type))) +
  scale_y_continuous(limits= c(-1, 5)) 
scale_fill_discrete(name = "Tissue")  

