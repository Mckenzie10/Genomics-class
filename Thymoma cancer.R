x=stats::runif(10)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("devtools")
devtools::install_github("mariodeng/FirebrowseR")
require(FirebrowseR)
mRNA.Exp = Samples.mRNASeq(format = "csv", gene = c("PTEN", "RUNX1"),tcga_participant_barcode = c("TCGA-GF-A4EO","TCGA-AC-A2FG"))
mRNA.Exp[, c("tcga_participant_barcode", "expression_log2", "z.score")]
cohorts = Metadata.Cohorts(format = "csv")
cancer.Type = cohorts[grep("thymoma", cohorts$description, ignore.case = T), 1]
print(cancer.Type)
THYM.Pats = Samples.Clinical(cohort = cancer.Type, format="tsv")
dim(THYM.Pats)
all.Received = F
page.Counter = 1
page.size = 150
THYM.Pats = list()
while(all.Received == F){
  THYM.Pats[[page.Counter]] = Samples.Clinical(format = "csv",cohort = cancer.Type, page_size = page.size,page = page.Counter)                                               
  if(page.Counter > 1)  
    colnames(THYM.Pats[[page.Counter]]) = colnames(THYM.Pats[[page.Counter-1]]    
                                                   if(nrow(THYM.Pats[[page.Counter]])<page.size){all.Received = T} else{page.Counter=page.Counter+1}   
                                                   THYM.Pats = do.call(rbind, THYM.Pats)
                                                   dim(THYM.Pats)
                                                   THYM.Pats = THYM.Pats[ which(THYM.Pats$vital_status == "dead"), ]
                                                   diff.Exp.Genes = c("ESR1", "GATA3", "XBP1", "FOXA1", "ERBB2", "GRB7", "EGFR","FOXC1", "MYC")                                                     
                                                   all.Found = F
                                                   page.Counter = 1
                                                   mRNA.Exp = list()
                                                   page.Size = 2000
                                                   while(all.Found == F){mRNA.Exp[[page.Counter]] = Samples.mRNASeq(format = "csv",mRNA.Exp[[page.Counter]] = Samples.mRNASeq(format = "csv", mRNA.Exp[[page.Counter]] = Samples.mRNASeq(format = "csv",gene = diff.Exp.Genes,cohort = "THYM",tcga_participant_barcode =tcga_participant_barcode =THYM.Pats$tcga_participant_barcode,page_size = page.Size,page_size = page.Size,page = page.Counter)