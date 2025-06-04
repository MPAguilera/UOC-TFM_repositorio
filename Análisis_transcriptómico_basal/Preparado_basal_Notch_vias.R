################################################################################
###############################- Look for Notch pathways and related - ########
################################################################################

################# - Upload libraries - #########################################

if (require(rstudioapi)==F){install.packages("rstudioapi")}
if (require(readxl)==F){install.packages("readx")}
if (require(openxlsx)==F){install.packages("openxlsx")}
if (require(gridExtra)==F){install.packages("gridExtra")}
if (require(ggplot2)==F){install.packages("ggplot2")}
if (require(paletteer)==F){install.packages("ggplot2")}

library(dbplyr)
library(msigdbr)
library(dplyr)

########################## - WORKING DIRECTORY and environment- ###############################

## Set WD to file location
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

## Create a directory for the results with the currentData
Directory_for_the_day <- paste0("NOTCH_Signalling", Sys.Date())
dir.create(Directory_for_the_day)

## Upload all the possible NOTCH pathways into the environment
msigdb_data <- msigdbr(species = "Homo sapiens")
notch_pathways <- msigdb_data %>%
  filter(grepl("Notch", gs_name, ignore.case = TRUE))
pthways<- unique(notch_pathways$gs_name)

## Look for the VSt file
files = list.files(recursive = T, full.names = T)
dat <- read.xlsx(files[grep ("Input*.*VST*.*xlsx", files)])
dat <- dat[,-1]; 

## Summary file
resumen <- read.xlsx(files[grep ("Input*.*resumen*.*xlsx", files)])

## Make sure that the names and order are correct
colnames(dat) <- gsub ("\\.", "-", colnames(dat)) %>% gsub ("OCI-LY", "OCI-Ly",.)
resumen$SampleID <- gsub ("\\.", "-", resumen$SampleID)
if (!all(sort (colnames(dat)[-(1:4)]) == sort(resumen$SampleID))){stop ("Los nombres no coinciden")}

## Reorder the names
resumen <- resumen[order (resumen$CB103),]
dat <- dat[, c(1:4, match (resumen$SampleID, colnames(dat)))]
if (!all(colnames(dat)[-(1:4)] == resumen$SampleID)){stop ("Los nombres no coinciden")}

########################## - Genes - ###############################

##Genes in correlation
corre <- data.frame(apply(dat[,-c(1:4)], 1, function(gen) cor(gen, resumen$CB103, method = "spearman")))
corre$Gene_symbol <- dat$external_gene_name
corre$ensembl_gene_id <- dat$ensembl_gene_id
corresig <- na.omit(corre[corre[,1] >= 0.65 |corre[,1] <= (-0.65 ), ])

## Look if any genes correlated are present in any notch pathways
genes <- unique(notch_pathways$ensembl_gene)
genes<-intersect(genes, corresig$ensembl_gene_id)
common <- corresig[na.omit(match(genes, corresig$ensembl_gene_id)),]

##Look to the pathways that hve those genes
colnames(common)[1] <- ("rho")
a <-notch_pathways[grep(paste0(common$ensembl_gene_id, collapse = "|"), notch_pathways$db_ensembl_gene),]
for (i in 1:nrow(common)){
  col<- colnames(common)
  if (length (a$gs_name[common$ensembl_gene_id[i] == a$ensembl_gene]) == 1){
    common$paths[i] <-a$gs_name[common$ensembl_gene_id[i] == a$ensembl_gene]
  } else {
   common$paths[i] <-a$gs_name[common$ensembl_gene_id[i] == a$ensembl_gene][1]
   z<- common[i, ]
   z$paths<- a$gs_name[common$ensembl_gene_id[i] == a$ensembl_gene][2]
   common <- rbind (common,z)
  }
}


## Wite the genes
write.xlsx(common, paste0("Cor_of_NOTCH_Pathway", Sys.Date(), ".xlsx"))

