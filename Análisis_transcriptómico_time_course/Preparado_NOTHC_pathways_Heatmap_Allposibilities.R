################################################################################
###############################- NOTCH pathways of interest- ###################


################# - Upload libraries - #########################################

if (require(rstudioapi)==F){install.packages("rstudioapi")}
if (require(readxl)==F){install.packages("readx")}
if (require(openxlsx)==F){install.packages("openxlsx")}
if (require(gridExtra)==F){install.packages("gridExtra")}
if (require(ggplot2)==F){install.packages("ggplot2")}
library(paletteer)
library(dbplyr)
library(pheatmap)
library(msigdbr)
library(dplyr)

########################## - WORKING DIRECTORY - ###############################

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


########################## - Genes of NOTCH in DEGS###############################
files = list.files(recursive = T, full.names = T)
files <-files[grep ("withreplicated*.*Significant*.*xlsx", files)]
lista <- list()
cell_lines <- c("WSU", "SUDHL", "U2932")
for (j in 1:length(cell_lines)){
  look <- grep (cell_lines[j], files)
  for (z in 1:length(look)){
    if ( nrow(read.xlsx(files[look[z]])) !=0){
    if (z ==1 ){b <- data.frame(Condition = gsub (paste0(".*",cell_lines[j]),"", files[look[z]]) %>% gsub ("Sig.*", "",.) %>% gsub (".*DESEQ2_","",.),
                                read.xlsx(files[look[z]]))} else { 
                                  b<-rbind (b, data.frame(Condition = gsub (paste0(".*",cell_lines[j]),"", files[look[z]]) %>% gsub ("Sig.*", "",.) %>% gsub (".*DESEQ2_","",.),
                                                          read.xlsx(files[look[z]])))}
  }
}
  
  ## Look for the info with cluster 
  files2 = list.files(recursive = T, full.names = T)
  files2 <-files2[grep (paste0("VST*.*",cell_lines[j] ,"*.*all_dif_genes*.*xlsx"), files2)]
  files2 <-files2[-grep("Sin replicas", files2)]
  genes_dif <-read.xlsx(files2)
  
  
  ## Get the pathways
  t <-notch_pathways[ notch_pathways$ensembl_gene %in%b$ensembl_gene_id,]
  
  genes <- unique(t$ensembl_gene)
  guardar <- data.frame (matrix(ncol = 6, nrow = 0))
  colnames(guardar) <- c("Ensembl_gene_id", "gene_name", "Condition","Significance", "Cluster","Pathways")
  for (i in 1:length(genes)){
   guardar <- rbind (guardar,  
    data.frame(Ensembl_gene_id = genes[i], 
               gene_name = b$external_gene_name[b$ensembl_gene_id == genes[i]][1],
               Condition = paste0(b$Condition[b$ensembl_gene_id == genes[i]], collapse  = "|"), 
                Significance = paste0(b$Significance[b$ensembl_gene_id == genes[i]], collapse = "|") ,
               Cluster = genes_dif$cluster_mfuzz[genes_dif$ensembl_gene_id == genes[i]],
               Pathways  = paste0(t$gs_name[t$ensembl_gene == genes[i]], collapse  = "|")))
    
    
  }

  guardar$Pathways <- 
gsub ("WP_APOPTOSISRELATED_NETWORK_DUE_TO_ALTERED_NOTCH3_IN_OVARIAN_CANCER", "", 	guardar$Pathways) %>%
  gsub ("VILIMAS_NOTCH1_TARGETS_UP", "",.) %>%
  gsub ("WP_NOTCH1_REGULATION_OF_ENDOTHELIAL_CELL_CALCIFICATION", "",.) %>%
  gsub ("BRAUNE_GEIST_20_GENE_NOTCH_SIG_BREAST_CANCER", "", .) %>%
  gsub ("GAVISH_3CA_METAPROGRAM_ENDOTHELIAL_NOTCH_SIGNALING", "",.) %>%
  gsub ("NGUYEN_NOTCH1_TARGETS_DN", "",.)

 guardar <- guardar[guardar$Pathways != "|" & guardar$Pathways != "",]
lista[[j]] <- guardar
}


## Wite the genes
write.xlsx(lista[[1]], paste0("Notch_DESEQ2_WSU", Sys.Date(), ".xlsx"))
## Wite the genes
write.xlsx(lista[[2]], paste0("Notch_DESEQ2_SUDHL", Sys.Date(), ".xlsx"))
## Wite the genes
write.xlsx(lista[[3]], paste0("Notch_DESEQ2_U2932", Sys.Date(), ".xlsx"))



