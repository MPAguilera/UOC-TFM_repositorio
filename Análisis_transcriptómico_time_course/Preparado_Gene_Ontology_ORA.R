###############################################################################
############################# - ORA GEO -####################################
################################################################################

################# - Automatized - ################# 

mensaje_error <- c("\n \n ............................ERROR...................................\n")

name <-args[[1]]
nombre <- args[[2]]
common<-  args[[3]]

################# - Upload libraries - #########################################
## Here package to be used from CRAN
packages<- c("rstudioapi","clusterProfiler","msigdbr","ggplot2", "enrichplot", "gridExtra","openxlsx","readxl","BiocManager", "corrplot", "pheatmap", "factoextra","RColorBrewer")

## Here package to be used from BiocManager
biocM <- c("biomaRt","clusterProfiler","pathview","wordcloud"); combine <- c(packages, biocM)

## If they are not already installed install them
for (i in packages){ if (!requireNamespace(i, quietly = TRUE)){install.packages(i)}}
for (i in biocM){if (!requireNamespace(i, quietly = TRUE)){BiocManager::install(i)}}

##Upload libraries
for (i in combine){suppressMessages(suppressWarnings(library(i, character.only = TRUE)))}
rm (packages, biocM, combine, i)

################# - Functions - ################################# 
draw <- function (edo, n, lista_resultados){

  if (!is.null(edo@result)){
  if (nrow (edo@result[edo@result$p.adjust< 0.05,])>=1){
  pdf (paste0(n, nombre, Sys.Date(), ".pdf"), width = 15, height = 15)
  #print (goplot(edo, showCategory = 25 ) + ggtitle (paste0("goplot", nombre)))
  print (barplot(edo, showCategory=20) + ggtitle (paste0("barplot", nombre)))
  print (dotplot(edo, showCategory=20) + ggtitle (paste0("dotplot", nombre)))
  print (cnetplot(edo,   color.params = list(foldChange = edo@result$p.adjust), node_label = "category", showCategory = 20) + ggtitle (paste0("cnetplot", nombre)))
  
    
    ## If there is enough pathways clusterize pathways
    if (nrow(edo@result[edo@result$p.adjust< 0.05, ]) >= 8){
         tryCatch({
          edox2 <- pairwise_termsim(edo)
          
          p1 <- treeplot(edox2) + ggtitle(nombre)
          p2 <- treeplot(edox2, hclust_method = "average") + ggtitle(paste0("treeplot", nombre))
          
          print(aplot::plot_list(p1, p2, tag_levels = 'A'))
        }, 
        error = function(e) {
          message("Ocurrió un error: ", e$message)
        })
      ## Emmaplot para ver relaciones biológicas
      tryCatch({ 
        edo2 <- pairwise_termsim(edo) 
        print (emapplot(edo2, layout="kk"))+ ggtitle(paste ("emmaplot", nombre))
      
      if (!is.null(dev.list())){dev.off()}
      if (!is.null(dev.list())){dev.off()}
        }, 
        error = function(e) {
          message("Ocurrió un error: ", e$message)
        })
      }
    
    if (length(lista_resultados) == 0){
      lista_resultados[1:2] <- list (edo@result, edo@result[edo@result$p.adjust< 0.05,])
      names(lista_resultados)[1:2] <- c(paste0("ORA_",n),paste0("ORA", n, "padjust"))
  
    } else {
        lista_resultados[(length(lista_resultados) +1 ) : (length(lista_resultados) +2) ] <- list (edo@result, edo@result[edo@result$p.adjust< 0.05,])
        names(lista_resultados)[(length(lista_resultados)-1): (length(lista_resultados)) ] <- c(paste0("ORA_",n), paste0("ORA_padjust",n))
      }
  } else {
     lista_resultados[(length(lista_resultados) +1 )] <- list (edo@result)
     names(lista_resultados)[(length(lista_resultados))] <- paste0("ORA_",n)
  }
  return (lista_resultados)
  }
}

################# - Working directory- #################################


## Create directory for the daty
Directory_for_day <- dir.create(paste0(name, Sys.Date()))

################# - Files - #################################

lista_resultados <- list ()
################# - GO- #################################

## Set the working directory
setwd(paste0(name, Sys.Date()))
Directory<- getwd()
dirnombre<- dir.create(paste0(nombre, Sys.Date()))
setwd(paste0(nombre, Sys.Date()))

## GO results
ORA_GO_BP<- enrichGO(gene = common, 
                keyType = "ENSEMBL",
                OrgDb =  "org.Hs.eg.db",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05, 
                pvalueCutoff  = 0.05,
                readable      = TRUE)

ORA_GO_BP <-clusterProfiler::simplify(ORA_GO_BP, cutoff = 0.7,by = "p.adjust", select_fun = min)

lista_resultados <- list()
lista_resultados <- draw (ORA_GO_BP, "GO_BP", lista_resultados); if (!is.null(dev.list())){dev.off()}; rm (ORA_GO_BP)

################# - KEGG - #################################

attributes <- c("ensembl_gene_id", 'entrezgene_id')

data_annotation <- getBM(attributes = attributes, 
                           filters = "ensembl_gene_id", 
                           values = common, 
                           mart = ensembl)

## Genvert gene ENSEMBL into ENTREZ ID
search_kegg_organism('Homo sapiens',  by='scientific_name')
ORA_KEGG <- enrichKEGG(gene= data_annotation$entrezgene_id,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05)
if(is.null(ORA_KEGG) == F){lista_resultados <- draw (ORA_KEGG, "ORA_KEGG", lista_resultados); if (!is.null(dev.list())){dev.off()}; rm (ORA_KEGG)}

## Wiki pathways
ORA_WIKI <- enrichWP(data_annotation$entrezgene_id, organism = "Homo sapiens") 
if(is.null(ORA_WIKI) == F){lista_resultados <-draw (ORA_WIKI, "ORA_WIKI", lista_resultados); if (!is.null(dev.list())){dev.off()}; rm (ORA_WIKI)}

################# - DISEASE - #################################
## Now the diseases, because we are focusing on an already specific disease, 
## better not use
## Disease Overrepresentation
diseases = F
if (diseases ==T){
  library(DOSE)
ORA_DO <- enrichDO(gene = data$entrezgene_id, 
                   organism = "hsa",
                   ont           =  "HDO",
                   pvalueCutoff  = 0.05,
                   pAdjustMethod = "BH",
                   #universe      = names(geneList),
                   minGSSize     = 5,
                   maxGSSize     = 500,
                   qvalueCutoff  = 0.05,
                   readable      = FALSE)


## Network of cancer gene
ORA_ncg <- enrichNCG(dat$entrezgene_id) 

## Disease gene network
ORA_DGN <- enrichDGN(dat$entrezgene_id) 
}

################# - msigdbr - #################################
library(msigdbr)
possible_genesets<- msigdbr_collections()
m_df <- msigdbr(species = "Homo sapiens")
H<- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
H<- enricher(data_annotation$entrezgene_id, TERM2GENE=H) 
if (!is.null(H)) {
  lista_resultados <-draw (H, "ORA_Hallmarks", lista_resultados);if (!is.null(dev.list())){dev.off()}; rm (H)
}

##  Biocarta
BIOCARTA<- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  filter(gs_subcollection== "CP:BIOCARTA") %>%
  dplyr::select(gs_name, entrez_gene)
BIOCARTA<- enricher(data_annotation$entrezgene_id, TERM2GENE=BIOCARTA) 

if (!is.null(BIOCARTA)) {
  lista_resultados <-draw (BIOCARTA, "BIOCARTA", lista_resultados); if (!is.null(dev.list())){dev.off()}; 
  rm (BIOCARTA)
}


##  WikiPath
WIKI<- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  filter(gs_subcollection== "CP:WIKIPATHWAYS") %>%
  dplyr::select(gs_name, entrez_gene)
WIKI<- enricher(data_annotation$entrezgene_id, TERM2GENE = WIKI) 
if (!is.null(WIKI)) {
  lista_resultados<-draw (WIKI, "msigdbr_WIkipathways", lista_resultados); if (!is.null(dev.list())){dev.off()};
  rm (WIKI)
}
##  Oncogenic_Signatures
C6<- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, entrez_gene)
Oncogenic_Signature<- enricher(data_annotation$entrezgene_id, TERM2GENE = C6) 
if (!is.null(Oncogenic_Signature)) {
lista_resultados<-draw (Oncogenic_Signature, "msigdbr_C6_Oncogenic_Signature",lista_resultados); if (!is.null(dev.list())){dev.off()}
rm (Oncogenic_Signature)
}

##  Inmune
REACTOME<- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  filter(gs_subcollection== "CP:REACTOME") %>%
  dplyr::select(gs_name, entrez_gene)
REACTOME<- enricher(data_annotation$entrezgene_id, TERM2GENE = REACTOME)
if (!is.null(REACTOME)) {
  lista_resultados<-draw (REACTOME, "REACTOME", lista_resultados); if (!is.null(dev.list())){dev.off()}; rm (REACTOME)
}
##  REACTOME
PID<- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  filter(gs_subcollection== "CP:PID") %>%
  dplyr::select(gs_name, entrez_gene)
PID<- enricher(data_annotation$entrezgene_id, TERM2GENE = PID) 
if (!is.null(PID)) { lista_resultados<-draw (PID, "PID", lista_resultados); if (!is.null(dev.list())){dev.off()}; rm (PID)}

## Cancer Gene Neighborhoods
CGN<- msigdbr(species = "Homo sapiens", category = "C4") %>% 
  filter(gs_subcollection== "CGN") %>%
  dplyr::select(gs_name, entrez_gene)
CGN<- enricher(data_annotation$entrezgene_id, TERM2GENE =CGN) 
if (!is.null(CGN)) { lista_resultados <-draw (CGN, "CGN", lista_resultados); if (!is.null(dev.list())){dev.off()}; rm (CGN)
}
cat ("...................... GUARDANDO TABLAS ............")

# Crear un archivo Excel
wb <- createWorkbook()
names(lista_resultados) <- gsub ("ORA_","",names(lista_resultados)) %>%
                              gsub ("_C6_", "", .) %>%
                                 gsub ("msigdbrOncogenic", "msigdbrONCO", .)
  
  
for (nom in names(lista_resultados)) {
  addWorksheet(wb, nom) 
  writeData(wb, sheet = nom, lista_resultados[[nom]])
}

saveWorkbook(wb, file = paste0(nombre, ".xlsx"), overwrite = TRUE)
capture.output(Sys.info(), file = paste0("R-SystemInfo_", Sys.Date(), ".txt"))
capture.output(sessionInfo(), file = paste0("R-SessionInfo_", Sys.Date(), ".txt"))


cat ("...................... TRABAJO TERMINADO ............")


## For plotting  the pathway
# hsa04110 <- pathview(gene.data  = data_annotation$entrezgene_id, gene.idtype = "entrez", pathway.id = "hsa04141", species = "hsa")