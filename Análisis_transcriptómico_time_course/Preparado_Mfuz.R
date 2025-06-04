##########################- Mfuzz- #############################################


################# - Upload libraries - #########################################

## Here package to be used from CRAN
packages<- c("rstudioapi","ggplot2","gridExtra","openxlsx","dplyr","readxl","BiocManager", "corrplot", "pheatmap", "factoextra","RColorBrewer", "ashr")

## Here package to be used from BiocManager
biocM <- c("sva", "DESeq2", "edgeR", "biomaRt"); combine <- c(packages, biocM)

## If they are not already installed install them
for (i in packages){ if (!requireNamespace(i, quietly = TRUE)){install.packages(i)}}
for (i in biocM){if (!requireNamespace(i, quietly = TRUE)){BiocManager::install(i)}}

##Upload libraries
for (i in combine){suppressMessages(suppressWarnings(library(i, character.only = TRUE)))}
rm (packages, biocM, combine, i)

################# - Scripts Instructions - #################################

## Mensaje error
mensaje_error <- c("...................... ERROR ............")

## Set WD to file location
workingD <- rstudioapi::getActiveDocumentContext()$pat
setwd(dirname(workingD))

################# - Working directory and data upload- ###############################
files <- list.files (recursive =T )
set.seed(11235)
resumen = read.xlsx(files[grep("Input*.*resumen", files)][1])

## Read significant genes for each cell comparisson
celllines <- c("WSU", "U2932","SUDHL")
lista_data <- list ()
for (j in 1:length(celllines)){
  for (i in 1:length(grep (paste0("withreplicates*.*",celllines[j],"*.*Sig*.*xlsx"), files))){
    a<-read.xlsx(files[grep (paste0("withreplicates*.*", celllines[j],"*.*Sig*.*xlsx"), files)][i])
    a <- a[, grep (paste0("ensembl_gene_id|external_gene_name|", celllines[j]), colnames(a))]
    if (i == 1){Sig <- a} else{ Sig <- rbind (Sig,a)}
    rm (a)
    Sig <- Sig[!duplicated(Sig$ensembl_gene_id),]
    lista_data[[j]] <- Sig
    names(lista_data)[j] <- celllines[j]
  }
}

## Create a directory for the results with the current Data
setwd(dirname(workingD))
Directory_for_the_day <- paste0("Mfuzz_VST_",Sys.Date())
dir.create(Directory_for_the_day)

################# - Perform Mfuzz  -#########################
VST <-read_xlsx(files[grep ("2ndround*.*VST", files)[2]])

for (i in 1:length(lista_data)){
  a <-data.frame(VST[match(lista_data[[i]]$ensembl_gene_id, VST$ensembl_gene_id),-c(1,3,4,5)])
  rownames(a) <- a$ensembl_gene_id; 
  a <- a[, -match ("ensembl_gene_id", colnames(a))]
  lista_data[[i]] <- a
  rm (a)
}

## Create data wqith only the expr of interest and no chr columns
setwd(Directory_for_the_day)
clusters <- c(3,5,3)
library(Mfuzz)

for (i in 1:length(lista_data)){
  Sig <- lista_data[[i]]
  Sig <- Sig[,grep (names(lista_data)[i], colnames(Sig))]
  Sig <- Sig [,c(grep ("DMSO", colnames(Sig)), grep ("_6H_", colnames(Sig)), grep ("_24H_", colnames(Sig)))]

  ## Generate the expressionset
  es<-ExpressionSet(assayData = as.matrix(Sig))

  ## Missing values
  m<- filter.NA(es, thres=0.25)
  m <- fill.NA(m,mode='knn')

  ## Filtering
  f<- filter.std(m,min.std=0)

  ## Standarization
  s<- standardise(f)
  
  ## Clsutering
  m1 <- mestimate(s)
  Dmin(es, m = m1, crange = (2:10), repeats = 3)

  cln <- clusters[i]
  cl <- mfuzz(s,c=cln,m=m1)
  cl$membership
  mfuzz.plot(s,cl=cl,mfrow=c(2,2), new.window = F, time.labels = c(0,0,0,6,6,6,24,24,24))
  O <- overlap(cl)
  write.xlsx(O, paste0("Overlap_clsuters_", celllines[j], "_",Sys.Date(),"_.xlsx"))
  print (O)
  Ptmp <- overlap.plot(cl,over=O,thres=0.05)


## Save results
Sig$cluster_mfuzz<- cl$cluster[match(rownames(Sig),names(cl$cluster))]
lista_data[[i]] <- Sig
write.xlsx(data.frame(ensembl_gene_id = rownames(Sig), Sig), paste0(names(lista_data)[i],"_all_dif_genes_", Sys.Date(), "_.xlsx"))
pdf(paste0(names(lista_data)[i],"_Mfuzz_Cluster", Sys.Date(), ".pdf"))
Dmin(es, m = m1, crange = (2:10), repeats = 3)
mfuzz.plot(s,cl=cl,mfrow=c(2,2), new.window = F, time.labels = c(0,0,0,6,6,6,24,24,24))
par (mfrow = c(1,1))
overlap.plot(cl,over=O,thres=0.05)
dev.off()
}

rm (m, f,es, cl, O, Ptmp, s, Sig)
############################# - ORA -####################################

## PerformORAin each cluster

## Set WD to file location
workingD <- rstudioapi::getActiveDocumentContext()$pat
setwd(Directory_for_the_day)

################# - List files - ################# 

cellines <- c("WSU", "U2932", "SUD")
lista <- list()
num <- 1
for (j in 1:length(lista_data)){
  numclusters <- max(unique(lista_data[[j]]$cluster_mfuzz))
  Sig <- lista_data[[j]]
  
  ## In each sublist fill with the genes for each cluster
  for (i in 1:numclusters){
    lista[[num]]<- rownames(Sig)[Sig$cluster_mfuzz == i]
    names(lista)[num] <- paste0("Cluster_", i, "_", names(lista_data)[j])
    num <- num +1
  }
}
rm (Sig)

##Aply the GEO 
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

for (is in 1:length(lista)) {
  setwd(dirname(workingD))
  print(names(lista)[is])
  args <- list()  # Resetear args
  args[[1]] <- paste0("ORA_", Sys.Date())
  args[[2]] <- paste0("ORA_", names(lista[is]))
  args[[3]] <- lista[[is]]
  print(args[[2]])
      
      source("Gene_Ontology_ORA.R")
}

