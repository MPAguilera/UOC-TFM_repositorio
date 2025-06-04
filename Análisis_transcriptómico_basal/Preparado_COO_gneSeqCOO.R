################################################################################
##########################- COO classification gnsSeqCOO- ##########################
################################################################################


################# - Upload libraries - #########################################

## Here package to be used from CRAN
packages<- c("rstudioapi","ggplot2","openxlsx","readxl","BiocManager", "ComplexHeatmap", "pheatmap", "dbplyr", "circlize","RColorBrewer", "remotes", "R.utils")

## Here package to be used from BiocManager
biocM <- c("sva", "DESeq2", "edgeR", "biomaRt"); combine <- c(packages, biocM)
library(circlize)
## Here the package that need to be downloaded from github
extra <- c("gneSeqCOO")

## If they are not already installed install them
for (i in packages){ if (!requireNamespace(i, quietly = TRUE)){install.packages(i)}}
for (i in biocM){if (!requireNamespace(i, quietly = TRUE)){BiocManager::install(i)}}
for (i in extra) {if (!requireNamespace(i, quietly = TRUE)){remotes::install_github("Genentech/gneSeqCOO")}}

##Upload libraries
for (i in combine){library(i, character.only = TRUE)}
rm (packages, biocM, combine, extra, i)


################# - Scripts Instructions - #################################

cat(" \n....................INSTRUCCIONES GENERALES DE USO .........................
\n Confirme que se ha introducido el directorio de trabajo en el que estará la carpeta Input de la siguiente forma Rscript analisis.R nombre_directorio_de_trabajo nombre_archivo_salida\n
    Esta carpeta Input debe contener: \n
    - Archivo batch que contenga el nombre contenga batch y únicamente tenga dos columnas, una con SampleID y otra con batch en xlsx 
    - Archivo con los counts que en el nombre contenga counts y que las columnas sean los samples y las filas los genes en xlsx o tsv
    - Archivo Resumen con los subtipos genéticos y demás información (debe tener dos columnas Grupo y SampleID con los nombres de los samples) en xlsx o tsv
\n") 

################# - Arguments to the automatization - #################################

# Before running the script for batch correction check that the arguments require have been provided
# Only run this code when you are making trials in r environment  args <- c("G:/Mi unidad/Linfomas/COO_gneSeqCOO","DLBCL")

mensaje_error <- c("\n \n ............................ERROR...................................\n")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {paste("Se ha introducido el directorio a emplear en el análisis", args[1])} else {
  stop("\n \n ............................ERROR...................................
       \n No se han proporcionado argumentos introduzca: Rscript analisis.R nombre_directorio_de_trabajo nombre_output", call. = FALSE)
}

## Set the working directory to the file indicated
setwd(args[1])

################# - Working directory and data upload- ###############################

## Check that the files to be used and the Input directory exists and read the files
if (any(grepl("/Input", list.dirs()))){
  setwd(list.dirs()[grep("/Input", list.dirs())])
} else {
  stop(paste0(mensaje_error, "El directorio Input no se encuentra en el directorio proporcionado"))
}

if (any(grepl("batch", list.files()))){
  if (grep ("batch*.*lnk", list.files ())){
    file <- R.utils::filePath(list.files ()[grep ("batch*.*lnk", list.files ())], expandLinks = "any")
  } 
  else{file = list.files()[grep("batch", list.files())]}
  if (grepl("xlsx", file)){
    batch = read.xlsx(file)
  } else { batch = read.delim(file) }
} else{
  stop(paste0(mensaje_error,"El archivo batch no se encuentra en el directorio Input"))}

if (any(grepl("counts", list.files()))){
  if (grep ("counts*.*lnk", list.files ())){
    file <- R.utils::filePath(list.files ()[grep ("counts*.*lnk", list.files ())], expandLinks = "any")
  }
  else {file = list.files()[grep("counts", list.files())] }
  if (grepl("xlsx", file)){
    counts = read.xlsx(file)
  } else { counts = read.delim(file) }
    
  rownames(counts) <- counts[,"ensembl_gene_id"]
  counts <- counts [,-1]
  if (any(grepl("biotype|symbol|description|external_gene_name|ensembl_gene_id", colnames(counts)))){
    counts <- counts[, -(grep("biotype|symbol|description|external_gene_name|ensembl_gene_id", colnames(counts)))]
  }
} else{
  stop(paste0(mensaje_error, "El archivo counts no se encuentra en el directorio Input"))}

if (any(grepl("resumen", list.files()))){
  if (grep ("resumen*.*lnk", list.files ())){
    file <- R.utils::filePath(list.files ()[grep ("resumen*.*lnk", list.files ())], expandLinks = "any")
  }
  else{ file = list.files()[grep("resumen", list.files())]}
  if (grepl("xlsx", file)){
    resumen = read.xlsx(file)
  } else { resumen = read.delim(file) }
} else{
  stop(paste0(mensaje_error, "El archivo resumen no se encuentra en el directorio Input"))}


## Print message that the files have been susccessful upload into the environment
dataframes <- sapply(ls(), function(x) is.data.frame(get(x)))
if (sum(dataframes) ==3){print ("Los tres dataframes se han cargado correctamente en el environment de R")} else {stop(paste0(mensaje_error,"Los dataframes no se han cargado correctamente"))}

## Check that the files have the correct columns
if (!all(c("SampleID", "batch") %in% colnames(batch))){stop("Las columnas de la tabla batch no son correctas")}
if (!all(c("SampleID", "Group") %in% colnames(resumen))){stop("Las columnas de la tabla resumen no son correctas")}

## Create a directory for the results with the current Data
setwd(args[1])
Directory_for_the_day <- paste0("gneSeqCOO", "_",args[2],"_",Sys.Date())
if(any(grepl(Directory_for_the_day, list.dirs()))){paste(Directory_for_the_day, "ya creado")} else {dir.create(Directory_for_the_day); paste ("Generado el directerio del día:", Directory_for_the_day)}

cat("\n \n ............................ADQUISICION DE DATOS FINALIZADA...................................\n")



################# - Check that the names are correct  -#########################

# Check that the number of samples are the same in the three files
if (((ncol(counts))== nrow (batch)) & (nrow (resumen) == ncol (counts))){print ("El número de samples en cada tabla es correcto")} else {stop ("El número de samples no coincide en una de las tablas")}

##Check the names from files are the same. First remove the mst common mistakes
colnames(counts) <- toupper(gsub ("\\.", "-", colnames(counts)))
batch$SampleID<-toupper(gsub ("\\.", "-", batch$SampleID))
resumen$SampleID<-toupper(gsub ("\\.", "-", resumen$SampleID))

if (all(sort(colnames(counts)) == sort(batch$SampleID))){print ("Los nombres de datos y batch son correctos y coinciden")} else {
  cat ("Los nombres que difieren entre el batch y los counts son:", setdiff(colnames(counts), batch$SampleID), "\n")
  stop("Los nombres de batch no cohinciden")}
if (all(sort(colnames(counts)) == sort(resumen$SampleID))){print ("Los nombres de datos y anotaciones son correctos y coinciden")} else {
  cat ("Los nombres que difieren entre el resumen y los counts son:", setdiff(colnames(counts), batch$SampleID), "\n")
  stop("Los nombres de batch no cohinciden")}

##Reorder the Annot and data to the batch
batch <- batch[match (colnames(counts), batch$SampleID),]
resumen <- resumen[match (colnames(counts), resumen$SampleID),]
if (all(colnames(counts) == batch$SampleID)){print ("Los nombres del batch están en el orden correcto")} else {
  cat ("Nombres no ordenados correctamente:", batch$SampleID[!colnames(counts) == batch$SampleID], "\n")
  stop ("Los nombres del batch no están en el orden correcto")}
if (all(colnames(counts) == resumen$SampleID)){print ("Los nombres del resumen están en el orden correcto")} else {
  cat ("Nombres no ordenados correctamente:",resumen$SampleID[!colnames(counts) == resumen$SampleID], "\n")
  stop ("Los nombres del resumen no están en el orden correcto")}
coldata <- data.frame(SampleID = colnames(counts),subtype = as.factor(resumen$Group), batch = as.factor(batch$batch))
rownames(coldata) <- coldata$SampleID

cat("\n \n ............................CHEQUEO DE NOMBRES Y ORDEN FINALIZADA...................................\n")


################# - gneSeqCOO  - ################################

cat("\n \n ............................Generando el gneSeqCOO...................................\n")
COMBAT_NO_NORMALIZED <-counts
if(!all(sort(colnames(COMBAT_NO_NORMALIZED)) == sort(rownames(coldata)))) {stop ("Los nombres de los samples del coldata y datos de COMBAT no coinciden")}
if(!all(colnames(COMBAT_NO_NORMALIZED) == rownames(coldata))) {stop ("Los nombres de los samples del coldata y datos de COMBAT no están en el orden correcto")}
dds <- DESeqDataSetFromMatrix(countData = round(COMBAT_NO_NORMALIZED),
                              colData = coldata,
                              design = ~ as.factor(subtype))

pred = gneSeqCOO::coo_rnaseq(dds)

