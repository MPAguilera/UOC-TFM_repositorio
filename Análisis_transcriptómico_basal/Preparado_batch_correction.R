
##########################- BATCH CORRECTION ##################################

################# - Upload libraries - #########################################

## Here package to be used from CRAN
packages<- c("RColorBrewer", "rstudioapi","ggplot2","gridExtra","openxlsx","readxl","BiocManager")

## Here package to be used from BiocManager
biocM <- c("sva", "DESeq2", "biomaRt"); combine <- c(packages, biocM)

## If they are not already installed install them
for (i in packages){ if (!requireNamespace(i, quietly = TRUE)){install.packages(i)}}
for (i in biocM){if (!requireNamespace(i, quietly = TRUE)){BiocManager::install(i)}}

##Upload libraries
for (i in combine){suppressMessages(suppressWarnings(library(i, character.only = TRUE)))}
rm (packages, biocM, combine, i)

###############################################################################
################ - WORKING DIRECTORY and checking of the Input files- #########
###############################################################################

## Instrucciones para correr el script
cat(" \n....................INSTRUCCIONES GENERALES DE USO .........................
\n Confirme que se ha introducido el directorio de trabajo en el que estará la carpeta Input de la siguiente forma Rscript analisis.R nombre_directorio_de_trabajo nombre_archivo_salida\n
    Esta carpeta Input debe contener: \n
    - Archivo batch que contenga en el nombre contenga batch y únicamente tenga dos columnas, una con SampleID y otra con batch en xlsx 
    - Archivo con los counts que en el nombre contenga counts y que las columnas sean los samples y las filas los genes en xlsx o tsv
    - Archivo Resumen con los subtipos genéticos y demás información (debe tener dos columnas Grupo y SampleID con los nombres de los samples) en xlsx o tsv ")

#### Before running the script for batch correction check that the arguments require have been provided
##Only run this code when you are making trials in r environment
mensaje_error <- c("\n \n ............................ERROR...................................\n")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {paste("\n Se ha introducido el directorio a emplear en el análisis", args[1])} else {
  stop("\n \n ............................ERROR...................................
       \n No se han proporcionado argumentos introduzca: Rscript analisis.R nombre_directorio_de_trabajo nombre_output", call. = FALSE)
}
setwd(args[1])


## Check that the files to be used and the Input directory exists and read the files
if (any(grepl("/Input", list.dirs()))){
  setwd(list.dirs()[grep("/Input", list.dirs())])
} else {
  stop(paste0(mensaje_error, "El directorio Input no se encuentra en el directorio proporcionado"))
}

if (any(grepl("batch", list.files()))){
  file = list.files()[grep("batch", list.files())]
  if (any(grepl("xlsx", file))){
    batch = read.xlsx(list.files()[grep("batch", list.files())])
  } else { batch = read.delim(list.files()[grep("batch", list.files())]) }
} else{
  stop(paste0(mensaje_error,"El archivo batch no se encuentra en el directorio Input"))}

if (any(grepl("counts", list.files()))){
  file = list.files()[grep("counts", list.files())]
  if (grepl("xlsx", file)){
    counts = read.xlsx(list.files()[grep("counts", list.files())])
  } else { counts = read.delim(list.files()[grep("counts", list.files())]) }
  row.names(counts) <- counts[,"ensembl_gene_id"]
  counts <- counts [,-1]
  if (any(grepl("biotype|symbol|description|external_gene_name|ensembl_gene_id", colnames(counts)))){
    counts <- counts[, -(grep("biotype|symbol|description|external_gene_name|ensembl_gene_id", colnames(counts)))]
  }
} else{
  stop(paste0(mensaje_error, "El archivo counts no se encuentra en el directorio Input"))}

if (any(grepl("resumen", list.files()))){
  file = list.files()[grep("resumen", list.files())]
  if (grepl("xlsx", file)){
    resumen = read.xlsx(list.files()[grep("resumen", list.files())])
  } else { resumen = read.delim(list.files()[grep("resumen", list.files())]) }
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
Directory_for_the_day <- paste0("Batchcorrection_Normalization_clustering", args[2],Sys.Date())
if(any(grepl(Directory_for_the_day, list.dirs()))){paste(Directory_for_the_day, "ya creado")} else {dir.create(Directory_for_the_day); paste ("Generado el directerio del día:", Directory_for_the_day)}

cat("\n \n ............................ADQUISICION DE DATOS FINALIZADA...................................\n")

########################## - CHECK THE CORRECT NAMES  -#########################

## Substitute the names so the . are substitute by - and all the letters are set to upper

## Check that the names are the same
if ((ncol(counts)== nrow (batch)) & (nrow (resumen) == ncol (counts))){print ("EL número de samples en cada tabla es correcto")} else {stop ("El número de samples no coincide en una de las tablas")}

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
coldata <- data.frame(SampleID = colnames(counts),subtype = resumen$Group, batch = batch$batch)
rownames(coldata) <- coldata$SampleID

cat("\n \n ............................CHEQUEO DE NOMBRES Y ORDEN FINALIZADA...................................\n")
########################## - COMBAT-SEQ -#####################################

cat("\n \n ............................INICIACION DE LA CORRECION DEL BATCH...................................\n")
COMBAT_NO_NORMALIZED = data.frame(ComBat_seq(as.matrix(counts), 
                                                     batch = as.factor(as.character(batch$batch)), 
                                                     group = as.factor(as.character(resumen$Group))))
colnames(COMBAT_NO_NORMALIZED) <- gsub ("\\.", "-", colnames(COMBAT_NO_NORMALIZED))

setwd(Directory_for_the_day)
write.xlsx(COMBAT_NO_NORMALIZED, paste0("Data_10_cell_lines_raw_counts_batch_corrected",Sys.Date(),".xlsx"))

cat("\n \n ............................FINALIZACION DE LA CORRECION DEL BATCH...................................\n")
