################################################################################
##########################- DEWSEQ and PCA ##########################
################################################################################


################# - Upload libraries - #########################################

## Here package to be used from CRAN
packages<- c("rstudioapi","ggplot2","gridExtra","openxlsx","readxl","BiocManager", "corrplot", "pheatmap", "factoextra","RColorBrewer", "ashr")

## Here package to be used from BiocManager
biocM <- c("sva", "DESeq2", "edgeR", "biomaRt"); combine <- c(packages, biocM)

## If they are not already installed install them
for (i in packages){ if (!requireNamespace(i, quietly = TRUE)){install.packages(i)}}
for (i in biocM){if (!requireNamespace(i, quietly = TRUE)){BiocManager::install(i)}}

##Upload libraries
for (i in combine){suppressMessages(suppressWarnings(library(i, character.only = TRUE)))}
rm (packages, biocM, combine, i)

################# - Scripts Instructions - #################################
################# - Arguments to the automatization - #################################

args <- c("Directorio","Nombre")
mensaje_error <- c("\n \n ............................ERROR...................................\n")
if (length(args) == 2) {paste("Se ha introducido el directorio a emplear en el análisis", args[1])} else {
  stop("\n \n ............................ERROR...................................
       \n No se han proporcionado argumentos introduzca: Rscript analisis.R nombre_directorio_de_trabajo nombre_output BOTH/NORM/COMP EXP CONTROL", call. = FALSE)
}
setwd(args[1])

################# - Working directory and data upload- ###############################

## Check that the files to be used and the Input directory exists and read the files
if (any(grepl("/Input", list.dirs()))){
  setwd(list.dirs()[grep("/Input", list.dirs())])
} else {
  stop(paste0(mensaje_error, "El directorio Input no se encuentra en el directorio proporcionado"))
}

if (any(grepl("batch", list.files()))){
  file = list.files()[grep("batch", list.files())]
  if (grepl("xlsx", file)){
    batch = read.xlsx(list.files()[grep("batch", list.files())])
  } else { batch = read.delim(list.files()[grep("batch", list.files())]) }
} else{
  stop(paste0(mensaje_error,"El archivo batch no se encuentra en el directorio Input"))}

if (any(grepl("counts", list.files()))){
  file = list.files()[grep("counts", list.files())]
  if (grepl("xlsx", file)){
    counts = read.xlsx(list.files()[grep("counts", list.files())])
  } else { counts = read.delim(list.files()[grep("counts", list.files())]) }
  rownames(counts) <- counts[,"ensembl_gene_id"]
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
coldata <- data.frame(SampleID = colnames(counts),subtype = resumen$Group, batch = batch$batch)
rownames(coldata) <- coldata$SampleID

cat("\n \n ............................CHEQUEO DE NOMBRES Y ORDEN FINALIZADA...................................\n")


################# - DESEQ2 NORMALIZATION  - ################################

cat("\n \n ............................INICIACION DE DESEQ2...................................\n")
COMBAT_NO_NORMALIZED <-counts
if(!all(sort(colnames(COMBAT_NO_NORMALIZED)) == sort(rownames(coldata)))) {stop ("Los nombres de los samples del coldata y datos de COMBAT no coinciden")}
if(!all(colnames(COMBAT_NO_NORMALIZED) == rownames(coldata))) {stop ("Los nombres de los samples del coldata y datos de COMBAT no están en el orden correcto")}
dds <- DESeqDataSetFromMatrix(countData = round(COMBAT_NO_NORMALIZED),
                              colData = coldata,
                              design = ~ subtype )#+ batch)

## Here the code to extract 100%% of counts
dds1 = DESeq(dds)
normalized_100 = counts (dds1, normalized = T)
VST_100 = vst(dds1)

## Filter data to select only the ones that ahve expression in 90% of data
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds = dds[keep,]

## Perform the DESEQ2
dds2 = DESeq(dds)
normalized_90 = counts (dds2,normalized = T)
VST_90 = vst (dds2)
cat("\n \n ............................FINALIZACION DE DESEQ2...................................\n")

################# - Anotar las tablas  - ################################

## Ask the query
cat("\n \n ............................INICIACION DE ANOTCACION...................................\n")

## Get the parameter and BioMart
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl = useEnsembl(biomart = "genes") ## Define the dataet to use (genes)
searchDatasets(mart =  ensembl, pattern = "hsapiens") ## Look fotr the sensembl organism anotation
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") ## Get the anotation
attributes <- c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype") ## Get the atributes to include in the data
relevant_info <- c(searchDatasets(mart =  ensembl, pattern = "hsapiens"))

## Execute the anotatio 
Annotate = function (data){
  original_len = ncol (data)
  data_annotation <- getBM(attributes = attributes, 
                           filters = "ensembl_gene_id", 
                           values = rownames(data), 
                           mart = ensembl)
  data$ensembl_gene_id = rownames(data)
  data=merge (data_annotation, data,  by = "ensembl_gene_id", all.x = TRUE) 
  return (data)
}
data_Anot_DESEQ_FILTER_COUNTS_100 = Annotate(data.frame(normalized_100))
data_Anot_DESEQ_FILTER_VST_100= Annotate(data.frame(assay(VST_100)))
data_Anot_DESEQ_FILTER_COUNTS_90= Annotate(data.frame(normalized_90))
data_Anot_DESEQ_FILTER_VST_90= Annotate(data.frame(assay(VST_90)))

cat("\n \n ............................FINALIZACION DE ANOTACION...................................\n")

## Look only for the protein coding
only_prot <-function (x){
  c <- x[x$gene_biotype =="protein_coding",]
  rownames(c) <- c[,1]; c<- c[,-c(1:4)]
  return (c)
}
prot_100 <- only_prot(data_Anot_DESEQ_FILTER_COUNTS_100)
prot_90 <- only_prot(data_Anot_DESEQ_FILTER_COUNTS_90)
prot_VST_100 <-  only_prot(data_Anot_DESEQ_FILTER_VST_100)
prot_VST_90 <-  only_prot(data_Anot_DESEQ_FILTER_VST_90)

## Now get the VST select the protein coding
VST100prot <-VST_100[data_Anot_DESEQ_FILTER_VST_100$ensembl_gene_id[data_Anot_DESEQ_FILTER_VST_100$gene_biotype == "protein_coding"]]
VST90prot <-VST_90[data_Anot_DESEQ_FILTER_VST_90$ensembl_gene_id[data_Anot_DESEQ_FILTER_VST_90$gene_biotype == "protein_coding"]]


## Set thw woorking directory for the day to botch types of settings
setwd(Directory_for_the_day)

######################### - Save the graphs #######################################

## PCA graphs
plots_PCA_DESEQ = function(dat, col, titulo){
  da <- dat
  vector1 <- c(nrow (assay(da)), 1000,500,100)
  i <- 1
  lista <- list()
    da <-  plotPCA(dat, intgroup = "subtype", ntop = vector1[i], returnData = T)
    da$titulo <- paste0("PCA VST 80% filtrado \n",  vector1[i], " genes codificadores de proteínas")
    lista[[i]] <- ggplot(da, aes(x = PC1, y = PC2, color = as.factor(subtype), label = col$SampleID)) +
      geom_point(size = 5) +  # Tamaño de los puntos
      labs(color = "Subtipo") + 
      xlab(paste0("PC1: ", round(100 * attr(da, "percentVar")[1], 2), "% varianza")) +
      ylab(paste0("PC2: ", round(100 * attr(da, "percentVar")[2], 2), "% varianza")) +
      geom_text(vjust = -1, size = 4) + 
      facet_grid(. ~ titulo) + 
      theme_bw()+
      xlim (-80,55)+
      ylim (-45,85)+
      scale_color_manual(values = c("BN2" = "purple", 
                                    "Other" = "#999999",
                                    "N1" = "#5fe75f", 
                                    "ST2"= "#D02F4B",
                                    "EZB" = "#d17944",  "EZB/MYC+" = "brown",
                                    "MCD" = "#93ccde")) +
      theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
            legend.text = element_text(size= 16), 
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.y = element_text(size = 16, face = "bold"),
            axis.title.x = element_text(size = 16, face = "bold"),
            legend.title = element_text(face = "bold", size= 18, ),
            strip.text = element_text(color = "black", face = "bold", size =18),
            strip.background = element_rect(color = "black", fill = "grey"))
  
  print (lista[[1]])
}
plots_PCA_DESEQ_IC50 = function(dat, col, titulo, IC, subtitulo){

  da <- dat
  vector1 <- c(nrow (assay(da)), 1000,500,100)
  i <- 1
  lista <- list()
  da <-  plotPCA(dat, intgroup = "subtype", ntop = vector1[i], returnData = T)
  da$titulo <- paste0("PCA ", subtitulo)
  lista[[i]] <- ggplot(da, aes(x = PC1, y = PC2, color = IC, label = col$SampleID)) +
    geom_point(size = 5) +  # Tamaño de los puntos
    labs(color = paste0("IC50 ", subtitulo)) + 
    xlab(paste0("PC1: ", round(100 * attr(da, "percentVar")[1], 2), "% varianza")) +
    ylab(paste0("PC2: ", round(100 * attr(da, "percentVar")[2], 2), "% varianza")) +
    geom_text(vjust = -1, size = 4) + 
    facet_grid(. ~ titulo) + 
    scale_color_gradient(low = "grey", high = "purple") +  
    theme_bw()+
    xlim (-80,55)+
    ylim (-45,85)+
    theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
          legend.text = element_text(size= 16), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.title.x = element_text(size = 16, face = "bold"),
          legend.title = element_text(face = "bold", size= 18, ),
          strip.text = element_text(color = "black", face = "bold", size =18),
          strip.background = element_rect(color = "black", fill = "grey"))
  
  print (lista[[1]])
}


setwd(Directory_for_the_day)
##Print graphs
tiff(paste0("Plots_DESEQ2_PCA_",args[2],Sys.Date(), ".tiff"), width = 20, height = 15, res = 150, unit = "cm")
plots_PCA_DESEQ(VST90prot, coldata, "")
dev.off()
resumen$Droga1

##Print graphs
if (!all (resumen$SampleID == colnames(VST90prot))){stop ("Los nombres no están ordenados")}
tiff(paste0("Plots_DESEQ2_PCA",args[2],Sys.Date(), ".tiff"), width = 20, height = 15, res = 150, unit = "cm")
plots_PCA_DESEQ_IC50(VST90prot, coldata, "", resumen$Droga1, "Droga1")
dev.off()

tiff(paste0("Plots_DESEQ2_PCA",args[2],Sys.Date(), ".tiff"), width = 20, height = 15, res = 150, unit = "cm")
plots_PCA_DESEQ_IC50(VST90prot, coldata, "", resumen$Droga2,"Droga2")
dev.off()


Write.xlsx(counts(dds, normalize = T), paste0("Nom_library_size_counts.xlsx"))
Write.xlsx(assay(vst(dds)), paste0("Nom_VST_size_counts.xlsx"))

cat("\n \n ............................GRAFICOS PCA GUARDADOS...................................\n")
