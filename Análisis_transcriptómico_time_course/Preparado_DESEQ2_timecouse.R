################################################################################
##########################- DESEQ2- ##########################


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
  if (any(grepl("biotype|symbol|description|external_gene_name|gene_id", colnames(counts)))){
    counts <- counts[, -(grep("biotype|symbol|description|external_gene_name|gene_id", colnames(counts)))]
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
setwd(dirname(workingD))
Directory_for_the_day <- paste0("DESEQ2_Cell_lines_treatment_model_~time",Sys.Date())
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
#resumen <- resumen[match (colnames(counts), resumen$SampleID),]
if (all(colnames(counts) == batch$SampleID)){print ("Los nombres del batch están en el orden correcto")} else {
  cat ("Nombres no ordenados correctamente:", batch$SampleID[!colnames(counts) == batch$SampleID], "\n")
  stop ("Los nombres del batch no están en el orden correcto")}

counts <- counts[, match (resumen$SampleID, colnames(counts))]
if (all(colnames(counts) == resumen$SampleID)){print ("Los nombres del resumen están en el orden correcto")} else {
  cat ("Nombres no ordenados correctamente:",resumen$SampleID[!colnames(counts) == resumen$SampleID], "\n")
  stop ("Los nombres del resumen no están en el orden correcto")}
coldata <- data.frame(SampleID = colnames(counts), 
                      replicate = as.factor(gsub(".*_", "",resumen$SampleID)), 
                      subtype = as.factor(resumen$Group), 
                      time = as.factor(resumen$Time))
rownames(coldata) <- coldata$SampleID

cat("\n \n ............................CHEQUEO DE NOMBRES Y ORDEN FINALIZADA...................................\n")


################# - DESEQ2 NORMALIZATION  - ################################

cat("\n \n ............................INICIACION DE DESEQ2...................................\n")

## Execute the normalization
COMBAT_NO_NORMALIZED <-counts
coldata$subtype=factor(coldata$subtype,levels=c("SUDHL6","U2932", "WSUNHL"))
coldata$time=factor(coldata$time,levels=c(0,6, 24))

## Perform the DESEQ2
if(!all(sort(colnames(COMBAT_NO_NORMALIZED)) == sort(rownames(coldata)))) {stop ("Los nombres de los samples del coldata y datos de COMBAT no coinciden")}
if(!all(colnames(COMBAT_NO_NORMALIZED) == rownames(coldata))) {stop ("Los nombres de los samples del coldata y datos de COMBAT no están en el orden correcto")}
ddsWSU <- DESeqDataSetFromMatrix(countData = round(COMBAT_NO_NORMALIZED[,-grep ("SUDHL6|U2932", colnames(COMBAT_NO_NORMALIZED))]),
                              colData =  coldata[-grep ("SUDHL6|U2932", coldata$SampleID),],
                              design = ~ replicate +time)

ddsSUD <- DESeqDataSetFromMatrix(countData = round(COMBAT_NO_NORMALIZED[,-grep ("WSUNHL|U2932", colnames(COMBAT_NO_NORMALIZED))]),
                                 colData =  coldata[-grep ("WSUNHL|U2932", coldata$SampleID),],
                                 design = ~ replicate+time)

ddsU2932 <- DESeqDataSetFromMatrix(countData = round(COMBAT_NO_NORMALIZED[,-grep ("SUDHL6|WSUNHL", colnames(COMBAT_NO_NORMALIZED))]),
                                 colData =  coldata[-grep ("SUDHL6|WSUNHL", coldata$SampleID),],
                                 design = ~ replicate +time)

ddsall<- DESeqDataSetFromMatrix(countData = round(COMBAT_NO_NORMALIZED),
                                   colData =  coldata, 
                                   design = ~ subtype + time)


## Here the code to extract 100%% of counts
ddsWSU = DESeq(ddsWSU)
ddsSUD = DESeq(ddsSUD)
ddsU2932 = DESeq(ddsU2932)
ddsall = DESeq(ddsall)

writeLines(capture.output(sessionInfo()), paste0("Session_info_",  Sys.Date(), ".tsv"))



######################### - Comparissons- ################################
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl = useEnsembl(biomart = "genes") ## Define the dataet to use (genes)
searchDatasets(mart =  ensembl, pattern = "hsapiens") ## Look for the sensembl organism anotation
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
plots_PCA_DESEQ = function(dat, col, titulo){
  da <- dat
  vector1 <- c(nrow (assay(da)), 1000,500,100)
  lista <- list ()
  for (i in 1:length(vector1)){
    da <-  plotPCA(dat, intgroup = "time", ntop = vector1[i], returnData = T)
    da$titulo <- paste0("PCA VST DESeq2 ", titulo, " top ", vector1[i], " genes")
    lista[[i]] <- ggplot(da, aes(x = PC1, y = PC2, color = time, label = col$SampleID, group = as.factor(time), shape = coldata$subtype)) +
      geom_point(size = 2) +  # Tamaño de los puntos
      labs(color = "Group") + 
      xlab(paste0("PC1: ", round(100 * attr(da, "percentVar")[1], 2), "% variance")) +
      ylab(paste0("PC2: ", round(100 * attr(da, "percentVar")[2], 2), "% variance")) +
      geom_text(vjust = -1, size = 2.5) + 
      facet_grid(. ~ titulo) + 
      theme(strip.background =   element_rect(color = "black", size = 1),
            panel.border =  element_rect(color = "black", fill = NA, size = 1),
            strip.text.x = element_text(colour = "black", face = "bold"),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),  # Cuadrícula mayor
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey90"))  # Cuadrícula menor
  }
  grid.arrange(grobs = lista[c(1:4)], ncol = 2, nrow = 2)
}
plots_PCA_DESEQ_spanish = function(dat, col, titulo){
  da <- dat
  vector1 <- c(nrow (assay(da)), 1000,500,100)
  lista <- list ()
  for (i in 1:length(vector1)){
    da <-  plotPCA(dat, intgroup = "time", ntop = vector1[i], returnData = T)
    da$titulo <- paste0("PCA VST DESeq2 ", titulo, " top ", vector1[i], " genes")
    lista[[i]] <- ggplot(da, aes(x = PC1, y = PC2, color = time, label = col$SampleID, group = as.factor(time), shape = coldata$subtype)) +
      geom_point(size = 2) +  # Tamaño de los puntos
      labs(color = "Línea celular") + 
      xlab(paste0("PC1: ", round(100 * attr(da, "percentVar")[1], 2), "% varianza")) +
      ylab(paste0("PC2: ", round(100 * attr(da, "percentVar")[2], 2), "% varianza")) +
      geom_text(vjust = -1, size = 2.5) + 
      facet_grid(. ~ titulo) + 
      theme_bw() + 
      theme(strip.background =   element_rect(color = "black", size = 1),
            panel.border =  element_rect(color = "black", fill = NA, size = 1),
            strip.text.x = element_text(colour = "black", face = "bold", size = 14),
            strip.text.y = element_text(colour = "black", face = "bold", size = 14),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),  # Cuadrícula mayor
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey90"))  # Cuadrícula menor
  }
  grid.arrange(grobs = lista[c(1:4)], ncol = 2, nrow = 2)
}
DEA = function (dds, grupo, treated, untreated, tipe, coldata){
  
  Directory_base <- getwd()
  
  ## Set the contrast  with a sig level of 0.05 and independent filtering
  res= results (dds, contrast = c(grupo, treated, untreated),
                               alpha = 0.05, independentFiltering = T)
  
  summary(res) ## Summary of statistics
  resultsa <- res
  res = data.frame(res) ## Save the results file
  
  ## Show if they are significance of not
  res$Significance = c()
  for (gen in 1:nrow(res)){
    if (is.na(res$pvalue[gen]) | is.na(res$padj[gen])){ res$Significance[gen] = NA}
    else if ((res$log2FoldChange[gen] >= 1 & res$padj[gen] <= 0.05) | (res$log2FoldChange[gen]<=-1 && res$padj[gen] <= 0.05)){
      if (res$log2FoldChange[gen] >= 1 & res$padj[gen] <= 0.05) {
        res$Significance[gen] = "Significant upregulated"
      }
      else if (res$log2FoldChange[gen] <=-1 & res$padj[gen] <= 0.05) {
        res$Significance[gen] = "Significant downregulated"
      }
    } 
    else {
      res$Significance[gen] = "No Significant"
    }
  }
 
  ## Info if they are considered significant or not
  table(res$Significance)
  
  ## Anotate the samples
  Annot_DESEQ2 <- Annotate(res)
  
  ## Get the counts from used samples to save together with the deseq2 Results
  
  ## Get the counts from used samples to save together with the deseq2 Results
  samples = paste0(coldata$SampleID[intersect(grep (paste0(treated,"|", untreated), coldata$time), 
                                              grep (tipe, coldata$subtype)) ], collapse = "|")
  normalized_filter = data.frame(counts (dds, normalize = T));
  colnames(normalized_filter) <- gsub ("\\.","-",colnames(normalized_filter))
  
  ## Combine the DESEQ2 results, and the normalize counts
  normalized_filter$ensembl_gene_id <- row.names(normalized_filter)
  Annot_DESEQ2= merge(Annot_DESEQ2, normalized_filter, by = "ensembl_gene_id", all.x = TRUE) 
  
  ## Save data into a directory
  dir.create(paste0(tipe, treated,"_vs_", untreated))
  setwd (paste0(tipe, treated,"_vs_", untreated))
  
  ## Save the table
  A = Annot_DESEQ2
  write.xlsx(A, paste0("DESEQ2_", treated, "_vs_", untreated, "_", Sys.Date(), ".xlsx")) ## All genes
  write.xlsx(subset(A, Significance %in% c("Significant downregulated", "Significant upregulated")),  ## Only significant
             paste0("DESEQ2_", treated, "_vs_", untreated, "_", "_Significant_", Sys.Date(), ".xlsx"))
  
  ## Save the summary
  writeLines(capture.output(summary(resultsa), "\n", coldata), paste0("Summary",treated, "_vs_",untreated, "_", ".txt"))
  
  ## Save the general plots
  par (mfrow = c(1,2))
  pdf(file = paste0("DESEQ2_plotMA_&_dispEST", treated, "_vs_", untreated, "_", Sys.Date(), ".pdf"));
  DESeq2::plotMA(resultsa, alpha = 0.05);
  plotDispEsts(dds); dev.off()
  par (mfrow = c(1,1))
  
  A$Significance = as.factor (A$Significance)
  pdf(file = paste0("DESEQ2_Volcano_plot", treated, "_vs_", untreated, "_", ".pdf"))
  a = ggplot(data = na.omit(A), 
             aes(x = log2FoldChange, y = -log10(padj), color = Significance)) + 
    geom_point() + 
    geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed', lwd = 0.75) + 
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed', lwd = 0.75) + 
    theme_light() + 
    scale_color_manual(values = c("Significant upregulated" = "red", 
                                  "Significant downregulated" = "navy", 
                                  "No Significant" = "grey")) + 
    ggtitle(paste(treated, "vs", untreated)) +  # Título del gráfico
    labs(color = "Regulation",
         x = expression("log"[2]*"FC"), 
         y = expression("-log"[10]*"padj"))
  print (a)
  dev.off()
  setwd(Directory_base )
  return (subset(A, Significance %in% c("Significant downregulated", "Significant upregulated")))
}

## Perform the comparisson for each cell line
setwd(Directory_for_the_day)
res_WSU_24h <-DEA (ddsWSU, "time", "24", "0", tipe = "WSUNHL", ddsWSU@colData)
res_WSU_6h <-DEA (ddsWSU, "time", "6", "0", tipe = "WSUNHL", ddsWSU@colData)
res_WSU_24_6h <-DEA (ddsWSU, "time", "24", "6", tipe = "WSUNHL", ddsWSU@colData)
lista_WSU<- list (res_WSU_24h, res_WSU_6h, res_WSU_24_6h)
names(lista_WSU)<- c("0 vs 24h","0 vs 6h","6 vs 24h")

res_SUDHL6_24h <-DEA (ddsSUD, "time", "24", "0", tipe = "SUDHL6", ddsSUD@colData)
res_SUDHL6_6h <-DEA (ddsSUD, "time", "6", "0", tipe = "SUDHL6", ddsSUD@colData)
res_SUDHL6_24_6h <-DEA (ddsSUD, "time", "24", "6" , tipe = "SUDHL6", ddsSUD@colData)
lista_SUDHL6<- list (res_SUDHL6_24h, res_SUDHL6_6h, res_SUDHL6_24_6h)
names(lista_SUDHL6)<- c("0 vs 24h","0 vs 6h","6 vs 24h")

res_U2932_24h <-DEA (ddsU2932, "time", "24", "0" , tipe = "U2932", ddsU2932@colData)
res_U2932_6h <-DEA (ddsU2932, "time", "6", "0" , tipe = "U2932", ddsU2932@colData)
res_U2932_24_6h <-DEA (ddsU2932, "time", "24", "6" , tipe = "U2932", ddsU2932@colData)
lista_U2932<- list (res_U2932_24h, res_U2932_6h, res_U2932_24_6h)
names(lista_U2932)<- c("0 vs 24h","0 vs 6h","6 vs 24h")

######################### - Complex Heatmaps- ################################

## Reorder the resumnen and data file 
lines <- c("WSU", "U2932", "SUDHL")
time <- c("DMSO", "_6H_", "_24H_")
orden <- c()
for (i in 1:length(lines)){
  for (j in 1:length(time)){
    orden <- append(orden, resumen$SampleID[grep (paste0(lines[i], "*.*", time[j]), resumen$SampleID)])
  }
}

resumen <- resumen[match(orden, resumen$SampleID),]

### Only three conditions 

pheatmap_only_3 <- function(lista, dds){
  for (i in 1:length(lista)){
    if (i ==1){a <- lista[[i]]
    } else {
      a <- rbind (a, lista[[i]])
    }
  }

  c <- data.frame(assay(vst(dds)))
  c<-c[a$ensembl_gene_id,]
  b <- data.frame(dds@colData)
  rownames(b) <- colnames(c)
  
  ## Rearrange the order
  c<-c[,c(grep ("DMSO", colnames(c)), grep ("6H", colnames(c)), grep ("24H", colnames(c)))]
  b <- b[c(grep ("DMSO", colnames(c)), grep ("6H", colnames(c)), grep ("24H", colnames(c))),]
  split <- factor(c(rep ("0",3), rep ("6",3), rep ("24",3)), levels = c("0", "6", "24"), labels = c("0", "6", "24"))
  
  a <- t(scale(t(c)))
  library(ComplexHeatmap)
  
  top_annotation  <- HeatmapAnnotation(
    simple_anno_size = unit(0.35, "cm"), 
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  
    gp = gpar(fontsize = 10), border = T, 
    Tiempo= split,
    col = list(
      Tiempo = c("0" = colorRampPalette(c("lightblue", "darkblue"))(3)[1], "6" = colorRampPalette(c("lightblue", "darkblue"))(3)[2],
                 "24" = colorRampPalette(c("lightblue", "darkblue"))(3)[3])
    )
  )
  Heatmap(a,
          column_title = NULL, column_names_rot = 60,
          column_names_gp = gpar (fontsize = 10), cluster_columns = F, 
          show_row_names = F,top_annotation = top_annotation,
          name = "Z-score", column_split = split, cluster_column_slices = F)
}
tiff(paste0("WSU-NHL_allsig", Sys.Date(), ".tiff"), res = 150, height = 15, width = 20, unit = "cm")
pheatmap_only_3 (lista_WSU, ddsWSU); dev.off()
tiff(paste0("SUDHL6L_allsig", Sys.Date(), ".tiff"), res = 150, height = 15, width = 20, unit = "cm")
pheatmap_only_3 (lista_SUDHL6, ddsSUD); dev.off()
tiff(paste0("U2932_allsig", Sys.Date(), ".tiff"), res = 150, height = 15, width = 20, unit = "cm")
pheatmap_only_3 (lista_U2932, ddsU2932); dev.off()

######################### - Print common genes between conditions- ################################

## Save common between 24 ans 6 hours in all conditions
save <-function (x, y, titulo){
common<-intersect(x$ensembl_gene_id, y$ensembl_gene_id)
save_results <-data.frame(matrix (ncol = 14, nrow = 1))
colnames(save_results)<- c("Comparisson", "Number_total_genes", "Number_up", "%Number_up", "Number_down", "%_Number_down", 
                           "Number common", "% Common genes from total", "n Same direction", " % Same direction","Number_up_common", "%_up_common", "Number_down_common", "% down common")

all <- data.frame(table (x$Significance[match (common, x$ensembl_gene_id)] == "Significant upregulated" & y$Significance[match (common, y$ensembl_gene_id)] == "Significant upregulated"))
npos <-  all$Freq[all$Var1 == T]
porpos <- npos /length(common)
all <- data.frame(table (x$Significance[match (common, x$ensembl_gene_id)] == "Significant downregulated" & y$Significance[match (common, y$ensembl_gene_id)] == "Significant downregulated"))
if (length(all$Freq[all$Var1 == F]) != length(all$Freq)){nneg<-  all$Freq[all$Var1 == T]; porneg <- nneg/length(common)} else {
  nneg<- 0
  porneg <- 0
}

save_results[1,] <- unlist(c(titulo, 
                             length(x$ensembl_gene_id), 
                             length(which(x$Significance == "Significant upregulated")), 
                             length(which(x$Significance == "Significant upregulated"))/(length(x$ensembl_gene_id)), 
                             length(which(x$Significance == "Significant downregulated")), 
                             length(which(x$Significance == "Significant downregulated"))/(length(x$ensembl_gene_id)), 
                             length(common), 
                             length(common) /(length(unique(c(x$ensembl_gene_id, y$ensembl_gene_id)))),
                             (npos + nneg), ((npos + nneg) /length (common)),
                             npos, porpos, 
                             nneg, porneg))
return (save_results)
}
common_24_6h <- rbind(
rbind(save (res_WSU_24h, res_WSU_6h, "res_WSU_24h"), save (res_WSU_6h, res_WSU_24h, "res_WSU_6H")),
rbind(save (res_U2932_24h, res_U2932_6h, "res_U2932_24h"), save (res_U2932_6h, res_U2932_24h, "res_U2932_6h")), 
rbind(save (res_SUDHL6_24h, res_SUDHL6_6h, "res_SUDHL6_24h"), save (res_SUDHL6_6h, res_SUDHL6_24h, "res_SUDHL6_6h"))
)
write.xlsx(common_24_6h, paste0("COmmon_0vs24_&_0_vs6h", Sys.Date(), ".xlsx"))

######################### - plotPCA- ################################

COMBAT_NO_NORMALIZED <-Annotate(COMBAT_NO_NORMALIZED)



## Svae graphs in pdf
pdf (file = paste0("PCA_plots_Base_CB103", Sys.Date(), ".pdf"), width  = 15, height = 12)
draw_function ("U-2932", ddsU2932,  ddsU2932@colData)
draw_function ("SU-DHL-6",ddsSUD,  ddsSUD@colData) 
draw_function ("WSU-NHL",ddsWSU,  ddsWSU@colData)  
draw_function ("ALL",ddsall,  ddsall@colData) 
dev.off()


## Plots en español
plots_PCA_DESEQ_Spanish= function(dat, col, titulo){
  da <- dat
  vector1 <- c(nrow (assay(da)))
  lista <- list ()
  for (i in 1:length(vector1)){
    da <-  plotPCA(dat, intgroup = "subtype", ntop = vector1[i], returnData = T)
    da$titulo <- paste0("PCA ", vector1[i], " genes","\n", titulo)
    lista[[i]] <- ggplot(da, aes(x = PC1, y = PC2, color = subtype, label = col$SampleID, group = as.factor(subtype))) +
      geom_point(size = 2) +  # Tamaño de los puntos
      labs(color = "Línea celular") + 
      scale_color_manual(values = c("SUDHL6" = "#A0522D", "U2932"= "#D02F4B", "WSUNHL" = "lightgreen")) +
      xlab(paste0("PC1: ", round(100 * attr(da, "percentVar")[1], 2), "% varianza")) +
      ylab(paste0("PC2: ", round(100 * attr(da, "percentVar")[2], 2), "% varianza")) +
      geom_text(vjust = -1, size = 2.5) + 
      facet_grid(. ~ titulo) + 
      theme_bw() + 
      ylim (-30,30)+
      xlim (-25,50)+
      theme(strip.background =   element_rect(color = "grey", size = 1),  
            legend.title = element_text(face = "bold"),
            axis.title.y = element_text(colour = "black", face = "bold", size = 10),
            axis.title.x = element_text(colour = "black", face = "bold", size = 10),
            strip.text.x = element_text(colour = "black", face = "bold", size = 10))
  }
  print (lista[[1]])
}
plots_PCA_DESEQ_Spanish_time= function(dat, col, titulo){
  da <- dat
  vector1 <- c(nrow (assay(da)))
  lista <- list ()
  for (i in 1:length(vector1)){
    da <-  plotPCA(dat, intgroup = "time", ntop = vector1[i], returnData = T)
    da$titulo <- paste0("PCA ", vector1[i], " genes","\n", titulo)
    lista[[i]] <- ggplot(da, aes(x = PC1, y = PC2, color = time, label = col$SampleID, group = as.factor(time))) +
      geom_point(size = 2) +  # Tamaño de los puntos
      labs(color = "Tiempo") + 
      xlab(paste0("PC1: ", round(100 * attr(da, "percentVar")[1], 2), "% varianza")) +
      ylab(paste0("PC2: ", round(100 * attr(da, "percentVar")[2], 2), "% varianza")) +
      geom_text(vjust = -1, size = 2.5) + 
      facet_grid(. ~ titulo) + 
      theme_bw() + 
      theme(strip.background =   element_rect(color = "grey", size = 1),  
            legend.title = element_text(face = "bold"),
            axis.title.y = element_text(colour = "black", face = "bold", size = 10),
            axis.title.x = element_text(colour = "black", face = "bold", size = 10),
            strip.text.x = element_text(colour = "black", face = "bold", size = 10))
  }
  print (lista[[1]])
}
draw <- function (ddsall, col, titulo, fun){
VST_U2932<- vst(ddsall)
VST_U2932_filer <- VST_U2932[rowSums(counts(ddsall)>= 10) >=2,]
protinfilter <-VST_U2932_filer[na.omit(match(COMBAT_NO_NORMALIZED$ensembl_gene_id[COMBAT_NO_NORMALIZED$gene_biotype == "protein_coding"], rownames(VST_U2932_filer)))]
tiff(paste0 (titulo, Sys.Date(), ".tiff"), res = 150, height = 12, width = 14, unit = "cm")
plots_PCA_DESEQ_Spanish (protinfilter,  protinfilter@colData,  "codificadores de proteínas filtrado un 20%")
dev.off()
}

draw(ddsall, ddsall@colData, "PCA_3_3cell_lines", plots_PCA_DESEQ_Spanish)
draw(ddsWSU, ddsWSU@colData, "PCA_3_WSU", plots_PCA_DESEQ_Spanish_time )
draw(ddsU2932, ddsU2932@colData, "PCA_3_U2932",plots_PCA_DESEQ_Spanish_time )
draw(ddsSUD, ddsSUD@colData, "PCA_3_SUDHL6",plots_PCA_DESEQ_Spanish_time )



######################### - Individual volcano plots- ################################
setwd(dirname(workingD))
setwd(Directory_for_the_day)

DEA = function (dds, grupo, treated, untreated, tipe, coldata, limx, limy){

  ## Set the contrast  with a sig level of 0.05 and independent filtering
  res= results (dds, contrast = c(grupo, treated, untreated),
                alpha = 0.05, independentFiltering = T)
  
  resultsa <- res
  res = data.frame(res) ## Save the results file
  
  ## Show if they are significance of not
  res$Significance = c()
  for (gen in 1:nrow(res)){
    if (is.na(res$pvalue[gen]) | is.na(res$padj[gen])){ res$Significance[gen] = NA}
    else if ((res$log2FoldChange[gen] >= 1 & res$padj[gen] <= 0.05) | (res$log2FoldChange[gen]<=-1 && res$padj[gen] <= 0.05)){
      if (res$log2FoldChange[gen] >= 1 & res$padj[gen] <= 0.05) {
        res$Significance[gen] = "Significant upregulated"
      }
      else if (res$log2FoldChange[gen] <=-1 & res$padj[gen] <= 0.05) {
        res$Significance[gen] = "Significant downregulated"
      }
    } 
    else {
      res$Significance[gen] = "No Significant"
    }
  }
  
  ## Anotate the samples
  Annot_DESEQ2 <- Annotate(res)
  library(ggrepel)
  
  ## Get the counts from used samples to save together with the deseq2 Results
  samples = paste0(coldata$SampleID[intersect(grep (paste0(treated,"|", untreated), coldata$time), 
                                              grep (tipe, coldata$subtype)) ], collapse = "|")
  normalized_filter = data.frame(counts (dds, normalize = T));
  colnames(normalized_filter) <- gsub ("\\.","-",colnames(normalized_filter))
  
  ## Combine the DESEQ2 results, and the normalize counts
  normalized_filter$ensembl_gene_id <- row.names(normalized_filter)
  Annot_DESEQ2= merge(Annot_DESEQ2, normalized_filter, by = "ensembl_gene_id", all.x = TRUE) 
  
  ## Save the table
  A = Annot_DESEQ2
  A$Significance = as.factor (A$Significance)
  
  tiff(file = paste0("DESEQ2_Volcano_plot", tipe,treated, "_vs_", untreated, "_", ".tiff"), res = 150, unit = "cm", width = 15, height = 14)
  a <- ggplot(data = na.omit(A), 
             aes(x = log2FoldChange, y = -log10(padj),color =  Significance)) + 
    geom_point(size = 1) + 
    xlim (-limx, limx) +
    geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed', lwd = 0.75) + 
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed', lwd = 0.75) + 
    theme_bw() +
    ylim (-1, limy) +
    scale_color_manual(values = c("Significant upregulated" = "red", 
                                  "Significant downregulated" = "navy", 
                                  "No Significant" = "grey"), 
                       labels = c("Significant upregulated" = "Significativo aumentado" , 
                                  "Significant downregulated" = "Significativo disminuido", 
                                  "No Significant" = "No significativo")) + 
    labs(color = "Significación",
         x = expression("log"[2]*"FC"), 
         y = expression("-log"[10]*"padj")) +
    #geom_text_repel(size = 2.7, max.overlaps =Inf, hjust = 1, vjust = 1) +
    theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1), 
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 11), 
          strip.text = element_text(color = "black", face = "bold", size =12), 
          strip.background = element_rect(color = "black", fill = "grey")) +
    facet_grid(. ~  paste (tipe, "\n",treated , "vs", untreated))

  print (a)
  dev.off()
}
## Perform the comparisson for each cell line
data <- rbind(res_WSU_24_6h, res_WSU_24h, res_WSU_6h)
limx <- max(data$log2FoldChange, max(abs(data$log2FoldChange), na.rm = T))
limy <- (-log10(min(data$padj, na.rm = T)))
DEA (ddsWSU, "time", "24", "0", tipe = "WSU-NHL", ddsWSU@colData, limx, limy)
DEA (ddsWSU, "time", "6", "0", tipe = "WSU-NHL", ddsWSU@colData, limx, limy)
DEA (ddsWSU, "time", "24", "6", tipe = "WSU-NHL", ddsWSU@colData, limx, limy)


data <- rbind(res_SUDHL6_24h, res_SUDHL6_24_6h, res_SUDHL6_6h)
limx <- max(data$log2FoldChange, max(abs(data$log2FoldChange))); limy <- -log10(min(data$pvalue))
DEA (ddsSUD, "time", "24", "0", tipe = "SU-DHL-6", ddsSUD@colData, limx, limy)
DEA (ddsSUD, "time", "6", "0", tipe = "SU-DHL-6", ddsSUD@colData, limx, limy)
DEA (ddsSUD, "time", "24", "6" , tipe = "SU-DHL-6", ddsSUD@colData, limx, limy)


data <- rbind(res_U2932_24_6h, res_U2932_24h, res_U2932_6h)
limx <- max(data$log2FoldChange, max(abs(data$log2FoldChange))); limy <- -log10(min(data$pvalue))
DEA (ddsU2932, "time", "24", "0" , tipe = "U2932", ddsU2932@colData, limx, limy)
DEA (ddsU2932, "time", "6", "0" , tipe = "U2932", ddsU2932@colData, limx, limy)
DEA (ddsU2932, "time", "24", "6" , tipe = "U2932", ddsU2932@colData, limx, limy)

