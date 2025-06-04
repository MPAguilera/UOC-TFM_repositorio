##########################- Correlacion Spearman IC50- ##########################

################# - Upload libraries - #########################################
set.seed(1299)
## Here package to be used from CRAN
packages<- c("dplyr","rstudioapi","ggplot2","gridExtra","openxlsx","readxl","BiocManager", "corrplot", "pheatmap", "factoextra","RColorBrewer", "ashr")

## Here package to be used from BiocManager
biocM <- c("sva", "DESeq2", "edgeR", "biomaRt"); combine <- c(packages, biocM)

## If they are not already installed install them
for (i in packages){ if (!requireNamespace(i, quietly = TRUE)){install.packages(i)}}
for (i in biocM){if (!requireNamespace(i, quietly = TRUE)){BiocManager::install(i)}}

##Upload libraries
for (i in combine){suppressMessages(suppressWarnings(library(i, character.only = TRUE)))}
rm (packages, biocM, combine, i)

################# - Upload and rearrange names - #################################

## Mensaje error
mensaje_error <- c("...................... ERROR ............")

## Set WD to file location
workingD <- rstudioapi::getActiveDocumentContext()$pat
setwd(dirname(workingD))

## Read files
files <- list.files (recursive = T, full.names = T)
IC50 <- read.xlsx(files[grep ("Input*.*IC50", files)])
datexprs <- read.xlsx(files[grep ("Input*.*VST_normalized", files)]); 
rownames(datexprs) <- datexprs$ensembl_gene_id
datexprs <- datexprs[,-1]

## Make sure than the names of samples are the same
colnames(datexprs) <- gsub ("\\.", "-", colnames(datexprs)) %>% gsub ("OCI-LY", "OCI-Ly", .)
if (!all(IC50$X1 %in% colnames(datexprs))){stop("Nombres diferentes entre IC50 y exprs")}

## Rearrange the IC50 and exprs names to be the same and in the same order
IC50 <-IC50[order(IC50$Droga1),]

original<- datexprs[, c(1:4, match(IC50$X1, colnames(datexprs)))]
datexprs <-datexprs[, c(2, match (IC50$X1,colnames(datexprs) ))]

## Create the daty WD
name_dir<- paste0("Correlation_IC50_&_gene_expression_ALL", Sys.Date())
dir.create(name_dir)

################# - Spearman correlation - #################################

## Correlation between IC50 value and gene expression
corre <- data.frame(matrix(ncol = 4, nrow = nrow(original)))
for (i in 1:nrow(original[,-c(1:4)])){
  gen <-c(as.numeric(original[i,-c (1:4)]))
  a<-cor.test(gen, IC50$Droga1, method = "spearman")
  corre[i, c(1:2)]<- c(a$p.value, a$estimate)
  rownames(corre)[i] <-original[i, 1]
  corre[i, (3:4)]<- original[i,c (1:2)]
}

colnames(corre) <- c("pvalue", "rho", "ensembl_gene_id", "external_gene_name")
corre <-data.frame(na.omit(corre))

## Graphic draw the expression and correlation of 0.85 and 0.65
resul_cor2 <- corre[(corre[,2]>=0.85 & corre[,1] <=0.05) | (corre[,2]<=(-0.85) & corre[,1] <=0.05),]
resul_cor2 <- resul_cor2[!duplicated(rownames(resul_cor2)), ]

draw<- original[match (resul_cor2$ensembl_gene_id, original$ensembl_gene_id),]
draw$external_gene_name[draw$external_gene_name == ""] <- draw$ensembl_gene_id[draw$external_gene_name == ""]
rownames(draw) <- draw$external_gene_name

resul_cor3 <- corre[(corre[,2]>=0.65 & corre[,1] <=0.05) | (corre[,2]<=(-0.65) & corre[,1] <=0.05),]
resul_cor3 <- resul_cor3[!duplicated(rownames(resul_cor2)), ]
draw2<- original[match (resul_cor3$ensembl_gene_id, original$ensembl_gene_id),]
draw2$external_gene_name[draw2$external_gene_name == ""] <- draw2$ensembl_gene_id[draw2$external_gene_name == ""]
rownames(draw2) <- draw2$external_gene_name


library(ComplexHeatmap)
library(circlize)
top_annotation  <- HeatmapAnnotation(
  simple_anno_size = unit(0.35, "cm"), 
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  # Cambiar el tamaño de la fuente del título de la anotación
  gp = gpar(fontsize = 10), border = T, 
  `IC50 CB-103`= IC50$Droga1,
  col = list(`IC50 CB-103` =  colorRamp2 (c(min (IC50$Droga1), max(IC50$Droga1)), c("white", "mediumpurple")))
  )
bottom_annotation  <- HeatmapAnnotation(
  simple_anno_size = unit(0.35, "cm"), 
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  # Cambiar el tamaño de la fuente del título de la anotación
  gp = gpar(fontsize = 10), border = T, 
  `Clasificador 2-S`= IC50$`2-S`,
  col = list( `Clasificador 2-S` = c ("EZB" = "#d17944", "EZB/MYC+"="brown","Other" = "grey", "ST2"= "#D02F4B", "BN2" = "purple", "MCD" = "lightblue", "N1" = "lightgreen"))
)

tiff(paste0("Spearman_0.85_cluster_Complex_Heatmap", Sys.Date(), ".tiff"), res = 150, height = 20, width = 15, unit = "cm")
scale_0.85 <- t(scale(t(draw[,-c(1:4)])))
Heatmap (scale_0.85, name = "Z-score", column_names_rot = 60, 
         row_names_gp = gpar (fontface = "italic", fontsize = 8),
         top_annotation = top_annotation, bottom_annotation = bottom_annotation,
         column_names_gp = gpar (fontface = "bold", fontsize = 8)
         ); dev.off()

tiff(paste0("Spearman_0.85_nocluster_Complex_Heatmap", Sys.Date(), ".tiff"), res = 150, height = 20, width = 15, unit = "cm")
Heatmap (scale_0.85, name = "Z-score", column_names_rot = 60, 
         row_names_gp = gpar (fontface = "italic", fontsize = 8),
         top_annotation = top_annotation,  bottom_annotation = bottom_annotation,
         column_names_gp = gpar (fontface = "bold", fontsize = 8), 
         cluster_columns = F); dev.off()

tiff(paste0("Spearman_0.65_cluster_Complex_Heatmap", Sys.Date(), ".tiff"), res = 150, height = 25, width = 15, unit = "cm")
scale_0.85 <- t(scale(t(draw2[,-c(1:4)])))
Heatmap (scale_0.85, name = "Z-score", column_names_rot = 60, show_row_names = F,
         row_names_gp = gpar (fontface = "italic", fontsize = 8),
         top_annotation = top_annotation,  bottom_annotation = bottom_annotation,
         column_names_gp = gpar (fontface = "bold", fontsize = 8)
); dev.off()

tiff(paste0("Spearman_0.65_nocluster_Complex_Heatmap", Sys.Date(), ".tiff"), res = 150, height = 25, width = 15, unit = "cm")
Heatmap (scale_0.85, name = "Z-score", column_names_rot = 60,  show_row_names = F,
         row_names_gp = gpar (fontface = "italic", fontsize = 8),
         top_annotation = top_annotation,  bottom_annotation = bottom_annotation,
         column_names_gp = gpar (fontface = "bold", fontsize = 8), 
         cluster_columns = F); dev.off()


## Save the ORA Results positive
args <- list()  # Resetear args
args[[1]] <- paste0("ORA_", Sys.Date())
args[[2]] <- paste0("ORA_", "Correlation_0.85pos")
args[[3]] <- rownames(resul_cor2[resul_cor2$rho>=0,])

## Run ORA
library(biomaRt)
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
setwd(dirname(workingD))
setwd(name_dir)
source("G:/Mi unidad/Linfomas/Scripts_automatizados/Gene_Ontology_ORA.R")


setwd(dirname(workingD))
setwd(name_dir)
## Save the ORA Results positive
args <- list()  # Resetear args
args[[1]] <- paste0("ORA_", Sys.Date())
args[[2]] <- paste0("ORA_", "Correlation_0.85neg")
args[[3]] <- rownames(resul_cor2[resul_cor2$rho<0,])

## Run ORA
library(biomaRt)
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

source("G:/Mi unidad/Linfomas/Scripts_automatizados/Gene_Ontology_ORA.R")

##Corelation up 0.65
setwd(dirname(workingD))
setwd(name_dir)
## Save the ORA Results
args <- list()  # Resetear args
args[[1]] <- paste0("ORA_", Sys.Date())
args[[2]] <- paste0("ORA_", "Correlation_0.65_pos")
args[[3]] <- rownames(corre)[corre$rho>=0.65 & corre$pvalue<0.05]
print(args[[2]])
source("Gene_Ontology_ORA.R")

## Correlation negative 0.65
setwd(dirname(workingD))
setwd(name_dir)

## Save the ORA Results
args <- list()  # Resetear args
args[[1]] <- paste0("ORA_", Sys.Date())
args[[2]] <- paste0("ORA_", "Correlation_0.65_neg")
args[[3]] <- rownames(corre)[corre$rho<=(-0.65) & corre$pvalue<0.05]
print(args[[2]])
source("Gene_Ontology_ORA.R")


################# - Save dataframe- #################################


## Save the correlation plots
resul_cor2$external_gene_name <- original$external_gene_name[match(rownames(resul_cor2), original$ensembl_gene_id)]
save_corre <- data.frame(corre)

wb = createWorkbook()
addWorksheet(wb, "Relevant_Spearman_Correlation") 
writeData(wb, sheet ="Relevant_Spearman_Correlation", data.frame(ensembl_gene_id = rownames(resul_cor2), resul_cor2))

addWorksheet(wb, "All_Spearman_Correlation") 
writeData(wb, sheet ="All_Spearman_Correlation", data.frame(ensembl_gene_id = rownames(save_corre), save_corre)) 

saveWorkbook(wb, file = paste0("Correlation_Droga1_IC50_gene_expression", ".xlsx"), overwrite = TRUE)
capture.output(Sys.info(), file = paste0("R-SystemInfo_", Sys.Date(), ".txt"))
capture.output(sessionInfo(), file = paste0("R-SessionInfo_", Sys.Date(), ".txt"))
