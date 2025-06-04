################################################################################
##########################- GSEA Graphs- #######################################
################################################################################


################################################################################
################# - Upload libraries - #########################################
################################################################################

if (require(rstudioapi)==F){install.packages("rstudioapi")}
if (require(readxl)==F){install.packages("readx")}
if (require(openxlsx)==F){install.packages("openxlsx")}
if (require(gridExtra)==F){install.packages("gridExtra")}
if (require(ggplot2)==F){install.packages("ggplot2")}
if (require(enrichplot)==F){BiocManager::install("enrichplot")}
if (require(clusterProfiler)==F){BiocManager::install("clusterProfiler")}
if (require(AnnotationDbi)==F){BiocManager::install("AnnotationDbi")}
if (require(DOSE)==F){iBiocManager::install("DOSE")}
if (require(dbplyr)==F){install.packages("dbplyr")}
library(ComplexHeatmap)
library(circlize)
###############################################################################
########################## - WORKING DIRECTORY - ###############################
###############################################################################

## Set WD to file location
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

## Create a directory for the results with the currentData
Directory_for_the_day <- paste0("ssGSEA_and_GSEA", Sys.Date())
dir.create(Directory_for_the_day)

################## -  Upload data into the environment - ######################

## Look for the GSEA results files
dirs = list.dirs(recursive = F)
dirs = dirs[grep ("GSEA_results_discrete", dirs)]
files = list.files (path = dirs, recursive = T)

## Get the table of GSEA
setwd(dirs)
files_gsea_result = files[grep ("gsea_report_for_*.*tsv", files)]
lista_gsea_results = list()
fil = 0

for (file in 1:length(files_gsea_result)){
  if (file%%2 == 0){
    fil = fil +1
    a =read.delim(files_gsea_result[file], header = T)
    b =read.delim(files_gsea_result[file -1], header = T) 
    lista_gsea_results[[fil]] = rbind (a, b)
    names(lista_gsea_results)[[fil]] = gsub (pattern = ".Gsea.*", "",files_gsea_result[file])
    rm (a,b)
  }
}

## Get the list of significative genes (p NOM value of 0.05 and FDRqvalue of 0.25)
lista_significativos = list(); max = 0;min = 0
for (file in 1:length(lista_gsea_results)){
  x = lista_gsea_results[[file]]
  x = x[x$FDR.q.val<=0.25  & x$NOM.p.val<=0.05,]
  lista_significativos[[file]] = x
  names(lista_significativos)[[file]] = names(lista_gsea_results)[[file]] 
  if (max<max(x$NES)){max = max(x$NES)}
  if (min>min(x$NES)){min = min(x$NES)}
  rm (x)
}

setwd(dirname(workingD))
clus <- c(4,4,5)
dirs <- list.dirs()
geneset <- c("Hallmarks","Lymphoma")
i = 1
setwd(dirname(workingD))
files <-list.files(path = dirs[grep ("ssGSEA_results", dirs)], full.names = T)
Hallmarks <- read.delim(files[grep (geneset[i], files)], skip = 2)

## Get the column with the time
Hallmarks <-data.frame(t(Hallmarks[,-2]))
colnames(Hallmarks) <- Hallmarks[1,]
Hallmarks <- Hallmarks[-1,]

Hallmarks$time <- rep(c(rep ("0", 3), rep ("6",3), rep ("24",3)),3)
Hallmarks$time <- factor (Hallmarks$time, levels = c("0", "6", "24"), labels =c("0", "6", "24") )
Hallmarks$Sample<- gsub ("_.*", "", rownames(Hallmarks))


#####################- ssGSEA- only ##########################
a <-data.frame(Sample =Hallmarks$Sample, time =Hallmarks$time, 
               apply(Hallmarks[, c(1:52)], MARGIN = 2, as.numeric))
rownames(a) <- rownames(Hallmarks)
rm(Hallmarks)

#####################- Combine sGSEA with GSEA results HEATMAP - ##########################

## Look for the significant genes
lista_guardar <- list()
pdf (paste0(geneset[i], "_individual_sig", Sys.Date(), ".pdf"), width = 10, height = 8)
for (j in 1:length(cell_lines)){
look <- grep (paste0(cell_lines[j],"*.*",geneset[i]), names(lista_significativos))

if(nrow(lista_significativos[[look[1]]]) !=0){d <-data.frame(time = "24vs0", lista_significativos[[look[1]]])}
if(nrow(lista_significativos[[look[2]]]) !=0){d <-rbind (d,data.frame(time = "24vs6", lista_significativos[[look[2]]]))}
if(nrow(lista_significativos[[look[3]]]) !=0){d <-rbind (d,data.frame(time = "6vs0", lista_significativos[[look[3]]]))}

## Draw a heatmap of all sig pathways
ssGSEA_all <- a[grep (cell_lines[j],rownames(a)),]
ssGSEA_all <- ssGSEA_all[,c(1:2, grep (paste0(unique(d$NAME), collapse = "|"), colnames(ssGSEA_all)))]


top_annotation  <- HeatmapAnnotation(
  simple_anno_size = unit(0.35, "cm"), 
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  # Cambiar el tamaño de la fuente del título de la anotación
  gp = gpar(fontsize = 10), border = T, 
  Tiempo= ssGSEA_all$time,
  col = list(
    Tiempo = c("0" = colorRampPalette(c("lightblue", "darkblue"))(3)[1], "6" = colorRampPalette(c("lightblue", "darkblue"))(3)[2],
               "24" = colorRampPalette(c("lightblue", "darkblue"))(3)[3])
  )
)

look <- grep (paste0(cell_lines[j],"*.*",geneset[i]), names(lista_gsea_results))
lista_plotting <- list(
  lista_gsea_results[[look[1]]][match(unique (d$NAME), lista_gsea_results[[look[1]]]$NAME),],
  lista_gsea_results[[look[2]]][match(unique (d$NAME), lista_gsea_results[[look[2]]]$NAME),],
  lista_gsea_results[[look[3]]][match(unique (d$NAME), lista_gsea_results[[look[3]]]$NAME),])

## Reorder for the anotation
for (path in 1:length(lista_plotting)){
  lista_plotting[[path]] <- lista_plotting[[path]][match (colnames(ssGSEA_all)[-c(1:2)], lista_plotting[[path]]$NAME),]
  lista_plotting[[path]]$sig <-lista_plotting[[path]]$FDR.q.val <=0.25 & lista_plotting[[path]]$NOM.p.val<= 0.05
  lista_plotting[[path]]$sig[lista_plotting[[path]]$sig == T]<- "*"
  lista_plotting[[path]]$sig[lista_plotting[[path]]$sig == F] <- NA
  }
path = 2
## Create the left annotation file
leftAno <- rowAnnotation( `6h vs 0h` = anno_simple(border = T,lista_plotting[[3]]$NES, pch= lista_plotting[[3]]$sig, 
                                               col = colorRamp2(c (-4,0,4), c("#06E9F9", "white", "#CF0178"))),
                         `24h vs 0h`= anno_simple(border = T,lista_plotting[[1]]$NES, pch= lista_plotting[[1]]$sig,
                                              col = colorRamp2(c (-4,0,4),  c("#06E9F9", "white", "#CF0178"))),
                         `24h vs 6h`= anno_simple(border = T, lista_plotting[[2]]$NES, pch= lista_plotting[[2]]$sig, 
                                              col = colorRamp2(c (-4,0,4),  c("#06E9F9", "white", "#CF0178"))), 
                         annotation_name_rot = 60, annotation_name_side = "bottom",
                         annotation_name_gp = gpar (fontface = "bold", fontsize = 10))

lgd_NES= Legend(title = "NES", col_fun = colorRamp2(c (-4,0,4), c("#06E9F9", "white", "#CF0178")), at = c(-4,-2, 0,2,4), 
                    labels = c("-4", "2", "0", "2", "4"))

   
## Print Heatmap
colnames(ssGSEA_all) <- gsub ("_", " ",colnames(ssGSEA_all)) %>% gsub ("HALLMARK","HM", .) %>% gsub ("BIOCARTA","BIO", .)
ssGSEA_all_draw <- t(scale(ssGSEA_all[,-c(1:2)]))
p <- Heatmap (ssGSEA_all_draw, name = "Z-score", 
         top_annotation = top_annotation, cluster_column_slices = F, cluster_columns = F,
         column_split = as.factor(ssGSEA_all$time), 
         column_names_rot = 60, column_names_gp = gpar(fontsize = 8),
         row_names_gp = gpar(fontsize = 8), 
         right_annotation = leftAno, column_title = paste(geneset[i], cell_lines[j])
         )

lista_guardar[[j]] <- draw(p, annotation_legend_list = list(lgd_NES))
}
dev.off()

for (j in 1:length(lista_guardar)){
  tiff(paste0(geneset[i], cell_lines[j],Sys.Date(), ".tiff"), res = 150, unit = "cm", width = 17, height = 15)
  print (lista_guardar[[j]])
  dev.off()
}

representar<-c()
mirar <- grep (geneset[i], names(lista_significativos))
for (s in mirar){
  representar <- append(representar,lista_significativos[[s]]$NAME)
}

ssGSEA_all <- a[,c(1:2, grep (paste0(unique(representar), collapse = "|"), colnames(a)))]

ssGSEA_all_draw <-t(scale(ssGSEA_all[,-c(1:2)] ))
Heatmap(ssGSEA_all_draw)

rm (p, ssGSEA_all, ssGSEA_all_draw, top_annotation, lgd_NES, leftAno,lista_plotting, lista_guardar)
