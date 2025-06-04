##########################- GSEA Graphs- #######################################
## Created by Marian Pérez on 15.05.2025

################# - Upload libraries - #########################################
if (require(rstudioapi)==F){install.packages("rstudioapi")}
if (require(readxl)==F){install.packages("readx")}
if (require(openxlsx)==F){install.packages("openxlsx")}
if (require(gridExtra)==F){install.packages("gridExtra")}
if (require(ggplot2)==F){install.packages("ggplot2")}
if (require(dbplyr)==F){install.packages("dbplyr")}
if (require(ComplexHeatmap)==F){install.packages("ComplexHeatmap")}
if (require(circlize)==F){install.packages("circlize")}


########################## - WORKING DIRECTORY - ###############################

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
lista_significativos = list()
nombres <-c()
for (file in 1:length(lista_gsea_results)){
  x = lista_gsea_results[[file]]
  x = x[x$FDR.q.val<=0.25  & x$NOM.p.val<=0.05,]
  lista_significativos[[file]]<- x
  names(lista_significativos)[file] <- names(lista_gsea_results)[file]
  nombres <- append (nombres, x$NAME)
  rm (x)
}

## Get the ssGSEA file
geneset <- c("Hallmarks")
i <-1
setwd(dirname(workingD))
dirs <- list.dirs()
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

## Make sure hallmarks is numeric
a <-data.frame(Sample =Hallmarks$Sample, time =Hallmarks$time, 
               apply(Hallmarks[, c(1:52)], MARGIN = 2, as.numeric))
rownames(a) <- rownames(Hallmarks)
rm(Hallmarks)

## Select the rows of interest
nombres  <- unique(nombres)
nombres <- nombres[grep ("HALLMARK", nombres)]
a <-a[,c(1:2,match (nombres, colnames(a)))]

## Draw a Heatmap of all Hallmarks pathways and all samples
annotation <- a[,c(1,2)]
rownames(annotation) <- rownames(a)
a$Sample<-rownames(a)

top_annotation  <- HeatmapAnnotation(
  simple_anno_size = unit(0.35, "cm"), 
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  # Cambiar el tamaño de la fuente del título de la anotación
  gp = gpar(fontsize = 10), border = T, 
  Tiempo= a$time,
  `Línea Celular` = gsub ("_DMSO.*", "",a$Sample) %>% gsub ("_CB103.*", "",.),
  col = list(
    `Línea Celular` = c("SU-DHL-6" = "#A0522D", "U-2932"= "#D02F4B", "WSU-NHL" = "lightgreen"),
    Tiempo = c("0" = colorRampPalette(c("lightblue", "darkblue"))(3)[1], "6" = colorRampPalette(c("lightblue", "darkblue"))(3)[2],
               "24" = colorRampPalette(c("lightblue", "darkblue"))(3)[3])
  )
)

## Look if they are sig or not 
WSUNHL <- c()
U2932 <- c()
SUDHL<- c()

look <- grep ("Hallmarks",names(lista_significativos))
lista_significativos <- lista_significativos[look]
n <- 1 
cell_lines <-c("WSU"  , "SUD"  , "U2932")
lista_rep <- list("WSU" = c(), "SUD" = c(), "U" = c())
for (i in 1:length(cell_lines)){
  t <-grep (cell_lines[i], names(lista_significativos))
  t <- rbind (lista_significativos[[t[1]]], lista_significativos[[t[2]]], lista_significativos[[t[3]]])
  
  for (j in 1:length(nombres)){
    if (any(grepl(nombres[j], t$NAME))){
      lista_rep[[i]][j] <- "*"
    } else {lista_rep[[i]][j] <- NA}
      
  }
}
a$Sample <-  gsub ("_DMSO.*", "",a$Sample) %>% gsub ("_CB103.*", "",.)
orden <- factor(paste0(a$Sample, a$time), levels = unique(paste0(a$Sample, a$time)), labels = unique(paste0(a$Sample, a$time)))

## Create the left annotation file
leftAno <- rowAnnotation( `WSU-NHL sig` = anno_simple(border = T,rep (0, length(lista_rep[[1]])), pch= lista_rep[[1]],  gp = gpar(col = "lightgreen")),
                         
                         `U-2932 sig`= anno_simple(border = T,rep (0, length(lista_rep[[1]])), pch= lista_rep[[3]],gp = gpar(col = "#D02F4B")), 
                         `SU-DHL-6 sig`= anno_simple(border = T,rep (0, length(lista_rep[[1]])), pch= lista_rep[[2]], gp = gpar(col = "#A0522D")),
                         annotation_name_rot = 60, annotation_name_side = "bottom",
                         annotation_name_gp = gpar (fontface = "bold", fontsize = 10))

## Print Heatmap
ssGSEA_all <- a[,-c(1:2)]
colnames(ssGSEA_all) <- gsub ("_", " ",colnames(ssGSEA_all)) %>% gsub ("HALLMARK","HM", .) %>% gsub ("BIOCARTA","BIO", .)
ssGSEA_all_draw <- t(scale(a[,-c(1:2)]))
Heatmap (ssGSEA_all_draw, name = "Z-score", 
         top_annotation = top_annotation, cluster_column_slices = F,
         column_names_rot = 60, column_names_gp = gpar(fontsize = 8),
         column_split = orden,  cluster_columns = F, cluster_row_slices = F,
         row_names_gp = gpar(fontsize = 8), right_annotation = leftAno)

#####################-  Heatmap NESS - ##########################
lista_gsea_results <-lista_gsea_results[grep ("Hallmarks", names(lista_gsea_results))]
representar <- data.frame (matrix(ncol = 18, nrow = length(nombres)))
colnames(representar) <- paste0(c("WSU 6vs0","WSU 24vs0","WSU 24vs6",
                           "SUD 6vs0","SUD 24vs0","SUD 24vs6",
                           "U2932 6vs0","U2932 24vs0","U2932 24vs6"), rep (c("NOM", "pvalue"), 9))
rownames(representar) <- nombres
tiempos <- c("6vs0", "24vs0", "24vs6")

ht <-Heatmap (ssGSEA_all_draw, name = "Z-score", column_title = NULL,
              top_annotation = top_annotation, cluster_column_slices = F,
              column_names_rot = 60, column_names_gp = gpar(fontsize = 8),
              column_split = c(rep (1,9), rep (2,9), rep (3,9)),  cluster_columns = F,
              cluster_row_slices = F, #row_split = orden_filas,
              row_names_gp = gpar(fontsize = 8), right_annotation = leftAno)

row_groups <- row_order(draw(ht))
cluster6_indices<- row_groups$`6`

#####################- Combine sGSEA with GSEA results ggplot - ##########################

## Make the table for ggplot
draw<- data.frame (matrix(ncol = 4, nrow = 0));
colnames(draw) <- c("Values", "Path", "Time", "cell_line")
look <- t(ssGSEA_all_draw[cluster6_indices,])
time <- gsub (".*CB103_","",rownames(look)) %>% gsub(".*DMSO.*", "0H", .) %>% gsub ("_.*","",.)
names <- gsub ("_.*","",rownames(look))
for (z in 1:ncol(look)){
  draw <- rbind (draw, data.frame (
    Path = rep (colnames(look)[z], nrow (look)), 
    Values = look[,z], 
    Time = time, 
    cell_line =names) )
  
}
draw$Time <- factor (draw$Time, labels = c("0H", "6H", "24H"), levels = c("0H", "6H", "24H"))

# cell line only clusters
settings <- list(labs(x= "Tiempo", y = "Valor ssGSEA normalizado", color = "Línea celular", fill = "Línea celular", group = "Línea celular"),
                 theme_bw(),
                 theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1), 
                       axis.title.x = element_text(size = 14, face = "bold"),
                       axis.title.y = element_text(size = 14, face = "bold"),
                       axis.text.x = element_text(size = 14, face = "bold"),
                       legend.title = element_text(face = "bold"),
                       strip.text = element_text(color = "black", face = "bold", size =14 ),
                       strip.background = element_rect(color = "black", fill = "grey"),
                       axis.text.y = element_text(size = 14, face = "bold")))


fig3 <-ggplot(draw, aes(x =Time, y = Values, color = as.factor(cell_line), group =as.factor(cell_line))) + 
  settings +
  geom_smooth(method = "loess", aes(fill = as.factor(cell_line)))+
  facet_grid(. ~  paste ("Modelo Loess",  ""))
ht


## Draw again heatmap but orden defined by me 
lista_rep[[1]] [is.na(lista_rep[[1]])] <- ""
lista_rep[[2]] [is.na(lista_rep[[2]])] <- ""
lista_rep[[3]] [is.na(lista_rep[[3]])] <- ""
orden_filas <- NA
orden_filas[lista_rep[[1]] == "*" &  lista_rep[[2]] == "*" & lista_rep[[3]] == "*"] <- 1
orden_filas[lista_rep[[2]] == "*"  & lista_rep[[3]] == "*" & is.na(orden_filas)] <- 2
orden_filas[lista_rep[[1]] == "*"  & lista_rep[[3]] == "*"  & is.na(orden_filas)] <- 3
orden_filas[lista_rep[[1]] == "*"  & lista_rep[[2]] == "*"  & is.na(orden_filas)] <- 4
orden_filas[lista_rep[[1]] == "*" & lista_rep[[2]] == ""  & lista_rep[[3]] == ""  & is.na(orden_filas)] <- 5
orden_filas[lista_rep[[1]] == ""  & lista_rep[[2]] == "*"  & lista_rep[[3]] == ""  & is.na(orden_filas)] <- 6
orden_filas[lista_rep[[1]] == ""  & lista_rep[[2]] == ""  & lista_rep[[3]] == "*"  & is.na(orden_filas)] <- 7




#####################- Combine sGSEA with GSEA results ggplot - ##########################

scale_WSU <- t(scale(ssGSEA_all[grep ("WSU",rownames(ssGSEA_all)),]))
scale_SUDHL<- t(scale(ssGSEA_all[grep ("SUDHL",rownames(ssGSEA_all)),]))
scale_U2932 <- t(scale(ssGSEA_all[grep ("U2932",rownames(ssGSEA_all)),]))

Heatmap (ssGSEA_all_draw, name = "Z-score", 
              top_annotation = top_annotation,
              cluster_column_slices = F,
              column_names_rot = 60, column_names_gp = gpar(fontsize = 8),
              column_split = c (rep(1,9), rep (2,9), rep (3,9)),  
              cluster_columns = F, row_title = NULL,
         row_split = orden_filas, cluster_row_slices = F,row_gap = unit(0.05, "cm"),
              row_names_gp = gpar(fontsize = 8) , column_title = NULL,
         
              right_annotation = leftAno
)

orden_filas <- factor (orden_filas, labels = c(1:7), levels = c(1:7))

h1 <-Heatmap (scale_WSU, name = "Z-score", 
         top_annotation = HeatmapAnnotation(
           simple_anno_size = unit(0.35, "cm"),
           annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  
           Tiempo= a$time[1:9],
           `Línea Celular` = gsub ("_DMSO.*", "",a$Sample[1:9]) %>% gsub ("_CB103.*", "",.),
           col = list(
             `Línea Celular` = c("SU-DHL-6" = "#A0522D", "U-2932"= "#D02F4B", "WSU-NHL" = "lightgreen"),
             Tiempo = c("0" = colorRampPalette(c("lightblue", "darkblue"))(3)[1], "6" = colorRampPalette(c("lightblue", "darkblue"))(3)[2],
                        "24" = colorRampPalette(c("lightblue", "darkblue"))(3)[3])
           )
         )
         ,
              column_names_rot = 60, column_names_gp = gpar(fontsize = 8),
         cluster_columns = F, cluster_rows = T,
              row_names_gp = gpar(fontsize = 8), 
         cluster_row_slices = F,
         #right_annotation = leftAno ,
         row_split = orden_filas
         )

h2 <-Heatmap (scale_U2932, name = "Z-score", 
              column_names_rot = 60, column_names_gp = gpar(fontsize = 8),
              top_annotation = HeatmapAnnotation(
                simple_anno_size = unit(0.35, "cm"), 
                annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  
                Tiempo= a$time[10:18],
                `Línea Celular` = gsub ("_DMSO.*", "",a$Sample[10:18]) %>% gsub ("_CB103.*", "",.),
                col = list(
                  `Línea Celular` = c("SU-DHL-6" = "#A0522D", "U-2932"= "#D02F4B", "WSU-NHL" = "lightgreen"),
                  Tiempo = c("0" = colorRampPalette(c("lightblue", "darkblue"))(3)[1], "6" = colorRampPalette(c("lightblue", "darkblue"))(3)[2],
                             "24" = colorRampPalette(c("lightblue", "darkblue"))(3)[3])
                )
              ),
              cluster_columns = F, cluster_rows = T,
              row_names_gp = gpar(fontsize = 8),
              cluster_row_slices = F,
              row_split = orden_filas
)

h3 <-Heatmap (scale_SUDHL, name = "Z-score", 
              top_annotation = HeatmapAnnotation(
                simple_anno_size = unit(0.35, "cm"), 
                annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  
                Tiempo= a$time[19:27],
                `Línea Celular` = gsub ("_DMSO.*", "",a$Sample[19:27]) %>% gsub ("_CB103.*", "",.),
                col = list(
                  `Línea Celular` = c("SU-DHL-6" = "#A0522D", "U-2932"= "#D02F4B", "WSU-NHL" = "lightgreen"),
                  Tiempo = c("0" = colorRampPalette(c("lightblue", "darkblue"))(3)[1], "6" = colorRampPalette(c("lightblue", "darkblue"))(3)[2],
                             "24" = colorRampPalette(c("lightblue", "darkblue"))(3)[3])
                )
              )
              ,
              column_names_rot = 60, column_names_gp = gpar(fontsize = 8),
              cluster_columns = F, cluster_rows = T,
              cluster_row_slices = F,
              row_names_gp = gpar(fontsize = 8) ,
              right_annotation = leftAno,
              row_split = orden_filas
)

draw(h1 + h2 + h3,  merge_legend = TRUE)

