###############################################################################
############################# - Complex Heatmap -####################################
################################################################################

################# - Libraries - ################# 
library(readxl)
library(dplyr)
library(biomaRt)
library(ComplexHeatmap)
library(rstudioapi)

################# - WD- ################# 

## Set WD to file location
workingD <- rstudioapi::getActiveDocumentContext()$pat
setwd(dirname(workingD))

## Create diurectory to save results
dir.create(paste0("Complex_heatmap", Sys.Date()))
Directory_for_the_Day <-paste0("Complex_heatmap", Sys.Date())

################# - Upload Data  - ################# 

## List files in the directory
files <- list.files(recursive = T)

## Look the DESEQ sig genes
read <- files[grep ("DESEQ2*.*Sig*.*xlsx", files)]
read<- read[-grep ("Sin replicas", read)]

## Compile all sig genes
genes <- c()
for (i in 1:length(read)){
  genes <- unique(append(genes, read_xlsx(read[i])$ensembl_gene_id))
}

## Upload Mfuzz data
mfuzz <- files[grep ("Mfuzz_VST_2025-05-23*.*all_dif*.*xlsx", files)]
mfuzz <- mfuzz[-grep ("Sin replicas", mfuzz)]

clusters <- list()
for (i in 1:length(mfuzz)){
  clusters[[i]] <- read_xlsx(mfuzz[i])
  names(clusters)[i]<- gsub (".*\\/", "", mfuzz[i]) %>% gsub ("_all.*", "", .)
}

## Upload the counts to draw
counts <-files[grep ("Normalization*.*VST*.*2ndround*.*xlsx",files)]
counts <- counts[-grep ("Filter", counts)]
dat <- data.frame(read_xlsx(counts))

## Rearange the dataframe
rownames (dat) <- dat$ensembl_gene_id
del <- grep ("ensembl_gene_id|external_gene_name|description|gene_biotype", colnames(dat))
dat <- dat[,-c(1:max(del))]


################# - Group the genes by clusters - ################# 
cluster_WSU <- clusters[[grep ("WSU",names(clusters))[1]]]
number_clust <- cluster_WSU$cluster_mfuzz
orden <-cluster_WSU[order(cluster_WSU$cluster_mfuzz),]

# Generate a anotacion file with the gene order 
anotacion <- data.frame (gene = genes, order = 0)
rownames(anotacion) <- anotacion$gene
anotacion$order <-orden$cluster_mfuzz[match(anotacion$gene, orden$ensembl_gene_id)]
rm(cluster_WSU, orden)

##If its NA then the cluster its 4
anotacion$order[is.na(anotacion$order)] <- 4


# Look what genes are different and in th clsuters of SU and U
match (anotacion$gene[is.na(anotacion$order)], clusters[["SUDHL"]]$ensembl_gene_id)
match (anotacion$gene[is.na(anotacion$order)], clusters[["WSU"]]$ensembl_gene_id)
match (anotacion$gene[is.na(anotacion$order)], clusters[["U2932"]]$ensembl_gene_id)

################# - Print heatmap - ################# 

## Look for the genes to represent
dat<- dat[match(genes, rownames(dat)),]

## Sort the anotation file 
anotacion <- anotacion[match (genes, anotacion$gene),]
anotacion$order <- factor (anotacion$order, labels = c(1: max(anotacion$order)), levels = c(1:max(anotacion$order)))

## Generate the Heatmap
scaled_dat <- t(scale(t(dat)))
orden_columnas <- rep (0, ncol(dat))
orden_columnas[grep ("SUDH*.*DMSO", colnames(scaled_dat))]<- "SU-DHL-6_0H"
orden_columnas[grep ("SUDH*.*6H", colnames(scaled_dat))]<- "SU-DHL-6_6H"
orden_columnas[grep ("SUDH*.*24", colnames(scaled_dat))]<- "SU-DHL-6_24H"
orden_columnas[grep ("U2932*.*DMSO", colnames(scaled_dat))]<- "U-2932_0H"
orden_columnas[grep ("U2932*.*6", colnames(scaled_dat))]<- "U-2932_6H"
orden_columnas[grep ("U2932*.*24", colnames(scaled_dat))]<- "U-2932_24H"
orden_columnas[grep ("WSU*.*DMSO", colnames(scaled_dat))]<- "WSU-NHL_0H"
orden_columnas[grep ("WSU*.*6", colnames(scaled_dat))]<- "WSU-NHL_6H"
orden_columnas[grep ("WSU*.*24", colnames(scaled_dat))]<- "WSU-NHL_24H"
orden_columnas <- factor(orden_columnas, labels = c("WSU-NHL_0H","WSU-NHL_6H", "WSU-NHL_24H","U-2932_0H","U-2932_6H", "U-2932_24H", "SU-DHL-6_0H", "SU-DHL-6_6H", "SU-DHL-6_24H"),
                         levels = c("WSU-NHL_0H","WSU-NHL_6H", "WSU-NHL_24H","U-2932_0H","U-2932_6H", "U-2932_24H",  "SU-DHL-6_0H", "SU-DHL-6_6H", "SU-DHL-6_24H"))

orden_columnas

## Anotation_columns
anotacion_columnas <- data.frame (Condition  =orden_columnas, Time = orden_columnas, Cell_lines = orden_columnas)
rownames (anotacion_columnas) <- colnames(scaled_dat)
anotacion_columnas$Time <- factor(gsub(".*_", "", anotacion_columnas$Time),levels = c("0H", "6H", "24H"),labels = c("0H", "6H", "24H"))
anotacion_columnas$Cell_lines <- factor(gsub("_.*", "", anotacion_columnas$Cell_lines), levels = c("WSU-NHL", "U-2932", "SU-DHL-6"), labels = c("WSU-NHL", "U-2932", "SU-DHL-6"))

## Reordenar el gráfico
cell_line<- c("WSUNHL", "U2932", "SUDHL6")
time <- c("DMSO", "_6H_", "_24H_")
replica <- c("_1", "_2", "_3")

or<- c()
for (i in 1:length(cell_line)){
  for (z in 1:length(time)){
    for (s in 1:length(replica)){
      or <- append(or, grep (paste0(cell_line[i], "*.*", time[z], "*.*", replica[s]), colnames(scaled_dat)))
    }
  }
}

scaled_dat <- scaled_dat[,or]
anotacion_columnas  <- anotacion_columnas[or,]


## Infom anotation columns heatmap


top_annotation  <- HeatmapAnnotation(
  simple_anno_size = unit(0.35, "cm"), 
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  # Cambiar el tamaño de la fuente del título de la anotación
  gp = gpar(fontsize = 10), 
  border = T, 
  `Línea Celular` = anotacion_columnas$Cell_lines,
  Tiempo = anotacion_columnas$Time, 
  col = list(
    `Línea Celular` = c("SU-DHL-6" = "#A0522D", "U-2932"= "#D02F4B", "WSU-NHL" = "lightgreen"),
    Tiempo = c("0H" = colorRampPalette(c("lightblue", "darkblue"))(3)[1], "6H" = colorRampPalette(c("lightblue", "darkblue"))(3)[2],
               "24H" = colorRampPalette(c("lightblue", "darkblue"))(3)[3])
  )
)



## Heatmap no scaled
Heatmap(scaled_dat, show_row_names = F)

Heatmap(scaled_dat,  show_row_names = F,  
        cluster_row_slices = F, cluster_rows = TRUE, row_split = anotacion$order, 
        column_split = anotacion_columnas$Cell_lines, cluster_column_slices = F
        )

setwd(Directory_for_the_Day)
tiff(paste0("heatmap_onlysig_sorted", Sys.Date(),".tiff"),res = 150, width = 30, height = 20, unit = "cm")
Heatmap(scaled_dat, name = "Z-score",
        show_row_names = F,
        column_dend_height = unit(0.2, "cm"), row_dend_width = unit(0.2, "cm"),
        top_annotation = top_annotation, column_names_rot = 45, 
        column_names_gp = gpar(fontsize = 8),
        left_annotation = rowAnnotation(cluster = anno_block(
          gp = gpar(fill =  c("1" = "#C0CAFB","2" = "#FBC0E8", "3" = "#FBF1C0","4" = "#C0FBD3")),
          labels = paste0("Cluster ", 1:4),
          labels_gp = gpar(col = "black", fontsize = 10, fontface = "bold"),
          labels_rot = 90,
          labels_just = "center"),width = unit(0.6, "cm")), 
        cluster_columns = F,
        border = T,
        cluster_row_slices = F, cluster_rows = TRUE, row_split = anotacion$order, 
        column_split = factor(c(rep ("WSU",9), rep ("U2932",9), rep ("SUD",9)),levels = c("WSU", "U2932", "SUD")),
        column_title = NULL, row_title = NULL, cluster_column_slices = F, column_gap = unit(2, "mm"), row_gap = unit(2, "mm"))
dev.off()



## Print the heatmpa with the clusters of WSU and SU-DHL-6

