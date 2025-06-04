################# - Oncoplots - ################################################

################# - Upload libraries - #########################################

if (require(rstudioapi)==F){install.packages("rstudioapi")}
if (require(readxl)==F){install.packages("readx")}
if (require(openxlsx)==F){install.packages("openxlsx")}
if (require(ggplot2)==F){install.packages("ggplot2")}
if (require(viridis)==F){install.packages("viridis")}
if (require(ComplexHeatmap)==F){install.packages("ComplexHeatmap")}
if (require(dplyr)==F){install.packages("dplyr")}
if (require(tidyr)==F){install.packages("tidyr")}

########################## - WORKING DIRECTORY - ###############################

## Set WD to file location
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

## Create a directory for the results with the currentData
Directory_for_the_day <- paste0("Complex_Map_Oncoplot_PLOTS_10celllines", Sys.Date())
dir.create(Directory_for_the_day)

## Get the input files from the input folder
files = list.files(recursive = T, full.names = T)
datamutations<-read.xlsx(files[grep("10*.*maquetada", files)])
dataclasificaciones<-read.xlsx(files[grep("Clasificaciones_Resumen", files)])

########################## - REARANGE THE  DATA FORMAT - ######################

## Make changes into the names of files
data = datamutations
data[data$Cell_Line== "SUDHL6",]$Cell_Line = "SU-DHL-6"
data[data$Cell_Line== "SUDHL4",]$Cell_Line = "SU-DHL-4"
data[data$Cell_Line== "Ocily19" ,]$Cell_Line = "OCI-Ly19"
data[data$Cell_Line== "Ocily10" ,]$Cell_Line = "OCI-Ly10"
data[data$Cell_Line== "Ocily3" ,]$Cell_Line = "OCI-Ly3"
data[data$Cell_Line== "WSUNHL" ,]$Cell_Line = "WSU-NHL"
data[data$Cell_Line== "U9232" ,]$Cell_Line = "U-2932"
data[data$Cell_Line== "HBL1" ,]$Cell_Line = "HBL-1"

datamutations = data
rm (data)

## Only get the samples of interest in all dataframes
dataclasificaciones <- dataclasificaciones[match(c("HBL-1", "OCI-LY10", "OCI-LY19", "OCI-LY3", "RIVA", "RL", "SU-DHL-4", "SU-DHL-6", "U2932", "WSU-NHL"), dataclasificaciones$Líneas.Celulares),]
datamutations <- datamutations[grep ("WSU-NHL|OCI-Ly3|OCI-Ly10|SU-DHL-6|SU-DHL-4|RL|OCI-Ly19|U-2932|RIVA|HBL-1", 
                    datamutations$Cell_Line),]

### Save the frequency of mutation within genes, the frequency of mutation within samples, and percentage of both
tabla_ancha <- as.data.frame(table(datamutations[,c(1,9)]))%>%
  pivot_wider(names_from = Gene.refGene, values_from = Freq, values_fill = list(Freq = 0))
tabla_ancha = data.frame(tabla_ancha)
rownames(tabla_ancha) = tabla_ancha[,1]; tabla_ancha = tabla_ancha[,-1]

## Convert into 0 and 1
tabla_ancha[tabla_ancha >= 1] = 1

########################## - Save file with frequencies absolute and relative - ######################
## Save the frequencies absolute and relative 
tabla_frecuencias_sample = data.frame(Frequency = apply (tabla_ancha, MARGIN = 1, sum), 
                     Proportion = apply (tabla_ancha, MARGIN = 1, function (x) sum(x)/ncol(tabla_ancha) ))
colnames(tabla_frecuencias_sample) = c("Freqency of mutated genes per sample", "Propotion of mutated genes per sample")
tabla_frecuencias_sample = rbind(tabla_frecuencias_sample, apply (tabla_frecuencias_sample, MARGIN = 2, mean))
rownames(tabla_frecuencias_sample) [nrow(tabla_frecuencias_sample)] = "Mean"

tabla_frecuencias_genes = data.frame(Frequency = apply (tabla_ancha, MARGIN = 2, sum), 
                               Proportion = apply (tabla_ancha, MARGIN = 2, function (x) sum(x)/nrow(tabla_ancha) ))
tabla_frecuencias_genes = rbind(tabla_frecuencias_genes, apply (tabla_frecuencias_genes, MARGIN = 2, mean))
rownames(tabla_frecuencias_genes) [nrow(tabla_frecuencias_genes)] = "Mean"
colnames(tabla_frecuencias_genes) = c("Freqency of mutated samples per gene", "Propotion of mutated samples per gene")

## Save the dataframes 
setwd(Directory_for_the_day)
wb <- createWorkbook()

# Add sheets to the file
addWorksheet(wb, "Samples")
writeData(wb, sheet = "Samples", data.frame(Cell_Line = rownames(tabla_frecuencias_sample),tabla_frecuencias_sample))

addWorksheet(wb, "Genes")
writeData(wb, sheet = "Genes", data.frame(Gene = rownames(tabla_frecuencias_genes),tabla_frecuencias_genes))

# Save excel file
saveWorkbook(wb, paste0("Tables_Oncoplot", Sys.Date(), ".xlsx"), overwrite = TRUE)


########################## - Change errors in anotation - ######################
## MYD88 265P
## PEST domains NOTCH1 y NOTCH2
## POU2AF1
##Spñlicing mutations

dat = datamutations[,c(1,9, 4,10,8)]
for (i in 1:nrow (dat)){
  if (dat$ExonicFunc.refGene[i] == "."){
    if (dat$Func.refGene[i] == "splicing") {dat$ExonicFunc.refGene[i]= "Splicing"}
    if (dat$Func.refGene[i] == "intronic") {dat$ExonicFunc.refGene[i]= "Splicing"}
  }
  if (dat$Gene.refGene[i] == "MYD88" & dat$Start[i] == 	38141150){ dat$ExonicFunc.refGene[i]= "MYD88 265P"}
  if (dat$Gene.refGene[i] == "NOTCH1" & dat$Start[i] <= 136497003){ dat$ExonicFunc.refGene[i]= "NOTCH PEST Dominio"}
  if (dat$Gene.refGene[i] == "NOTCH2" & dat$Start[i] <= 119916527){ dat$ExonicFunc.refGene[i]= "NOTCH PEST Dominio"}    
  if (dat$Gene.refGene[i] == "EZH2" & (dat$Start[i] == 148811635 |dat$Start[i] == 148811636 )){ dat$ExonicFunc.refGene[i]= "EZH2 646"}  
}

dat$ExonicFunc.refGene[dat$ExonicFunc.refGene == "nonsynonymous SNV" | dat$ExonicFunc.refGene == "nonframeshift substitution" ] = "Mutación de cambio de sentido" ## Missense varinat
dat$ExonicFunc.refGene[dat$ExonicFunc.refGene == "frameshift deletion" | dat$ExonicFunc.refGene == "frameshift insertion" ] = "Inserción/delección con cambio del marco de lectura" ## Frameshift Indel
dat$ExonicFunc.refGene[dat$ExonicFunc.refGene == "nonframeshift insertion" | dat$ExonicFunc.refGene == "nonframeshift deletion" ]  = "Insercción/delección sin cambio de marco de lectura" #"NonFrameshift Indel" 
dat$ExonicFunc.refGene[dat$ExonicFunc.refGene == "stoploss" ] = "Pérdida codón final"
dat$ExonicFunc.refGene[dat$ExonicFunc.refGene == "startloss" ] = "Ganancia codón final"
dat$ExonicFunc.refGene[dat$ExonicFunc.refGene == "stopgain" ] = "Ganancia codón inicio"


## Rename the MEF2B gene
dat$Gene.refGene[dat$Gene.refGene == "BORCS8-MEF2B;MEF2B"] = "MEF2B"

## Create the matrix 
dat = dat[,c(1,2,4)]
matriz = matrix (ncol = length(unique(dat[,1])), nrow = length(unique(dat[,2])), "")
colnames(matriz) = unique(dat[,1]); rownames(matriz) = unique(dat[,2])

##If there is two mutations put them together in the matrix file
for (i in 1:nrow(dat)){
  for (j in 1:ncol (matriz)){
     if (colnames(matriz)[j] == dat$Cell_Line[i]){
        for (gen1 in 1:nrow (matriz)){
          if (dat$Gene.refGene[i] == rownames(matriz)[gen1]){
            if (matriz[gen1, j] == ""){matriz[gen1, j] = dat$ExonicFunc.refGene[i]}
            else {matriz[gen1, j] = matriz[gen1, j] = "Multi Hit"}   
            #else {matriz[gen1, j] = matriz[gen1, j] = paste0(matriz[gen1, j],";",dat$ExonicFunc.refGene[i])}            
           }
        }
      }
  }
}




########################## - Create Heatmap - ######################

## Create a alter_fun
alter_fun = list(background = function(x, y, w, h) {grid.rect(x, y, w-unit(1, "pt"), h-unit(1, "pt"), 
                                                              gp = gpar(fill = "#CCCCCC", col = NA))},  
                 `EZH2 646`= function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                  gp = gpar(fill = col["EZH2 646"], col = NA)),
                 `Mutación de cambio de sentido` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                       gp = gpar(fill = col["Mutación de cambio de sentido"], col =NA)),
                 `Multi Hit` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                       gp = gpar(fill = col["Multi Hit"], col = NA)),
                 `MYD88 265P` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                              gp = gpar(fill = col["MYD88 265P"], col = NA)),
                 `Inserción/delección con cambio del marco de lectura` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                              gp = gpar(fill = col["Inserción/delección con cambio del marco de lectura"], col = NA)),
                 `Ganancia codón inicio`= function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                             gp = gpar(fill = col["Ganancia codón inicio"], col = NA)),

                 `Splicing`= function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                              gp = gpar(fill = col["Splicing"], col = NA)),
                 `NOTCH PEST Dominio`= function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                            gp = gpar(fill = col["NOTCH PEST Dominio"], col = NA)),
                 `Insercción/delección sin cambio de marco de lectura`= function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                               gp = gpar(fill = col["Insercción/delección sin cambio de marco de lectura"], col = NA)),
                 
                 `Pérdida codón final`= function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,
                                            gp = gpar(fill = col["Pérdida codón final"], col = NA)),
                 `Ganancia codón final`= function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,
                                                                       gp = gpar(fill = col["Ganancia codón final"], col = NA))
                
                 
                 
                 )

col =  c("#F8766D","#343a40","#00BADE","#8ac926","#00BD5C", 
         "#8893FF",
         "#457b9d",  "#EF67EB","#e9c46a", "#DB8E00", "#a4133c")

names(col) = c("EZH2 646",  "MYD88 265P", "Multi Hit",  "Inserción/delección con cambio del marco de lectura" ,"Mutación de cambio de sentido",
               "Ganancia codón final", 
               "Splicing", "NOTCH PEST Dominio", 
                "Insercción/delección sin cambio de marco de lectura", "Pérdida codón final", "Ganancia codón inicio")



## Rearrange the clinic data to include 
if (!all(sort(dataclasificaciones$Líneas.Celulares) == sort(colnames(matriz)))){
  dataclasificaciones$Líneas.Celulares = gsub ("DOHH2", "DoHH2", dataclasificaciones$Líneas.Celulares)  %>%
                                        gsub ("FARAGE", "Farage",.) %>%
                                        gsub ("HBL1", "HBL-1", .)%>%
                                        gsub ("KARPAS", "Karpas", .) %>%
                                        gsub ("OCI-LY", "OCI-Ly", .) %>%
                                        gsub ("PFEIFFER", "Pfeiffer", .) %>%
                                        gsub ("TOLEDO", "Toledo", .) %>%
                                        gsub ("U2932", "U-2932", .) %>%
                                        gsub ("WSU-DLBCL2", "WSU-DLCL2", .)
}
                                        
if (!all(sort(dataclasificaciones$Líneas.Celulares) == sort(colnames(matriz)))){stop("Los nombres no coinciden entre antoaciones y matriz de mtuaciones")}
  

dataclasificaciones = dataclasificaciones[match(colnames(matriz),dataclasificaciones$Líneas.Celulares),]
if (!all(colnames(matriz) == dataclasificaciones$Líneas.Celulares)){stop ("Los nombres de las matrices no están en el orden correcto")}


transformar <- c ("Ampl", "NO T", "ganancia", "No T ")
for (i in 1:length(dataclasificaciones$fBCL2)){
  for (j in 1:length(transformar)){    if (grepl("T, ganancia ", dataclasificaciones$fBCL2[i])){ 
      dataclasificaciones$fBCL2[i]<- "T"  }
    
    if (grepl(transformar[j], dataclasificaciones$fBCL2[i])){ 
      dataclasificaciones$fBCL2[i]<- "No T"
    }  }}

## BCL6
transformar <- c ("Ampl", "NO T", "ganancia", "Patrón anómalo", "No T ")
for (i in 1:length(dataclasificaciones$fBCL6)){
  for (j in 1:length(transformar)){    if (grepl("T, ganancia ", dataclasificaciones$fBCL6[i])){ 
      dataclasificaciones$fBCL6[i]<- "T"    }
    if (grepl("T y amplificada", dataclasificaciones$fBCL6[i])){ 
      dataclasificaciones$fBCL6[i]<- "T"    }
    if (grepl(transformar[j], dataclasificaciones$fBCL6[i])){ 
      dataclasificaciones$fBCL6[i]<- "No T"
    }  }}

## BfMYC
for (i in 1:length(dataclasificaciones$fMYC)){
  for (j in 1:length(transformar)){ if (grepl("T, ganancia ", dataclasificaciones$fMYC[i])){ 
      dataclasificaciones$fMYC[i]<- "T" }
    if (grepl("T, patrón", dataclasificaciones$fMYC[i])){ 
      dataclasificaciones$fMYC[i]<- "T"}
    if (grepl(transformar[j], dataclasificaciones$fMYC[i])){ 
      dataclasificaciones$fMYC[i]<- "No T"
    }} }


## Rearrange the anmes
dataclasificaciones = dataclasificaciones[,c(1:ncol(dataclasificaciones)),]

dataclasificaciones$fBCL2[dataclasificaciones$fBCL2 == "No T"]= "FISH negativo"; dataclasificaciones$fBCL2[dataclasificaciones$fBCL2 == "T"]= "FISH positivo"
dataclasificaciones$fBCL6[dataclasificaciones$fBCL6 == "No T"]= "FISH negativo"; dataclasificaciones$fBCL6[dataclasificaciones$fBCL6 == "T"]= "FISH positivo"
dataclasificaciones$fMYC[dataclasificaciones$fMYC == "No T"]= "FISH negativo"; dataclasificaciones$fMYC[dataclasificaciones$fMYC == "T"]= "FISH positivo"
dataclasificaciones$gneSeqCOO[dataclasificaciones$gneSeqCOO == "ABC"] = "ABC"
dataclasificaciones$gneSeqCOO[dataclasificaciones$gneSeqCOO ==  "Unclassified"] = "No clasificados"


## Bottom annotation
bottom_annotation  <- HeatmapAnnotation(
  simple_anno_size = unit(0.35, "cm"), 
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  
  gp = gpar(fontsize = 10), 
  border = T, 
  `Clasificador 2-S` = c(dataclasificaciones$`subtype.2-S`), 

  `FISH MYC`= c(dataclasificaciones$fMYC),
  `FISH BCL2`= c(dataclasificaciones$fBCL2),
  `FISH BCL6`= c(dataclasificaciones$fBCL6),
   gneSeqCOO= c(dataclasificaciones$gneSeqCOO),

  
  
  col = list(
    `Clasificador 2-S` = c ("EZB" = "#d17944", "EZB/MYC+"="brown","Other" = "grey", "ST2"= "#D02F4B", "BN2" = "purple", "MCD" = "lightblue", "N1" = "lightgreen"),
        `FISH MYC` = c("FISH positivo" = "#000000", "FISH negativo" = "#9c9c9c"), 
    `FISH BCL2` = c("FISH positivo" ="#000000", "FISH negativo" = "#9c9c9c"),  
    `FISH BCL6` = c("FISH positivo" = "#000000", "FISH negativo" = "#9c9c9c"),  
    gneSeqCOO = c("GC" = "orange", "ABC" = "#29d6cb", "No clasificados" = "grey")
  

    )
)


Subtype = c ("EZB" = "#d17944","EZB/MYC+"="brown","Other" = "grey","ST2"= "#D02F4B", "BN2" = "purple", "MCD" = "lightblue", "N1" = "lightgreen")

## Order of genes
orden_genes_split <- rownames(matriz)
orden_genes_split <- gsub (paste0(c("CREBBP", "BCL2", "KMT2D", "EZH2", "EP300","IRF8", "TNFRSF14"), collapse = "|"), "EZB",   orden_genes_split) %>%
                    gsub (paste0(c("TET2", "SOCS1", "SGK1",  "STAT3"), collapse = "|"), "ST2",   .) %>%
                    gsub (paste0(c("CCND3", "NOTCH2", "CD70", "DTX1"), collapse = "|"), "BN2",   .) %>%
                    gsub (paste0(c("CD79B","MYD88", "PIM1","PRDM1","BTG1" ,"PIM2"), collapse = "|"), "MCD",   .) %>%
                    gsub ("NOTCH1", "N1",   .)
                    
orden_genes_split[!grepl(c("EZB|ST2|BN2|MCD|N1"), orden_genes_split)] <- "Otros"

orden_genes_split = factor (orden_genes_split, levels = c("EZB", "ST2", "BN2", "MCD", "N1", "Otros"))      

left_annotation =  rowAnnotation( simple_anno_size = unit(0.4, "cm"),
                                  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  
                                  gp = gpar(fontsize = 3), 
  Bars = anno_block(gp = gpar(fill = Subtype[levels(orden_genes_split)]),
                    height = unit(0.4, "cm"), 
                    width = unit(0.4, "cm"), 
                   
                    labels = NULL) )


## Orden of columns
orden1 = c("SU-DHL-6" , "OCI-Ly19", "RL",  "SU-DHL-4" ,"U-2932" , "RIVA", "OCI-Ly10", "HBL-1","OCI-Ly3", "WSU-NHL") 
orden2<-c("RL", "SU-DHL-4" ,"HBL-1",  "OCI-Ly10", "OCI-Ly19", "OCI-Ly3","SU-DHL-6" , "U-2932" ,  "RIVA","WSU-NHL")
orden_columnas_split = dataclasificaciones$`subtype.2-S`
orden_columnas_split  <- factor(orden_columnas_split, levels = c("EZB", "ST2","BN2", "MCD","N1", "Other"))
orden_genes2 = c ("CREBBP", "BCL2", "KMT2D", "EZH2", "EP300","IRF8", "TNFRSF14", ## EZB 
                  "TET2", "SOCS1", "SGK1",  "STAT3",  ## ST2
                  "CCND3", "NOTCH2", "CD70", "DTX1",  # BN2
                  "CD79B","MYD88", "PIM1","PRDM1","BTG1" ,"PIM2", # MCD
                  "NOTCH1")## Print the oncoplot 

########################## - DRAW THE MAPS #################################

column_title = paste("Oncoplot 29 Cell lines","Genes and samples order by 2Step")
plot3 = ComplexHeatmap::oncoPrint(matriz, alter_fun = alter_fun, alter_fun_is_vectorized = FALSE,   
                          col = col, 
                          heatmap_legend_param =  list(title = "Mutaciones"), 
                          column_order = orden1, 
                          column_title = paste("Oncoplot 29 Cell lines","Genes Ordered by frequency"),
                          show_column_names = T, 
                          column_names_rot = 45, 
                          pct_side = "right", row_names_side = "left", 
                          row_names_gp = gpar(fontsize = 10, fontface = "italic" ), pct_gp =gpar(fontsize = 7), column_names_gp = gpar(fontsize = 10,fontface = "bold"),
                          row_title_gp = gpar(fontsize =10, fontface = "bold"),
                          row_gap = unit(1, "mm"),  
                          column_gap = unit(1, "mm"), 
                          
                          bottom_annotation =  bottom_annotation, 
                          right_annotation = rowAnnotation(row_barplot = anno_oncoprint_barplot( names(col), ylim = c(0,20), border = F, height = unit(4, "cm"),  axis_param = list(side = "top", labels_rot = 0, at = c(0,10, 20), labels = c("0", "10", "20")), )),
                          top_annotation = columnAnnotation(column_barplot = anno_oncoprint_barplot(border = FALSE, 
                                                                                    height = unit(3, "cm"), 
                                                                                    ylim = c(0,20),
                                                                                    axis_param = list(labels_rot = 0,
                                                                                                      at = c(0,5, 10,15,20), 
                                                                                         labels = c("0", "5","10", "15","20"))
                                            )))

print (plot3)

plot4 = ComplexHeatmap::oncoPrint(matriz, alter_fun = alter_fun, alter_fun_is_vectorized = FALSE, 
                                  col = col,
                                  heatmap_legend_param =  list(title = "Mutations"), 
                                  column_order = orden1, 
                                  column_title = paste("Oncoplot 29 Cell lines","Genes Ordered by frequency"),
                                  show_column_names = T, 
                                  column_names_rot = 45, 
                                  pct_side = "right", row_names_side = "left", 
                                  row_names_gp = gpar(fontsize = 10, fontface = "italic" ), 
                                  pct_gp =gpar(fontsize = 10), 
                                  column_names_gp = gpar(fontsize = 10,fontface = "bold"),
                                  row_title_gp = gpar(fontsize =10, fontface = "bold"),
                                  row_gap = unit(1, "mm"), 
                                  row_split = orden_genes_split, 
                                  column_gap = unit(1, "mm"),  
                                  
                                  bottom_annotation =  bottom_annotation, 
                                  right_annotation = rowAnnotation(row_barplot = anno_oncoprint_barplot( names(col), ylim = c(0,20), border = F, height = unit(4, "cm"),  axis_param = list(side = "top", labels_rot = 0, at = c(0,10, 20), labels = c("0", "10", "20")), )),
                                  top_annotation = columnAnnotation(column_barplot = anno_oncoprint_barplot(border = FALSE, 
                                                                                                            height = unit(3, "cm"), 
                                                                                                            ylim = c(0,20),
                                                                                                            axis_param = list(labels_rot = 0,
                                                                                                                              at = c(0,5, 10,15,20), 
                                                                                                                              labels = c("0", "5","10", "15","20"))
                                  )))

print (plot4)

capture.output(sessionInfo(), file = paste0("SessionInfoR_Oncoplot_ComplexHeatmap", Sys.Date(), ".txt"))



## Save into tiff files 
tiff(filename = paste0("Oncoplot_10_cell_lines_order_by_frequency", Sys.Date(),".tiff"), width = 30, height = 20, units = "cm", res = 100)
plot3
dev.off()
tiff(filename = paste0("Oncoplot_10_cell_lines_order_by_2S", Sys.Date(),".tiff"),  width = 30, height = 20,  units = "cm", res = 100)
plot4
dev.off()
