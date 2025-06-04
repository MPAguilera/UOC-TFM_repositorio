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
if (require(dbplyr)==F){install.packages("dbplyr")}

###############################################################################
########################## - WORKING DIRECTORY - ###############################
###############################################################################

## Set WD to file location
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

## Create a directory for the results with the currentData
Directory_for_the_day <- paste0("Filtering_and_plots_GSEA", Sys.Date())
dir.create(Directory_for_the_day)


################## -  Upload data into the environment - ######################

## Look for the GSEA results files
dirs = list.dirs(recursive = F)
dirs = dirs[grep ("GSEA_results", dirs)]
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


########################## - DOTPLOTS- ########################################

## Rearrange the names in Hallmarks
library(dbplyr)
lista_index <- grep ("Hallmarks", names(lista_gsea_results))
for (i in 1:length(lista_index)){
  lista_gsea_results[[lista_index[i]]]$NAME <- 
    gsub ("HALLMARK", "HM", lista_gsea_results[[lista_index[i]]]$NAME)  %>% 
    gsub ("_", " ", . )
}

## Create the function that will create the dotplot graphs, two per each comparisson
## ALL and only significant
plot_GSEA<-function (x, nombre) { 
  
  ## Get the gene count altered by pathway
  library(stringr)
  data_percent_altered_genes <- c()
  for (i in 1:length(x$LEADING.EDGE)){
    cadena <- c(str_split(x$LEADING.EDGE[i], "\\s+"))[[1]]
    data_percent_altered_genes<- append (x = data_percent_altered_genes, (as.numeric(substr(cadena[1], start = 6, stop = nchar(cadena[1])-2)))/100)
  }
  x$`Percentage of altered genes per pathway`<-data_percent_altered_genes*100 
  x$Count <- round(data_percent_altered_genes * x$SIZE) 
  x = x[!(x$NAME == "ENNISHI_DHITSIG_POS" | x$NAME == "ENNISHI_DHITSIG_NEG" | x$NAME ==  "REDDY_COO_ABC"  | x$NAME ==  "REDDY_COO_GC"  |x$NAME ==  "LYMPH2CX_COO_GC" |x$NAME ==   "LYMPH2CX_COO_ABC" | x$NAME == "SHA_MHG"),] ## Delete the pathways of no interest
  
  data_sig<-x[(x$FDR.q.val<=0.25 & x$NOM.p.val<= 0.05),] ## Filter only data of interst
  x$NOM.p.val[(x$FDR.q.val>0.25 | x$NOM.p.val>0.05)] = NA
  
  for (i in 1:2){
    if (i == 1){da = x; range = 8; sizes = 15; type = "no significant"} else {da = data_sig; range = 12; type = "significant"; sizes = 12}
    
    plot <-ggplot(da, aes(x= NES, y= reorder(NAME, NES), size=  `Percentage of altered genes per pathway`, fill = NOM.p.val)) +
      geom_point(shape = 21, color = "black",  stroke = 0.5) +
      geom_vline(xintercept = 0, linetype = 2, color = "black") + 
      labs(x = "NES", y = "", size = "Percentage of Genes") + 
      scale_fill_gradient(low = "red", high = "blue", na.value = "grey", limits = c(0, 0.05)) + 
      scale_size_continuous(range = c(2, range), breaks = seq(10, 100, by = 10), 
                            guide = guide_legend(override.aes = list(fill = "black", color = "black")), 
                            limits = c(0, 100)) +
      theme_light()+
      theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1), 
            axis.title.x = element_text(size = 20, face = "bold"), 
            axis.text.y = element_text(size = sizes, face = "bold"),
            axis.text.x = element_text(size = sizes, face = "bold"),
            legend.title = element_text(size=18), 
            legend.text = element_text(size=16),
            strip.text = element_text(color = "black", face = "bold", size =15 ),
            strip.background = element_rect(color = "black", fill = "grey")) + xlim (min, max) + 
    
      facet_grid(. ~  paste (gsub ("_", " ", nombre)))
    
    print (plot)}
}

## Now draw the dotplot save them into a pdf file
## Save the actual direcotry
setwd(dirname(workingD))
setwd(Directory_for_the_day)

## Print the graphs
pdf (file = paste0("PLOTs_Dotplot_GSEA_29_cell_lines", Sys.Date(), ".pdf"),  width = 16, height = 12)
for (i in 1: length(lista_gsea_results)){
  x = lista_gsea_results[[i]]
  nombre = names(lista_gsea_results)[i]
  plot_GSEA(x,nombre)
}
## Make sure to close the pdf 
dev.off() 

## Save the files in tiff format
for (i in 1: length(lista_gsea_results)){
  x = lista_gsea_results[[i]]
  nombre = names(lista_gsea_results)[i]
  tiff(filename = paste0(nombre, ".tiff"), width = 37, height = 29, res = 600, units =  "cm")
  plot_GSEA(x,nombre)
  dev.off()
}
 dev.off()

## Save the significant genes into a excel
saveall = function (x, titulo){

## Names of files to save
type = gsub (".*geneset_","" , names(lista_gsea_results)) 
wb <- createWorkbook()

for (z in 1:length(type)){
  ## Create a new sheet
  addWorksheet(wb, type[z])
  look <- grep (type[z], names(x))
  writeData(wb, sheet = type[z], x[[look]] )
  
}

saveWorkbook(wb, file = paste0(titulo, ".xlsx"), overwrite = TRUE)

}

saveall (lista_gsea_results, "All_GSEA_results")
saveall (lista_significativos, "only_sig_GSEA_results")
writeLines(capture.output(sessionInfo()), paste0("Session_Info_dotplots", Sys.Date(), ".txt"))


