##########################- Ora Plots- ##########################


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

## Create the daty WD
name_dir<- paste0("ORA_Draw_ggplots", Sys.Date())
dir.create(name_dir)

## Read files
files <- list.files (recursive = T, full.names = T)
files <-files[grep ("ORA*.*Cluster*.*xlsx", files)]
files <- files[-grep ("Sin replicas", files)]

cell_lines <- c("WSU", "SUDHL", "U2932")
d = 2
for (d in 1:length(cell_lines)){
files2 <-grep (cell_lines[d], files)

for (i in 1:length(files2)){
name <-files[files2[i]]
if (exists("study") == F & any(grepl("Representar", excel_sheets(name))) ){
  study <-data.frame(read.xlsx(name, sheet = "Representar"), 
                               cluster = gsub(".*Cluster_", "",name)%>% gsub (paste0("_",cell_lines[d], ".*"),"",.))} 
if (exists("study") == F & any(grepl("Representar", excel_sheets(name))) == F){print ("No incluir cluster")}
if (exists("study") == T & any(grepl("Representar", excel_sheets(name))) == F){print ("No incluir cluster")}
if (exists("study") == T & any(grepl("Representar", excel_sheets(name))) ){
    study <-rbind (study,data.frame(read.xlsx(name, sheet = "Representar"), 
                   cluster = gsub(".*Cluster_", "",name)%>% gsub (paste0("_",cell_lines[d], ".*"),"",.)))

}


}



################# - Ggplots - #################################
library(circlize)
library(dplyr)
setwd(name_dir)
study$Description <-gsub ("HALLMARK","HM " , study$Description) %>% gsub ("_", " ",.)
study$Description<-tolower(study$Description)
cell_lines[d]<-gsub ("WSU", "WSU-NHL",cell_lines[d]) %>% gsub ("SUDHL", "SU-DHL-6",.)
study$Description <-paste0(cell_lines[d],"-",study$cluster, " ", study$Description)
study <- study[!duplicated(study$Description),]
tiff(paste0(gsub (".*2025-05-24/", "", name) %>% gsub (".xlsx", "",.),".tiff"),res = 150, width = 23, height =13, unit = "cm")
settings <- list(labs(x= "nº de genes", y = "", fill = "p valor \n ajustado"),
                 theme_bw(),
                 theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1), 
                       axis.title.x = element_text(size = 12, face = "bold"),
                       axis.text.x = element_text(size = 12, face = "bold"),
                       legend.title = element_text(face = "bold"),
                       axis.text.y =  element_text(size = 8),
                       strip.text = element_text(color = "black", face = "bold", size =10 ),
                       strip.background = element_rect(color = "black", fill = "grey")))


print(ggplot(study, aes(x = Count, y = reorder(Description, -Orden), fill = p.adjust)) +
    geom_bar(stat = "identity", color = "black", width = 0.5) + 
  scale_fill_gradient(low = "red",high = "blue") +
  settings + 
  scale_fill_gradient(
  low = "red", 
    high = "blue", 
      limits = c(0, 0.05)) +
  
  facet_grid(. ~  paste ("ORA \n", cell_lines[d])))
dev.off()
setwd(dirname(workingD))
rm (study)
}


