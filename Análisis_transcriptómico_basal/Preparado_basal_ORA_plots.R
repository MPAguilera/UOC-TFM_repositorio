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

## Read files
files <- list.files (recursive = T, full.names = T)
dat<- read.xlsx(files[grep ("ORA_for_draw", files)])

## Create the daty WD
name_dir<- paste0("ORA_Draw_ggplots", Sys.Date())
dir.create(name_dir)

################# - Ggplots - #################################
library(circlize)
dat$Count[dat$Correlation == "Neg"]<- dat$Count[dat$Correlation == "Neg"]*-1
dat$Orden_num <-as.numeric(dat$Orden_num)

library(dplyr)
dat$rename <-paste0(gsub ("Hallmarks","HM " , dat$Database)  %>% 
         gsub ("Reactome","R ", .)  %>% gsub ("Gene Ontology","GO ", .) %>% gsub ("Wikipathways","WIKI ", .),
       gsub ("_"," " ,tolower(dat$Description))  %>% gsub ("reactome ","", .))

tiff(paste0("ORA_0.65_all.tiff"),res = 150, width = 20, height =  20, unit = "cm")
settings <- list(labs(x= "nº de genes", y = "", fill = "p valor \n ajustado"),
                 theme_bw(),
                 theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1), 
                       axis.title.x = element_text(size = 12, face = "bold"),
                       axis.text.x = element_text(size = 12, face = "bold"),
                       legend.title = element_text(face = "bold"),
                       axis.text.y =  element_text(size = 6),
                       strip.text = element_text(color = "black", face = "bold", size =14 ),
                       strip.background = element_rect(color = "black", fill = "grey")))


ggplot(dat, aes(x = Count, y = reorder(rename, -Orden_num), fill = p.adjust)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_gradient(low = "red",high = "blue") +
  settings + 
  facet_grid(. ~  paste ("ORA \ngenes correlacionados"))
dev.off()

