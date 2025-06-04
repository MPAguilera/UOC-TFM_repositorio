################################################################################
##########################- Combine and filter Annovar results -################
################################################################################


################# - Upload libraries - #########################################
if (require(rstudioapi)==F){install.packages("rstudioapi")}
if (require(readxl)==F){install.packages("readx")}
if (require(openxlsx)==F){install.packages("openxlsx")}
if (require(vcfR)==F){install.packages("vcfR")}

########################## - Working Directory - ###############################

## Set WD to file location
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

## Create a directory for the results with the currentData
Directory_for_the_day <- paste0("Cell_lines_Combined_dataframes_NfcoreSarek", Sys.Date())
dir.create(Directory_for_the_day)

########################## - Upload files into environment- ###############################


## get the directories where the Annovar files are located
## Get the directories
dirs = list.dirs()
dirs<-dirs[grep ("./Nfcore_results_mutect2/Annotation_Annovar_*.*/Mutec2_filtered_prepared",  dirs)] ## introducir el nombre del directorio con archivos annovar

## Get the files
files = list.files(path = dirs, full.names = T)
files = files[grep ("vcf*.*Multianno.hg38_multianno.txt", files)] ## Introducir el final de los archivos annotar

## Create a list which contains the annovar results into R 
lista_multianno = list()
for (file in 1:length(files)){
  lista_multianno[[file]] = read.delim(files[file])

  ## Rename columns
  colnames(lista_multianno[[file]])[89:98] = c("chr", "pos", "ID", "ref.1", "alt.1", "QUAL", "FILTER", "INFO", "FORMAT", "FORMAT_DATA")
  
  ## Get the relevant info ans split into different columns
  GT = c(); AD = c(); AF = c(); DP = c(); F1R2 = c(); F2R2 = c(); FAD = c(); PGT = c(); PID = c()  ; PS = c(); SB = c()
  
  for (i in 1:dim(lista_multianno[[file]])[1]){
    A = strsplit(lista_multianno[[file]][i,98], split = ":")[[1]]
    if (length(A) == 8){ GT = append (GT, A[1]); AD = append (AD, A[2]); AF = append (AF, A[3]); DP = append (DP, A[4]); F1R2 = append (F1R2, A[5]); F2R2 = append (F2R2, A[6]);
    FAD = append (FAD, A[7]);  PGT = append (PGT, NA);  PID = append (PID, NA);  PS= append (PS, NA);  SB = append (SB, A[8]); 
    
    } else{
      GT = append (GT, A[1]); AD = append (AD, A[2]); AF = append (AF, A[3]); DP = append (DP, A[4]); F1R2 = append (F1R2, A[5]); F2R2 = append (F2R2, A[6]);
      FAD = append (FAD, A[7]);  PGT = append (PGT, A[8]);  PID = append (PID, A[9]);  PS= append (PS, A[10]);  SB = append (SB, A[11]); 
      
    }
  }
  lista_multianno[[file]][,98:108] = data.frame(GT,AD, AF, DP, F1R2, F2R2, FAD, PGT, PID, PS, SB)
  
  nombre = gsub (".mutect2.filtered.vcf*.*Multianno.hg38_multianno.txt", "", files[file])
  nombre = gsub ("./Nfcore_results_mutect2/Annotation_Annovar_*.*/Mutec2_filtered_prepared_*.*/", "", nombre)
  names(lista_multianno)[file] = nombre
}

## Bind the Annovar info into only one dataframe
tabla_annovar = data.frame()
for (file in 1:length(lista_multianno)){
  lista_multianno[[file]]$Cell_line = names(lista_multianno)[file]
  tabla_annovar = rbind (tabla_annovar, lista_multianno[[file]])
}

########################## - Rearrange data- ###############################

##Rearrange the allele frequency
count = 0
alleleanterior = 0
tabladef = tabla_annovar
colnames(tabladef)[colnames(tabladef)== "AF"] =   "Allele_Frequency"
colnames(tabladef)[colnames(tabladef)== "DP"] =   "Deep_coverage"
colnames(tabladef)[colnames(tabladef)== "ref"] =   "ref_mutect2"
colnames(tabladef)[colnames(tabladef)== "alt"] =   "alt_mutect2"
colnames(tabladef)[colnames(tabladef)== "pos"] =   "pos_mutect2"
colnames(tabladef)[colnames(tabladef)== "FILTER"] =   "FILTER_Mutect"

## Convert the allele frequency into numeric (sometimes there is more than one data)
for (i in 1:dim(tabladef)[1]){
  allele = tabladef$Allele_Frequency[i]
  if (allele == alleleanterior){
    alleles = strsplit(allele,",")[[1]]
    count = count +1
    if (length(alleles) >=2){
      tabladef$Allele_Frequency[i] = alleles[count]
    }
  } else {
    alleles = strsplit(allele,",")[[1]]
    if (length(alleles) >=2){
      count = 1
      tabladef$Allele_Frequency[i] = alleles[count]
    } else { count = 0}
  }
  alleleanterior =allele
}


## Convert  columns into numeric variables
tabladef$Deep_coverage = as.numeric(tabladef$Deep_coverage)
tabladef$Allele_Frequency = as.numeric(tabladef$Allele_Frequency)
tabladef$ExAC_ALL = as.numeric(tabladef$ExAC_ALL)
tabladef$X1000g2015aug_all = as.numeric(tabladef$X1000g2015aug_all)
tabladef$X1000g2015aug_all = as.numeric(tabladef$X1000g2015aug_eur)
tabladef$X1000g2015aug_all = as.numeric(tabladef$X1000g2015aug_eas)
tabladef$X1000g2015aug_all = as.numeric(tabladef$X1000g2015aug_afr)

## Generate one column per sample
dimension <- dim(tabladef)[2]
ncol1 <- dimension +1
ncol2 <- dimension + 10
tabladef[,c(ncol1:ncol2)] <- 0
colnames(tabladef)[ncol1:ncol2] <- names(lista_multianno)

## Look if the mutation found is ppresen in any other line, and save in into the correct column
for (line in 1:dim(tabladef)[1]){
  num = grep(paste0(tabladef$Start[line],tabladef$End[line], tabladef$Ref[line],  tabladef$Alt[line]), paste0(tabladef$Start, tabladef$End, tabladef$Ref,  tabladef$Alt))
  if (length(num)>1){
    for (i in 1:length(num)){
      linea <- num[i]
      columna<- dimension + (grep(paste0(tabladef$Cell_line[num[i]],"_AF"),  paste0(names(lista_multianno),"_AF")))
      tabladef[line, columna] <- tabladef$Allele_Frequency[linea]
    }
  } else {
    for (i in 1:length(num)){
      linea <- num[i]
      columna<- dimension + (grep(tabladef$Cell_line[num[i]], names(lista_multianno)))
      tabladef[line, columna] <- tabladef$Allele_Frequency[linea]
    }
  }
}

## Filter the data
tabladeffiltrada = tabladef[tabladef$Deep_coverage>=50,]
tabladeffiltrada = tabladeffiltrada[tabladeffiltrada$Allele_Frequency>=0.1,]
tabladeffiltrada = tabladeffiltrada[grepl("weak_evidence", tabladeffiltrada$FILTER_Mutect) ==F, ]
tabladeffiltrada = tabladeffiltrada[grepl("orientation", tabladeffiltrada$FILTER_Mutect) == F,]
tabladeffiltrada = tabladeffiltrada[grepl("base_qual", tabladeffiltrada$FILTER_Mutect) == F,]
tabladeffiltrada = tabladeffiltrada[grepl("contamination", tabladeffiltrada$FILTER_Mutect) == F,]
tabladeffiltrada = tabladeffiltrada[grepl("strand_bias", tabladeffiltrada$FILTER_Mutect) == F,]
tabladeffiltrada = tabladeffiltrada[grepl("map_qual", tabladeffiltrada$FILTER_Mutect) == F,]
tabladeffiltrada = tabladeffiltrada[grepl("slippage", tabladeffiltrada$FILTER_Mutect) == F,]
tabladeffiltrada = tabladeffiltrada[grepl("fragment", tabladeffiltrada$FILTER_Mutect) == F,]
tabladeffiltrada = tabladeffiltrada[grepl("downstream", tabladeffiltrada$Func.refGene) == F,]
tabladeffiltrada = tabladeffiltrada[grepl("upstream", tabladeffiltrada$Func.refGene) == F,]
tabladeffiltrada = tabladeffiltrada[grepl("intergenic", tabladeffiltrada$Func.refGene) == F,]
tabladeffiltrada = tabladeffiltrada[grepl("ncRNA_intronic", tabladeffiltrada$Func.refGene) == F,]
tabladeffiltrada = tabladeffiltrada[grepl("panel_of_normals", tabladeffiltrada$FILTER_Mutect) == F,]
tabladeffiltrada = tabladeffiltrada[(tabladeffiltrada$ExAC_ALL <= 0.01 | is.na(tabladeffiltrada$ExAC_ALL)), ]
tabladeffiltrada = tabladeffiltrada[(tabladeffiltrada$X1000g2015aug_all <= 0.01 | is.n0a(tabladeffiltrada$X1000g2015aug_all)), ]

## Filter the AD
delete <- c()
for (i in 1:ncol (tabladeffiltrada)){
  a <-  as.numeric(unlist(strsplit(tabladeffiltrada$AD[i], ",")))
  if (all(a[2:length(a)]<5)){delete <- append (delete, i)}
}
if (is.null(delete) == F){tabladeffiltrada <- tabladeffiltrada[-delete, ]}

########################## - Save data- ###############################
write.xlsx(tabladeffiltrada, paste0("tabla_mutect2filter_ANNOVAR&VEP_filtros_laxos", ".xlsx"))










