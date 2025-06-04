###### Save the definitive dataframe o filter variables ################


################# - Upload libraries - #########################################
if (require(rstudioapi)==F){install.packages("rstudioapi")}
if (require(readxl)==F){install.packages("readx")}
if (require(openxlsx)==F){install.packages("openxlsx")}
if (require(vcfR)==F){install.packages("vcfR")}

########################## - WORKING DIRECTORY - ###############################

## Set WD to file location
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

## Create a directory for the results with the currentData
Directory_for_the_day <- paste0("Tabla_ya_filtrada_maquetada", Sys.Date())
dir.create(Directory_for_the_day)

##GET THE ANNOVAR ANOTATION
## Get the directories
files = list.files(recursive = T) 
files = files[grep ("filtros_duros", files)] ## Introducir el final de los archivos annotar
data = read.xlsx(files)

#### Change the colnames
colnames(data)
colnames(data)[c(2,15, 18, 19,20, 21:24, 27:31)] = c("Cell_Line", "Allele_Depth", "F1R2(Forward 1st,reverse 2nd)", "F2R1(Reverse 1st, foward 2nd)", "Filtered_Allelic_Depths",
                                            "Phased Genotype", 
                                                  "Phasing ID", 
                                                  "Phase Set", 
                                                  "Strand Bias", 
                                                "ClinVar Allele ID", 
                                                "ClinVar Disease Name", 
                                                "ClinVar Disease Database", 
                                                "ClinVar Review Status", 
                                                "ClinVar Clinical Significance")
dat = data

## Save the principal refAllele into only one column
for (i in 1:nrow(data)){
  info = strsplit(data$GeneDetail.refGene[i], ";")[[1]]
  data$GeneDetail.refGene[i] = info[1]
  
  infoprot = strsplit(data$AAChange.refGene[i], ",")[[1]]
  data$AAChange.refGene[i] = infoprot[1]
  
  if (length(info) != 1){
    data$GeneDetail.otherRef[i] = paste(info[2:length(info)], collapse = ";") } else { data$GeneDetail.otherRef[i] = "." }
  
  if (length(infoprot) != 1){
    data$AAChange.otherRef[i] = paste(infoprot[2:length(infoprot)], collapse = ";") } else { data$AAChange.otherRef[i] = "." }
  }
  
## Split the AAChange and GeneDetail into separate
for (i in 1:nrow(data)){
  if (data$AAChange.refGene[i] == "." & data$GeneDetail.refGene[i] != "."){
    data$AAChange.refGene[i] = paste0(data$Gene.refGene[i], ":", data$GeneDetail.refGene[i], ":", ".")
  }
  if(data$AAChange.refGene[i] == "." & data$GeneDetail.refGene[i] == "."){
    data$AAChange.refGene[i] = ".:.:.:.:."
  }
}

for (i in 1:nrow(data)){
  if (data$AAChange.refGene[i] != "."){
    info3 = strsplit(data$AAChange.refGene[i], ":")[[1]]
    data$refGene[i] =info3[2]
    data$exon[i] = gsub (x = info3[3], pattern = "exon", "")
    data$c.Change[i] = info3[4]
    data$p.Change[i] = info3[5]
  }
}

## Rearrange the columns so its easier to read files
dat = data[, c(2:10, 12, 153:156, 14:98, 109:139, 147:149, 151,152)]


## REARRANGE THE COSMIC70 INFO into two different columns
for (i in 1:nrow(dat)){
  info = strsplit(dat$cosmic70mod[i], ";")[[1]]
  if (length(info) != 1){
    ocurrence = gsub("OCCURENCE=", "",info[2])
    numers = (strsplit(gsub ("\\(.*?\\)", "", ocurrence), ","))[[1]]
    dat$cosmic70mod_ID[i] = gsub ("ID=","", info[1])
    dat$COSMIC_n_samples_mutated[i] = sum(as.numeric(numers))
    dat$COSMIC_samples_mutated[i] = ocurrence
    
  }
   else {
     dat$cosmic70mod_ID[i] = "."
     dat$COSMIC_n_samples_mutated[i] = "."
     dat$COSMIC_samples_mutated[i] = "."
   }
}

dat = dat [, c(1:27, 136:138, 28:135)]

## Get the zigosity
dat$Genotype  = unlist(lapply (dat$Format_INFO, FUN = function (row) strsplit(row, ":")[[1]][1]))
dat$Heterozygosity = dat$Genotype
dat$Heterozygosity[dat$Heterozygosity == "0/1" | dat$Heterozygosity == "1|0" | dat$Heterozygosity == "0|1" | dat$Heterozygosity == "0/1/2" ] = "YES"

## Rearrange the ocily1, 10, and 19
dat$Ocily10[dat$Ocily1 == dat$Ocily10 & dat$Ocily1 != 0] = 0
dat$Ocily19[dat$Ocily1 == dat$Ocily19 & dat$Ocily1 != 0] = 0

### Rearrange the MYD88 annotation
dat$ExonicFunc.refGene[dat$Start == 38141150] = "nonsynonymous SNV"
dat$c.Change[dat$Start == 38141150] = "c.818T>C"
dat$p.Change[dat$Start == 38141150] = "p.L273P"
dat$exon[dat$Start == 38141150] = 4

## Save the data
dat_def = dat[, c(1:25,139,140,28:138)]
setwd("./Cell_lines_Combined_dataframes_NfcoreSarek2024-09-29")
write.xlsx (dat_def, paste0("29_cell_lines_mutations_tabla_maquetada", Sys.Date(), ".xlsx"))
write.csv (dat_def, paste0("29_cell_lines_mutations_tabla_maquetada", Sys.Date(), ".csv"))
write.table(dat_def, paste0("29_cell_lines_mutations_tabla_maquetada", Sys.Date(), ".tsv"))









