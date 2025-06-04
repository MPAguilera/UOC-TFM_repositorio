### Upload libraries
if (require(rstudioapi)==F){install.packages("rstudioapi")}
if (require(readxl)==F){install.packages("readx")}
if (require(openxlsx)==F){install.packages("openxlsx")}

## Set WD to file location
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

## Create a directory for the 2S if there is none already
if (any(grepl("./2S", list.dirs(recursive = F))) == F){dir.create("2S")} else{ print ("El directorio 2S ya está creado, se creará una carpeta para almacenar los datos en el mismo")}

Directory_for_the_day <- paste0("2Step/Data_for_2S_", Sys.Date())
dir.create(Directory_for_the_day)

##Check that the neccesary data is avaible in the folder
files <- list.files(path= getwd(), full.names = TRUE,recursive = TRUE) ## List all files in the folder
if (any(grepl("Traslocaciones", files))){print ("El archivo de traslocaciones ha sido encontrado correctamente")} else {stop("El archivo de traslocaciones no ha podido ser hallado")} ## Look for the traslocation data
if (any(grepl("filtros_duros", files))){ print ("El archivo de mutaciones ha sido encontrado correctamente")} else {stop("El archivo de mutaciones no ha podido ser hallado")} ## Look for the traslocation data


## Load the mutatiosn data
data<-read_xlsx(files[grep("filtros_duros", files)]) ## Select and read the data

##Inspect the mutations data

##First look at data
head(data)
str(data)
names(data)
dim(data)

## Make changes into the names of files
data[data$Cell_line.x== "SUDHL6",]$Cell_line.x = "SU-DHL-6"
data[data$Cell_line.x== "SUDHL4",]$Cell_line.x = "SU-DHL-4"
data[data$Cell_line.x== "Ocily19" ,]$Cell_line.x = "OCI-LY19"
data[data$Cell_line.x== "Ocily10" ,]$Cell_line.x = "OCI-LY10"
data[data$Cell_line.x== "Ocily3" ,]$Cell_line.x = "OCI-LY3"
data[data$Cell_line.x== "WSUNHL" ,]$Cell_line.x = "WSU-NHL"
data[data$Cell_line.x== "U9232" ,]$Cell_line.x = "U2932"

## General vision of the data
cat("Mutation found", dim(data)[1])
cat("Cell lines and samples:")
table (data$Cell_line.x) ## Cell lines used and  number of samples per line
cat ("Type of Mutec")
table (data$FILTER_Mutect)

## Generate a blank dataframe with as many rows as cellLines samples and as many columns as mutations (including the genes that are necessary for the data.frame and are not in this dataset)
mut_genes <- unique(data$Gene.refGene)
gene <- c("BCL2.fusion", "BCL6.fusion", "fMYC", "BCL6", "TNFAIP3", "CD58", "BCL10", "CD79B", "KMT2D", "PIM2", "CCND3","EP300", "PIM1",
          "EZH2","BTG1", "DTX1", "STAT3", "CD70", "UBE2A", "PRDM1", "CREBBP", "TNFRSF14", "SGK1", "NOTCH2", "TET2", "BCL2", "IRF8","SOCS1","NOTCH1","MYD88")
mut_genes <- c(mut_genes, c(gene[gene %in% mut_genes == F])) ## Combine the genes in the mutations data with the genes in the 2S


cellines<- unique(data$Cell_line.x)
Data_for_2Step<-data.frame(matrix(rep (0, ), nrow = length(cellines), ncol = length(mut_genes)))
colnames(Data_for_2Step)<- mut_genes
row.names(Data_for_2Step)<- cellines


## Change the names of columns: sampleID is the name of the Cell_line and `Gene_ensGene` the gen name

colnames(data)[colnames(data) == "Cell_line.x"] = "SampleID"
colnames(data)[colnames(data) == "Gene.refGene"] = "Gene ensGene"

##Delete the notch1 and NOTCH2 non relevant mutation (nt in the PEST domain)
paste("Mutaciones de NOTCH1 encontradas", length(grep("NOTCH1", data$`Gene ensGene`)), "de las cuales ", dim(data[(data$`Gene ensGene` == "NOTCH1" & data$Start <= 136497003),])[1] ,"son en el dominio PEST")
paste("Mutaciones de NOTCH2 encontradas", length(grep("NOTCH2", data$`Gene ensGene`)), "de las cuales ", dim(data[(data$`Gene ensGene` == "NOTCH2" & data$Start <=  117696872),])[1],"son en el dominio PEST")

data <- data[!(data$`Gene ensGene` == "NOTCH1" & data$Start >= 136497003 & data$ExonicFunc.refGene == "nonsynonymous SNV"), ]
data <-  data[!(data$`Gene ensGene` == "NOTCH2" & data$Start >= 119916694 & data$ExonicFunc.refGene == "nonsynonymous SNV"),]

## Generate the dataframe counting the mutated genes per cellLine
for (sample in 1:dim(data)[1]){  ## Looks into each row from the original data (each row represents a mutation)
  for (celline in 1:length(cellines)){ ## Look into each row from the final datafram (each row represents a cellLine)
    if (data$SampleID[sample] == cellines[celline]){ ## If the cellLine is the same in both dataframes
      mutated_gene<-(data$`Gene ensGene`[sample]) ## Save the name of the mutated gene
      for (gene in 1:length(mut_genes)){ 
        if (mutated_gene == mut_genes[gene]){
            Data_for_2Step[celline,gene] = Data_for_2Step[celline,gene] + 1
        }
      }
    }
  }
}

## Preparation of the dataframe
data_traslocation<-data.frame(read_xlsx(files[grepl("Traslocaciones", files)]))
row.names(data_traslocation) <-data_traslocation$Cell_lines
colnames(data_traslocation)<- c("celline", "BCL2.fusion","BCL6.fusion", "fMYC")
data_traslocation<-data_traslocation[-1]

##Change here the names
data_traslocation<- data_traslocation[rownames(Data_for_2Step), ]

## Check than the names from samples for each table are correct
if (all(rownames(data_traslocation) == rownames(Data_for_2Step))){print ("Los nombres y orden de los datos en la tabla Data_for_2S y traslocaciones son los mismos")} else {stop("El orden o nombre de los samples es distinto, ordene o cambie los nombres antes de continuar")}

## Function to fill the data for 2S with the traslocation information
for (i in 1:dim(data_traslocation)[1]){
  if (data_traslocation$BCL2.fusion[i] =="T" | data_traslocation$BCL2.fusion[i] =="T, ganancia señal roja"  ){data_traslocation$BCL2.fusion[i] <- 1} else {data_traslocation$BCL2.fusion[i] <- 0}
  if (data_traslocation$BCL6.fusion[i] =="T" | data_traslocation$BCL6.fusion[i] == "T y amplificada" ){data_traslocation$BCL6.fusion[i] <- 1} else {data_traslocation$BCL6.fusion[i] <- 0}
  if (data_traslocation$fMYC[i] =="T" | data_traslocation$fMYC[i] =="T, patrón variado" | data_traslocation$fMYC[i] =="T, ganancia de señal roja (1F 2-3R)" 
      | data_traslocation$fMYC[i] =="T, ganancia de señal verde"  | data_traslocation$fMYC[i] =="T, ganancia de señal roja" 
      ){data_traslocation$fMYC[i] <- 1} else {data_traslocation$fMYC[i] <- 0}
}


Data_for_2Step$BCL2.fusion <- as.numeric(data_traslocation$BCL2.fusion)
Data_for_2Step$BCL6.fusion <- as.numeric(data_traslocation$BCL6.fusion)
Data_for_2Step$fMYC <- as.numeric(data_traslocation$fMYC)

## Preparation of the dataframe to print
Data_for_2Step$SampleID <- c(row.names(Data_for_2Step))
Data_for_2Step <- Data_for_2Step[, c(dim(Data_for_2Step)[2], rep (1:(dim(Data_for_2Step)[2]-1)))]
row.names(Data_for_2Step)<- NULL

## Save the data with all mutations quantification
setwd (paste0("2Step/", "Data_for_2S_", Sys.Date()))
data_name<- paste0(Sys.Date(), "_Data_", "for_2Classification_all_quantification")
write.table(Data_for_2Step, paste0(data_name,".file"))
write.xlsx(Data_for_2Step, paste0(data_name,".xlsx"))

### Generate the dataframe with only 1 and 0.
Data_for_2Step_0_1<- Data_for_2Step
for (column in 2:length(Data_for_2Step_0_1)){
  for (row in 1:dim(Data_for_2Step_0_1)[1]){
    if (Data_for_2Step_0_1 [row, column] >1 ){
      Data_for_2Step_0_1[row, column] = 1 
    }
  }
}


## Save files
data_name<- paste0(Sys.Date(), "_Data_", "for_2Classification_only_0_&_1")
write.xlsx(Data_for_2Step_0_1, paste0(data_name, ".xlsx"), colNames = TRUE, rowNames = FALSE)
write.table(Data_for_2Step_0_1, paste0(data_name,".file"))
