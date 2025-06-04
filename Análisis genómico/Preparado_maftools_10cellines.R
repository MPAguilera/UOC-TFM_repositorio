################# - Maftools - ################################################
# Based on  Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP. 2018. Maftools: efficient and comprehensive analysis of somatic variants in cancer. Genome Research. http://dx.doi.org/10.1101/gr.239244.118

################# - Upload libraries - #########################################

if (require(rstudioapi)==F){install.packages("rstudioapi")}
if (require(readxl)==F){install.packages("readx")}
if (require(openxlsx)==F){install.packages("openxlsx")}
if (require(ggplot2)==F){install.packages("ggplot2")}
if (require(maftools)==F){BiocManager::install("maftools")}
if (require(viridis)==F){install.packages("viridis")}
if (require(dplyr)==F){install.packages("dplyr")}

########################## - Working directory - ###############################

## Set WD to file location
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

## Create a directory for the results with the current Data
Directory_for_the_day <- paste0("Mafplots_10_cellines", Sys.Date())
dir.create(Directory_for_the_day)

## Get the input files from the input folder
files = list.files(recursive = T, full.names = T)
datamutations<-read.xlsx(files[grep("10*.*filtros_duros", files)])
dataclasificaciones<-read.xlsx(files[grep("Clasificaciones_Resumen", files)])

########################## - Rearrange data format - ######################

##Rearrange the names of the sampels
## Make changes into the names of files
change_names <- function (a){
                  a <- gsub ("SUDHL","SU-DHL-", a) %>%
                  gsub ("Ocily", "OCI-Ly", .) %>%
                  gsub ("WSUNHL" , "WSU-NHL",.) %>%
                  gsub ("U9232",  "U-2932",.) %>%
                  gsub ("U2932",  "U-2932",.) %>%
                  gsub ( "HBL1", "HBL-1",.)%>%
                  gsub ( "OCI-LY", "OCI-Ly",.)%>%
                    gsub ( "RIVA", "Riva",.)%>%
                  gsub ( "\\.", "-",.)
                return (a)
                }
datamutations$Cell_line.x <- change_names (datamutations$Cell_line.x)
dataclasificaciones$Líneas.Celulares <- change_names(dataclasificaciones$Líneas.Celulares)

## Get only the relevant data from the annotation and mutation file and reorder samples
dataclasificaciones <- dataclasificaciones[match(c("RL", "SU-DHL-4", "OCI-Ly3", "OCI-Ly10", "OCI-Ly19", "HBL-1","Riva" ,"WSU-NHL", "U-2932", "SU-DHL-6"), 
                          dataclasificaciones$Líneas.Celulares),]

## Create a column with the protein annotation
datamutations$AAChange <- sapply(strsplit(datamutations$AAChange.refGene, ":"), function (x) x[5])

## Arrange data for elimination of unnecesary columns
df<-datamutations[, c("Gene.refGene", "Chr", "Start", "End", 
                      "Ref", "Alt", "Cell_line.x",
                      "ExonicFunc.refGene","ExonicFunc.refGene","Func.refGene", "AAChange")]

colnames (df) <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", 
                   "Reference_Allele", "Tumor_Seq_Allele2", 
                   "Tumor_Sample_Barcode",
                   "Variant_Classification", "Variant_Type", "Biotype", "AAChange")

## Change the chromosomes to only 1-22 or X,Y
for (i in 1:length(df$Chromosome)){
  df$Chromosome[i]<-gsub("chr", "",df$Chromosome[i])
} 


## Change the Variant Classification to the maf file
for (i in 1:length(df$Hugo_Symbol)){
  ## For exonic
  if (grepl("exonic" , df$Biotype[i])){
    if (grepl("missense variant", df$Variant_Classification[i]) | grepl("nonsynonymous", df$Variant_Classification[i]) | grepl("nonframeshift substitution", df$Variant_Classification[i]) ){
      df$Variant_Classification[i] <- "Missense_Mutation"
        if (length(strsplit(df$Tumor_Seq_Allele2[i],"")[[1]])== 1){
          df$Variant_Type[i] <- "SNP"
        } else {
          if (length(strsplit(df$Tumor_Seq_Allele2[i],"")[[1]])== 2) {
            df$Variant_Type[i] <- "TNP"
          }
          if (length(strsplit(df$Tumor_Seq_Allele2[i],"")[[1]])> 2) {
            df$Variant_Type[i] <- "ONP"
          }
        }
      
    } else {
        
      if (grepl("nonframeshift deletion", df$Variant_Classification[i]) | grepl("inframe deletion", df$Variant_Classification[i])){
        df$Variant_Classification[i] <- "In_Frame_Del"
        df$Variant_Type[i] <- "DEL"
        }
        
      if (grepl("frameshift deletion", df$Variant_Classification[i])){
        df$Variant_Classification[i] <- "Frame_Shift_Del"
        df$Variant_Type[i] <- "DEL"
      }
      if (grepl("frameshift insertion", df$Variant_Classification[i])) {
        df$Variant_Classification[i] <- "Frame_Shift_Ins"
        df$Variant_Type[i] <- "INS"
      }
      if (grepl("synonymous SNV", df$Variant_Classification[i])) {
        df$Variant_Classification[i] <- "Silent"
        df$Variant_Type[i] <- "SNP"
      }
      
      
      if (grepl("stop gained", df$Variant_Classification[i])| grepl("stopgai", df$Variant_Classification[i])) {
        df$Variant_Classification[i] <- "Nonsense_Mutation"
        df$Variant_Type[i] <- "SNP"
      }
      if (grepl("stop_loss", df$Variant_Classification[i]) | grepl("stoploss", df$Variant_Classification[i])) {
        df$Variant_Classification[i] <- "Nonstop_Mutation"
        df$Variant_Type[i] <- "SNP"
      }
      
      if (grepl("start lost", df$Variant_Classification[i]) | grepl("startloss", df$Variant_Classification[i])) {
        df$Variant_Classification[i] <- "Translation_Start_Site"
        df$Variant_Type[i] <- "SNP"
      }
      if (grepl("splice_donor_var", df$Variant_Classification[i]) | grepl("splice_acceptor_", df$Variant_Classification[i])  ) {
        df$Variant_Classification[i] <- "Splice_Site"
        df$Variant_Type[i] <- "SNP"
      }
    
      if (grepl("frameshift varia", df$Variant_Classification[i])){
        if(grepl("-", df$Tumor_Seq_Allele2[i])) {
          df$Variant_Classification[i] <- "Frame_Shift_Del"
          df$Variant_Type[i] <- "DEL"
        }
        if(grepl("-", df$Tumor_Seq_Allele2[i]) == F) {
          df$Variant_Classification[i] <- "Frame_Shift_Ins"
          df$Variant_Type[i] <- "INS"
        }
      }
      
      
      
      } 
    
  }
  
  ## For the rest
  else {
    if (grepl("UTR3" , df$Biotype[i])){
      df$Variant_Classification[i] <- "3'UTR"
    } 
    
    if (grepl("UTR5" , df$Biotype[i])){
        df$Variant_Classification[i] <- "5'UTR"
    } 
    
    if (grepl("intron" , df$Biotype[i])){
      df$Variant_Classification[i] <- "Splice_Site"
    }     
    
    if (grepl("downstream" , df$Biotype[i])){
      df$Variant_Classification[i] <- "3'Flank"
    }   
    if (grepl("upstream" , df$Biotype[i])){
      df$Variant_Classification[i] <- "5'Flank"
    } 
    
    if (grepl("splicing" , df$Biotype[i])){
      df$Variant_Classification[i] <- "Splice_Site"
    } 
    
    if (length(strsplit(df$Tumor_Seq_Allele2[i],"")[[1]])== 1){
      df$Variant_Type[i] <- "SNP"
    } else {
      if (length(strsplit(df$Tumor_Seq_Allele2[i],"")[[1]])== 2) {
        df$Variant_Type[i] <- "TNP"
      }
      if (length(strsplit(df$Tumor_Seq_Allele2[i],"")[[1]])> 2) {
        df$Variant_Type[i] <- "ONP"
      }
    }
  }
}

df <- df[, - c(10)] ## delete Biotype column

## Make sure to change the SNP correctly
df$Variant_Type[df$Tumor_Seq_Allele2 == "-"] <- "DEL"

## Rearange  the clinical data to have the same lnames as the cell lines
if (all (sort(dataclasificaciones$Líneas.Celulares) == sort(unique(df$Tumor_Sample_Barcode)))){print ("Los nombres de los samples son correctos")} else {stop("Los nombres no coinciden")}
colnames(dataclasificaciones)[1] = "Tumor_Sample_Barcode"

## Correct the format for clinical data
## BCL2
transformar <- c ("Ampl", "NO T", "ganancia", "No T ")
for (i in 1:length(dataclasificaciones$fBCL2)){
  for (j in 1:length(transformar)){
    
    if (grepl("T, ganancia ", dataclasificaciones$fBCL2[i])){ 
      dataclasificaciones$fBCL2[i]<- "T"
    }
    
    if (grepl(transformar[j], dataclasificaciones$fBCL2[i])){ 
      dataclasificaciones$fBCL2[i]<- "No T"
    }
  }

}

## BCL6
transformar <- c ("Ampl", "NO T", "ganancia", "Patrón anómalo", "No T ")
for (i in 1:length(dataclasificaciones$fBCL6)){
  for (j in 1:length(transformar)){
    
    if (grepl("T, ganancia ", dataclasificaciones$fBCL6[i])){ 
      dataclasificaciones$fBCL6[i]<- "T"
    }
    if (grepl("T y amplificada", dataclasificaciones$fBCL6[i])){ 
      dataclasificaciones$fBCL6[i]<- "T"
    }
    
    if (grepl(transformar[j], dataclasificaciones$fBCL6[i])){ 
      dataclasificaciones$fBCL6[i]<- "No T"
    }
  }
}

## BfMYC
for (i in 1:length(dataclasificaciones$fMYC)){
  for (j in 1:length(transformar)){
    
    if (grepl("T, ganancia ", dataclasificaciones$fMYC[i])){ 
      dataclasificaciones$fMYC[i]<- "T"
    }
    if (grepl("T, patrón", dataclasificaciones$fMYC[i])){ 
      dataclasificaciones$fMYC[i]<- "T"
    }
    
    if (grepl(transformar[j], dataclasificaciones$fMYC[i])){ 
      dataclasificaciones$fMYC[i]<- "No T"
    }
  }
}
dataclasificaciones = dataclasificaciones[,c(1:9),]
colnames(dataclasificaciones)[2] = "TwoStep"

########################## - Generate maf file - ###############################

## Write the maf file with the classification

maf_header <- c(
  "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode\tVariant_Classification\tVariant_Type\tAAChange")
maf_data <- paste(
  df$Hugo_Symbol, df$Chromosome, df$Start_Position, df$End_Position, 
  df$Reference_Allele, df$Tumor_Seq_Allele2, df$Tumor_Sample_Barcode,
  df$Variant_Classification,df$Variant_Type, df$AAChange,
  sep = "\t"
)

output_file <- "output_1.maf"
writeLines(c(maf_header, maf_data), con = output_file)
m<-read.maf ("output_1.maf", clinicalData = dataclasificaciones,useAll = TRUE)

########################## - Generate the plots of interest- ############################
setwd(Directory_for_the_day)
tiff(paste0("Summary_mutations_10_celllines", Sys.Date(), ".tiff"), width = 20, height =  15,  units =  "cm", res = 150)
plotmafSummary(maf =m, )
dev.off()

tiff(paste0("NOTCH1_mutations_", Sys.Date(), ".tiff"), width = 20, height =  15,  units =  "cm", res = 150)
lollipopPlot(
  maf = m, showDomainLabel = T, defaultYaxis = T, 
  pointSize = 3,bgBorderCol = "black", domainLabelSize = 0.8, 
  gene = 'NOTCH1', labPosAngle = 0,
  showMutationRate = TRUE, AACol = "AAChange"
) 
dev.off()

tiff(paste0("NOTCH2_mutations_", Sys.Date(), ".tiff"), width = 20, height =  15,  units =  "cm", res = 150)
lollipopPlot(
  maf = m, showDomainLabel = T, defaultYaxis = T, 
  pointSize = 3,bgBorderCol = "black", domainLabelSize = 0.8, 
  gene = 'NOTCH2', refSeqID = "NM_024408", labPosAngle = 0,
  showMutationRate = TRUE, AACol = "AAChange"
) 
dev.off()

