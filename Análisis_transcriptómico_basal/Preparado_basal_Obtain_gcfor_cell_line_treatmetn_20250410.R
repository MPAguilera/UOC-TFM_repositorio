###############################################################################
#################### - Obtain GCT for GSEA - ##############
###############################################################################

### Upload libraries
if (require(rstudioapi)==F){install.packages("rstudioapi")}
if (require(readxl)==F){install.packages("readxl")}
if (require(dplyr)==F){install.packages("dplyr")}
if (require(openxlsx)==F){install.packages("openxlsx")}

## Set WD to file location
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
directory = getwd()

#Get the  count file and the annot file
dat = read.xlsx(list.files(recursive = T)[grep ("Input*.*counts", list.files(recursive = T))])
rownames(dat) <- dat$ensembl_gene_id; dat <- dat[,-c(1:2)]
dat <- dat[,-match (c("external_gene_name", "description", "gene_biotype"), colnames(dat))]
resumen <-  read.xlsx(list.files(recursive = T)[grep ("Input*.*resumen", list.files(recursive = T))])

## Check taht the names are the same
colnames(dat) <- gsub ("\\.", "-", colnames(dat)) %>% gsub ("OCI-LY","OCI-Ly", .)
if (!all(sort (colnames(dat)) == sort (resumen$SampleID))){stop ("Los nombres no son los mismos")}

## chek that the names are in the corrcect order
resumen <- resumen[order(resumen$CB103), ] ## Make sure that the names are in the order by CB1103 IC50
dat <- dat[, match(resumen$SampleID, colnames(dat))]
if (!all (resumen$SampleID == colnames(dat))){stop ("Los nombres no están ordenados correctamente")}

## Change the names from ENSEMBL to HCGN 
library(biomaRt)
dat$ensembl_gene_id <- rownames(dat)
dat <-  dat[, c(length(dat), (1:(length(dat)-1)))]

mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") ## Get the anotation
attributes <- c("ensembl_gene_id" , "external_gene_name", "hgnc_symbol") ## Get the atributes to include in the data
data_annotation <- getBM(attributes = attributes,  filters = "ensembl_gene_id", 
                         values = dat$ensembl_gene_id, 
                         mart = ensembl)
da=merge (data_annotation, dat, by = "ensembl_gene_id", all.x = TRUE) 

for (i in 1:nrow(da)){
  
  if (da$hgnc_symbol[i] == ""){
    da$hgnc_symbol[i] = da$ensembl_gene_id[i]
  }

}

da = da[,-c(1)]; colnames(da)[1] = "NAME"; colnames(da)[2] = "Description"

## If the name is null then change it to the Ensembl Gene_id
da$NAME[da$NAME == ""] <- da$Description[da$NAME == ""]

## Create the directory for the GSEA input files
dir.create("Input_for_GSEA_files", showWarnings = FALSE)
setwd("Input_for_GSEA_files")

# Write the cls file 
archivo = paste0("10_cell_lines_order_by_continuous",Sys.Date(), ".cls")
write(c("#numeric", "#IC50", paste(resumen$CB103, collapse = "\t")),  file = archivo)

# Write the gct file 
archivo = paste0("10_cell_lines_", Sys.Date(), ".gct")
num_rows <- nrow(da)
num_cols <- ncol(da) - 2 
encabezado <- sprintf("#1.2\n%d\t%d", num_rows, num_cols)

writeLines(encabezado, archivo)
write.table(da, file = archivo, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE)

