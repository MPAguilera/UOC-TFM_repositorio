################################################################################
##########################- Comparisson IC50 - #####################
################################################################################

## Here introduce the Drug to be used and if any samples should be excluded from the cluster analysis
## If you do not want to exclude any sample, just introduce NA in variable excluded
set.seed(1299)
study <-b ## Here introduce the study
excluded <- #"RIVA" 
  c("SU-DHL-4", "RL")

## Here te language that the graphs should be produced 
lenguage <- "Spanish" # English

ylim <- 90000 
  #75000

################# - Upload libraries - #########################################

## Here package to be used from CRAN
packages<- c("rstudioapi","mclust" ,"FSA", "car", "ggpubr","nortest","ggplot2","openxlsx","readxl","BiocManager", "ComplexHeatmap", "pheatmap", "dbplyr", "circlize","RColorBrewer", "remotes", "R.utils")
library(gridExtra)
## If they are not already installed install them
for (i in packages){ if (!requireNamespace(i, quietly = TRUE)){install.packages(i)}}

##Upload libraries
for (i in packages){suppressMessages(suppressWarnings(library(i, character.only = TRUE)))}
rm (packages, i)

################# - WD and files - #########################################

## Set WD to file location
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

## Create a directory for the results with the currentData
Directory_for_the_day <- paste0("Analysis_IC50_and_mutations_",study, "_" ,Sys.Date())
dir.create(Directory_for_the_day)

## Look for the files of interest
dir <- list.dirs()[grep ("Input", list.dirs())]
files <- list.files(recursive = T, full.names = T, path = dir)
dat <- read.xlsx(files[grep ("2Classification", files)])
rownames(dat) <- dat[,1];  dat <- dat[,-1]
IC50 <- read.xlsx(files[grep ("IC50", files)])

## Make sure the names are the same and rearrange the columns
IC50$X1<- gsub ( "SU-SHL-4","SU-DHL-4", IC50$X1) %>%
            gsub ("Riva", "RIVA", .)
variable <- IC50[, grep (paste0(study, "*.*IC50"), colnames(IC50))]
IC50 <- IC50[order(variable),] 

## Make sure that both tables have the information in the same order 
if (!all(IC50$X1 %in%  rownames(dat))){stop ("los nombres no son iguales")}
dat <- dat[match (IC50$X1, rownames(dat)), ]
if (!all(IC50$X1 ==rownames(dat))){stop ("los nombres no son tiene  el mismo orden")}

## Select the IC50 for study
variable <- IC50[, grep (paste0(study, "*.*IC50"), colnames(IC50))] 

################# - Delete non relevant columns- #########################################

##Delete mutations that all 1 or that are all 0, so there is more than 2 values per gene
dat2 <- dat
dat2$NOTCH <- dat2$NOTCH1 + dat2$NOTCH2
dat2<- dat2[, -match (c("NOTCH1", "NOTCH2"), colnames(dat2))]
keep <- sapply(dat2, function(x) length(unique(x))) >=2
dat2<-dat2[,keep]

################# - Perform a non supervised clusterization- #########################################

################# - Compare IC50 between groups #########################################

## Create blank dataframe
Subtypes <- data.frame(matrix (ncol = 6, nrow = (length(unique(IC50$`2-S`)) +3)))
colnames(Subtypes) <- c("Comparisson", "Test", "Pvalue", "Sample statistic difference", "Post-HOC", "pvalueposthoc")

## Compare all subtypes as they are (BN2 and N1 by separate)
KW <- kruskal.test(variable ~ IC50$`2-S`)
if (KW$p.value <= 0.05) {res <- "Some"} else {res <- ("None")}
if (any(KW$p.value <= 0.1)) {
  post <- paste0(dunn$res$Comparison[dunn$res$P.adj <=0.1], collapse = "|")
  p<- paste0(dunn$res$P.adj[dunn$res$P.adj <=0.1], collapse = "|")} else {post <- ""; p <- ""}
Subtypes[1,] <- c("Comparisson all subtypes", "KW test", round(KW$p.value,3),
                  res, post, p)
rm (res, KW,  post, p)

## Compare KW but with BN2 and N1 as one unique group
IC50$`2-Smod1` <- IC50$`2-S`
IC50$`2-Smod1`[IC50$`2-Smod1` == "N1" | IC50$`2-Smod1` == "BN2"] <- "N1/BN2"
IC50$`2-SmodEZB`<- IC50$`2-Smod1`
IC50$`2-SmodEZB`[IC50$`2-Smod1` == "EZB/MYC+" | IC50$`2-Smod1` == "EZB"] <- "EZBall"
for (sub in 1:2){
if (sub ==1){z <-IC50$`2-Smod1`; n <- 2; title <- "All subtypes with BN2/N1"} else {z <- IC50$`2-SmodEZB`; n <- 3; title <- "All subtypes with BN2/N1 and EZB together"}
KW <- kruskal.test(variable ~ z)
dunn <- dunnTest(variable ~ z, data = IC50, method = "bonferroni")
  if (any (dunn$res$P.adj <=0.1)){
    post <- paste0(dunn$res$Comparison[dunn$res$P.adj <=0.1], collapse = "|")
    p<- paste0(dunn$res$P.adj[dunn$res$P.adj <=0.1], collapse = "|")
  } else {post = ""; p = ""}
if (KW$p.value <= 0.1) {res <- "Some"} else {res <- ("None")}
Subtypes[n,] <- c(title, "KW test", round(KW$p.value,3), 
                  res,post, p)

rm (res, KW, post, p, title)
}


## Compare each subtype vs the REST
n <- 3
subtypes <- unique (IC50$`2-Smod1`)
for (i in 1:(length(subtypes) +1)){
  n <- n+1
  var <- IC50$`2-Smod1`
  if (n <=8){
    var[var != subtypes[i]]<- "RESTO"
  } else if (n == 9) {
    subtypes[i] <- "EZBall"
    var[var != "EZB" & var != "EZB/MYC+" ]<- "RESTO"
    var[var == "EZB" | var == "EZB/MYC+" ]<- subtypes[i]
  }

  ## Look the normality 
  if (length (var[var == "RESTO"]) <= 2 | length (var[var != "RESTO"]) <= 2){
    normalidad <- F
  }  else {
    M <- shapiro.test(variable[var == "RESTO"])
    NM<- shapiro.test(variable[var != "RESTO"])
    if (M$p.value>=0.05 & NM$p.value>=0.05){normalidad = T} else{normalidad = F}
  }
  
  ## Look homocedasticity
  homoce <- leveneTest(variable, as.factor(var)  )
  if (homoce$`Pr(>F)`[1]>=0.05){
    homoce  <- T } else {homoce <- F}

  ## Decide the test to use
  if(normalidad == T & homoce == T) {
    resultado <- t.test(x = variable[var == "RESTO"],y = variable[var != "RESTO"])
    Subtypes[n, 2] <- "T.test"
  } 
  else {
    resultado <- wilcox.test( variable[var == "RESTO"],variable[var != "RESTO"])
    Subtypes[n, 2] <- "U-Mann-Whitney"
  }
  
  Subtypes[n, c(1,3)] <-c(paste(subtypes[i], "vs REST"), round(resultado$p.value,3))
  IC50[, ncol(IC50) +1] <- factor(var, levels = c("RESTO", subtypes[i]), labels = c("RESTO", subtypes[i]))
  colnames(IC50)[ncol(IC50)] <- subtypes[i]
 rm (var, resultado, normalidad, homoce)
}


# Print the Graphs
settings <-  list(ylab("IC50"),  xlab (""), 
                  theme_bw(), 
                  theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1), 
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        axis.title.y = element_text(size = 10, face = "bold"),
                        legend.title = element_text(face = "bold"),
                        strip.text = element_text(color = "black", face = "bold", size =10),
                        strip.background = element_rect(color = "black", fill = "grey")), 
                 labs(color = "Subtipo"), ylim(0,ylim),
                 scale_color_manual(values = c("BN2" = "purple", 
                                               "N1" = "#5fe75f", 
                                               "N1" = "#5fe75f",
                                               "ST2"= "#D02F4B",
                                               "EZB" = "#d17944",  
                                               "EZB/MYC+" = "brown",
                                               "MCD" = "#93ccde",
                                               "N1/BN2" = "#572364", 
                                               "RESTO"= "grey50",
                                               "EZBall" = "brown")),
                 facet_grid(. ~  paste0(study)))




        

pdf (paste0("Graphs_IC50_comparisson_subtypes", study, Sys.Date(), ".pdf"), width = 17, height = 8)

IC50$Groups <- factor(IC50$`2-S`, labels = c("N1", "BN2", "ST2", "EZB/MYC+", "MCD", "EZB"), levels = c("N1", "BN2", "ST2", "EZB/MYC+", "MCD", "EZB"))
IC50$Groups2 <- factor(IC50$`2-Smod1`, labels = c("N1/BN2", "ST2", "EZB/MYC+", "MCD", "EZB"), levels =  c("N1/BN2", "ST2", "EZB/MYC+", "MCD", "EZB"))

a<- ggboxplot(IC50, x = "Groups", y = "variable", color = "Groups", add = "jitter") +
  stat_compare_means(method ="kruskal.test", label.y = ylim-10000, label.x = 1.5) + 
  settings  

b <- ggboxplot(IC50, x = "Groups2", y = "variable", color = "Groups2", add = "jitter") +
  settings + 
  stat_compare_means(method = "kruskal.test",  label.y = ylim-10000, label.x = 1.5)   
# + facet_grid(. ~  paste("Comparación IC50 subtipos genéticos \n agrupando BN2/N1", study))
library(cowplot)
print (plot_grid(a, b, labels = c("A","B"), nrow = 2, ncol = 2))


colnames(IC50)[colnames(IC50) == "EZB/MYC+" ] <- "EZBMYC"
subtypes[subtypes == "EZB/MYC+"] <- "EZBMYC"
colnames(IC50)[colnames(IC50) == "N1/BN2" ] <- "N1BN2"
subtypes[subtypes == "N1/BN2"] <- "N1BN2"


lista <- list()
for (i in 1:length(subtypes)){
c <- ggboxplot(IC50, x = subtypes[i], y = "variable", color =  subtypes[i], add = "jitter") +
  settings + 
  stat_compare_means(method = "wilcox.test", label.y = ylim-10000, label.x = 0.75) 
  # + facet_grid(. ~ paste("Comparación IC50 subtipos genéticos \n agrupando BN2/N1 vs Resto", study))
lista[[i]] <- c
}

do.call(grid.arrange, c(lista, nrow = 2, ncol = 3))
dev.off()

## Only fot the TFM; BN2/N1 and EZB 
library(cowplot)
if (study == "CB103") {labels <-c("A", "B") } else {labels <- c( "C", "D")}
tiff(paste0("Graphsall", study, Sys.Date(), ".tiff"), res = 150, height = 20, width = 30, units = "cm")
print (plot_grid(a, b, labels = labels, nrow = 2, ncol = 2)); dev.off()

tiff(paste0("Graphsseparate", study, Sys.Date(), ".tiff"), res = 150, height = 20, width = 30, units = "cm")
if (study == "CB103") {labels <-c("E", "F") } else {labels <- c("G", "H")}
c <-lista[[match("N1BN2", subtypes)]]
d<-lista[[match("EZB", subtypes)]]
print (plot_grid(c,d, labels = labels, nrow = 2, ncol = 2))
dev.off()

write.xlsx(Subtypes, paste0("IC50_by2Stepsubtype", study , Sys.Date(), ".xlsx"))


################# - Compare IC50 between mutated no mutated samples #########################################

## Blank dataframe
test_per_mut <- data.frame(matrix(ncol = 4, nrow = ncol(dat3)))
colnames(test_per_mut) <- c("Gene", "Normality", "Homocedaticity", "pvalue")
for (i in 1:ncol (dat3)){
  
  ##Generate the groups to compare
  Mut<- IC50[dat3[,i] == "Mutado", look]
  NoMut<- IC50[dat3[,i] == "No Mutado", look]
  test_per_mut[i, 1] <- colnames(dat3)[i]
  
  ## Perform the normality test
  if (length (Mut) <= 2 | length (NoMut) <= 2  ){
    normalidad <- F
  } else {
    M <- shapiro.test(Mut)
    NM<- shapiro.test(NoMut)
    if (M$p.value>=0.05 & NM$p.value>=0.05){normalidad = T} else{normalidad = F}
  }
  test_per_mut [i, 2] <- normalidad
  
  ## Perform the homocedasticity test
  homoce <- leveneTest(IC50[, look], as.factor(dat3[,i])  )
  if (homoce$`Pr(>F)`[1]>=0.05){
    homoce  <- T } else {homoce <- F}
  test_per_mut [i, 3] <-homoce
  
  ## Decide the test to use
  if((test_per_mut [i, 2] == T) & (test_per_mut [i, 3] == T) ){
    resultado <- t.test(x = Mut,y = NoMut)
    
  } else {
    resultado <- wilcox.test(Mut, NoMut)
  }
  
  test_per_mut [i, 4] <-resultado$p.value
  
}


pdf (paste0("Mutated vs no mutated IC50", "_",study, "_",  Sys.Date(),".pdf"))
par (mfrow  = c(3,3))
for (i in 1:ncol (dat3)){

  ##Generate the groups to compare
  Mut<- IC50[dat3[,i] == "Mutado", look]
  NoMut<- IC50[dat3[,i] == "No Mutado", look]
  boxplot (Mut, NoMut, names = c("Mutation", "No Mutation"),
                         main = paste0("Boxplot ", colnames(dat3) [i]), 
                         ylab = "IC50 value", xlab = "Mutation",
                         col = c("lightblue", "coral"))
  
  
}

dev.off()


settings <-  list(ylab("IC50"),  xlab (""), 
                  theme_bw(), 
                  theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1), 
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        axis.title.y = element_text(size = 10, face = "bold"),
                        legend.title = element_text(face = "bold"),
                        strip.text = element_text(color = "black", face = "bold", size =10),
                        strip.background = element_rect(color = "black", fill = "grey")), 
                  labs(color = "Subtipo"), ylim(0,ylim))


dat3 <- as.data.frame(dat3)
dat3[,"NOTCH"] <- as.factor (dat3[,"NOTCH"])
a <- ggboxplot(dat3, x = "NOTCH", y = "variable", color = "NOTCH", add = "jitter") + 
  settings + ylim(0, ylim) + 
  stat_compare_means(method = "wilcox.test",  label.y = ylim-10000, label.x = 0.75) +
  facet_grid(. ~  paste0(study, "\n Mutación NOTCH dominio PEST"))

dat3[,"NOTCH1noPEST"] <- as.factor (dat3[,"NOTCH1noPEST"])
b <- ggboxplot(dat3, x = "NOTCH1noPEST", y = "variable", color = "NOTCH1noPEST", add = "jitter") + 
  settings + ylim(0, ylim) +
  stat_compare_means(method = "wilcox.test",  label.y = ylim-10000, label.x = 0.75) +
  facet_grid(. ~  paste0(study, "\n Mutación NOTCH fuera dominio PEST"))

if (study == "CB103") {labels <-c("A", "B") } else {labels <- c( "C", "D")}
tiff(paste0("Graphs_M;utation_NOTCH", study, Sys.Date(), ".tiff"), res = 150, height = 20, width = 30, units = "cm")
print (plot_grid(a, b, labels = labels, nrow = 2, ncol = 2)); dev.off()

write.xlsx( test_per_mut, paste0("Testcompare mutNoMut IC50", "_", study,"_" , Sys.Date(),".xlsx") )


## Queeda por rutas ( mitar si las mutaciones por rutas tienen sentido)

################# - Compare the COO #########################################

## Blank dataframe
test_per_COO <- data.frame(matrix(ncol = 4, nrow = 2))
colnames(test_per_COO) <- c("COO", "Normality", "Homocedaticity", "pvalue")
test_per_COO[,1] <- colnames(IC50)[grep ("COO", colnames(IC50))]

n <- 1
for (i in grep ("COO", colnames(IC50))){
## Generate the groups to compare
  ABC<- IC50[IC50[,i] == "ABC", look]
  noABC<- IC50[IC50[,i] == "GC", look ]

## Perform the normality test
  if (length (ABC) <= 2 | length (noABC) <= 2 ){
    normalidad <- F
  } else {
    M <- shapiro.test(ABC)
    NM<- shapiro.test(noABC)
    if (M$p.value>=0.05 & NM$p.value>=0.05){normalidad = T} else{normalidad = F}
  }
  test_per_COO [n, 2] <- normalidad
  
## Perform the homocedasticity test
  homoce <- leveneTest(IC50[,look], as.factor(IC50[,i])  )
  if (homoce$`Pr(>F)`[1]>=0.05){
    homoce  <- T } else {homoce <- F}
  test_per_COO [n, 3] <-homoce
  
  ## Decide the test to use
  if((test_per_COO [n, 2] == T) & (test_per_mut [i, 3] == T) ){
    resultado <- t.test(x = ABC,y = noABC)
    
  } else {
    resultado <- wilcox.test(ABC, noABC)
  }
  
  test_per_COO [n, 4] <-resultado$p.value
  n <- n +1
}

## Save graphs


# Print the Graphs
settings <-  list(ylab("IC50"),  xlab (""), 
                  theme_bw(), 
                  theme(axis.title.y = element_text(size = 10, face = "bold"),
                         panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1), 
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        strip.text = element_text(color = "black", face = "bold", size =8.5),
                        strip.background = element_rect(color = "black", fill = "grey")), 
                  labs(color = "Subtipo"), ylim(0,ylim),
                  scale_color_manual(values = c("GC" = "orange", 
                                                "ABC" = "#29d6cb")))

IC50$COO_hans<- factor (IC50$COO_hans, levels = c("ABC", "GC"), labels = c("ABC", "GC"))

pdf (paste0("Graphs_IC50_COO", "_", study, "_",Sys.Date(), ".pdf"), width = 17, height = 8)
IC502<-IC50[IC50$COO_gneSeqCOO != "Unclass",]
variable2 <-variable[IC50$COO_gneSeqCOO != "Unclass"]
IC502$COO_gneSeqCOO<- factor (IC502$COO_gneSeqCOO, levels = c("ABC", "GC"), labels = c("ABC", "GC"))
a<- ggboxplot(IC502, x = colnames(IC502)[grep ("COO", colnames(IC502))[1]], 
              y = "variable2", color = colnames(IC502)[grep ("COO", colnames(IC502))[1]], add = "jitter") +
  scale_color_discrete(guide = "none") +
  stat_compare_means(method ="wilcox.test", label.y = ylim-10000, label.x = 1)+ 
  facet_grid(. ~  paste0(study, "\n COO gneSeqCOO")) +
  settings  
b<- ggboxplot(IC50, x = colnames(IC50)[grep ("COO", colnames(IC50))[2]], y = "variable", color = colnames(IC50)[grep ("COO", colnames(IC50))[2]], add = "jitter") +
  stat_compare_means(method ="wilcox.test", label.y = ylim-10000, label.x = 1) + 
  scale_color_discrete(guide = "none") +
  facet_grid(. ~  paste0(study, "\n COO Hans")) +
  settings  

library(cowplot)

print (plot_grid(a, b, labels = c("A","B"), nrow = 2, ncol = 2))
dev.off()
tiff(paste0("COO_",study,Sys.Date() ,".tiff"),, res = 150, height = 20, width = 30, units = "cm")
print (plot_grid(a, b, labels = c("A","B"), nrow = 2, ncol = 2))
dev.off()

write.xlsx( test_per_COO, paste0("Testcompare_COO_IC50", "_", study,"_" , Sys.Date(),".xlsx"))
