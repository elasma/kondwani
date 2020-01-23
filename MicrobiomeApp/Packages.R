check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

check.packages.bioc <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    BiocManager::install(new.pkg)
  sapply(pkg, require, character.only = TRUE)
}
check.packages(c("devtools","BiocManager"))
PkgsBioC=c("dada2","DESeq2","phyloseq", "phangorn", "edgeR", "limma", "metagenomeSeq","DECIPHER",
           "DirichletMultinomial","coin")
check.packages.bioc(PkgsBioC)
packs= c ("TMB","ggplot2" ,"svDialogs", "lme4" , "zCompositions" , "lattice" , "igraph" , "reshape2" ,"vegan",
          "knitr","tidyverse","foreign","dplyr","OTUtable","data.table","ape","microbiome","plyr",
          "ggpubr","ggforce","ggplus","nnet","digest","miLineage","HMP","doParallel","DT",
          "exactRankTests","foreach","Rcpp","shiny","qwraps2","magrittr",
          "snakecase","stringr")

check.packages(packs)






