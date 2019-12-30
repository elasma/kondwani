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

packs= c ("ANCOM","TMB","ggplot2" ,"svDialogs", "lme4" , "zCompositions" , "lattice" , "igraph" , "reshape2" ,"vegan",
          "knitr","tidyverse","foreign","dplyr","OTUtable","data.table","ape","microbiome","plyr",
          "ggpubr","ggforce","ggplus","nnet","digest","miLineage","HMP","doParallel","DT",
          "exactRankTests","foreach","Rcpp","shiny","devtools","BiocManager","qwraps2","magrittr",
          "sjPlot","sjmisc","sjlabelled","snakecase","stringr")

check.packages(packs)

PkgsBioC=c("phyloseq","dada2", "phangorn", "edgeR", "limma","DESeq2", "metagenomeSeq","DECIPHER",
           "DirichletMultinomial","microbiomeX","coin")
check.packages.bioc(PkgsBioC)