# gsea_hallmark_ism1_lung.R
# Gene Set Enrichment Analysis (GSEA) of Ism1-/- vs Wild Type RNA-seq data
# Dataset: Mus musculus lung tissue, postnatal day 7 (GSE262122)
# Gene sets: MSigDB Hallmark mouse collection (Mm.H)
# Pre-ranked by Wald statistic from DESeq2 results
#
# Author: Mario Antonio Rodriguez Diaz
# MSc Bioinformatics, University of Edinburgh, 2025-2026
# Functional Genomic Technologies coursework
# Adapted from class script: gsea_script.R
# (Tomlinson, University of Edinburgh, 2026)
#
# Input:  result_treatment_sorted_not_ann.csv  (output of rnaseq_deseq2_ism1_lung.Rmd)
# Output: fgsea results table + enrichment plots

library(fgsea)
library(dplyr)
#Loading DEG results
results <-read.csv("result_treatment_sorted_not_ann.csv",header =T)

head(results)

#grab the hallmark gene set from bioinf.wehi.edu.au
#generally it is better to cache download results, not download each time, but this file is small!
#You may also wish to check for updates to the annotation.
download.file("https://bioinf.wehi.edu.au/software/MSigDB/mouse_H_v5p2.rdata", destfile = "mouse_H_v5p2.rdata")
load("mouse_H_v5p2.rdata")

#We need to map to EntrezIDs- the annotation uses EntrezID
head(Mm.H)

library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
#listAttributes(mart)
ens2entrez <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"), mart=mart)
head(ens2entrez)

#Map gene names to entrezgene_id
#if NA remove (we only need to consider mappable genes)

results2 <- inner_join(results, ens2entrez, by=join_by("X"=="ensembl_gene_id"))

head(results2,n=40)

#remove any that are NA
results2 <-results2[is.na(results2$entrezgene_id)==FALSE,]
dim(results2)

head(results2)

#This analysis takes a ranked list of numbers with entrez IDs as names
#We need to choose how to order the hits
#Here I order by stat value, but I could have use bidirectional ordering using FC 
#or some combination of FC and FDR
results2 <-results2[order(results2$stat,decreasing=TRUE),]

results2 <-results2[,c("stat","entrezgene_id")]

#remove non-unique row names- we don't mind deleting these in this case
#it is wise to check the genes that are duplicated to see if they should be excluded
results2 <-results2[!duplicated(results2$entrezgene_id),]


rownames(results2)<-results2[,"entrezgene_id"]
results2["entrezgene_id"] <-NULL


#fgsea requires a named numeric vector (Entrez IDs as names)
rn <-rownames(results2)
colnames(results2) <-NULL
results2 <-results2[,1]
names(results2)<-rn

head(results2)

#Run pre-ranked GSEA against MSigDB Hallmark mouse gene sets
results2  <- results2[!(is.na(names(results2)))]
results2  <- results2[!(is.na(results2))]

fgseaRes <- fgsea(Mm.H, results2, minSize=25, maxSize = 500)

#Display downregulated pathways first
fgseaRes_down <-fgseaRes[order(fgseaRes$NES,decreasing=FALSE),]
head(fgseaRes_down, n=10)

#Display upregulated pathways first
fgseaRes_up <-fgseaRes[order(fgseaRes$NES,decreasing=TRUE),]
head(fgseaRes_up, n=10)


#Plot enrichment curves for top upregulated and downregulated pathways
library(patchwork)
library(ggplot2)

p1 <- plotEnrichment(Mm.H[["HALLMARK_HEME_METABOLISM"]],results2)+
  labs(title = "Heme Metabolism Pathway")
p2 <- plotEnrichment(Mm.H[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],results2)+
  labs(title = "IF-gamma response Pathway")
p1 + p2


