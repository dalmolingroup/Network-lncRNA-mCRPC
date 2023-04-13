# #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
#
#                     
#            Rafaella Sousa Ferraz <rafaellaferraz.16@hotmail.com>
#
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

#Packages used in this analysis:
library(tidyverse)
library(DESeq2)
library(readxl)

################################################################################
####################  DIFFERENTIAL EXPRESSION ANALYSIS (DEA) ###################
################################################################################

#CRPC primario dataset

load("clinical_CRPC.RData") #clinical_CRPC
load("counts_all_elements_CPRCprimario.RData") #matrix_crpc

all(rownames(clinical_CRPC) == colnames(matrix_CRPC))

#Select groups
clinical_CRPC <- clinical_CRPC[clinical_CRPC$sample_type != "Primary Tumor",]
matrix_CRPC <- matrix_CRPC %>%
  select(clinical_CRPC$sample)

sapply(clinical_CRPC, class)
clinical_CRPC$sample_type <- str_replace_all(clinical_CRPC$sample_type, " ", "_")
clinical_CRPC$sample_type <- as.factor(clinical_CRPC$sample_type)
clinical_CRPC %>%
  dplyr::count(sample_type)



dds_crpc <- DESeqDataSetFromMatrix(countData = as.matrix(matrix_CRPC),
                                   colData = clinical_CRPC,
                                   design = ~ sample_type)
dds_crpc <- DESeq(dds_crpc)

resultsNames(dds_crpc)

plotPCA(vst(dds_crpc), intgroup="sample_type")

stand <- vst(dds_crpc)

# Computar PCA
tt <- pca(assay(stand), metadata = as.matrix(colData(dds_crpc)))

# Plotar PC1 e PC2
PCAtools::biplot(tt, colby = "sample_type", shape = "sample_type", legendPosition = "left")

# Proporção de variabilidade de cada componente
PCAtools::screeplot(tt)

PCAtools::eigencorplot(tt, components = getComponents(tt, 1:10), metavars = c("sizeFactor", "sample_type"))

#Metastatic x normal
res_crpc_norm_met <- results(dds_crpc, contrast = c("sample_type", "Metastatic", "Solid_Tissue_Normal"))
#Metastatic x primary
res_crpc_prim_met <- results(dds_crpc, contrast = c("sample_type", "Metastatic", "CRPC"))
#Primary x normal
res_crpc_prim_norm <- results(dds_crpc, contrast = c("sample_type", "CRPC", "Solid_Tissue_Normal"))

################################################################################
################################  DEA TO TNA ###################################
################################################################################

#METASTATIC X NORMAL
res_crpc_norm_met <- res_crpc_norm_met %>% as.data.frame() %>% drop_na()

res_crpc_norm_met <- merge(res_crpc_norm_met, gff, by.x = "row.names", by.y = "gene_id", all.x = T) 
res_crpc_norm_met$ensg <- sapply(strsplit(as.character(res_crpc_norm_met$Row.names),'\\.'), "[", 1)

for (i in 1:nrow(res_crpc_norm_met)) {
  if (is.na(res_crpc_norm_met[i,"gene_name"])) {
    res_crpc_norm_met[i, "gene_name"] <- res_crpc_norm_met[i,"ensg"]
  }
}

#duplicated genes
res_crpc_norm_met[duplicated(res_crpc_norm_met$gene_name),c("ensg","gene_name")]

dup11 <- unique(res_crpc_norm_met[duplicated(res_crpc_norm_met$gene_name),c("gene_name")])
local_dup11 <- row.names(res_crpc_norm_met[res_crpc_norm_met$gene_name %in% dup11,])
for (x in local_dup11){
  res_crpc_norm_met[x, "gene_name"] <- res_crpc_norm_met[x, "ensg"]
}


#foldchange 
fold_crpc_norm_met <- res_crpc_norm_met$log2FoldChange
names(fold_crpc_norm_met) <- res_crpc_norm_met$gene_name

#name of differentially expressed genes
DEG_crpc_norm_met <- res_crpc_norm_met %>%
  filter(padj < 0.01, log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>%
  arrange(desc(log2FoldChange))

symbol_crpc_norm_met <- DEG_crpc_norm_met$gene_name


#METASTATIC X PRIMARY
res_crpc_prim_met <- res_crpc_prim_met %>% as.data.frame() %>% drop_na()

res_crpc_prim_met <- merge(res_crpc_prim_met, gff, by.x = "row.names", by.y = "gene_id", all.x = T) 
res_crpc_prim_met$ensg <- sapply(strsplit(as.character(res_crpc_prim_met$Row.names),'\\.'), "[", 1)

for (i in 1:nrow(res_crpc_prim_met)) {
  if (is.na(res_crpc_prim_met[i,"gene_name"])) {
    res_crpc_prim_met[i, "gene_name"] <- res_crpc_prim_met[i,"ensg"]
  }
}

#duplicated genes
res_crpc_prim_met[duplicated(res_crpc_prim_met$gene_name),c("ensg","gene_name")]

dup12 <- unique(res_crpc_prim_met[duplicated(res_crpc_prim_met$gene_name),c("gene_name")])
local_dup12 <- row.names(res_crpc_prim_met[res_crpc_prim_met$gene_name %in% dup12,])
for (x in local_dup12){
  res_crpc_prim_met[x, "gene_name"] <- res_crpc_prim_met[x, "ensg"]
}

#foldchange
fold_crpc_prim_met <- res_crpc_prim_met$log2FoldChange
names(fold_crpc_prim_met) <- res_crpc_prim_met$gene_name

#name of differentially expressed genes
DEG_crpc_prim_met <- res_crpc_prim_met %>%
  filter(padj < 0.01, log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>%
  arrange(desc(log2FoldChange))

symbol_crpc_prim_met <- DEG_crpc_prim_met$gene_name

#PRIMARY X NORMAL
res_crpc_prim_norm <- res_crpc_prim_norm %>% as.data.frame() %>% drop_na()

res_crpc_prim_norm <- merge(res_crpc_prim_norm, gff, by.x = "row.names", by.y = "gene_id", all.x = T) 
res_crpc_prim_norm$ensg <- sapply(strsplit(as.character(res_crpc_prim_norm$Row.names),'\\.'), "[", 1)

for (i in 1:nrow(res_crpc_prim_norm)) {
  if (is.na(res_crpc_prim_norm[i,"gene_name"])) {
    res_crpc_prim_norm[i, "gene_name"] <- res_crpc_prim_norm[i,"ensg"]
  }
}

#duplicated genes
res_crpc_prim_norm[duplicated(res_crpc_prim_norm$gene_name),c("ensg","gene_name")]

dup13 <- unique(res_crpc_prim_norm[duplicated(res_crpc_prim_norm$gene_name),c("gene_name")])
local_dup13 <- row.names(res_crpc_prim_norm[res_crpc_prim_norm$gene_name %in% dup13,])
for (x in local_dup13){
  res_crpc_prim_norm[x, "gene_name"] <- res_crpc_prim_norm[x, "ensg"]
}

#foldchange
fold_crpc_prim_norm <- res_crpc_prim_norm$log2FoldChange
names(fold_crpc_prim_norm) <- res_crpc_prim_norm$gene_name

#name of differentially expressed genes
DEG_crpc_prim_norm <- res_crpc_prim_norm %>%
  filter(padj < 0.01, log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>%
  arrange(desc(log2FoldChange))

symbol_crpc_prim_norm <- DEG_crpc_prim_norm$gene_name

################################################################################
################################  EMAT-SIGNATURE ###############################
################################################################################

library(readxl)
signature_emat <- read_excel("./singature_MAET.xlsx")
signature_emat <- merge(signature_emat, gff, by.x = "EMAT_genes", "gene_name", all.x = T)
signature_emat$ensg <- sapply(strsplit(as.character(signature_emat$gene_id),'\\.'), "[", 1)

signature_emat[10,1] <- ""
signature_emat[11,1] <- "KANK2"
signature_emat[31,1] <- "DEPP1"
signature_emat[32,1] <- "ADIRF"
signature_emat[33,1] <- "ZCCHC24"
signature_emat[34,1] <- "SHFL"
signature_emat[35,1] <- "C1orf116"
signature_emat[37,1] <- "KIZ"
signature_emat[39,1] <- "NREP"
signature_emat[40,1] <- "ADTRP"
signature_emat[78,1] <- "CCN2"
signature_emat[81,1] <- "CTSV"
signature_emat[115,1] <- "EVA1A"
signature_emat[120,1] <- "TENT5C"
signature_emat[132,1] <- "MICU1"
signature_emat[133,1] <- "SYBU"
signature_emat[151,1] <- "H2BC12"
signature_emat[152,1] <- "H4C15"
signature_emat[179,1] <- "FAM169A"
signature_emat[180,1] <- "ERMP1"
signature_emat[195,1] <- "P3H2"
signature_emat[200,1] <- ""
signature_emat[201,1] <- ""
signature_emat[202,1] <- ""
signature_emat[203,1] <- ""
signature_emat[204,1] <- ""
signature_emat[224,1] <- ""
signature_emat[252,1] <- "FERMT2"
signature_emat[258,1] <- "PLPP3"
signature_emat[270,1] <- "NECTIN3"
signature_emat[273,1] <- "PLAAT4"
signature_emat[275,1] <- "ESRP1"
signature_emat[276,1] <- "ESRP2"
signature_emat[298,1] <- "SCARNA9"
signature_emat[336,1] <- "SYNC"
signature_emat[338,1] <- "EPCAM"
signature_emat[343,1] <- "ELOA3BP"
signature_emat[351,1] <- "CEMIP2"
signature_emat[361,1] <- "TP63"

################################################################################
################################  PCNA-SIGNATURE ###############################
################################################################################
pcna <- read_excel("./meta_PCNA.xlsx")
pcna <- merge(pcna, gff, by.x = "SYMBOL", by.y = "gene_name", all.x = T)

pcna[11,1] <- "TSPO2"
pcna[12,1] <- "MIS18A"
pcna[16,1] <- "CDK1"
pcna[18,1] <- "CDC45"
pcna[29,1] <- "DDX39A"
pcna[30,1] <- "HJURP"
pcna[34,1] <- "AHSP"
pcna[48,1] <- "H3-3A"
pcna[53,1] <- "PCLAF"
pcna[63,1] <- "KIF18B"
pcna[76,1] <- "CENPU"
pcna[86,1] <- "ORC6"
pcna[109,1] <- "SRSF2"
pcna[130,1] <- "NSD2"