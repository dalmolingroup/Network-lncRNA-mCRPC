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
library(TCGAbiolinks)
library(tidyverse)
library(readxl)
library(rtracklayer)

#TCGA: primary and normal samples
query_PRAD <- GDCquery(project = "TCGA-PRAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       experimental.strategy = "RNA-Seq",
                       workflow.type = "STAR - Counts",
                       sample.type = c("Primary Tumor", "Solid Tissue Normal"))
results_query <- getResults(query_PRAD)
GDCdownload(query = query_PRAD,
            method = "api")
PRAD <- GDCprepare(query = query_PRAD,
                   summarizedExperiment = T,
                   directory = "./GDCdata/")

#TCGA: metastatic sample
query_met <- GDCquery(project = "WCDT-MCRPC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      experimental.strategy = "RNA-Seq",
                      workflow.type = "STAR - Counts")

GDCdownload(query = query_met,
            method = "api")

MET <- GDCprepare(query = query_met,
                  summarizedExperiment = T,
                  directory = "./GDCdata/")

#GEO: primary CRPC
meta_CRPC <- read.table("./Projeto_TFs_mama_prostata/Artigo_Tfs_Mama_Prostata/Metadata_SRP073789.txt",
                        header = T, sep = ",")

counts_CRPC <- read.table("./Projeto_TFs_mama_prostata/Artigo_Tfs_Mama_Prostata/raw_counts_SRP073789_CRPCprim.txt",
                          header = T, sep = "\t")
match_CRPC <- match(meta_CRPC[meta_CRPC$progression_step == "CRPC","Run"],colnames(counts_CRPC))

counts_CRPC_filter <- counts_CRPC[,match_CRPC]

#GFF
library(rtracklayer)
gff <- import.gff3("./Projeto_TFs_mama_prostata/Artigo_Tfs_Mama_Prostata/gencode.v36.annotation.gff3")
gff <- as.data.frame(gff@elementMetadata)
gff <- gff[,c("gene_id" ,"gene_type", "gene_name")]
gff <- distinct(gff)

sum(gff$gene_type %in% c("lncRNA", "protein_coding"))

gff40 <- import.gff3("./Projeto_TFs_mama_prostata/Artigo_Tfs_Mama_Prostata/gencode.v40.annotation.gff3")
gff40 <- as.data.frame(gff40@elementMetadata)
gff40 <- gff40[,c("gene_id" ,"gene_type", "gene_name")]
gff40 <- distinct(gff40)

sum(gff40$gene_type %in% c("lncRNA", "protein_coding"))


#Merge dataframes
library(SummarizedExperiment)

all(rownames(assay(PRAD)) == rownames(assay(MET)))

matrix_tcga <- merge(assay(PRAD), assay(MET), by = "row.names")
matrix_tcga <- matrix_tcga %>%
  column_to_rownames(var = "Row.names")

sum(rownames(matrix_tcga) %in% counts_CRPC$gene_id)

row.names(counts_CRPC_filter) <- counts_CRPC$gene_id

matrix_CRPC <- merge(matrix_tcga, counts_CRPC_filter, by = "row.names")

matrix_CRPC <- matrix_CRPC %>%
  column_to_rownames(var = "Row.names")

match_ln_prot <- match(gff40[gff40$gene_type %in% c("protein_coding" ,"lncRNA"),"gene_id"],row.names(matrix_CRPC))

matrix_CRPC_filter <- matrix_CRPC[match_ln_prot,]
matrix_CRPC_filter <- matrix_CRPC_filter[order(row.names(matrix_CRPC_filter)) , ]


save(matrix_CRPC, file = "counts_all_elements_CPRCprimario.RData")
save(matrix_CRPC_filter, file = "counts_ln_prot_CRPCprimario.RData")



#Clinical Data
all(rownames(colData(PRAD)) %in% colnames(assay(PRAD)))
all(rownames(colData(MET)) %in% colnames(assay(MET)))

a <- as.data.frame(colData(PRAD))
a <- a[,c("barcode", "sample_type")]
colnames(a)[1] <- "sample"
b <- as.data.frame(colData(MET))
b <- b[,c("sample", "sample_type")]
c <- meta_CRPC[meta_CRPC$progression_step == "CRPC", c("Run", "progression_step")]
colnames(c) <- c("sample", "sample_type")
rownames(c) <- c$sample
clinical_CRPC <- rbind(a,b,c)

save(clinical_CRPC, file = "clinical_CRPC.RData")
