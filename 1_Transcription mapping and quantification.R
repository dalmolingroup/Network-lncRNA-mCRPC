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

#Data specified: GSE126078

#Packages used in this analysis:
library(tximport)
library(tidyverse)
library(GenomicFeatures)
library(ensembldb)
library(readxl)



################################################################################
########################### IMPORT SAMPLE METADATA  ############################
################################################################################
metadata <- read.table("metadata.txt", header = T, stringsAsFactors = F, sep = ",")
#Filter only samples of metastatic tissue
metadata <- metadata[grep("CRPC", metadata$source_name),]
#Exclude unnecessary columns
metadata <- metadata[,c("Run", "Cell_type", "source_name", "molecular_phenotype", "PATIENT", "tumor_site")]
#Put rownames
row.names(metadata) <- metadata$Run
#Clean metadata
metadata$tumor_site <- str_replace_all(metadata$tumor_site, " ", "_")
metadata$AR <- sapply(strsplit(as.character(metadata$molecular_phenotype),'_'), "[", 1)
metadata$NE <- sapply(strsplit(as.character(metadata$molecular_phenotype),'_'), "[", 2)
metadata$tumor_site <- as.factor(metadata$tumor_site)
metadata$PATIENT <- as.factor(metadata$PATIENT)

#Exclude individuals with repeated metastatic site

`%ni%` <- Negate(`%in%`)

exclusao <- c("SRR8528519", "SRR8528523","SRR8528525", "SRR8528527", "SRR8528529",
              "SRR8528532", "SRR8528536", "SRR8528537", "SRR8528540", "SRR8528543",
              "SRR8528545", "SRR8528547", "SRR8528549", "SRR8528550", "SRR8528552",
              "SRR8528553", "SRR8528555", "SRR8528557", "SRR8528561",
              "SRR8528563", "SRR8528564", "SRR8528566", "SRR8528570", "SRR8528571",
              "SRR8528573", "SRR8528575", "SRR8528577", "SRR8528584", "SRR8528585",
              "SRR8528586", "SRR8528589", "SRR8528591", "SRR8528593", "SRR8528594",
              "SRR8528596", "SRR8528598", "SRR8528601", "SRR8528604", "SRR8528606",
              "SRR8528609", "SRR8528610", "SRR8528613", "SRR8528614")


metadata_55 <- metadata[metadata$Run %ni% exclusao,]

metadata_55 %>%
  dplyr::count(tumor_site)


################################################################################
########################### IMPORT SALMON SAMPLES   ############################
################################################################################

#Open list with file names
lista <- read.table("./lista_arquivos_sf.txt", header = F, stringsAsFactors = F)
lista <- lista %>%
  mutate("ID" = str_extract(V1, "SRR.*[0-9]"))

#tx2gene object
txdb_transcript <- GenomicFeatures::makeTxDbFromGFF(file = "./gencode.v40.annotation.gff3")
ke_transcript <- keys(txdb_transcript, keytype = "TXNAME")
tx2gene_transcript <- ensembldb::select(txdb_transcript, ke_transcript, "GENEID", "TXNAME")
#clean tx2gene object
tx2gene_transcript$TXNAME <- sapply(strsplit(as.character(tx2gene_transcript$TXNAME),'\\.'), "[", 1)

#Import files
txi_transcript <- tximport(paste("./quant_sf_transcript/",lista$V1,sep=""), type = "salmon", 
                           tx2gene = tx2gene_transcript, 
                           ignoreTxVersion = T, txOut = F, txIn = T, ignoreAfterBar = T)
save(txi_transcript, file = "txi_transcript.RData")

#Get count data
count_data <- as.data.frame(txi_transcript$counts)
colnames(count_data) <- lista$ID

#Filter only lncRNAs and mRNA coding_protein using annotation of Salmon
sf <- read.delim("./quant_sf_transcript/SRR8528517_quant.sf", header=T, comment.char="#")
sf <- data.frame(do.call("rbind", strsplit(as.character(sf$Name), "|", fixed = TRUE)))
sf_filter <- sf[sf$X8 %in% c("lncRNA", "protein_coding"),]

unique(sf_filter$X8)
sum(sf_filter$X8 == "lncRNA") #50746 lncRNAs
sum(sf_filter$X8 == "protein_coding") #87591 protein-coding

#Count data:lncRNA and protein-coding
sum(row.names(count_data) %in% unique(sf_filter$X2))
count_data_filter <- count_data[row.names(count_data) %in% unique(sf_filter$X2),]
#Count data: lncRNAs
sum(sf_filter$X8 == "lncRNA")
lncRNAs_all <- sf_filter[sf_filter$X8 == "lncRNA", "X2"]
lncRNA_filter <- row.names(count_data_filter[row.names(count_data_filter) %in% lncRNAs_all,])
#Count data: protein_coding
sum(sf_filter$X8 == "protein_coding")
coding_all <- sf_filter[sf_filter$X8 == "protein_coding", "X2"]
coding_filter <- row.names(count_data_filter[row.names(count_data_filter) %in% coding_all,])

#Filter lncRNAs with 90% variance
count_lncRNA <- count_data_filter %>%
  dplyr::filter(row.names(count_data_filter) %in% lncRNA_filter)

count_lncRNA$variance  <- apply(count_lncRNA, 1, var)
count_lncRNA <- count_lncRNA[count_lncRNA$variance >= quantile(count_lncRNA$variance, c(0.90)),]
count_lncRNA$variance <- NULL

count_coding <- count_data_filter %>%
  dplyr::filter(row.names(count_data_filter) %in% coding_filter )
all(colnames(count_lncRNA) == colnames(count_coding))
count_network <- rbind(count_lncRNA, count_coding)

save(count_network, file = "count_network.RData")