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
library(rtracklayer)
library(parallel)
library(snow)
library(RTN)

################################################################################
###################### TRANCRIPTIONAL NETWORK INFERENCE ########################
################################################################################
load("count_network.RData") #count_network

#Import gff file
library(rtracklayer)
gff <- import.gff3("./gencode.v40.annotation.gff3")
gff <- as.data.frame(gff@elementMetadata) 
gff <- gff[,c("gene_id" ,"gene_type", "gene_name")]
gff <- distinct(gff)
gff$id <- sapply(strsplit(as.character(gff$gene_id),'\\.'), "[", 1)

#Convert ENSEMBL ID to SYMBOL
matrix_name <- merge(count_network, gff, by.x = "row.names", by.y = "gene_id", all.x = T)
#duplicated genes
matrix_name[duplicated(matrix_name$gene_name),c("gene_name")]
dup <- unique(matrix_name[duplicated(matrix_name$gene_name),c("gene_name")])
local_dup <- row.names(matrix_name[matrix_name$gene_name %in% dup,])
matrix_name$Row.names <- sapply(strsplit(as.character(matrix_name$Row.names),'\\.'), "[", 1)
for (x in local_dup){
  matrix_name[x, "gene_name"] <- matrix_name[x, "Row.names"]
}
#Put rownames
row.names(matrix_name) <- matrix_name$gene_name

#Rowannotation
RowAnnot <- matrix_name[,c("Row.names", "gene_name", "gene_type")]
row.names(RowAnnot) <- matrix_name$gene_name
names(RowAnnot) <- c("ENSEMBL","SYMBOL", "TYPE")

#Converte dataframe to matrix
matrix_name <- data.matrix(matrix_name[,-c(1,100,101)], rownames.force = T)

#Regulatory elements
lncRNAs_symbol <- unique(sf_filter[sf_filter$X8 == "lncRNA", "X6"])

#rtni object
rtni_constructor_90lnc_55sample <- tni.constructor(expData = matrix_name,
                                                   regulatoryElements = lncRNAs_symbol,
                                                   colAnnotation = metadata,
                                                   rowAnnotation = RowAnnot)

save(rtni_constructor_90lnc_55sample, file = "rtni_constructor_90lnc_98sample.RData")


library("parallel")
options(cluster=snow::makeCluster(spec=32, "SOCK"))

#permutation
rtni_55_sample <- tni.permutation(rtni_constructor_90lnc_55sample, pValueCutoff = 0.01)
save(rtni_55_sample, file = "rtni_90_perm_55_sample.RData")
#bootstrap
rtni_55_sample <- tni.bootstrap(rtni_55_sample)
save(rtni_98_sample, file = "rtni_90_boot_55_sample.RData")
#ARACNe
rtni_55_sample <- tni.dpi.filter(rtni_55_sample, eps = 0)
save(rtni_55_sample, file = "rtni_90_final_55_sample.RData")

stopCluster(getOption("cluster"))


######### Summary

tni.regulon.summary(rtni_55_sample)
regulons_55 <- tni.get(rtni_55_sample, what = "regulons.and.mode")
regs <- tni.get(rtni_55_sample, what = "regulonSize", idkey = "SYMBOL")

#count target
count <- 0
for (i in 1:length(regulons_55)){
  if (length(regulons_55[[i]]) >= 15){
    count <- count + 1
  }
}