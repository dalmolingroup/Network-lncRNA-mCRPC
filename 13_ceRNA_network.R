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

#Libraries
devtools::install_github("RafaellaFerraz/encoriR")
library(encoriR)
library(RedeR)
library(igraph)

#ceRNA Investigation
#Using ENCORI

ceRNA <- encoriR::mirna_target(assembly = "hg38", 
                               geneType = "lncRNA", 
                               miRNA = "all", 
                               target = c("HELLPAR", "SNHG18"))
#Filter All cell line
ceRNA_complete <- ceRNA[,c("miRNAname", "geneName", "cellline.tissue")]
ceRNA_complete <- dplyr::distinct(ceRNA_complete)

#Filter to 22rv1 (prostate cancer cell line)
ceRNA_filter <- ceRNA_complete[grep(pattern = "22RV1", ceRNA_complete[["cellline.tissue"]]),]
rownames(ceRNA_filter) <- NULL

#####################

load("./regulons_55.RData")
hellpar_targets <- names(regulons_55[["HELLPAR"]])

#GFF
gff <- rtracklayer::import.gff3("./gencode.v40.annotation.gff3")
gff <- as.data.frame(gff@elementMetadata) 
gff <- gff[,c("gene_id" ,"gene_type", "gene_name")]
gff <- dplyr::distinct(gff)
gff$id <- sapply(strsplit(as.character(gff$gene_id),'\\.'), "[", 1)

############ HELLPAR

hellpar_names <- merge(as.data.frame(hellpar_targets), gff,by.x = "hellpar_targets", by.y = "gene_name")

ceRNA_hellpar_mRNA <- encoriR::mirna_target(assembly = "hg38", 
                                            geneType = "mRNA", 
                                            miRNA = "all", 
                                            target = dplyr::filter(hellpar_names, gene_type == "protein_coding")[["hellpar_targets"]])


ceRNA_hellpar_lncRNA <- encoriR::mirna_target(assembly = "hg38", 
                                              geneType = "lncRNA", 
                                              miRNA = "all", 
                                              target = dplyr::filter(hellpar_names, gene_type == "lncRNA")[["hellpar_targets"]])


cell_prost <- unique(ceRNA_hellpar_mRNA[["cellline.tissue"]][grep("22RV", ceRNA_hellpar_mRNA[["cellline.tissue"]])])
miRNA_mRNA_HELLPAR <- unique(dplyr::filter(ceRNA_hellpar_mRNA, cellline.tissue %in% cell_prost)[["miRNAname"]])

cell_prost2 <- unique(ceRNA_hellpar_lncRNA[["cellline.tissue"]][grep("22RV", ceRNA_hellpar_lncRNA[["cellline.tissue"]])])
miRNA_lncRNA_HELLPAR <- unique(dplyr::filter(ceRNA_hellpar_lncRNA, cellline.tissue %in% cell_prost2)[["miRNAname"]])

miRNA_targets_HELLPAR <- unique(c(miRNA_mRNA_HELLPAR, miRNA_lncRNA_HELLPAR))

intersect_HELPAR_targets <- intersect(ceRNA_filter[["miRNAname"]], miRNA_targets_HELLPAR)

#ceRNA network - 22rv1
gr_22rv1 <- igraph::graph_from_edgelist(as.matrix(ceRNA_filter[,c("miRNAname", "geneName")]), 
                                        directed = F)

rdp <- RedeR::RedPort()
RedeR::calld(rdp, checkcalls = T)
RedeR::addGraph(rdp, gr_22rv1, layout = NULL)
RedeR::selectNodes(rdp, intersect_HELPAR_targets)


#Save table
a <- dplyr::filter(ceRNA_hellpar_mRNA, cellline.tissue %in% cell_prost)[,c("miRNAname", "geneName", "cellline.tissue")]
b <- dplyr::filter(ceRNA_hellpar_lncRNA, cellline.tissue %in% cell_prost2)[,c("miRNAname", "geneName", "cellline.tissue")]
ceRNA_hell <- distinct(rbind(ceRNA_filter, a, b))
ceRNA_hell$cellline.tissue <- gsub(",", "/", ceRNA_hell$cellline.tissue)
write.csv(ceRNA_hell, file = "ceRNA_network.csv", quote = F, row.names = F, col.names = T)
