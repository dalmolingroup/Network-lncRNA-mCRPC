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

################################################################################
############################## SYMBOLIC REGRESSION MATRIX ######################
################################################################################

#Symbolic Regression 
tab_reg <- matrix_CRPC
rownames(tab_reg) <- sapply(strsplit(as.character(rownames(tab_reg)),'\\.'), "[", 1)
tab_reg <- merge(tab_reg, gff, by.x = "row.names", by.y = "id", all.x = T)

#not ENSG00000289080, ENSG00000289194,ENSG00000289443
all_mras_crpc <- unique(c(mra_crpc_norm_met$Regulon,
                          mra_crpc_prim_met$Regulon,
                          mra_crpc_prim_norm$Regulon))

tab_reg <- tab_reg[tab_reg$gene_name %in% all_mras_crpc,]
rownames(tab_reg) <- tab_reg$gene_name
tab_reg <- tab_reg[,-c(1,165:167)]

#Normalizacao
tab_reg <- as.data.frame(edgeR::cpm(tab_reg))

tab_reg <- t(tab_reg)
all(row.names(clinical_CRPC) == row.names(tab_reg))


tab_reg <- merge(tab_reg, clinical_CRPC, by = "row.names", all.x = T )
colnames(tab_reg)[1] <- "run"
colnames(tab_reg)[30] <- "phenotype"
tab_reg$phenotype_reg <- ifelse(tab_reg$phenotype == "Solid_Tissue_Normal", 0,
                                ifelse(tab_reg$phenotype == "CRPC", 1, 2))

tab_reg <- tab_reg[,-29]

write.csv(tab_reg, file = "crpc_run2_MRs_prim.csv", row.names = F)
