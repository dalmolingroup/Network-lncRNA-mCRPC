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
library(RTN)

#------------------------------ CREATE RTNA-OBJECT ----------------------------#

#EMAT-signature

rtna_emat_norm_met <- tni2tna.preprocess(object = rtni_55_sample,
                                         phenotype = fold_crpc_norm_met,
                                         hits = signature_emat$EMAT_genes)

rtna_emat_prim_met <- tni2tna.preprocess(object = rtni_55_sample,
                                         phenotype = fold_crpc_prim_met,
                                         hits = signature_emat$EMAT_genes)

rtna_emat_prim_norm <- tni2tna.preprocess(object = rtni_55_sample,
                                          phenotype = fold_crpc_prim_norm,
                                          hits = signature_emat$EMAT_genes)


#PCNA signature
rtna_pcna_norm_met <- tni2tna.preprocess(object = rtni_55_sample,
                                         phenotype = fold_crpc_norm_met,
                                         hits = pcna$SYMBOL)

rtna_pcna_prim_met <- tni2tna.preprocess(object = rtni_55_sample,
                                         phenotype = fold_crpc_prim_met,
                                         hits = pcna$SYMBOL)

rtna_pcna_prim_norm <- tni2tna.preprocess(object = rtni_55_sample,
                                          phenotype = fold_crpc_prim_norm,
                                          hits = pcna$SYMBOL)
#CRPC-signature

rtna_crpc_norm_met <- tni2tna.preprocess(object = rtni_55_sample,
                                         phenotype = fold_crpc_norm_met,
                                         hits = symbol_crpc_norm_met)

rtna_crpc_prim_met <- tni2tna.preprocess(object = rtni_55_sample,
                                         phenotype = fold_crpc_prim_met,
                                         hits = symbol_crpc_prim_met)

rtna_crpc_prim_norm <- tni2tna.preprocess(object = rtni_55_sample,
                                          phenotype = fold_crpc_prim_norm,
                                          hits = symbol_crpc_prim_norm)

#--------------------------- MASTER REGULATOR ANALYSIS ------------------------#

#EMAT-SIGNATURE
#metastatic x normal
rtna_emat_norm_met <- tna.mra(rtna_emat_norm_met, pValueCutoff = 0.01, pAdjustMethod = "BH",
                              minRegulonSize = 15, tnet = "dpi")

mra_emat_norm_met <- tna.get(rtna_emat_norm_met, what = "mra")

#metastatic x primary
rtna_emat_prim_met <- tna.mra(rtna_emat_prim_met, pValueCutoff = 0.01, pAdjustMethod = "BH",
                              minRegulonSize = 15, tnet = "dpi")

mra_emat_prim_met <- tna.get(rtna_emat_prim_met, what = "mra")

#primary x normal
rtna_emat_prim_norm <- tna.mra(rtna_emat_prim_norm, pValueCutoff = 0.01, pAdjustMethod = "BH",
                               minRegulonSize = 15, tnet = "dpi")

mra_emat_prim_norm <- tna.get(rtna_emat_prim_norm, what = "mra")


#PCNA-SIGNATURE
#metastatic x normal
rtna_pcna_norm_met <- tna.mra(rtna_pcna_norm_met, pValueCutoff = 0.01, pAdjustMethod = "BH",
                              minRegulonSize = 15, tnet = "dpi")

mra_pcna_norm_met <- tna.get(rtna_pcna_norm_met, what = "mra")

#metastatic x primary
rtna_pcna_prim_met <- tna.mra(rtna_pcna_prim_met, pValueCutoff = 0.01, pAdjustMethod = "BH",
                              minRegulonSize = 15, tnet = "dpi")

mra_pcna_prim_met <- tna.get(rtna_pcna_prim_met, what = "mra")

#primary x normal
rtna_pcna_prim_norm <- tna.mra(rtna_pcna_prim_norm, pValueCutoff = 0.01, pAdjustMethod = "BH",
                               minRegulonSize = 15, tnet = "dpi")

mra_pcna_prim_norm <- tna.get(rtna_pcna_prim_norm, what = "mra")

#CRPC-signature
#metastatic x normal
rtna_crpc_norm_met <- tna.mra(rtna_crpc_norm_met, pValueCutoff = 0.01, pAdjustMethod = "BH",
                              minRegulonSize = 15, tnet = "dpi")

mra_crpc_norm_met <- tna.get(rtna_crpc_norm_met, what = "mra")

#metastatic x primary
rtna_crpc_prim_met <- tna.mra(rtna_crpc_prim_met, pValueCutoff = 0.01, pAdjustMethod = "BH",
                              minRegulonSize = 15, tnet = "dpi")

mra_crpc_prim_met <- tna.get(rtna_crpc_prim_met, what = "mra")

#primary x normal
rtna_crpc_prim_norm<- tna.mra(rtna_crpc_prim_norm, pValueCutoff = 0.01, pAdjustMethod = "BH",
                              minRegulonSize = 15, tnet = "dpi")

mra_crpc_prim_norm <- tna.get(rtna_crpc_prim_norm, what = "mra")