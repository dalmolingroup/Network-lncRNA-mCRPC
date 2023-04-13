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

################################################################################
################################## lncRNA-lncRNA  ##############################
################################################################################

count_long <- NULL
regula_long <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(regula_long) <- c("lncRNAs", "target_lncRNAs")
linha <- 1
for (long in all_met_crpc){
  alvos <- names(regulons_55[[long]])
  if (! rlang::is_empty(alvos) & 
      any(alvos %in% gff[gff$gene_type == "lncRNA","gene_name"]) &
      length(alvos) >= 15){
    count_long <- append(count_long, long)
    ids <- which(alvos %in% gff[gff$gene_type == "lncRNA","gene_name"])
    alvos_ll <- NULL
    for (x in ids){
      alvos_ll <- append(alvos_ll, alvos[x]) 
    }
    alvos_ll <- paste(alvos_ll, collapse = ",")
    regula_long[linha, "lncRNAs"] <- long
    regula_long[linha, "target_lncRNAs"] <- alvos_ll
    linha <- linha + 1
  }
}

rm(alvos_ll, alvos, linha, ids, long)

unlist(strsplit(regula_long[17,2], split = ",")) %in% all_met_crpc

write.csv(regula_long, file = "met_lncRNA_target.csv")