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
library(UpSetR)

#List of all MRs
lista_crpc <- list("Primary x\nNormal" = row.names(mra_crpc_prim_norm), 
                   "Metastatic x\nPrimary" = row.names(mra_crpc_prim_met),
                   "Metastatic x\nNormal" = row.names(mra_crpc_norm_met))

#List of MRs from MetPri and MetNorm
names(lista_crpc) <- c("PriNorm", "MetPri", "MetNorm")
lista_met <- lista_crpc[c(2,3)]
lista_met[["MetPri"]] <- lista_met[["MetPri"]][c(2,4,5,6,9,10)]
lista_met[["MetNorm"]] <- lista_met[["MetNorm"]][2:16]

#Upset of all MRs
upset(fromList(lista_crpc), order.by = "freq",
      mainbar.y.label = "Intersection Size", sets.x.label = "MRs by Signature",
      text.scale = 1.5,
      #      sets.bar.color = c("turquoise4","dodgerblue2","hotpink3"),
      main.bar.color = c("turquoise4","dodgerblue2","hotpink3","sienna3", "gold","red4"),
      line.size = 1,
      shade.alpha = 0.25)

#Upset of MRs from MetPri and MetNorm
upset(fromList(lista_met), order.by = "freq",
      mainbar.y.label = "Intersection Size", sets.x.label = "MRs by Signature",
      text.scale = 1.5,
      #      sets.bar.color = c("turquoise4","#c83153"),
      main.bar.color = c("turquoise4","#c83153", "gold"),
      line.size = 1,
      shade.alpha = 0.25)