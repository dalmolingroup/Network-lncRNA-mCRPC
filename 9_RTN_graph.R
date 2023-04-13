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
library(RedeR)
library(igraph)

################################################################################
################################## RTN GRAPH  ##################################
################################################################################

### ----------------------------- Network ---------------------------------- ###

### Netowrk of MRs from MetPri and Metnorm
g_55 <- tni.graph(rtni_55_sample, tnet = "dpi", minRegulonSize = 15, gtype = "rmap",
                  regulatoryElements = all_met_crpc)


rdp <- RedPort()
calld(rdp, checkcalls = T)
addGraph(rdp, g_55, layout = NULL)
addLegend.color(rdp, g_55, type="edge")
addLegend.shape(rdp, g_55)
relax(rdp, ps = TRUE)

### HELLPAR network
g_55 <- tni.graph(rtni_55_sample, tnet = "dpi", minRegulonSize = 15, gtype = "rmap",
                  regulatoryElements = "HELLPAR")
map_HELLPAR <- data.frame(target = c("HELLPAR", names(regulons_55[["HELLPAR"]])),
                          MI = c(0, regulons_55[["HELLPAR"]]))
rownames(map_HELLPAR) <- map_HELLPAR$target
all(igraph::as_data_frame(x = g_55, what ="vertices")$SYMBOL == map_HELLPAR$target)
g_55 <- att.mapv(g_55, dat = map_HELLPAR, refcol = 1)
color_code <- RColorBrewer::brewer.pal(n = 6, name = "Reds")
breaks <- seq(0, 1, by = 0.2)

g_55 <- att.setv(g = g_55, from = "MI", to = "nodeColor", breaks = breaks, cols = color_code)

calld(rdp, checkcalls = T)
addGraph(rdp, g_55, layout = layout.kamada.kawai(g_55))
addLegend.color(rdp, 
                colvec=g_55$legNodeColor$scale, 
                labvec=g_55$legNodeColor$legend, 
                title="Node Color (Mutual Information)")

### SNHG18 netowrk

g_55 <- tni.graph(rtni_55_sample, tnet = "dpi", minRegulonSize = 15, gtype = "rmap",
                  regulatoryElements = "SNHG18")
map_SNHG18 <- data.frame(target = c("SNHG18", names(regulons_55[["SNHG18"]])),
                         MI = c(0, regulons_55[["SNHG18"]]))
rownames(map_SNHG18) <- map_SNHG18$target
all(igraph::as_data_frame(x = g_55, what ="vertices")$SYMBOL == map_SNHG18$target)
g_55 <- att.mapv(g_55, dat = map_SNHG18, refcol = 1)
color_code <- RColorBrewer::brewer.pal(n = 6, name = "Reds")
breaks <- seq(0, 0.5, by = 0.1)
g_55 <- att.setv(g = g_55, from = "MI", to = "nodeColor", breaks = breaks, cols = color_code)

calld(rdp, checkcalls = T)
addGraph(rdp, g_55, layout = layout.kamada.kawai(g_55))
addLegend.color(rdp, 
                colvec=g_55$legNodeColor$scale, 
                labvec=g_55$legNodeColor$legend, 
                title="Node Color (Mutual Information)")

#### --------------------------- Tree-and-Leaf ----------------------------- ###

load("./tree_55_symbol_15_size.RData")
tree_55_symbol_15_size <- tni.graph(rtni_55_sample, tnet = "dpi", minRegulonSize = 15, gtype = "amapDend")


rdp <- RedPort()
calld(rdp, checkcalls = T)
addGraph(rdp, tree_55_symbol_15_size$g, layout.kamada.kawai(tree_55_symbol_15_size$g))
selectNodes(rdp,all_met_crpc)

#### --------------------------- Association-map --------------------------- ###

g_r <- tni.graph(rtni_55_sample, tnet = "dpi", minRegulonSize = 15,
                 gtype = "amap",
                 regulatoryElements = all_met_crpc, amapFilter = "phyper",
                 amapCutoff = 1)

rdp <- RedPort()
calld(rdp, checkcalls = T)
addGraph(rdp, g_r)
relax(rdp, p8=100, p4=1, p3=30, p1=40)
selectEdges(rdp, c("BMPR1B-DT", "ENSG00000268230"))
