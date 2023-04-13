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
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
library(ggplot2)

#Select exclusive MRs

shared_metastatic <- c("SNHG18", "MEG3")
MN_crpc <- c("FAM99A","ENSG00000289080","ENSG00000261012","HELLPAR","HAND2-AS1","LINC01485","LINC01595",
             "LINC00261","ENSG00000289194","BMPR1B-DT","ENSG00000240801","ENSG00000270504","ENSG00000268230")
MP_crpc <- c("ENSG00000289443","GTF3C2-AS2","ENSG00000261098","ENSG00000226332")
NP_crpc <- c("ENSG00000238260","ENSG00000289154","H1-10-AS1","MMP25-AS1","ENSG00000286223","ENSG00000275185","LINC01341" )


cluster_exclusive <- list()

for (lnc in NP_crpc){
  long = names(regulons_55[[lnc]])
  
  eg_lnc = bitr(long, fromType= "SYMBOL", toType = c("SYMBOL","ENTREZID"),
                OrgDb = org.Hs.eg.db, drop = T)
  
  cluster_exclusive[lnc] <- list(eg_lnc$ENTREZID)
}

rm(long, eg_lnc, lnc)

#-------------------------------- Gene Ontology -------------------------------#

CompareGO_NP = compareCluster(cluster_exclusive, fun="enrichGO", pvalueCutoff=0.01,
                              pAdjustMethod="BH", OrgDb=org.Hs.eg.db, ont="BP",
                              readable=T)
CompareGO_NP_table <- as_tibble(CompareGO_NP)


#Table with all enrichment
all_enrich <- rbind(CompareGO_MP_table, CompareGO_MN_table, CompareGO_NP_table, CompareGO_shared_table)
all_enrich$Group <- NA
for (i in 1:nrow(all_enrich)){
  if(all_enrich$Cluster[i] %in% MN_crpc){
    all_enrich$Group[i] <- "MetNorm" 
  } else if (all_enrich$Cluster[i] %in% MP_crpc){
    all_enrich$Group[i] <- "MetPri"
  } else if (all_enrich$Cluster[i] %in% NP_crpc){
    all_enrich$Group[i] <- "PriNorm"
  } else {
    all_enrich$Group[i] <- "MetNorm_∩_MetPri"
  }
}
all_enrich <- all_enrich[,c(11,1:10)]
write.csv(all_enrich, file = "all_enrich.csv")


#Function to Count Genes in GO
countGO <- function(pathways){
  go_ <- matrix(nrow = length(pathways), ncol = 2)
  colnames(go_) <- c("Description","Count")
  linha <- 0
  for(i in pathways){
    linha <- linha + 1
    tabela <- all_enrich[all_enrich$Description == i, ]
    genes_ <- NULL
    for (j in 1:nrow(tabela)){
      genes_ <- append(genes_, unlist(strsplit(as.character(tabela[j,"geneID"]), "/")))
    }
    go_[linha,"Description"] <- i
    go_[linha,"Count"] <- length(unique(genes_))
  }
  go_ <- as.data.frame(go_)
  go_$Count <- as.numeric(go_$Count)
  return(go_)
  rm(tabela, genes_, linha, go_)
}


#Barplot of GO enrichment

#GO_MetPri
path_metpri <- all_enrich$Description[c(2,10,12)]
go_metpri <- countGO(pathways = path_metpri)
ggplot(data=go_metpri, aes(x=Count,y=Description )) +
  geom_bar(position="dodge",stat="identity", fill = "#c83153") +
  ylab("")+
  theme_classic()


#GO_MetNorm
path_metnorm <- all_enrich$Description[c(22,23,30,45,47,78,86,310,349,350,624,632,640,643,677)]

go_metnorm <- countGO(pathways = path_metnorm)

lim <- go_metnorm %>%
  arrange(Count) %>%
  pull(Description)

ggplot(data=go_metnorm, aes(x=Count,y=Description)) +
  geom_bar(position="dodge",stat="identity", fill = "turquoise4") +
  ylab("")+
  scale_y_discrete(name = " ", limits = lim)+
  theme_classic()



#GO_PriNorm
path_prinorm <- all_enrich$Description[c(850,854,855,860,879,899)]
go_prinorm <- countGO(pathways = path_prinorm)

lim <- go_prinorm %>%
  arrange(Count) %>%
  pull(Description)

ggplot(data=go_prinorm, aes(x=Count,y=Description )) +
  geom_bar(position="dodge",stat="identity", fill = "steelblue") +
  ylab("")+
  scale_y_discrete(name = " ", limits = lim)+
  theme_classic()

#Shared

path_shared <- all_enrich$Description[c(905,917,930,935,948)]
go_shared <- countGO(pathways = path_shared)

lim <- go_shared %>%
  arrange(Count) %>%
  pull(Description)

ggplot(data=go_shared, aes(x=Count,y=Description )) +
  geom_bar(position="dodge",stat="identity", fill = "gold") +
  ylab("")+
  scale_y_discrete(name = " ", limits = lim)+
  theme_classic()