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
library(DESeq2)
library(ggpubr)

################################################################################
############################### BOXPLOT OF MRS #################################
################################################################################
# Reference:https://rpubs.com/turnersd/plot-deseq-results-multipage-pdf

#MRs of  prim x met | norm x met | and both
#Link one by one and match in dds-tcga
#Are not in dds_tcga: ENSG00000289443.1, ENSG00000289194.1

met <- c("BMPR1B-DT","ENSG00000226332","ENSG00000261012","ENSG00000261098",
         "ENSG00000268230","ENSG00000270504","FAM99A","GTF3C2-AS2","HAND2-AS1",
         "HELLPAR","LINC00261","LINC01485","LINC01595","MEG3", "SNHG18")

#All MRs
#Are not in dds_tcga: ENSG00000289443, ENSG00000289194.1, ENSG00000289080.1
ensg_MRs <- gff[gff$gene_name %in% met,"gene_id"]
ensg_MRs <- sapply(strsplit(as.character(ensg_MRs),'\\.'), "[", 1)

ensg_MRs <- substr(paste0(ensg_MRs, sep = "|", collapse = ""),1,
                   nchar(paste0(ensg_MRs, sep = "|", collapse = ""))-1)

dds_crpc_MRs <- dds_crpc[ grep(ensg_MRs, row.names(dds_crpc)) , ]


#Normalization

tcounts <- t(log2((DESeq2::counts(dds_crpc_MRs, normalized=TRUE, replaced=FALSE)+0.5))) %>%
  merge(SummarizedExperiment::colData(dds_crpc_MRs), ., by="row.names") %>%
  tidyr::gather(gene, expression, (ncol(.)-length(rownames(dds_crpc_MRs))+1):ncol(.))


#Gene symbol
tcounts$gene_id <- sapply(strsplit(as.character(tcounts$gene),'\\.'), "[", 1)
tcounts <- merge(tcounts, gff, by.x = "gene_id", by.y = "id", all.x = T)
tcounts$sample_type <- as.character(tcounts$sample_type)

tcounts$sample_type[tcounts$sample_type == 'CRPC'] <- 'Primary_Tumor'

#modified: ENSG00000237125.11 (HAND2-AS1), ENSG00000250786.2 (SNGH18),
#ENSG00000245526.12(LINC00461), ENSG00000235142.11 (LINC02532),
#ENSG00000130600.20 (H19), ENSG00000259974.4 (LINC00261)
#ENSG00000240801.1 without counts

my_compare <- list(c("Metastatic", "Primary_Tumor"), c("Metastatic", "Solid_Tissue_Normal"), c("Solid_Tissue_Normal", "Primary_Tumor"))
ggboxplot(tcounts, x = "sample_type", y = "expression",
          color = "sample_type", palette = "jco") + 
  facet_wrap(~gene_name, scales="free_y") + 
  scale_x_discrete(labels = c("Met", "Pri", "Norm")) + 
  scale_color_discrete(name = "Groups",
                       labels = c("Metastatic", "Primary Tumor", "Normal Tissue")) +
  stat_compare_means( comparisons = my_compare, label = "p.signif", p.adjust.method = "bonferroni") + 
  labs(x="", y="Expression (log normalized counts)")

#Option 2

library(rstatix)
stat_test <- tcounts %>%
  dplyr::group_by(gene_name) %>%
  rstatix::wilcox_test(expression ~ sample_type) %>%
  rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_significance("p.adj")

save(tcounts, file = "tcounts.RData")
stat_test <- stat_test %>%
  add_xy_position(x = "gene_name", dodge = 0.8)

stat_test$y.position <- rep(c(8,9,10), 15)


#Ordenar
tcounts_2 <- tcounts %>% dplyr::arrange(factor(gene_name, levels = unique(stat_test$gene_name)))

#Boxplot
bxp <- ggboxplot(tcounts_2, x = "gene_name", y = "expression", 
                 color = "sample_type")

bxp + 
  stat_pvalue_manual(stat_test, label = "p.adj.signif", 
                     tip.length = 0.01,bracket.nudge.y = 6, vjust = 0) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1, size = 10)) +
  labs(x="", y="Expression (log normalized counts)") +
  scale_color_discrete(name = "Groups",
                       labels = c("Metastatic", "Primary Tumor", "Normal Tissue"))
