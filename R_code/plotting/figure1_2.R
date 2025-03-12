####################################################################################
# Script making the PDFs with Var(county) diagrams #################################
# Figures 1 and 2; Supplementary Figures with EA4 PGS and with unrelated subsample #
####################################################################################


library(data.table)
library(ggplot2)
library(reshape2)
library(scales)
library(ggnewscale)
library(grid)
library(gridExtra)

library(Rcpp)
library(optparse)



option_list = list(
  make_option(c("-c", "--pc"), type="character", default=NULL, 
              help="file with F-statistics and Var(county) info for PCs", metavar="character"),
  make_option(c("-s", "--pgs"), type="character", default=NULL, 
              help="file with F-statistics and Var(county) info for PGSs", metavar="character"),
  make_option(c("-a", "--pgsadj"), type="character", default=NULL, 
              help="file with F-statistics and Var(county) info for PGSs adjusted for PGSea", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output figure path", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



# Upload a file with Fstat for PCs
Fstat_pca <- fread(opt$pc, header = T)
Fstat_pca$Trait <- factor(Fstat_pca$Trait, levels = paste0("PC",1:100))

# PGSs UK Biobank
Fstat <- fread(opt$pgs, header = T)

# PGSs UK Biobank adjusted for PGS for EA
Fstat_pgsEAadj <- fread(opt$pgsadj)

# Get min and max p-values for the line coloring
Fstat[, pmax := max(pval), by = c("Trait", "Nat", "Sex", "Age", "YoA")]
Fstat[, pmin := min(pval), by = c("Trait", "Nat", "Sex", "Age", "YoA")]

# Get the order of PGSs according to the difference between Var(county) for POR and POB
Fstat_wide <- as.data.table(dcast(Fstat[Nat == "E" & Sex == "All" & Age == "All" & YoA == "All" , .(Trait, BorR, prop_var)], 
                                  formula = Trait ~ BorR, value.var = "prop_var"))
Fstat_wide[, delta_prop_var := PoR - PoB]
Fstat_wide <- Fstat_wide[, .(Trait, delta_prop_var)]
Fstat_wide <- Fstat_wide[, n := rank(delta_prop_var, ties.method = "first")]
Fstat <- merge(Fstat, Fstat_wide, by = "Trait")
# write.table(unique(Fstat[, .(Trait, n)]), "VarCounty_ranks.tsv", row.names = F, quote = F, sep = "\t")

# Prepare table with PGSs adjusted for PGS for EA
Fstat_pgsEAadj[, pmax := max(pval), by = c("Trait", "Nat", "Sex", "Age", "YoA")]
Fstat_pgsEAadj[, pmin := min(pval), by = c("Trait", "Nat", "Sex", "Age", "YoA")]
Fstat_wide[, Trait := paste0(Trait, "_adjPGSea")]
Fstat_pgsEAadj <- merge(Fstat_pgsEAadj, Fstat_wide, by = "Trait")




pdf(opt$out, width=7, height=6)
for(nat in c("E")){
  
  # Plot PCs
  Fstat_pca_tmp <- Fstat_pca[Nat==nat & Sex=="All" & Age=="All" & YoA == "All", ]
  
  plot_pca_fstat <- ggplot(data = Fstat_pca_tmp, aes(y = Trait, x = prop_var*100)) +
    geom_line(aes(group = Trait, color=(P<0.05/100)))+#, linewidth=1) +
    geom_point(aes(fill = BorR, color=(pval<0.05/100)),
               stroke = 0.7, pch=21, size=1.5) +
    scale_fill_manual(values = c("#CC0000", "#009900"), name = "By place of",
                      labels = c("POB", "POR")) +
    scale_color_manual(values = c("orange", "grey"), name = "FWER",
                       labels = c("Significant", "Not significant"),
                       breaks=c(T, F)) +
    coord_flip() + labs(x = bquote(paste(Var[county], ", %")),
                        y = "Principal components") +
    theme_bw() + scale_shape_manual(values=c(21:25)) +
    theme(legend.position="none") +
    scale_y_discrete(breaks = paste0("PC", c(1, seq(10, 100, 10))),
                     labels = c(1, seq(10, 100, 10))) +
    scale_x_continuous(trans='log10', labels = comma)
  
  # Plot PGSs
  Fstat_tmp <- Fstat[Nat==nat & Sex=="All" & Age == "All" & YoA == "All", ]
  Fstat_tmp <- Fstat_tmp[grep("_adj2$", Trait)]
  plot_pgs_order_fstat <- ggplot(data = Fstat_tmp[Trait != "EA4_adj2",], aes(x = n, y = prop_var*100)) + 
    geom_line(aes(group = Trait, color=pmax<0.05/169 & P<0.05/169 | 
                    pmax>0.05/169 & pmin<0.05/169), linewidth=0.7) +
    geom_point(aes(fill = BorR, color=(pval<0.05/169)),
               stroke = 0.7, pch=21, size=1.5) +
    scale_color_manual(values = c("orange", "grey"), name = "FWER",
                       labels = c("Significant", "Not significant"),
                       breaks=c(T, F), guide = "none") +
    new_scale_color() +
    scale_color_manual(values = c("#CC0000", "#009900"), name = "By place of",
                       labels = c("POB", "POR"), guide = "none") +
    scale_fill_manual(values = c("#CC0000", "#009900"), name = "By place of",
                      labels = c("POB", "POR")) +
    labs(y = bquote(paste(Var[county], ", %")), 
         x = bquote(paste(rank[PGS], "(", Var[county], "[POR] - ", Var[county], "[POB])"))) +
    theme_bw() + scale_shape_manual(values=c(21:25)) +
    theme(text = element_text(size = 10),
          legend.position="none"#,

    ) + 
    scale_y_continuous(breaks = seq(0, 2.25, 0.25), limits = c(0, 1.2)) #+
  
  # Plot PGSs adjusted for PGE for EA
  Fstat_tmp <- Fstat_pgsEAadj[Nat==nat & Sex=="All" & Age == "All" & YoA == "All", ]
  Fstat_tmp <- Fstat_tmp[grep("_adj2_adjPGSea$", Trait)]
  plot_pgs_adj_order_fstat <- ggplot(data = Fstat_tmp, aes(x = n, y = prop_var*100)) +
    scale_y_continuous(breaks = seq(0, 2.25, 0.25), limits = c(0, 1.2)) +
    geom_line(aes(group = Trait, color=pmax<0.05/169 & P<0.05/169 | 
                    pmax>0.05/169 & pmin<0.05/169), linewidth=0.7) +
    geom_point(aes(fill = BorR, color=(pval<0.05/169)),
               stroke = 0.7, pch=21, size=1.5) +
    scale_color_manual(values = c("orange", "grey"), name = "FWER",
                       labels = c("Significant", "Not significant"),
                       breaks=c(T, F)) +
    new_scale_color() +
    scale_color_manual(values = c("#CC0000", "#009900"), name = "By place of",
                       labels = c("POB", "POR")) +
    scale_fill_manual(values = c("#CC0000", "#009900"), name = "By place of",
                      labels = c("POB", "POR")) +
    labs(y = bquote(paste(Var[county], ", %")), 
         x = bquote(paste(rank[PGS], "(", Var[county], "[POR] - ", Var[county], "[POB])"))) +
    theme_bw() + scale_shape_manual(values=c(21:25)) +
    theme(text = element_text(size = 10),
          legend.direction = "horizontal", legend.box = "horizontal",
          legend.position = "bottom",
          legend.title=element_blank(),
          legend.spacing.y = unit(0, "pt")
    )
    
  
  # Make a single plot
  plots <- list(plot_pca_fstat, plot_pgs_order_fstat, plot_pgs_adj_order_fstat)
  
  print(
    grid.arrange(
      grobs = plots,
      heights = c(1, 1, 1.3),
      layout_matrix = matrix(1:3, ncol = 1)
    )
  )
  
  
  grid.text("A", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
  grid.text("B", x = 0.02, y = 0.68, gp = gpar(fontsize=14, fontface = "bold"))
  grid.text("C", x = 0.02, y = 0.38, gp = gpar(fontsize=14, fontface = "bold"))
  
}

dev.off()

