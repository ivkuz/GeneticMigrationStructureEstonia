#############################################################################
# Script making the PDFs with maps with mean PGS per county #################
# and forest plots for the difference POR-POB in population and in siblings #
# Figure 3 ##################################################################
#############################################################################


library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(sf)
library(geos)
library(scales) 

library(Rcpp)
library(optparse)


# Transform EA catigories to the years of education
transformEAtoEduYears <- function(EA_vector){
  
  ea_names <- EA_vector
  ea_names[!(ea_names %in% 1:8)] <- NA
  ea_names[ea_names == 8] <- 22
  ea_names[ea_names == 7] <- 20
  ea_names[ea_names == 6] <- 18
  ea_names[ea_names == 5] <- 15
  ea_names[ea_names == 4] <- 13
  ea_names[ea_names %in% c(2, 3)] <- 10
  ea_names[ea_names == 1] <- 7
  ea_names[ea_names == 0] <- 1
  ea_names <- as.numeric(ea_names)
  return(ea_names)
}


# Adjust PGS for demography, PCs and EA
adjustPRSandEA <- function(ebb){
  
  lm_prs2 <- lm(paste0("PRS ~ Sex + Age + I(Age^2) + I(Sex*Age) + ", paste("PC", 1:100, sep = "", collapse = " + ")),
                data = ebb)
  ebb$PRS_adj2 <- ebb$PRS - predict(object = lm_prs2, newdata = ebb)
  ebb[, PRS_adj2 := (PRS_adj2-mean(PRS_adj2))/sd(PRS_adj2)]
  
  lm_prs3 <- lm(paste0("PRS ~ Sex + Age + I(Age^2) + I(Sex*Age) + EA +", paste("PC", 1:100, sep = "", collapse = " + ")),
                data = ebb)
  ebb$PRS_adj3 <- ebb$PRS - predict(object = lm_prs3, newdata = ebb)
  ebb[, PRS_adj3 := (PRS_adj3-mean(PRS_adj3, na.rm = T))/sd(PRS_adj3, na.rm = T)]
  
  lm_prs4 <- lm(paste0("PRS ~ Sex + Age + I(Age^2) + I(Sex*Age) + EduYears +", paste("PC", 1:100, sep = "", collapse = " + ")),
                data = ebb)
  ebb$PRS_adj4 <- ebb$PRS - predict(object = lm_prs4, newdata = ebb)
  ebb[, PRS_adj4 := (PRS_adj4-mean(PRS_adj4, na.rm = T))/sd(PRS_adj4, na.rm = T)]
  
  lm_ey2 <- lm(paste0("EduYears ~ Sex + Age + I(Age^2) + I(Sex*Age) +", paste("PC", 1:100, sep = "", collapse = " + ")),
               data = ebb)
  ebb$EY_adj2 <- ebb$EduYears - predict(object = lm_ey2, newdata = ebb)
  ebb[, EY_adj2 := (EY_adj2-mean(EY_adj2, na.rm = T))/sd(EY_adj2, na.rm = T)]
  
  return(ebb)
}


# Plot maps
plotMaps_pgs <- function(ebb, title, maakonds, sib=FALSE){
  
  ebb <- adjustPRSandEA(ebb)
  ebb_summary_PoB <- ebb[, .(mean(EY_adj2, na.rm = T), mean(PRS_adj2, na.rm = T), mean(PRS_adj3, na.rm = T), mean(PRS_adj4, na.rm = T),
                             sd(EY_adj2, na.rm = T), sd(PRS_adj2, na.rm = T), sd(PRS_adj3, na.rm = T), sd(PRS_adj4, na.rm = T)), 
                         by=PoB]
  ebb_N_PoB <- ebb[, .N, 
                   by=PoB]
  ebb_summary_PoB <- merge(ebb_summary_PoB, ebb_N_PoB)
  
  # Get asterisks for POB
  colnames(ebb_summary_PoB) <- c("MNIMI", "mean_EY", "mean_PRS_adj2", "mean_PRS_adj3", "mean_PRS_adj4",
                                 "sd_EY", "sd_PRS_adj2", "sd_PRS_adj3", "sd_PRS_adj4", "N")
  ebb_summary_PoB$p_EY <- p.adjust(2*pt(-abs(ebb_summary_PoB$mean_EY)/(ebb_summary_PoB$sd_EY/sqrt(ebb_summary_PoB$N)), df=ebb_summary_PoB$N - 1), method = 'BH')
  ebb_summary_PoB$p_PRS_adj2 <- p.adjust(2*pt(-abs(ebb_summary_PoB$mean_PRS_adj2)/(ebb_summary_PoB$sd_PRS_adj2/sqrt(ebb_summary_PoB$N)), df=ebb_summary_PoB$N - 1), method = 'BH')
  ebb_summary_PoB$p_PRS_adj3 <- p.adjust(2*pt(-abs(ebb_summary_PoB$mean_PRS_adj3)/(ebb_summary_PoB$sd_PRS_adj3/sqrt(ebb_summary_PoB$N)), df=ebb_summary_PoB$N - 1), method = 'BH')
  ebb_summary_PoB$p_PRS_adj4 <- p.adjust(2*pt(-abs(ebb_summary_PoB$mean_PRS_adj4)/(ebb_summary_PoB$sd_PRS_adj4/sqrt(ebb_summary_PoB$N)), df=ebb_summary_PoB$N - 1), method = 'BH')
  ebb_summary_PoB[p_EY < 0.05, sign_EY := "*"]
  ebb_summary_PoB[p_EY >= 0.05, sign_EY := ""]
  ebb_summary_PoB[p_PRS_adj2 <  0.05, sign_PRS_adj2 := "*"]
  ebb_summary_PoB[p_PRS_adj2 >= 0.05, sign_PRS_adj2 := ""]
  ebb_summary_PoB[p_PRS_adj3 <  0.05, sign_PRS_adj3 := "*"]
  ebb_summary_PoB[p_PRS_adj3 >= 0.05, sign_PRS_adj3 := ""]
  ebb_summary_PoB[p_PRS_adj4 <  0.05, sign_PRS_adj4 := "*"]
  ebb_summary_PoB[p_PRS_adj4 >= 0.05, sign_PRS_adj4 := ""]
  
  
  ebb_summary_PoR <- ebb[, .(mean(EY_adj2, na.rm = T), mean(PRS_adj2, na.rm = T), mean(PRS_adj3, na.rm = T), mean(PRS_adj4, na.rm = T),
                             sd(EY_adj2, na.rm = T), sd(PRS_adj2, na.rm = T), sd(PRS_adj3, na.rm = T), sd(PRS_adj4, na.rm = T)), 
                         by=PoR]
  ebb_N_PoR <- ebb[, .N, 
                   by=PoR]
  ebb_summary_PoR <- merge(ebb_summary_PoR, ebb_N_PoR)
  
  # Get asterisks for POR
  colnames(ebb_summary_PoR) <- c("MNIMI", "mean_EY", "mean_PRS_adj2", "mean_PRS_adj3", "mean_PRS_adj4",
                                 "sd_EY", "sd_PRS_adj2", "sd_PRS_adj3", "sd_PRS_adj4", "N")
  ebb_summary_PoR$p_EY <- p.adjust(2*pt(-abs(ebb_summary_PoR$mean_EY)/(ebb_summary_PoR$sd_EY/sqrt(ebb_summary_PoR$N)), df=ebb_summary_PoR$N - 1), method = 'BH')
  ebb_summary_PoR$p_PRS_adj2 <- p.adjust(2*pt(-abs(ebb_summary_PoR$mean_PRS_adj2)/(ebb_summary_PoR$sd_PRS_adj2/sqrt(ebb_summary_PoR$N)), df=ebb_summary_PoR$N - 1), method = 'BH')
  ebb_summary_PoR$p_PRS_adj3 <- p.adjust(2*pt(-abs(ebb_summary_PoR$mean_PRS_adj3)/(ebb_summary_PoR$sd_PRS_adj3/sqrt(ebb_summary_PoR$N)), df=ebb_summary_PoR$N - 1), method = 'BH')
  ebb_summary_PoR$p_PRS_adj4 <- p.adjust(2*pt(-abs(ebb_summary_PoR$mean_PRS_adj4)/(ebb_summary_PoR$sd_PRS_adj4/sqrt(ebb_summary_PoR$N)), df=ebb_summary_PoR$N - 1), method = 'BH')
  ebb_summary_PoR[p_EY < 0.05, sign_EY := "*"]
  ebb_summary_PoR[p_EY >= 0.05, sign_EY := ""]
  ebb_summary_PoR[p_PRS_adj2 <  0.05, sign_PRS_adj2 := "*"]
  ebb_summary_PoR[p_PRS_adj2 >= 0.05, sign_PRS_adj2 := ""]
  ebb_summary_PoR[p_PRS_adj3 <  0.05, sign_PRS_adj3 := "*"]
  ebb_summary_PoR[p_PRS_adj3 >= 0.05, sign_PRS_adj3 := ""]
  ebb_summary_PoR[p_PRS_adj4 <  0.05, sign_PRS_adj4 := "*"]
  ebb_summary_PoR[p_PRS_adj4 >= 0.05, sign_PRS_adj4 := ""]
  
  
  # Get asterisks for POR - POB
  ebb_summary_dif <- ebb_summary_PoR[, 1:10]
  ebb_summary_dif[, 2:5] <- ebb_summary_dif[, 2:5] - ebb_summary_PoB[, 2:5]
  ebb_summary_dif[, 6:9] <- sqrt(ebb_summary_dif[, 6:9]^2/ebb_summary_dif$N + ebb_summary_PoB[, 6:9]^2/ebb_summary_PoB$N)
  ebb_summary_dif[, N := NULL]
  ebb_summary_dif[, dgfr := ebb_summary_PoB$N + ebb_summary_PoR$N - 2]
  ebb_summary_dif$p_EY <- p.adjust(2*pt(-abs(ebb_summary_dif$mean_EY)/(ebb_summary_dif$sd_EY), df = ebb_summary_dif$dgfr), method = 'BH')
  ebb_summary_dif$p_PRS_adj2 <- p.adjust(2*pt(-abs(ebb_summary_dif$mean_PRS_adj2)/(ebb_summary_dif$sd_PRS_adj2), df = ebb_summary_dif$dgfr), method = 'BH')
  ebb_summary_dif$p_PRS_adj3 <- p.adjust(2*pt(-abs(ebb_summary_dif$mean_PRS_adj3)/(ebb_summary_dif$sd_PRS_adj3), df = ebb_summary_dif$dgfr), method = 'BH')
  ebb_summary_dif$p_PRS_adj4 <- p.adjust(2*pt(-abs(ebb_summary_dif$mean_PRS_adj4)/(ebb_summary_dif$sd_PRS_adj4), df = ebb_summary_dif$dgfr), method = 'BH')
  ebb_summary_dif[p_EY < 0.05, sign_EY := "*"]
  ebb_summary_dif[p_EY >= 0.05, sign_EY := ""]
  ebb_summary_dif[p_PRS_adj2 <  0.05, sign_PRS_adj2 := "*"]
  ebb_summary_dif[p_PRS_adj2 >= 0.05, sign_PRS_adj2 := ""]
  ebb_summary_dif[p_PRS_adj3 <  0.05, sign_PRS_adj3 := "*"]
  ebb_summary_dif[p_PRS_adj3 >= 0.05, sign_PRS_adj3 := ""]
  ebb_summary_dif[p_PRS_adj4 <  0.05, sign_PRS_adj4 := "*"]
  ebb_summary_dif[p_PRS_adj4 >= 0.05, sign_PRS_adj4 := ""]
  
  # Locate asterisks
  nudge_x_stars <- c(0, -5, 0, -5, 0, -5, 10, -10, 0, 5, 10, -15, 0, 5, 0)*1000
  nudge_y_stars <- c(-10, -5, -10, -10, -10, -10, -10, 10, -15, -5, 2, -20, -3, -5, -15)*1000
  lim_min_prs <- min(ebb_summary_PoB$mean_PRS_adj2, ebb_summary_PoR$mean_PRS_adj2)
  lim_max_prs <- max(ebb_summary_PoB$mean_PRS_adj2, ebb_summary_PoR$mean_PRS_adj2)
  lim_min_ea <- min(ebb_summary_PoB$mean_EY, ebb_summary_PoR$mean_EY)
  lim_max_ea <- max(ebb_summary_PoB$mean_EY, ebb_summary_PoR$mean_EY)
  
  
  
  # Maps for POB
  plot_map <- merge(maakonds, ebb_summary_PoB, by = "MNIMI")
  
  pl_2 <- ggplot(data = plot_map)+
    geom_sf(aes(geometry = geometry, fill = mean_PRS_adj2))+
    geom_sf_text(aes(label = sign_PRS_adj2), nudge_x = nudge_x_stars,
                 nudge_y = nudge_y_stars, size = 5, color = "black") +
    ggtitle(bquote(paste(PGS[" EA, POB"]))) +
    scale_fill_gradient2(name="", limits=c(lim_min_prs-0.02, lim_max_prs), 
                         breaks=seq(-1, 1, by=round((lim_max_prs-lim_min_prs+0.02)/4, 1)), 
                         labels = comma)
  
  
  # Maps for POR
  plot_map <- merge(maakonds, ebb_summary_PoR, by = "MNIMI")
  
  pl_4 <- ggplot(data = plot_map) +
    geom_sf(aes(geometry = geometry, fill = mean_PRS_adj2))+
    geom_sf_text(aes(label = sign_PRS_adj2), nudge_x = nudge_x_stars,
                 nudge_y = nudge_y_stars, size = 5, color = "black") +
    ggtitle(bquote(paste(PGS[" EA, POR"]))) +
    scale_fill_gradient2(name="", limits=c(lim_min_prs-0.02, lim_max_prs), 
                         breaks=seq(-1, 1, by=round((lim_max_prs-lim_min_prs+0.02)/4, 1)), 
                         labels = comma)
  
  
  # Forest-plot for POR-POB
  ebb_summary_dif$MNIMI <- factor(ebb_summary_dif$MNIMI, levels = ebb_summary_dif$MNIMI[order(ebb_summary_dif$mean_PRS_adj2)])
  
  pl_8 <- ggplot(ebb_summary_dif, aes(x = mean_PRS_adj2, y = MNIMI, label = sign_PRS_adj2)) + #, fill = mean_PRS_adj2
    geom_vline(aes(xintercept = 0), linewidth = .25, linetype = "dashed") +
    geom_errorbarh(aes(xmax = mean_PRS_adj2 + 1.96*sd_PRS_adj2, xmin = mean_PRS_adj2 - 1.96*sd_PRS_adj2), linewidth = .3, height = .2, color = "gray50") +
    geom_point(size = 2, shape=21, fill = "gray") +
    coord_trans(x = scales:::exp_trans(10)) +
    scale_x_continuous(breaks = round(seq(-0.15, 0.05, 0.05), digits = 2), labels = round(seq(-0.15, 0.05, 0.05), digits = 2),
                       limits = c(-0.15, 0.05)) + # limits = c(-0.16, 0.05)) +
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          legend.position="none",
          text = element_text(size = 7)) +
    ylab("County") +
    xlab(bquote(paste(PGS[" EA, POR"], " - ", PGS[" EA, POB"]))) +
    geom_text(aes(x = mean_PRS_adj2 + 1.96*sd_PRS_adj2), hjust = -0.6, vjust = 0.8, size = 4, color = "black") +
    theme()
  
  
  plots <- list(pl_2, pl_4, pl_8)

  # Add the coordinates of the cities on the map
  city <- data.frame(
    name = c("Tartu", "Tallinn"),
    x_p = c(660578, 541078), 
    y_p = c(6471700, 6588700),
    x_t = c(689578, 500078), 
    y_t = c(6490700, 6601700)
  )
  
  # Add common features to the plots
  for(i in 1:2){
    plots[[i]] <- plots[[i]] + 
      theme_classic() + 
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size=8, hjust = 0.2, vjust = -0.8), 
            legend.key.width = unit(1, 'cm'),
            legend.key.height = unit(0.1, 'cm'),
            legend.position = "bottom",
            text = element_text(size = 7)) +
      geom_point(data=city, aes(x=x_p, y=y_p), shape=21, stroke=1.2) +
      geom_label(data=city, aes(x=x_t, y=y_t, label=name), 
                 size=2.5, label.padding=unit(0.5, "mm"))
  }
  
  return(plots)
  
  
}




option_list = list(
  make_option(c("-e", "--pca_est"), type="character", default=NULL, 
              help="file with 100 PCs of estonians", metavar="character"),
  make_option(c("-r", "--pca_rus"), type="character", default=NULL, 
              help="file with 100 PCs of russians", metavar="character"),
  make_option(c("-f", "--pheno"), type="character", default=NULL, 
              help="file with phenotypes", metavar="character"),
  make_option(c("-s", "--score"), type="character", default=NULL, 
              help="file with vkoods and PRS", metavar="character"),
  make_option(c("-m", "--map"), type="character", default=NULL, 
              help="file with geographic map in .shp format", metavar="character"),
  make_option(c("-s", "--sib"), type="character", default=NULL, 
              help="file with results on distribution of within-sibship deviations of PGSs in the counties", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output figure path", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



# File with PRS
prs <- fread(opt$score)
colnames(prs) <- c("vkood", "PRS")
ebb <- merge(ebb, prs)

# File with PCs
PCA_est <- fread(opt$pca_est)
PCA_rus <- fread(opt$pca_rus)
PCA <- rbind(PCA_est, PCA_rus)
PCA <- PCA[, FID:=NULL]
colnames(PCA)[1] <- "vkood"

# File with phenotypes
ebb <- fread(opt$pheno)

# Combine PRS, PCA, Phenotypes
ebb <- mergeData(PRS=PRS, PCA=PCA, Phenotypes=ebb)


# Upload geographical data
maakonds <-  st_read(opt$map)
maakonds$MNIMI <- gsub(x=maakonds$MNIMI, pattern = " maakond", replacement = "")


pdf(opt$out, width=7*2/3, height = 5)

# plots for the main Estonian sample
plots <- plotMaps_pgs(ebb[Nat=="Eestlane",], maakonds = maakonds)

# Forest-plot for POR-POB for siblings
results <- fread(opt$sib)
results <- results[Variable == "EA4_adj2" & Residence == "PoR" & County != "All", ]
results[, p_adj := p.adjust(p, method = "BH")]
results[p_adj < 0.05, sign := "*"]
results[p_adj >= 0.05, sign := ""]

results$County <- factor(results$County, levels = c("Järva", "Lääne-Viru", "Valga", 
                                                    "Jõgeva", "Pärnu", "Võru",
                                                    "Rapla", "Viljandi", "Saare",
                                                    "Ida-Viru", "Hiiu", "Põlva",
                                                    "Lääne", "Tartu", "Harju"))
psib <- ggplot(results, aes(x = mean, y = County, label = sign)) + # , fill = mean
  geom_vline(aes(xintercept = 0), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = mean + 1.96*sd/sqrt(N), xmin = mean - 1.96*sd/sqrt(N)), linewidth = .3, height = 0.2, color = "gray50") +
  geom_point(size = 2, shape=21, fill = "gray") +
  coord_trans(x = scales:::exp_trans(10)) +
  scale_x_continuous(breaks = round(seq(-0.15, 0.05, 0.05), digits = 2), labels = round(seq(-0.15, 0.05, 0.05), digits = 2),
                     limits = c(-0.15, 0.05)) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position="none",
        text = element_text(size = 7)) +
  ylab("County") +
  xlab(bquote(paste("Within-family ", PGS[" EA4, POR"]))) +
  geom_text(aes(x = mean + 1.96*sd/sqrt(N)), hjust = -0.6, vjust = 0.8, size = 4, color = "black")

# Plot all the panels together
plots[[4]] <- psib
print(grid.arrange(grobs = plots, layout_matrix = rbind(c(1, 2), c(3, 4))))
grid.text("A", x = 0.03, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.53, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.03, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.53, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


