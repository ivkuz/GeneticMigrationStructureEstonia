#####################################################################
# Script making the PDFs with maps - mean PGS and EA in the regions #
# Figure 3; Supplementary Figures 14-25 #############################
#####################################################################


library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(sf)
library(geos)
library(scales) 
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

# Adjustment PRS and EA for demography, PCs
adjustPRSandEA <- function(ebb){
  
  lm_prs2 <- lm(paste0("PRS ~ Sex + Age + I(Age^2) + I(Sex*Age) + ", paste("PC", 1:100, sep = "", collapse = " + ")),
                data = ebb)
  ebb$PRS_adj2 <- ebb$PRS - predict(object = lm_prs2, newdata = ebb)
  ebb[, PRS_adj2 := (PRS_adj2-mean(PRS_adj2))/sd(PRS_adj2)]
  
  lm_ey2 <- lm(paste0("EduYears ~ Sex + Age + I(Age^2) + I(Sex*Age) +", paste("PC", 1:100, sep = "", collapse = " + ")),
               data = ebb)
  ebb$EY_adj2 <- ebb$EduYears - predict(object = lm_ey2, newdata = ebb)
  ebb[, EY_adj2 := (EY_adj2-mean(EY_adj2, na.rm = T))/sd(EY_adj2, na.rm = T)]
  
  return(ebb)
}

# Function for plotting maps
plotMaps_pgs <- function(ebb, title, maakonds, sib=FALSE){
  
  # Adjust PRS
  ebb <- adjustPRSandEA(ebb)
  
  # Calculate mean values, SD and N for counties of birth (PoB)
  ebb_summary_PoB <- ebb[, .(mean(EY_adj2, na.rm = T), mean(PRS_adj2, na.rm = T),
                             sd(EY_adj2, na.rm = T), sd(PRS_adj2, na.rm = T)), 
                         by=PoB]
  ebb_N_PoB <- ebb[, .N, by=PoB]
  ebb_summary_PoB <- merge(ebb_summary_PoB, ebb_N_PoB)
  colnames(ebb_summary_PoB) <- c("MNIMI", "mean_EY", "mean_PRS_adj2",
                                 "sd_EY", "sd_PRS_adj2", "N")
  
  # Mark whether the mean values are significantly different from zero after BH adjustment
  ebb_summary_PoB$p_EY <- p.adjust(2*pt(-abs(ebb_summary_PoB$mean_EY)/(ebb_summary_PoB$sd_EY/sqrt(ebb_summary_PoB$N)), df=ebb_summary_PoB$N - 1), method = 'BH')
  ebb_summary_PoB$p_PRS_adj2 <- p.adjust(2*pt(-abs(ebb_summary_PoB$mean_PRS_adj2)/(ebb_summary_PoB$sd_PRS_adj2/sqrt(ebb_summary_PoB$N)), df=ebb_summary_PoB$N - 1), method = 'BH')
  ebb_summary_PoB[p_EY < 0.05, sign_EY := "*"]
  ebb_summary_PoB[p_EY >= 0.05, sign_EY := ""]
  ebb_summary_PoB[p_PRS_adj2 <  0.05, sign_PRS_adj2 := "*"]
  ebb_summary_PoB[p_PRS_adj2 >= 0.05, sign_PRS_adj2 := ""]

  # Repeat the same for PoR (counties of residence)
  ebb_summary_PoR <- ebb[, .(mean(EY_adj2, na.rm = T), mean(PRS_adj2, na.rm = T),
                             sd(EY_adj2, na.rm = T), sd(PRS_adj2, na.rm = T)), 
                         by=PoR]
  ebb_N_PoR <- ebb[, .N, 
                   by=PoR]
  ebb_summary_PoR <- merge(ebb_summary_PoR, ebb_N_PoR)
  
  colnames(ebb_summary_PoR) <- c("MNIMI", "mean_EY", "mean_PRS_adj2",
                                 "sd_EY", "sd_PRS_adj2", "N")
  ebb_summary_PoR$p_EY <- p.adjust(2*pt(-abs(ebb_summary_PoR$mean_EY)/(ebb_summary_PoR$sd_EY/sqrt(ebb_summary_PoR$N)), df=ebb_summary_PoR$N - 1), method = 'BH')
  ebb_summary_PoR$p_PRS_adj2 <- p.adjust(2*pt(-abs(ebb_summary_PoR$mean_PRS_adj2)/(ebb_summary_PoR$sd_PRS_adj2/sqrt(ebb_summary_PoR$N)), df=ebb_summary_PoR$N - 1), method = 'BH')
  ebb_summary_PoR[p_EY < 0.05, sign_EY := "*"]
  ebb_summary_PoR[p_EY >= 0.05, sign_EY := ""]
  ebb_summary_PoR[p_PRS_adj2 <  0.05, sign_PRS_adj2 := "*"]
  ebb_summary_PoR[p_PRS_adj2 >= 0.05, sign_PRS_adj2 := ""]

  # Repeat the same for the difference between PoR and PoB
  ebb_summary_dif <- ebb_summary_PoR[, 1:10]
  ebb_summary_dif[, 2:5] <- ebb_summary_dif[, 2:5] - ebb_summary_PoB[, 2:5]
  ebb_summary_dif[, 6:9] <- sqrt(ebb_summary_dif[, 6:9]^2/ebb_summary_dif$N + ebb_summary_PoB[, 6:9]^2/ebb_summary_PoB$N)
  ebb_summary_dif[, N := NULL]
  ebb_summary_dif[, dgfr := ebb_summary_PoB$N + ebb_summary_PoR$N - 2]
  ebb_summary_dif$p_EY <- p.adjust(2*pt(-abs(ebb_summary_dif$mean_EY)/(ebb_summary_dif$sd_EY), df = ebb_summary_dif$dgfr), method = 'BH')
  ebb_summary_dif$p_PRS_adj2 <- p.adjust(2*pt(-abs(ebb_summary_dif$mean_PRS_adj2)/(ebb_summary_dif$sd_PRS_adj2), df = ebb_summary_dif$dgfr), method = 'BH')
  ebb_summary_dif[p_EY < 0.05, sign_EY := "*"]
  ebb_summary_dif[p_EY >= 0.05, sign_EY := ""]
  ebb_summary_dif[p_PRS_adj2 <  0.05, sign_PRS_adj2 := "*"]
  ebb_summary_dif[p_PRS_adj2 >= 0.05, sign_PRS_adj2 := ""]

  # Parameters for plotting
  nudge_x_stars <- c(0, -5, 0, -5, 0, -5, 10, -10, 0, 5, 10, -15, 0, 5, 0)*1000
  nudge_y_stars <- c(-10, -5, -10, -10, -10, -10, -10, 10, -15, -5, 2, -20, -3, -5, -15)*1000
  lim_min_prs <- min(ebb_summary_PoB$mean_PRS_adj2, ebb_summary_PoR$mean_PRS_adj2)
  lim_max_prs <- max(ebb_summary_PoB$mean_PRS_adj2, ebb_summary_PoR$mean_PRS_adj2)
  lim_min_ea <- min(ebb_summary_PoB$mean_EY, ebb_summary_PoR$mean_EY)
  lim_max_ea <- max(ebb_summary_PoB$mean_EY, ebb_summary_PoR$mean_EY)
  
  # Make plot
  # Plot only sPGS maps for sibship-based GWAS as the phenotype maps would be the same
  
  plot_map <- merge(maakonds, ebb_summary_PoB, by = "MNIMI")
  
  # Map for PRS PoB
  pl_2 <- ggplot(data = plot_map)+
    geom_sf(aes(geometry = geometry, fill = mean_PRS_adj2))+
    geom_sf_text(aes(label = sign_PRS_adj2), nudge_x = nudge_x_stars,
                 nudge_y = nudge_y_stars, size = 5, color = "black") +
    scale_fill_gradient2(name="", limits=c(lim_min_prs-0.02, lim_max_prs), 
                         breaks=seq(-1, 1, by=round((lim_max_prs-lim_min_prs+0.02)/4, 1)), 
                         labels = comma)
  
  
  plot_map <- merge(maakonds, ebb_summary_PoR, by = "MNIMI")
  
  # Map for PRS PoB
  pl_4 <- ggplot(data = plot_map) +
    geom_sf(aes(geometry = geometry, fill = mean_PRS_adj2))+
    geom_sf_text(aes(label = sign_PRS_adj2), nudge_x = nudge_x_stars,
                 nudge_y = nudge_y_stars, size = 5, color = "black") +
    scale_fill_gradient2(name="", limits=c(lim_min_prs-0.02, lim_max_prs), 
                         breaks=seq(-1, 1, by=round((lim_max_prs-lim_min_prs+0.02)/4, 1)), 
                         labels = comma)
  
  
  plot_map <- merge(maakonds, ebb_summary_dif, by = "MNIMI")
  
  # Map for PRS PoR-PoB
  pl_6 <- ggplot(data = plot_map)+
    geom_sf(aes(geometry = geometry, fill = mean_PRS_adj2))+
    geom_sf_text(aes(label = sign_PRS_adj2),
                 nudge_x = nudge_x_stars, nudge_y = nudge_y_stars, size = 5, color = "black") + 
    scale_fill_gradient2(name="", limits=c(min(plot_map$mean_PRS_adj2), max(plot_map$mean_PRS_adj2)), 
                         breaks=seq(-0.9, 1.5, by=round((max(plot_map$mean_PRS_adj2) - min(plot_map$mean_PRS_adj2))/5, 2)), 
                         labels = comma)
  
  if(!sib){
    
    plot_map <- merge(maakonds, ebb_summary_PoB, by = "MNIMI")
    
    # Map for EA PoB
    pl_1 <- ggplot(data = plot_map)+
      geom_sf(aes(geometry = geometry, fill = mean_EY))+
      geom_sf_text(aes(label = sign_EY), nudge_x = nudge_x_stars,
                   nudge_y = nudge_y_stars, size = 5, color = "black") +
      ggtitle(bquote(paste(EA[" POB"]))) +
      scale_fill_gradient2(name="", limits=c(lim_min_ea-0.02, lim_max_ea), 
                           breaks=seq(-1, 1, by=round((lim_max_ea-lim_min_ea+0.02)/4, 1)), 
                           labels = comma)
    
    plot_map <- merge(maakonds, ebb_summary_PoR, by = "MNIMI")
    
    # Map for EA PoR
    pl_3 <- ggplot(data = plot_map) +
      geom_sf(aes(geometry = geometry, fill = mean_EY))+
      geom_sf_text(aes(label = sign_EY), nudge_x = nudge_x_stars,
                   nudge_y = nudge_y_stars, size = 5, color = "black") +
      ggtitle(bquote(paste(EA[" POR"]))) +
      scale_fill_gradient2(name="", limits=c(lim_min_ea-0.02, lim_max_ea), 
                           breaks=seq(-1, 1, by=round((lim_max_ea-lim_min_ea+0.02)/4, 1)), 
                           labels = comma)
    
    plot_map <- merge(maakonds, ebb_summary_dif, by = "MNIMI")
    
    # Map for EA PoR-PoB
    pl_5 <- ggplot(data = plot_map)+
      geom_sf(aes(geometry = geometry, fill = mean_EY))+
      geom_sf_text(aes(label = sign_EY), nudge_x = nudge_x_stars, 
                   nudge_y = nudge_y_stars, size = 5, color = "black") + 
      ggtitle(bquote(paste(EA[" POR"], " - ", EA["POB"]))) +
      scale_fill_gradient2(name="", limits=c(min(plot_map$mean_EY), max(plot_map$mean_EY)), 
                           breaks=seq(-0.9, 1.5, by=round((max(plot_map$mean_EY) - min(plot_map$mean_EY))/5, 2)), 
                           labels = comma)
    
    # Add titles
    pl_2 <- pl_2 + ggtitle(bquote(paste(PGS[" EA, POB"])))
    pl_4 <- pl_4 + ggtitle(bquote(paste(PGS[" EA, POR"])))
    pl_6 <- pl_6 + ggtitle(bquote(paste(PGS[" EA, POR"], " - ", PGS[" EA, POB"])))
    
    plots <- list(pl_1, pl_2, pl_3, pl_4, pl_5, pl_6)

  } else{
    
    # Add titles
    pl_2 <- pl_2 + ggtitle(bquote(paste(sPGS[" EA, POB"])))
    pl_4 <- pl_4 + ggtitle(bquote(paste(sPGS[" EA, POR"])))
    pl_6 <- pl_6 + ggtitle(bquote(paste(sPGS[" EA, POR"], " - ", sPGS[" EA, POB"])))
    
    plots <- list(pl_2, pl_4, pl_6)

  }
  
  # Make cirles for the cities
  city <- data.frame(
    name = c("Tartu", "Tallinn"),
    x_p = c(660578, 541078), 
    y_p = c(6471700, 6588700),
    x_t = c(689578, 500078), 
    y_t = c(6490700, 6601700)
  )
  
  # Add universal parameters to all the plots
  for(i in 1:length(plots)){
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
  make_option(c("-f", "--pheno"), type="character", default=NULL, 
              help="file with basic phenotypes for filtered sample", metavar="character"),
  make_option(c("-g", "--geo"), type="character", default=NULL, 
              help="file with info on places of birth and residence", metavar="character"),
  make_option(c("-l", "--lastedu"), type="character", default=NULL, 
              help="file with info on last Education", metavar="character"),
  make_option(c("-e", "--pca_est"), type="character", default=NULL, 
              help="file with 100 PCs of estonians", metavar="character"),
  make_option(c("-r", "--pca_rus"), type="character", default=NULL, 
              help="file with 100 PCs of russians", metavar="character"),
  make_option(c("-p", "--pop_prs"), type="character", default=NULL, 
              help="Population GWAS polygenic score for EA", metavar="character"),
  make_option(c("-s", "--sib_prs"), type="character", default=NULL, 
              help="Sibship-based GWAS polygenic score for EA", metavar="character"),
  make_option(c("-u", "--nonrel"), type="character", default="non_relatives.tsv", 
              help="list of non-related individuals [default= %default]", metavar="character"),
  make_option(c("-m", "--map"), type="character", default="maakond_shp/maakond_20230201.shp", 
              help="map with counties' borders (.shp object) [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="results.tsv", 
              help="output directory name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# Read the main table with basic data for the filtered sample
estbb_filtered <- fread(opt$pheno)
# Add data on places of birth and residence
ebb <- fread(opt$geo)
ebb <- ebb[, c("Person skood", "PersonLocation birthParishName", "PersonLocation residencyParishName")]
colnames(ebb) <- c("skood", "ParishBirth", "ParishRes")
ebb <- merge(estbb_filtered, ebb, by="skood")
# Make age2 and age*sex covariates
ebb[, Age2 := Age^2]
ebb[, SxA := Age*Sex]
# Add education data
ebb2 <- fread(opt$lastedu)
ebb2 <- ebb2[, c("Person skood", "PersonPortrait lastEducation code")]
colnames(ebb2) <- c("skood", "EA")
ebb <- merge(ebb, ebb2, by="skood")
# Transform education level to Years of Education (EduYears) and binary "University degree" (EA)
ebb[, EduYears := ransformEAtoEduYears(EA)]
ebb[EA<6, EA := 0]
ebb[EA>=6, EA := 1]
# Add PC coordinates for Estonian and Russian subsamples separately
pca_est <- fread(opt$pca_est)
pca_est <- pca_est[,c("IID", paste0("PC", 1:100))]
colnames(pca_est)[1] <- "vkood"
pca_rus <- fread(opt$pca_rus)
pca_rus <- pca_rus[,c("IID", paste0("PC", 1:100))]
colnames(pca_rus)[1] <- "vkood"
pca <- rbind(pca_est, pca_rus)
ebb <- merge(ebb, pca, by="vkood")
# Add sibship-based GWAS PRS
ebb_sib <- fread(opt$sib_prs)
colnames(ebb_sib) <- c("vkood", "PRS")
ebb_sib <- merge(ebb, ebb_sib)
# Add population-based GWAS PRS
prs <- fread(opt$pop_prs)
colnames(prs) <- c("vkood", "PRS")
ebb <- merge(ebb, prs)
# List of unrelated Estonian individuals
non_rel <- fread(opt$nonrel)

# Read Shp file with the spatial data on the counties' borders
maakonds <-  st_read(opt$map)
maakonds$MNIMI <- gsub(x=maakonds$MNIMI, pattern = " maakond", replacement = "")


# Make PDFs with plots for the subgroups and population/sibship PRS

pdf(paste0(opt$out, "/Maps_E_V.pdf"), width=7, height=2.5)
plots <- plotMaps_pgs(ebb[Nat=="Eestlane",], maakonds = maakonds)
print(grid.arrange(grobs = plots[c(2,4,6)], layout_matrix = rbind(c(1, 2, 3))))
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.351, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.684, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
print(grid.arrange(grobs = plots[c(1,3,5)], layout_matrix = rbind(c(1, 2, 3))))
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.351, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.684, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

plots <- plotMaps_pgs(ebb[Nat=="Venelane",], maakonds = maakonds)
print(grid.arrange(grobs = plots[c(2,4,6)], layout_matrix = rbind(c(1, 2, 3))))
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.351, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.684, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
print(grid.arrange(grobs = plots[c(1,3,5)], layout_matrix = rbind(c(1, 2, 3))))
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.351, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.684, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf(paste0(opt$out, "/Maps_E_nonrel.pdf"), width=7, height=2.5)
plots <- plotMaps_pgs(ebb[vkood %in% non_rel$vkood,], maakonds = maakonds)
print(grid.arrange(grobs = plots[c(2,4,6)], layout_matrix = rbind(c(1, 2, 3))))
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.351, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.684, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
print(grid.arrange(grobs = plots[c(1,3,5)], layout_matrix = rbind(c(1, 2, 3))))
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.351, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.684, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf(paste0(opt$out, "/Maps_E_sib.pdf"), width=7, height=2.5)
plots <- plotMaps_pgs(ebb_sib[Nat=="Eestlane",], maakonds = maakonds, sib=T)
print(grid.arrange(grobs = plots, layout_matrix = rbind(c(1, 2, 3))))
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.351, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.684, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf(paste0(opt$out, "/Maps_sex.pdf"), width=7, height=5)
plots_m <- plotMaps_pgs(ebb[Nat=="Eestlane" & Sex=="1",], maakonds = maakonds)
plots_f <- plotMaps_pgs(ebb[Nat=="Eestlane" & Sex=="2",], maakonds = maakonds)
print(grid.arrange(grobs = c(plots_m[c(2,4,6)], plots_f[c(2,4,6)]), 
                   layout_matrix = rbind(c(1, 2, 3),
                                         c(4, 5, 6))))
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.351, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.684, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("E", x = 0.351, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("F", x = 0.684, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))

print(grid.arrange(grobs = c(plots_m[c(1,3,5)], plots_f[c(1,3,5)]),
                   layout_matrix = rbind(c(1, 2, 3),
                                         c(4, 5, 6))))
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.351, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.684, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("E", x = 0.351, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("F", x = 0.684, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf(paste0(opt$out, "/Maps_yoa.pdf"), width=7, height=5)
plots_pre2016 <- plotMaps_pgs(ebb[Nat=="Eestlane" & YoA<=2016,], maakonds = maakonds)
plots_post2016 <- plotMaps_pgs(ebb[Nat=="Eestlane" & YoA>2016,], maakonds = maakonds)
print(grid.arrange(grobs = c(plots_pre2016[c(2,4,6)], plots_post2016[c(2,4,6)]), 
                   layout_matrix = rbind(c(1, 2, 3),
                                         c(4, 5, 6))))
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.351, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.684, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("E", x = 0.351, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("F", x = 0.684, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))

print(grid.arrange(grobs = c(plots_pre2016[c(1,3,5)], plots_post2016[c(1,3,5)]),
                   layout_matrix = rbind(c(1, 2, 3),
                                         c(4, 5, 6))))
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.351, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.684, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("E", x = 0.351, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("F", x = 0.684, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
dev.off()


pdf(paste0(opt$out, "/Maps_age.pdf"), width=7, height=10)
plots_age1 <- plotMaps_pgs(ebb[Nat=="Eestlane" & Age>=18 & Age<=24,], maakonds = maakonds)
plots_age2 <- plotMaps_pgs(ebb[Nat=="Eestlane" & Age>=25 & Age<=48,], maakonds = maakonds)
plots_age3 <- plotMaps_pgs(ebb[Nat=="Eestlane" & Age>=49 & Age<=64,], maakonds = maakonds)
plots_age4 <- plotMaps_pgs(ebb[Nat=="Eestlane" & Age>=65,], maakonds = maakonds)
print(grid.arrange(grobs = c(plots_age1[c(2,4,6)], plots_age2[c(2,4,6)],
                             plots_age3[c(2,4,6)], plots_age4[c(2,4,6)]), 
                   layout_matrix = rbind(c(1, 2, 3),
                                         c(4, 5, 6),
                                         c(7, 8, 9),
                                         c(10, 11, 12))))
grid.text("A", x = 0.018, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.351, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.684, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.018, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("E", x = 0.351, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("F", x = 0.684, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("G", x = 0.018, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("H", x = 0.351, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("I", x = 0.684, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("J", x = 0.018, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("K", x = 0.351, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("L", x = 0.684, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))

print(grid.arrange(grobs = c(plots_age1[c(1,3,5)], plots_age2[c(1,3,5)],
                             plots_age3[c(1,3,5)], plots_age4[c(1,3,5)]), 
                   layout_matrix = rbind(c(1, 2, 3),
                                         c(4, 5, 6),
                                         c(7, 8, 9),
                                         c(10, 11, 12))))
grid.text("A", x = 0.018, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.351, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.684, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.018, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("E", x = 0.351, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("F", x = 0.684, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("G", x = 0.018, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("H", x = 0.351, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("I", x = 0.684, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("J", x = 0.018, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("K", x = 0.351, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("L", x = 0.684, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))
dev.off()
