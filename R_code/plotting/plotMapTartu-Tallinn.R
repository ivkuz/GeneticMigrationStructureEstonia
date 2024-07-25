#########################################################################################################################
# Script making the PDFs with maps - differences in mean PGS and EA between migrants to Tallinn or Tartu in the regions #
# Supplementary Figures 32-37 ###########################################################################################
#########################################################################################################################


library(tidyverse)
library(sf)
library(data.table)


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
plotMaps <- function(ebb, title, maakonds, s=FALSE){

  # Keep individuals residing in Tartu or Tallinn and who wasn't born there
  ebb <- ebb[ParishRes %in% c("Tartu linn", "Tallinn") & (!(PoB %in% c("Tartu", "Harju")) | !(ParishBirth %in% c("", "Tartu linn", "Tallinn"))), ]
  
  # Adjust PRS
  ebb <- adjustPRSandEA(ebb)

  # Calculate mean values, SD and N for counties of birth (PoB)
  ebb_summary <- ebb[, .(mean(EY_adj2, na.rm = T), mean(PRS_adj2, na.rm = T),
                             sd(EY_adj2, na.rm = T), sd(PRS_adj2, na.rm = T)), 
                         by=c("PoB", "PoR")]
  ebb_N <- ebb[, .N, 
                   by=c("PoB", "PoR")]
  ebb_summary <- merge(ebb_summary, ebb_N)
  colnames(ebb_summary) <- c("PoB", "PoR", "mean_EY", "mean_PRS_adj2", "mean_PRS_adj3", "mean_PRS_adj4",
                                 "sd_EY", "sd_PRS_adj2", "sd_PRS_adj3", "sd_PRS_adj4", "N")
  
  # Calculate differences between migrants to Tallinn and Tartu for every county of birth
  ebb_diff <- data.frame()
  for(county in unique(ebb$PoB)){
    ebb_tmp <- ebb_summary[PoB == county,]
    if(nrow(ebb_tmp) < 2 | 1 %in% ebb_tmp$N){
      next
    }
    dif <- ebb_tmp[PoR == "Harju", .(mean_EY, mean_PRS_adj2)] - 
      ebb_tmp[PoR == "Tartu", .(mean_EY, mean_PRS_adj2)]
    dif_se <- sqrt( ebb_tmp[PoR == "Harju", .(sd_EY, sd_PRS_adj2)]^2/ c(ebb_tmp[PoR == "Harju", .(N)]) +
                      ebb_tmp[PoR == "Tartu", .(sd_EY, sd_PRS_adj2)]^2/c(ebb_tmp[PoR == "Tartu", .(N)]) )
    dif_n <- ebb_tmp[, min(N)]
    ebb_diff <- rbind(ebb_diff, c(county, as.numeric(dif), as.numeric(dif_se), dif_n))
  }
  colnames(ebb_diff) <- c("MNIMI", "mean_EY", "mean_PRS_adj2",
                          "se_EY", "se_PRS_adj2", "N")
  ebb_diff[, 2:10] <- lapply(ebb_diff[, 2:ncol(ebb_diff)], as.numeric)
  ebb_diff <- as.data.table(ebb_diff)
  
  # Mark whether the mean values are significantly different from zero after BH adjustment
  ebb_diff$p_EY <- p.adjust(2*pt(-abs(ebb_diff[, mean_EY/se_EY]), df=ebb_diff$N - 1), method = 'BH')
  ebb_diff$p_PRS_adj2 <- p.adjust(2*pt(-abs(ebb_diff[, mean_PRS_adj2/se_PRS_adj2]), df=ebb_diff$N - 1), method = 'BH')
  ebb_diff[p_EY < 0.05, sign_EY := "*"]
  ebb_diff[p_EY >= 0.05, sign_EY := ""]
  ebb_diff[p_PRS_adj2 <  0.05, sign_PRS_adj2 := "*"]
  ebb_diff[p_PRS_adj2 >= 0.05, sign_PRS_adj2 := ""]

  # Parameters for plotting
  nudge_x_stars <- c(0, -5, 0, -5, 0, -5, 10, -10, 0, 5, 10, -15, 0, 5, 0)*1000
  nudge_y_stars <- c(-10, -5, -10, -10, -10, -10, -10, 10, -15, -5, 2, -20, -3, -5, -15)*1000
  lim_min_prs <- min(ebb_diff$mean_PRS_adj2)
  lim_max_prs <- max(ebb_diff$mean_PRS_adj2)
  lim_min_ea <- min(ebb_diff$mean_EY)
  lim_max_ea <- max(ebb_diff$mean_EY)
  
  plot_map <- merge(maakonds, ebb_diff, by = "MNIMI", all.x = TRUE)
  
  # Map for PRS/sPRS
  pl_1 <- ggplot(data = plot_map)+
    geom_sf(aes(geometry = geometry, fill = mean_PRS_adj2))+
    geom_sf_text(aes(label = sign_PRS_adj2), nudge_x = nudge_x_stars,
                 nudge_y = nudge_y_stars, size = 5, color = "black") +
    scale_fill_gradient2(name="", limits=c(lim_min_prs-0.02, lim_max_prs), 
                         breaks=seq(-1, 1, by=round((lim_max_prs-lim_min_prs+0.02)/4, 1)), 
                         labels = comma)
  if(s){
    pl_1 <- pl_1 + ggtitle(bquote(paste(sPGS[" EA"])))
  } else{
    pl_1 <- pl_1 + ggtitle(bquote(paste(PGS[" EA"])))
  }
  
  # Map for EA
  pl_2 <- ggplot(data = plot_map)+
    geom_sf(aes(geometry = geometry, fill = mean_EY))+
    geom_sf_text(aes(label = sign_EY), nudge_x = nudge_x_stars,
                 nudge_y = nudge_y_stars, size = 5, color = "black") +
    ggtitle("Years of Education") +
    scale_fill_gradient2(name="", limits=c(lim_min_ea-0.02, lim_max_ea), 
                         breaks=seq(-1, 1, by=round((lim_max_ea-lim_min_ea+0.02)/4, 1)), 
                         labels = comma)
  
  plots <- list(pl_1, pl_2)
  
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



estbb_filtered <- fread(opt$pheno)
ebb <- fread(opt$geo)
ebb <- ebb[, c("Person skood", "PersonLocation birthParishName", "PersonLocation residencyParishName")]
colnames(ebb) <- c("skood", "ParishBirth", "ParishRes")
ebb <- merge(estbb_filtered, ebb, by="skood")
ebb[, Age2 := Age^2]
ebb[, SxA := Age*Sex]

ebb2 <- fread(opt$lastedu)
ebb2 <- ebb2[, c("Person skood", "PersonPortrait lastEducation code")]
colnames(ebb2) <- c("skood", "EA")
ebb <- merge(ebb, ebb2, by="skood")
ebb[, EduYears := ransformEAtoEduYears(EA)]
ebb[EA<6, EA := 0]
ebb[EA>=6, EA := 1]

pca_est <- fread(opt$pca_est)
pca_est <- pca_est[,c("IID", paste0("PC", 1:100))]
colnames(pca_est)[1] <- "vkood"
pca_rus <- fread(opt$pca_rus)
pca_rus <- pca_rus[,c("IID", paste0("PC", 1:100))]
colnames(pca_rus)[1] <- "vkood"
pca <- rbind(pca_est, pca_rus)
ebb <- merge(ebb, pca, by="vkood")

ebb_sib <- fread(opt$sib_prs)
colnames(ebb_sib) <- c("vkood", "PRS")
ebb_sib <- merge(ebb, ebb_sib)

prs <- fread(opt$pop_prs)
colnames(prs) <- c("vkood", "PRS")
ebb <- merge(ebb, prs)

maakonds <-  st_read(opt$map)
maakonds$MNIMI <- gsub(x=maakonds$MNIMI, pattern = " maakond", replacement = "")

non_rel <- fread(opt$nonrel)


# Make PDFs with plots for the subgroups and population/sibship PRS

pdf(paste0(opt$out, "/MapsTT_E_R.pdf"), width=7*2/3, height=2.5)
plots <- plotMaps(ebb[Nat=="Eestlane",], maakonds = maakonds)
print(grid.arrange(grobs = plots, layout_matrix = rbind(c(1, 2))))
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.518, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

plots <- plotMaps(ebb[Nat=="Venelane",], maakonds = maakonds)
print(grid.arrange(grobs = plots, layout_matrix = rbind(c(1, 2))))
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.518, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf(paste0(opt$out, "/MapsTT_E_nonrel.pdf"), width=7*2/3, height=2.5)
plots <- plotMaps(ebb[vkood %in% non_rel$vkood,], maakonds = maakonds)
print(grid.arrange(grobs = plots, layout_matrix = rbind(c(1, 2))))
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.518, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf(paste0(opt$out, "/MapsTT_E_sib.pdf"), width=7/3, height=2.5)
plots <- plotMaps(ebb_sib[Nat=="Eestlane",], maakonds = maakonds, s=T)
print(grid.arrange(grobs = plots[1], layout_matrix = rbind(c(1))))

dev.off()


pdf(paste0(opt$out, "/MapsTT_sex.pdf"), width=7*2/3, height=5)
plots_m <- plotMaps(ebb[Nat=="Eestlane" & Sex=="1",], maakonds = maakonds)
plots_f <- plotMaps(ebb[Nat=="Eestlane" & Sex=="2",], maakonds = maakonds)
print(grid.arrange(grobs = c(plots_m, plots_f), 
                   layout_matrix = rbind(c(1, 2),
                                         c(3, 4))))
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.518, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.518, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf(paste0(opt$out, "/MapsTT_yoa.pdf"), width=7*2/3, height=5)
plots_pre2016 <- plotMaps(ebb[Nat=="Eestlane" & YoA<=2016,], maakonds = maakonds)
plots_post2016 <- plotMaps(ebb[Nat=="Eestlane" & YoA>2016,], maakonds = maakonds)
print(grid.arrange(grobs = c(plots_pre2016, plots_post2016), 
                   layout_matrix = rbind(c(1, 2),
                                         c(3, 4))))
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.518, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.518, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf(paste0(opt$out, "/MapTT_age.pdf"), width=7*2/3, height=10)
plots_age1 <- plotMaps(ebb[Nat=="Eestlane" & Age>=18 & Age<=24,], maakonds = maakonds)
plots_age2 <- plotMaps(ebb[Nat=="Eestlane" & Age>=25 & Age<=48,], maakonds = maakonds)
plots_age3 <- plotMaps(ebb[Nat=="Eestlane" & Age>=49 & Age<=64,], maakonds = maakonds)
plots_age4 <- plotMaps(ebb[Nat=="Eestlane" & Age>=65,], maakonds = maakonds)
print(grid.arrange(grobs = c(plots_age1, plots_age2,
                             plots_age3, plots_age4), 
                   layout_matrix = rbind(c(1, 2),
                                         c(3, 4),
                                         c(5, 6),
                                         c(7, 8))))
grid.text("A", x = 0.018, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.518, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.518, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("E", x = 0.018, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("F", x = 0.518, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("G", x = 0.018, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("H", x = 0.518, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()
