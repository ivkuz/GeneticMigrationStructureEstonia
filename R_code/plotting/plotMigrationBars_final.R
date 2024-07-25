#############################################################################################
# Script making the PDFs with barplots - mean PGS and EA in the migration groups ############
# Figure 4; Supplementary Figures 26-31, 38-63 ##############################################
#############################################################################################
# Also makes PDF with plot with mean PRS in the groups of born in or migrated to the cities #
#############################################################################################
# Also makes PDF with plot with PRS distributions in different migration groups for Q&A #####
#############################################################################################


library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(optparse)


# Transform EA catigories to the years of education
ransformEAtoEduYears <- function(EA_vector){
  
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


# Assign individuals to one of nine migration groups by PoB and PoR
# Harju and Tartu counties and ORE as meta-regions
annotatteMigrGroups <- function(ebb){
  
  ebb[, migr_group_counties := ""]
  
  ebb[!(PoB %in% c("Tartu", "Harju")) & 
        !(PoR %in% c("Tartu", "Harju")),
      migr_group_counties := "far2far"]
  
  ebb[!(PoB %in% c("Tartu", "Harju")) & 
        (PoR == "Tartu"),
      migr_group_counties := "far2Tartu"]
  
  ebb[!(PoB %in% c("Tartu", "Harju")) & 
        (PoR == "Harju"),
      migr_group_counties := "far2Harju"]
  
  ebb[(PoB == "Tartu") & 
        !(PoR %in% c("Tartu", "Harju")),
      migr_group_counties := "Tartu2far"]
  
  ebb[(PoB == "Tartu") & 
        (PoR == "Tartu"),
      migr_group_counties := "Tartu2Tartu"]
  
  ebb[(PoB == "Tartu") & 
        (PoR == "Harju"),
      migr_group_counties := "Tartu2Harju"]
  
  ebb[(PoB == "Harju") & 
        !(PoR %in% c("Tartu", "Harju")),
      migr_group_counties := "Harju2far"]
  
  ebb[(PoB == "Harju") & 
        (PoR == "Tartu"),
      migr_group_counties := "Harju2Tartu"]
  
  ebb[(PoB == "Harju") & 
        (PoR == "Harju"),
      migr_group_counties := "Harju2Harju"]
  
  return(ebb)
}


# Tallinn and Tartu linn and ORE as meta-regions
annotatteMigrGroups2 <- function(ebb){
  ebb[, migr_group_cities := ""]
  
  ebb[!(PoB %in% c("Tartu", "Harju")) & 
        !(PoR %in% c("Tartu", "Harju")),
      migr_group_cities := "far2far"]
  
  ebb[!(PoB %in% c("Tartu", "Harju")) & 
        ParishRes == "Tartu linn",
      migr_group_cities := "far2Tartu"]
  
  ebb[!(PoB %in% c("Tartu", "Harju")) & 
        ParishRes == "Tallinn",
      migr_group_cities := "far2Tallinn"]
  
  ebb[ParishBirth == "Tartu linn" & 
        !(PoR %in% c("Tartu", "Harju")),
      migr_group_cities := "Tartu2far"]
  
  ebb[ParishBirth == "Tartu linn" & 
        ParishRes == "Tartu linn",
      migr_group_cities := "Tartu2Tartu"]
  
  ebb[ParishBirth == "Tartu linn" & 
        ParishRes == "Tallinn",
      migr_group_cities := "Tartu2Tallinn"]
  
  ebb[ParishBirth == "Tallinn" & 
        !(PoR %in% c("Tartu", "Harju")),
      migr_group_cities := "Tallinn2far"]
  
  ebb[ParishBirth == "Tallinn" & 
        ParishRes == "Tartu linn",
      migr_group_cities := "Tallinn2Tartu"]
  
  ebb[ParishBirth == "Tallinn" & 
        ParishRes == "Tallinn",
      migr_group_cities := "Tallinn2Tallinn"]
  
  return(ebb)
}


# Rough assignment for Q&A
annotatteMigrGroups3 <- function(ebb){
  
  ebb[, migr_group_E_C := ""]
  
  ebb[!(PoB %in% c("Tartu", "Harju")) & 
        !(PoR %in% c("Tartu", "Harju")),
      migr_group_E_C := "far2far"]
  
  ebb[!(PoB %in% c("Tartu", "Harju")) & 
        ParishRes %in% c("Tartu linn", "Tallinn"),
      migr_group_E_C := "far2city"]
  
  ebb[ParishBirth %in% c("Tartu linn", "Tallinn") & 
        ParishRes  %in% c("Tartu linn", "Tallinn"),# &
      migr_group_E_C := "city2city"]
  
  return(ebb)
  
}


# Adjustment PRS and EA for demography, PCs (and EA)
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


# Make barplots
plotBars <- function(ebb, variable = "PRS_adj2", s=FALSE){
  
  # Adjust PRS and EA
  ebb <- adjustPRSandEA(ebb)
  
  # Rename the variable of interest
  colnames(ebb)[which(colnames(ebb)) == variable] <- "variable"
  
  # Calculate mean values, SD and N for the migration groups (Tallinn and Taru cities as meta-regions)
  ebb_cities_summary <- ebb[, .(mean(variable, na.rm = T), sd(variable, na.rm = T)), 
                            by=migr_group_cities]
  ebb_cities_N <- ebb[, .N, 
                      by=migr_group_cities]
  ebb_cities_summary <- merge(ebb_cities_summary, ebb_cities_N)
  colnames(ebb_cities_summary) <- c("migr_group", "mean_var", "mean_var", "N")
  ebb_cities_summary$migr_group <- factor(ebb_cities_summary$migr_group, levels = c("far2far", "far2Tartu", "far2Tallinn", 
                                                                                    "Tartu2far", "Tartu2Tartu", "Tartu2Tallinn", 
                                                                                    "Tallinn2far", "Tallinn2Tartu", "Tallinn2Tallinn"))
  
  # Calculate mean values, SD and N for the migration groups (Harju and Taru counties as meta-regions)
  ebb_counties_summary <- ebb[, .(mean(variable, na.rm = T), sd(variable, na.rm = T)), 
                              by=migr_group_counties]
  ebb_counties_N <- ebb[, .N, 
                        by=migr_group_counties]
  ebb_counties_summary <- merge(ebb_counties_summary, ebb_counties_N)
  colnames(ebb_counties_summary) <- c("migr_group", "mean_var", "mean_var", "N")
  ebb_counties_summary$migr_group <- factor(ebb_counties_summary$migr_group, levels = c("far2far", "far2Tartu", "far2Harju", 
                                                                                        "Tartu2far", "Tartu2Tartu", "Tartu2Harju", 
                                                                                        "Harju2far", "Harju2Tartu", "Harju2Harju"))
  

  # Parameters for the future plots
  y_min <- min(ebb_counties_summary[,mean_var-1.96*sd_var/sqrt(N)], ebb_cities_summary[,mean_var-1.96*sd_var/sqrt(N)])
  y_max <- max(ebb_counties_summary[,mean_var+1.96*sd_var/sqrt(N)], ebb_cities_summary[,mean_var+1.96*sd_var/sqrt(N)])
  y_lim = c(y_min, y_max)
  # Define the scale of the vertical axis
  if(y_max - y_min < 0.2){
    if(variable == "PRS_adj2"){
      sc <- 0.05
    } else if(variable == "PRS_adj4"){
      sc <- 0.04
    }
  } else{
    sc <- 0.2
  }
  
  # Plot with Harju and Tartu counties
  pl_counties <- ggplot(ebb_counties_summary, aes(x=migr_group, y=mean_var, fill=migr_group)) +
    geom_bar(stat = "identity") + 
    geom_errorbar(aes(ymin=(mean_PRS_adj2-1.96*sd_var/sqrt(N)), 
                      ymax=(mean_PRS_adj2+1.96*sd_var/sqrt(N))), 
                  width=.2, linewidth = 0.1) + 
    theme_bw() +
    theme(text = element_text(size=10),
          legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank(),
          panel.border = element_rect(size=0.1)
    ) +
    scale_y_continuous(breaks = seq(-1, 1, by=sc), limits = y_lim) +
    scale_fill_manual(values=c("#FFD3D5", "#FA4E5A", "#830717",
                               "#D8F1DC", "#40AA5F", "#005923",
                               "#D7E8FF", "#5898D6", "#00467D"))
  
  # Plot with Tallinn and Tartu cities
  pl_cities <- ggplot(ebb_cities_summary[migr_group != "",], aes(x=migr_group, y=mean_var, fill=migr_group)) +
    geom_bar(stat = "identity") + 
    geom_errorbar(aes(ymin=(mean_var-1.96*sd_var/sqrt(N)), 
                      ymax=(mean_var+1.96*sd_var/sqrt(N))), 
                  width=.2, linewidth = 0.1) +
    theme_bw() +
    theme(text = element_text(size=10),
          axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank(),
          panel.border = element_rect(size=0.1),
          legend.key.width = unit(0.5, 'cm'),
          legend.key.height = unit(0.5, 'cm')
    ) +
    scale_y_continuous(breaks = seq(-1, 1, by=sc), limits = y_lim) +
    scale_x_discrete(labels=c("Estonia -> Estonia",
                              "Estonia -> Tartu",
                              "Estonia -> Tallinn",
                              "Tartu -> Estonia",
                              "Tartu -> Tartu",
                              "Tartu -> Tallinn",
                              "Tallinn -> Estonia",
                              "Tallinn -> Tartu",
                              "Tallinn -> Tallinn")
    ) +
    scale_fill_manual(name="POB -> POR",
                      labels=c("ORE -> ORE",
                               "ORE -> Tartu",
                               "ORE -> Tallinn",
                               "Tartu -> ORE",
                               "Tartu -> Tartu",
                               "Tartu -> Tallinn",
                               "Tallinn -> ORE",
                               "Tallinn -> Tartu",
                               "Tallinn -> Tallinn"),
                      values=c("#FFD3D5", "#FA4E5A", "#830717",
                               "#D8F1DC", "#40AA5F", "#005923",
                               "#D7E8FF", "#5898D6", "#00467D"))
  
  # Add the respective y labels
  if(variable == "prs_adj2"){
    
    if(s){
      pl_counties <- pl_counties + ylab(bquote(paste(sPGS[" EA"])))
      pl_cities <- pl_cities + ylab(bquote(paste(sPGS[" EA"])))
    } else{
      pl_counties <- pl_counties + ylab(bquote(paste(PGS[" EA"])))
      pl_cities <- pl_cities + ylab(bquote(paste(PGS[" EA"])))
    }
    
  } else if(variable == "prs_adj3"){
    
    if(s){
      pl_counties <- pl_counties + ylab(bquote(paste(sPGS[" EA"], " adj. University Deg.")))
      pl_cities <- pl_cities + ylab(bquote(paste(sPGS[" EA"], " adj. University Deg.")))
    } else{
      pl_counties <- pl_counties + ylab(bquote(paste(PGS[" EA"], " adj. University Deg.")))
      pl_cities <- pl_cities + ylab(bquote(paste(PGS[" EA"], " adj. University Deg.")))
    }
  
  } else if(var == "prs_adj4"){
  
    if(s){
      pl_counties <- pl_counties + ylab(bquote(paste(sPGS[" EA"], " adj. Years of Edu.")))
      pl_cities <- pl_cities + ylab(bquote(paste(sPGS[" EA"], " adj. Years of Edu.")))
    } else{
      pl_counties <- pl_counties + ylab(bquote(paste(PGS[" EA"], " adj. Years of Edu.")))
      pl_cities <- pl_cities + ylab(bquote(paste(PGS[" EA"], " adj. Years of Edu.")))
    }
    
  } else if(var == "EY_adj2"){
    
    pl_counties <- pl_counties + ylab("Years of Education")
    pl_cities <- pl_cities + ylab("Years of Education")
    
  } else if(var == "EA"){
    
    pl_counties <- pl_counties + ylab("University Degree")
    pl_cities <- pl_cities + ylab("University Degree")
    
  }

  
  plot_list <- list(pl_prs2_counties, pl_prs2_cities)
  
  return(plot_list)
  
  
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

# Assign individuals to one of nine migration groups by PoB and PoR
ebb <- annotatteMigrGroups(ebb)
ebb <- annotatteMigrGroups2(ebb)


############################################################################
# Make PDFs with plots for the subgroups and population/sibship PRS and EA #
############################################################################

# PGS and sPGS adjusted for the demographic covariates and 100 PCs

pdf("figures/MigrationBars_E_R.pdf", width=7, height = 2.8)

plot_list <- plotBars(ebb[Nat=="Eestlane",], variable = "prs_adj2", s=F)
print(
  grid.arrange(
    grobs = plot_list,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

plot_list <- plotBars(ebb[Nat=="Venelane",], variable = "prs_adj2", s=F)
print(
  grid.arrange(
    grobs = plot_list,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf("figures/MigrationBars_E_nonrel.pdf", width=7, height = 2.8)

plot_list <- plotBars(ebb[vkood %in% non_rel$vkood,], variable = "prs_adj2", s=F)
print(
  grid.arrange(
    grobs = plot_list,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf("figures/MigrationBars_E_sib.pdf", width=7, height = 2.8)

plot_list <- plotBars(ebb_sib[Nat=="Eestlane",], variable = "prs_adj2", s=T)
print(
  grid.arrange(
    grobs = plot_list,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf("figures/MigrationBars_Sex_MF.pdf", width=7, height = 5.6)
plot_list1 <- plotBars(ebb[Nat=="Eestlane" & Sex==1,], variable = "prs_adj2", s=F) #Men
plot_list2 <- plotBars(ebb[Nat=="Eestlane" & Sex==2,], variable = "prs_adj2", s=F) #Women
plot_list <- c(plot_list1, plot_list2)

print(
  grid.arrange(
    grobs = plot_list,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4))
  )
)
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.418, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf("figures/MigrationBars_Age_young-old.pdf", width=7, height = 11)

plot_list1 <- plotBars(ebb[Nat=="Eestlane" & Age>=18 & Age<=24,], variable = "prs_adj2", s=F)
plot_list2 <- plotBars(ebb[Nat=="Eestlane" & Age>=25 & Age<=48,], variable = "prs_adj2", s=F)
plot_list3 <- plotBars(ebb[Nat=="Eestlane" & Age>=49 & Age<=64,], variable = "prs_adj2", s=F)
plot_list4 <- plotBars(ebb[Nat=="Eestlane" & Age>=65,], variable = "prs_adj2", s=F)
plot_list <- c(plot_list1, plot_list2, plot_list3, plot_list4)

print(
  grid.arrange(
    grobs = plot_list,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4),
                          c(5, 6),
                          c(7, 8))
  )
)
grid.text("A", x = 0.018, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.418, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("E", x = 0.018, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("F", x = 0.418, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("G", x = 0.018, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("H", x = 0.418, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf("figures/MigrationBars_YoA.pdf", width=7, height = 5.6)
plot_list1 <- plotBars(ebb[Nat=="Eestlane" & YoA<=2016,], variable = "prs_adj2", s=F)
plot_list2 <- plotBars(ebb[Nat=="Eestlane" & YoA>=2017,], variable = "prs_adj2", s=F)

plot_list <- c(plot_list1, plot_list2)

print(
  grid.arrange(
    grobs = plot_list,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4))
  )
)
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.418, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()




# EduYears adjusted for the demographic covariates and 100 PCs and binary EA

pdf("figures/MigrationBars_suppl_E_R.pdf", width=7, height = 2.8)

plot_list1 <- plotBars_suppl(ebb[Nat=="Eestlane",], variable = "EY_adj2")
plot_list2 <- plotBars_suppl(ebb[Nat=="Eestlane",], variable = "EA")
print(
  grid.arrange(
    grobs = plot_list2,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
print(
  grid.arrange(
    grobs = plot_list2,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

plot_list1 <- plotBars_suppl(ebb[Nat=="Venelane",], variable = "EY_adj2")
plot_list2 <- plotBars_suppl(ebb[Nat=="Venelane",], variable = "EA")
print(
  grid.arrange(
    grobs = plot_list1,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
print(
  grid.arrange(
    grobs = plot_list2,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf("figures/MigrationBars_suppl_E_nonrel.pdf", width=7, height = 2.8)

plot_list1 <- plotBars_suppl(ebb[vkood %in% non_rel$vkood,], variable = "EY_adj2")
plot_list2 <- plotBars_suppl(ebb[vkood %in% non_rel$vkood,], variable = "EA")
print(
  grid.arrange(
    grobs = plot_list1,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
print(
  grid.arrange(
    grobs = plot_list2,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))


dev.off()


pdf("figures/MigrationBars_suppl_Sex_MF.pdf", width=7, height = 5.6)
plot_list1.1 <- plotBars_suppl(ebb[Nat=="Eestlane" & Sex==1,], variable = "EY_adj2") #Men
plot_list1.2 <- plotBars_suppl(ebb[Nat=="Eestlane" & Sex==1,], variable = "EA") #Men
plot_list2.1 <- plotBars_suppl(ebb[Nat=="Eestlane" & Sex==2,], variable = "EY_adj2") #Women
plot_list2.2 <- plotBars_suppl(ebb[Nat=="Eestlane" & Sex==2,], variable = "EA") #Women
plot_list_ey <- c(plot_list1.1, plot_list2.1)
plot_list_ea <- c(plot_list1.2, plot_list2.2)

print(
  grid.arrange(
    grobs = plot_list_ey,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4))
  )
)
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.418, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))

print(
  grid.arrange(
    grobs = plot_list_ea,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4))
  )
)
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.418, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf("figures/MigrationBars_suppl_Age_young-old.pdf", width=7, height = 11)

plot_list1.1 <- plotBars_suppl(ebb[Nat=="Eestlane" & Age>=18 & Age<=24,], variable = "EY_adj2")
plot_list2.1 <- plotBars_suppl(ebb[Nat=="Eestlane" & Age>=25 & Age<=48,], variable = "EY_adj2")
plot_list3.1 <- plotBars_suppl(ebb[Nat=="Eestlane" & Age>=49 & Age<=64,], variable = "EY_adj2")
plot_list4.1 <- plotBars_suppl(ebb[Nat=="Eestlane" & Age>=65,], variable = "EY_adj2")
plot_list1.2 <- plotBars_suppl(ebb[Nat=="Eestlane" & Age>=18 & Age<=24,], variable = "EA")
plot_list2.2 <- plotBars_suppl(ebb[Nat=="Eestlane" & Age>=25 & Age<=48,], variable = "EA")
plot_list3.2 <- plotBars_suppl(ebb[Nat=="Eestlane" & Age>=49 & Age<=64,], variable = "EA")
plot_list4.2 <- plotBars_suppl(ebb[Nat=="Eestlane" & Age>=65,], variable = "EA")
plot_list_ey <- c(plot_list1.1, plot_list2.1, plot_list3.1, plot_list4.1)
plot_list_ea <- c(plot_list1.2, plot_list2.2, plot_list3.2, plot_list4.2)

print(
  grid.arrange(
    grobs = plot_list_ey,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4),
                          c(5, 6),
                          c(7, 8))
  )
)
grid.text("A", x = 0.018, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.418, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("E", x = 0.018, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("F", x = 0.418, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("G", x = 0.018, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("H", x = 0.418, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))

print(
  grid.arrange(
    grobs = plot_list_ea,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4),
                          c(5, 6),
                          c(7, 8))
  )
)
grid.text("A", x = 0.018, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.418, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("E", x = 0.018, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("F", x = 0.418, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("G", x = 0.018, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("H", x = 0.418, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf("figures/MigrationBars_suppl_YoA.pdf", width=7, height = 5.6)
plot_list1.1 <- plotBars_suppl(ebb[Nat=="Eestlane" & YoA<=2016,], variable = "EY_adj2")
plot_list2.1 <- plotBars_suppl(ebb[Nat=="Eestlane" & YoA>=2017,], variable = "EY_adj2")
plot_list1.2 <- plotBars_suppl(ebb[Nat=="Eestlane" & YoA<=2016,], variable = "EA")
plot_list2.2 <- plotBars_suppl(ebb[Nat=="Eestlane" & YoA>=2017,], variable = "EA")
plot_list_ey <- c(plot_list1.1, plot_list2.1)
plot_list_ea <- c(plot_list1.2, plot_list2.2)

print(
  grid.arrange(
    grobs = plot_list_ey,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4))
  )
)
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.418, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))

print(
  grid.arrange(
    grobs = plot_list_ea,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4))
  )
)
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.418, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()




# PGS adjusted for the demographic covariates, 100 PCs and EdyYears or binary EA

pdf("figures/MigrationBars_suppl2_E_R.pdf", width=7, height = 2.8)

plot_list1 <- plotBars_suppl2(ebb[Nat=="Eestlane",], variable = "prs_adj4", s=F)
plot_list2 <- plotBars_suppl2(ebb[Nat=="Eestlane",], variable = "prs_adj3", s=F)
print(
  grid.arrange(
    grobs = plot_list1,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
print(
  grid.arrange(
    grobs = plot_list2,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

plot_list1 <- plotBars_suppl2(ebb[Nat=="Venelane",], variable = "prs_adj4", s=F)
plot_list2 <- plotBars_suppl2(ebb[Nat=="Venelane",], variable = "prs_adj3", s=F)
print(
  grid.arrange(
    grobs = plot_list1,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
print(
  grid.arrange(
    grobs = plot_list2,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf("figures/MigrationBars_suppl2_E_nonrel.pdf", width=7, height = 2.8)

plot_list1 <- plotBars_suppl2(ebb[vkood %in% non_rel$vkood,], variable = "prs_adj4", s=F)
plot_list2 <- plotBars_suppl2(ebb[vkood %in% non_rel$vkood,], variable = "prs_adj3", s=F)
print(
  grid.arrange(
    grobs = plot_list1,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
print(
  grid.arrange(
    grobs = plot_list2,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf("figures/MigrationBars_suppl2_E_sib.pdf", width=7, height = 2.8)

plot_list1 <- plotBars_suppl2(ebb_sib[Nat=="Eestlane",], variable = "prs_adj4", s=T)
plot_list2 <- plotBars_suppl2(ebb_sib[Nat=="Eestlane",], variable = "prs_adj3", s=T)
print(
  grid.arrange(
    grobs = plot_list1,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
print(
  grid.arrange(
    grobs = plot_list2,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2))
  )
)
grid.text("A", x = 0.018, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.95, gp = gpar(fontsize=14, fontface = "bold"))


dev.off()


pdf("figures/MigrationBars_suppl2_Sex_MF.pdf", width=7, height = 5.6)
plot_list1.1 <- plotBars_suppl2(ebb[Nat=="Eestlane" & Sex==1,], variable = "prs_adj4", s=F) #Men
plot_list2.1 <- plotBars_suppl2(ebb[Nat=="Eestlane" & Sex==2,], variable = "prs_adj4", s=F) #Women
plot_list1.2 <- plotBars_suppl2(ebb[Nat=="Eestlane" & Sex==1,], variable = "prs_adj3", s=F) #Men
plot_list2.2 <- plotBars_suppl2(ebb[Nat=="Eestlane" & Sex==2,], variable = "prs_adj3", s=F) #Women
plot_list_ey <- c(plot_list1.1, plot_list2.1)
plot_list_ea <- c(plot_list1.2, plot_list2.2)

print(
  grid.arrange(
    grobs = plot_list_ey,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4))
  )
)
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.418, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))

print(
  grid.arrange(
    grobs = plot_list_ea,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4))
  )
)
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.418, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf("figures/MigrationBars_suppl2_Age_young-old.pdf", width=7, height = 11)

plot_list1.1 <- plotBars_suppl2(ebb[Nat=="Eestlane" & Age>=18 & Age<=24,], variable = "prs_adj4", s=F)
plot_list2.1 <- plotBars_suppl2(ebb[Nat=="Eestlane" & Age>=25 & Age<=48,], variable = "prs_adj4", s=F)
plot_list3.1 <- plotBars_suppl2(ebb[Nat=="Eestlane" & Age>=49 & Age<=64,], variable = "prs_adj4", s=F)
plot_list4.1 <- plotBars_suppl2(ebb[Nat=="Eestlane" & Age>=65,], variable = "prs_adj4", s=F)
plot_list1.2 <- plotBars_suppl2(ebb[Nat=="Eestlane" & Age>=18 & Age<=24,], variable = "prs_adj3", s=F)
plot_list2.2 <- plotBars_suppl2(ebb[Nat=="Eestlane" & Age>=25 & Age<=48,], variable = "prs_adj3", s=F)
plot_list3.2 <- plotBars_suppl2(ebb[Nat=="Eestlane" & Age>=49 & Age<=64,], variable = "prs_adj3", s=F)
plot_list4.2 <- plotBars_suppl2(ebb[Nat=="Eestlane" & Age>=65,], variable = "prs_adj3", s=F)
plot_list_ey <- c(plot_list1.1, plot_list2.1, plot_list3.1, plot_list4.1)
plot_list_ea <- c(plot_list1.2, plot_list2.2, plot_list3.2, plot_list4.2)

print(
  grid.arrange(
    grobs = plot_list_ey,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4),
                          c(5, 6),
                          c(7, 8))
  )
)
grid.text("A", x = 0.018, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.418, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("E", x = 0.018, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("F", x = 0.418, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("G", x = 0.018, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("H", x = 0.418, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))

print(
  grid.arrange(
    grobs = plot_list_ea,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4),
                          c(5, 6),
                          c(7, 8))
  )
)
grid.text("A", x = 0.018, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.9875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.418, y = 0.7375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("E", x = 0.018, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("F", x = 0.418, y = 0.4875, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("G", x = 0.018, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("H", x = 0.418, y = 0.2375, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


pdf("figures/MigrationBars_suppl2_YoA.pdf", width=7, height = 5.6)
plot_list1.1 <- plotBars_suppl2(ebb[Nat=="Eestlane" & YoA<=2016,], variable = "prs_adj4", s=F)
plot_list2.1 <- plotBars_suppl2(ebb[Nat=="Eestlane" & YoA>=2017,], variable = "prs_adj4", s=F)
plot_list1.2 <- plotBars_suppl2(ebb[Nat=="Eestlane" & YoA<=2016,], variable = "prs_adj3", s=F)
plot_list2.2 <- plotBars_suppl2(ebb[Nat=="Eestlane" & YoA>=2017,], variable = "prs_adj3", s=F)
plot_list_ey <- c(plot_list1.1, plot_list2.1)
plot_list_ea <- c(plot_list1.2, plot_list2.2)

print(
  grid.arrange(
    grobs = plot_list_ey,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4))
  )
)
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.418, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))

print(
  grid.arrange(
    grobs = plot_list_ea,
    widths = c(1, 1.5),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4))
  )
)
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.418, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("C", x = 0.018, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.418, y = 0.475, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()




#######################################################################################
# Make PDF with plot with mean PRS in the groups of born in or migrated to the cities #
#######################################################################################

# Keep only Estonians
ebb_est <- ebb[Nat=="Eestlane",]
ebb_est <- adjustPRSandEA(ebb_est)

# Make table for the groups by the age of birth
ore_city_comp <- data.frame()
pval <- c()
X_list <- list(c(1905, 1930), c(1931, 1940), c(1941, 1950), c(1951, 1960), c(1961, 1970), c(1971, 1980), c(1981, 1990), c(1991, 2005))

# For each group get summary data on PRS
for(X in X_list){
  
  filtered_ebb_city <- ebb_est[YoB >= X[1] & YoB <= X[2] & ParishBirth %in% c("Tartu linn", "Tallinn"), ]
  filtered_ebb_ore <- ebb_est[YoB >= X[1] & YoB <= X[2] & !(PoB %in% c("Tartu", "Harju")), ]
  filtered_ebb_ore_moved <- ebb_est[YoB >= X[1] & YoB <= X[2] & !(PoB %in% c("Tartu", "Harju"))& ParishRes %in% c("Tartu linn", "Tallinn"), ]
  filtered_ebb_ore_resid <- ebb_est[YoB >= X[1] & YoB <= X[2] & !(PoB %in% c("Tartu", "Harju")) & !(PoR %in% c("Tartu", "Harju")), ]

  ore_city_comp <- rbind(ore_city_comp, c("Tallinn/Tartu", paste(X[1], X[2], sep="-"), filtered_ebb_city[, mean(PRS_adj2)], filtered_ebb_city[, sd(PRS_adj2)]/sqrt(nrow(filtered_ebb_city))))
  ore_city_comp <- rbind(ore_city_comp, c("ORE", paste(X[1], X[2], sep="-"), filtered_ebb_ore[, mean(PRS_adj2)], filtered_ebb_ore[, sd(PRS_adj2)]/sqrt(nrow(filtered_ebb_ore))))
  ore_city_comp <- rbind(ore_city_comp, c("ORE_moved", paste(X[1], X[2], sep="-"), filtered_ebb_ore_moved[, mean(PRS_adj2)], filtered_ebb_ore_moved[, sd(PRS_adj2)]/sqrt(nrow(filtered_ebb_ore_moved))))
  ore_city_comp <- rbind(ore_city_comp, c("ORE_resid", paste(X[1], X[2], sep="-"), filtered_ebb_ore_resid[, mean(PRS_adj2)], filtered_ebb_ore_resid[, sd(PRS_adj2)]/sqrt(nrow(filtered_ebb_ore_resid))))
  
}

colnames(ore_city_comp) <- c("PoB", "YoB", "mean_PRS_adj2", "se_PRS_adj2")
ore_city_comp$mean_PRS_adj2 <- as.numeric(ore_city_comp$mean_PRS_adj2)
ore_city_comp$se_PRS_adj2 <- as.numeric(ore_city_comp$se_PRS_adj2)
ore_city_comp$PoB <- factor(ore_city_comp$PoB, levels = c("ORE", "ORE_resid", "ORE_moved", "Tallinn/Tartu"))
ore_city_comp <- as.data.table(ore_city_comp)

# Make plot for individuals born in ORE or the cities
pl_1 <- ggplot(ore_city_comp[PoB %in% c("ORE", "Tallinn/Tartu"),], aes(x = YoB, y = mean_PRS_adj2, fill = PoB)) + 
  geom_bar(stat = 'identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=(mean_PRS_adj2-1.96*se_PRS_adj2), 
                    ymax=(mean_PRS_adj2+1.96*se_PRS_adj2)),
                width=.2, position=position_dodge(0.9)) + 
  theme_bw() + ylab(bquote(paste(PGS[" EA"]))) +
  scale_fill_manual(name="POB", labels = c("ORE", "Tartu/Tallinn"), values=c("#D6878C", "#329494")) +
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_blank()) + ylim(c(-0.22, 0.33))

# Make plot for individuals residing in ORE or moved to the cities
pl_2 <- ggplot(ore_city_comp[PoB %in% c("ORE_resid", "ORE_moved"),], aes(x = YoB, y = mean_PRS_adj2, fill = PoB)) + 
  geom_bar(stat = 'identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=(mean_PRS_adj2-1.96*se_PRS_adj2), 
                    ymax=(mean_PRS_adj2+1.96*se_PRS_adj2)),
                width=.2, position=position_dodge(0.9)) + 
  theme_bw() + ylab(bquote(paste(PGS[" EA"]))) +
  scale_fill_manual(name="POB -> POR", labels = c("ORE -> ORE", "ORE -> Tallinn/Tartu"), values=c("#FFD3D5", "#AD3B43")) +
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(angle = 45, hjust = 1)) + ylim(c(-0.22, 0.33))

plots <- list(pl_1, pl_2)

# Make the final plot
pdf("figures/MigrationBars_refEstonia_YoB.pdf", width=7, height = 5.6)
print(
  grid.arrange(
    grobs = plots,
    widths = c(1, 0.065),
    heights = c(1, 1.2),
    layout_matrix = rbind(c(1, NA),
                          c(2, 2))
  )
)
grid.text("A", x = 0.018, y = 0.975, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.018, y = 0.520, gp = gpar(fontsize=14, fontface = "bold"))
dev.off()




###################################################################################
# Make PDF with plot with PRS distributions in different migration groups for Q&A #
###################################################################################


ebb_est <- annotatteMigrGroups3(ebb_est)
ebb_est <- ebb_est[migr_group_E_C %in% c("city2city", "far2far"), ]
ebb_est$migr_group_E_C <- factor(ebb_est$migr_group_E_C, levels = c("far2far", "city2city"))

pdf("figures/distribution_comp_QA.pdf", width=7, height = 5.6)
ggplot(ebb_est, aes(x=PRS_adj2, fill=migr_group_E_C)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  theme(text = element_text(size=15)) +
  xlab(label = "Polygenic score") +
  ylab(label = "Distribution density")+
  scale_fill_manual(name = "Born and living in", labels = c("other regions", "the cities"),
                    values=c("#FFD3D5", "#00467D"))
dev.off()



