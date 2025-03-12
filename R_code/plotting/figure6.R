#########################################################################
# Script making the PDFs with forest plots for PGS effects on migration #
# Figure 6 ##############################################################
#########################################################################


library(data.table)
library(ggplot2)
library(see)
library(grid)
library(gridExtra)

library(Rcpp)
library(optparse)


# Get PGS effects on the migration phenotype from a prepared file
getBetas <- function(lm_results, n_line = 2, effect_level){

  traits <- names(lm_results)
  b <- lapply(lm_results, function(x){return(x[n_line, c(1,2,4)])})
  b <- as.data.table(t(as.data.table(b)))
  colnames(b) <- c("b", "se", "p")
  b[, Trait := traits]
  b[, effect := effect_level]
  
}

# Get PGS effects on the migration phenotype from a prepared file
getAllBetas <- function(prefix, traits = "PRSs", model = "fixed"){
  
  if(model == "fixed"){
    
    # PRSs population pop effects
    unrellist <- readRDS(paste0(opt$prefix, "population2.0_PRS_effects_0pc_cities_", traits, ".rds"))
    unrel_b <- getBetas(lm_results = unrellist, n_line = 2, effect_level = "Population")

    m <- "LM"
    
  } else if(model == "LMM"){
    
    m <- "LMM"
    
  }
  
  # PRSs siblings pop effects fixed effects
  siblist <- readRDS(paste0(opt$prefix, "sibship_PRS2.0_pop_effects_0pc_", m, "_cities_", traits, ".rds"))
  pop_b <- getBetas(lm_results = siblist, n_line = 2, effect_level = "Population\nsibs sample")
  
  # PRSs siblings within-between family effects fixed effects
  siblist <- readRDS(paste0(opt$prefix, "sibship_PRS2.0_fam_effects_0pc_", m, "_cities_", traits, ".rds"))
  between_b <- getBetas(lm_results = siblist, n_line = 2, effect_level = "Between-sibship")
  within_b <- getBetas(lm_results = siblist, n_line = 3, effect_level = "Within-sibship")
  
  
  if(model == "fixed"){
    
    b <- rbind(unrel_b, pop_b, between_b, within_b)
    b$effect <- factor(b$effect, levels = rev(c("Population", "Population\nsibs sample", "Within-sibship", "Between-sibship")))
    
  } else if(model == "LMM"){
    
    b <- rbind(pop_b, between_b, within_b)
    b$effect <- factor(b$effect, levels = rev(c("Population\nsibs sample", "Within-sibship", "Between-sibship")))
    
  }

  if(traits == "PRSs_adjEAadj"){
    b[, Trait := gsub("_adj2_adjPGSea", "", Trait)]
    ea4_row <- data.table(b = NA_real_, se = NA_real_, p = NA_real_, Trait = "EA4", effect = "Population\nsibs sample")
    b <- rbind(b, ea4_row)
  } else {
    if(traits == "PRSs"){
      b[, Trait := gsub("_adj2", "", Trait)]
    } else if(traits == "PCs"){
      b[, Trait := gsub("PC", "", Trait)]
    } 
    
    b_list <- b[effect == "Within-sibship",]$Trait[order(b[effect == "Within-sibship",]$b)]
    b$Trait <- factor(b$Trait, levels = b_list)
    b_list_sign <- b[effect == "Within-sibship" & p < 0.05/169, Trait]
    b[Trait %in% b_list_sign, sign := T]
    b[!Trait %in% b_list_sign, sign := F]
    
  }
  
  

  return(b)
    
}




option_list = list(
  make_option(c("-", "--input"), type="character", default=NULL, 
              help="directory with input effect files", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output directory", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)





# Forest plots for PRSs with effect estimates from fixed effects model

# Effects of PRSs adjusted for demography and PCs
b_PRSs_fixed <- getAllBetas(prefix = opt$input, traits = "PRSs", model = "fixed")
pl1 <- ggplot(b_PRSs_fixed[sign == T,],
              aes(x = exp(b), y = Trait, color = effect)) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = exp(b + 1.96*se), xmin = exp(b - 1.96*se)), 
                 position = position_dodge(width = 0.7), linewidth = .5, height = 0) +
  geom_point(position = position_dodge(width = 0.7), size = 1.5, shape=19) + # , color = "orange"
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        legend.position = "none") +
  scale_y_discrete(labels=rev(c("EA4", "University degree", "A/AS levels", "Age completed\nfull time education",
                                "O levels/GCSEs", "Fluid intelligence", "Household income", "Other proffecional\nqualification",
                                "Cheese intake", "Time watching TV", "Time outdoors", "No qualification"))) +
  xlim(c(0.7, 1.5)) +
  labs(y = "PGS", x = "Odds Ratio") +
  guides(colour = guide_legend(reverse=T)) + scale_color_okabeito(order = c(1, 5, 7, 9))

# Effects of PRSs adjusted for demography and PCs
b_PRSs_EAadj_fixed <- getAllBetas(prefix = opt$input, traits = "PRSs_adjEAadj", model = "fixed")
b_PRSs_EAadj_fixed$Trait <- factor(b_PRSs_EAadj_fixed$Trait, levels = levels(b_PRSs_fixed$Trait))
pl2 <- ggplot(b_PRSs_EAadj_fixed[Trait %in% b_PRSs_fixed[sign == T, Trait], ],
              aes(x = exp(b), y = Trait, color = effect, order = effect)) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = exp(b + 1.96*se), xmin = exp(b - 1.96*se)), 
                 position = position_dodge(width = 0.7), linewidth = .5, height = 0) +
  geom_point(position = position_dodge(width = 0.7), size = 1.5, shape=19) + # , color = "orange"
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        legend.spacing.y = unit(0, "pt")) +
  scale_y_discrete(labels= rep("     ", 13)) +
  xlim(c(0.7, 1.5)) +
  labs(y = NULL, x = "Odds Ratio", color = "Effects") +
  guides(colour = guide_legend(reverse=T)) + scale_color_okabeito(order = c(1, 5, 7, 9))



plot_list <- list(pl1, pl2)

pdf(paste0(opt$out, "Figure6.2_fixed.pdf"), width=7, height = 5)
print(
  grid.arrange(
    grobs = plot_list,
    widths = c(1, 1.13),
    layout_matrix = matrix(c(1, 2), nrow = 1, byrow = TRUE)
  )
)
grid.text("A", x = 0.02, y = 0.97, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.49, y = 0.97, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()






# Forest plots for PRSs with effect estimates from mixed effects model

# Effects of PRSs adjusted for demography and PCs
b_PRSs_LMM <- getAllBetas(prefix = opt$input, traits = "PRSs", model = "LMM")
pl1 <- ggplot(b_PRSs_LMM[sign == T,],
              aes(x = exp(b), y = Trait, color = effect)) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = exp(b + 1.96*se), xmin = exp(b - 1.96*se)), 
                 position = position_dodge(width = 0.7), linewidth = .5, height = 0) +
  geom_point(position = position_dodge(width = 0.7), size = 1.5, shape=19) + # , color = "orange"
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        legend.position = "none") +
  scale_y_discrete(labels=rev(c("EA4", "University degree", "A/AS levels", "Age completed\nfull time education",
                                "O levels/GCSEs", "Fluid intelligence", "Other proffecional\nqualification", "Household income",
                                "Cheese intake", "Cereal type - Muesli", "Time watching TV", "Time outdoors", "No qualification"))) +
  xlim(c(0.67, 1.65)) +
  labs(y = "PGS", x = "Odds Ratio") +
  guides(colour = guide_legend(reverse=T)) + scale_color_okabeito(order = c(1, 5, 7, 9))

# Effects of PRSs adjusted for demography and PCs
b_PRSs_EAadj_LMM <- getAllBetas(prefix = opt$input, traits = "PRSs_adjEAadj", model = "LMM")
b_PRSs_EAadj_LMM$Trait <- factor(b_PRSs_EAadj_LMM$Trait, levels = levels(b_PRSs_LMM$Trait))
pl2 <- ggplot(b_PRSs_EAadj_LMM[Trait %in% b_PRSs_LMM[sign == T, Trait], ],
              aes(x = exp(b), y = Trait, color = effect, order = effect)) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = exp(b + 1.96*se), xmin = exp(b - 1.96*se)), 
                 position = position_dodge(width = 0.7), linewidth = .5, height = 0) +
  geom_point(position = position_dodge(width = 0.7), size = 1.5, shape=19) + # , color = "orange"
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        legend.spacing.y = unit(0, "pt")) +
  scale_y_discrete(labels= rep("     ", 13)) +
  xlim(c(0.67, 1.65)) +
  labs(y = NULL, x = "Odds Ratio", color = "Effects") +
  guides(colour = guide_legend(reverse=T)) + scale_color_okabeito(order = c(1, 5, 7, 9))



plot_list <- list(pl1, pl2)

pdf(paste0(opt$out, "Figure6.2_LMM.pdf"), width=7, height = 4)
print(
  grid.arrange(
    grobs = plot_list,
    widths = c(1, 1.13),
    layout_matrix = matrix(c(1, 2), nrow = 1, byrow = TRUE)
  )
)
grid.text("A", x = 0.02, y = 0.97, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.49, y = 0.97, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()





# Supplementary figures with all the PRSs
pl1 <- ggplot(b_PRSs_LMM[Trait %in% b_PRSs_LMM[effect == "Population\nsibs sample" & p < 0.05/169, Trait], ],
              aes(x = exp(b), y = Trait, color = effect)) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = exp(b + 1.96*se), xmin = exp(b - 1.96*se)), 
                 position = position_dodge(width = 0.7), linewidth = .5, height = 0) +
  geom_point(position = position_dodge(width = 0.7), size = 1.5, shape=19) + # , color = "orange"
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        legend.spacing.y = unit(0, "pt")) +
  labs(y = "PGS", x = "Odds Ratio") +
  guides(colour = guide_legend(reverse=T)) + scale_color_okabeito(order = c(1, 5, 7, 9))

pl2 <- ggplot(b_PRSs_fixed[Trait %in% b_PRSs_fixed[effect == "Population" & p < 0.05/169, Trait], ],
              aes(x = exp(b), y = Trait, color = effect)) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = exp(b + 1.96*se), xmin = exp(b - 1.96*se)), 
                 position = position_dodge(width = 0.7), linewidth = .5, height = 0) +
  geom_point(position = position_dodge(width = 0.7), size = 1.5, shape=19) + # , color = "orange"
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        legend.spacing.y = unit(0, "pt")) +
  labs(y = "PGS", x = "Odds Ratio") +
  guides(colour = guide_legend(reverse=T)) + scale_color_okabeito(order = c(1, 5, 7, 9))

pdf(paste0(opt$out, "Figure6.SM_PRS.pdf"), width=7, height = 10)
print(pl1)
print(pl2)
dev.off()




# Supplementary figures with PCs
b_PCs_fixed <- getAllBetas(prefix = opt$input, traits = "PCs", model = "fixed")
b_PCs_fixed$Trait <- factor(b_PCs_fixed$Trait, levels = 100:1)
b_PCs_LMM <- getAllBetas(prefix = opt$input, traits = "PCs", model = "LMM")
b_PCs_LMM$Trait <- factor(b_PCs_LMM$Trait, levels = 100:1)

pl1.1 <- ggplot(b_PCs_LMM[Trait %in% 1:50, ],
              aes(x = exp(b), y = Trait, color = effect)) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = exp(b + 1.96*se), xmin = exp(b - 1.96*se)), 
                 position = position_dodge(width = 0.7), linewidth = .5, height = 0) +
  geom_point(position = position_dodge(width = 0.7), size = 1.5, shape=19) + # , color = "orange"
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        legend.position = "none") +
  labs(y = "PC", x = "Odds Ratio") +
  xlim(c(0.67, 1.41)) +
  guides(colour = guide_legend(reverse=T)) + scale_color_okabeito(order = c(1, 5, 7, 9))

pl1.2 <- ggplot(b_PCs_LMM[Trait %in% 51:100, ],
                aes(x = exp(b), y = Trait, color = effect)) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = exp(b + 1.96*se), xmin = exp(b - 1.96*se)), 
                 position = position_dodge(width = 0.7), linewidth = .5, height = 0) +
  geom_point(position = position_dodge(width = 0.7), size = 1.5, shape=19) + # , color = "orange"
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        legend.spacing.y = unit(0, "pt")) +
  labs(y = "PC", x = "Odds Ratio") +
  xlim(c(0.67, 1.41)) +
  guides(colour = guide_legend(reverse=T)) + scale_color_okabeito(order = c(1, 5, 7, 9))

pl2.1 <- ggplot(b_PCs_fixed[Trait %in% 1:50, ],
                aes(x = exp(b), y = Trait, color = effect)) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = exp(b + 1.96*se), xmin = exp(b - 1.96*se)), 
                 position = position_dodge(width = 0.7), linewidth = .5, height = 0) +
  geom_point(position = position_dodge(width = 0.7), size = 1.5, shape=19) + # , color = "orange"
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        legend.position = "none") +
  labs(y = "PC", x = "Odds Ratio") +
  xlim(c(0.71, 1.32)) +
  guides(colour = guide_legend(reverse=T)) + scale_color_okabeito(order = c(1, 5, 7, 9))

pl2.2 <- ggplot(b_PCs_fixed[Trait %in% 51:100, ],
                aes(x = exp(b), y = Trait, color = effect)) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = exp(b + 1.96*se), xmin = exp(b - 1.96*se)), 
                 position = position_dodge(width = 0.7), linewidth = .5, height = 0) +
  geom_point(position = position_dodge(width = 0.7), size = 1.5, shape=19) + # , color = "orange"
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        legend.spacing.y = unit(0, "pt")) +
  labs(y = "PC", x = "Odds Ratio") +
  xlim(c(0.71, 1.32)) +
  guides(colour = guide_legend(reverse=T)) + scale_color_okabeito(order = c(1, 5, 7, 9))


plot_list1 <- list(pl1.1, pl1.2)
plot_list2 <- list(pl2.1, pl2.2)

pdf(paste0(opt$out, "Figure6.SM_PCs.pdf"), width=7, height = 10)

for(plot_list in list(plot_list1, plot_list2)){
  print(
    grid.arrange(
      grobs = plot_list,
      widths = c(1, 1.55),
      layout_matrix = matrix(c(1, 2), nrow = 1, byrow = TRUE)
    )
  )
}

dev.off()
