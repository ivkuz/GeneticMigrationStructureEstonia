##############################################################
# Script making the PDF with barplots for sibling design PGS #
##############################################################


library(data.table)
library(pROC)
library(Rcpp)


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



# Read the sibling table
k <- fread(k, "~/tmp/siblings.tsv")

# table with future results
results <- data.frame()

# Work with EA4 PRS
prs <- "EA4"

# Calculate deviations from sibship means
names(k)[names(k) == paste0(prs, ".x")] <- "prs.x"
names(k)[names(k) == paste0(prs, ".y")] <- "prs.y"

k_m <- k[, .(mean(prs.x), mean(prs.y)), by=FID]
k_m[, prs := rowMeans(k_m[, c("V1", "V2")])]
k <- merge(k, k_m[, .(FID, prs)], by="FID")

k[, prs1 := prs.x - prs]
k[, prs2 := prs.y - prs]

ksib1 <- k[, .(ID1, PoB.x, PoR.x, prs1)]
ksib2 <- k[, .(ID2, PoB.y, PoR.y, prs2)]
colnames(ksib1) <- c("vkood", "PoB", "PoR", "prs")
colnames(ksib2) <- c("vkood", "PoB", "PoR", "prs")
ksib <- rbind(ksib1, ksib2)
ksib <- unique(ksib)


# The data table is ready
ebb <- copy(ksib)
# Assign migration groups
ebb <- annotatteMigrGroups(ebb)


# Calculate subbary statistics for the migration groups
ebb_counties_summary <- ebb[, .(mean(prs), sd(prs)), by=migr_group_counties]
ebb_counties_N <- ebb[, .N, by=migr_group_counties]
ebb_counties_summary <- merge(ebb_counties_summary, ebb_counties_N)
colnames(ebb_counties_summary) <- c("migr_group", "mean_PRS", "sd_PRS", "N")
ebb_counties_summary$migr_group <- factor(ebb_counties_summary$migr_group, levels = c("far2far", "far2Tartu", "far2Harju", 
                                                                                      "Tartu2far", "Tartu2Tartu", "Tartu2Harju", 
                                                                                      "Harju2far", "Harju2Tartu", "Harju2Harju"))

# Make the barplot
pl_prs2_counties <- ggplot(ebb_counties_summary, aes(x=migr_group, y=mean_PRS, fill=migr_group)) +
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin=(mean_PRS-1.96*sd_PRS/sqrt(N)), 
                    ymax=(mean_PRS+1.96*sd_PRS/sqrt(N))), 
                width=.2, linewidth = 0.1) + #ylim(y_lim) + #ylim(c(-0.30, 0.55)) + 
  theme_bw() +
  ylab(bquote(paste(PGS[" EA4"], " - ", PGS[" EA4 sibMean"]))) +
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_rect(linewidth=0.1),
        legend.key.width = unit(0.5, 'cm'),
        legend.key.height = unit(0.5, 'cm')
  ) +

  scale_x_discrete(labels=c("far2far",
                            "far2Tartu",
                            "far2Harju",
                            "Tartu2far",
                            "Tartu2Tartu",
                            "Tartu2Harju",
                            "Harju2far",
                            "Harju2Tartu",
                            "Harju2Harju")
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


# Save the plot
pdf("~/figures/MigrBars_sib_family.pdf", width=4, height = 2.8)
print(pl_prs2_counties)
dev.off()