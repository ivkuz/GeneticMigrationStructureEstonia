#################################################################
# Script prepares .rds file with Figure 4C, input for figure4.R #
#################################################################


library(data.table)
library(pROC)
library(Rcpp)
library(optparse)
library(ggplot2)

library(Rcpp)
library(optparse)



# recode PoB from codes to particular counties. If multiple -> NA
transformPoB <- function(PoB_vector, countycode){
  pob <- strsplit(as.character(PoB_vector), "(;| )+")
  pob_names <- c()
  for(i in pob){
    counties <- c()
    for(ccode in i){
      counties <- c(counties, countycode[code == ccode, county])
    }
    if(length(unique(counties)) == 1){
      pob_names <- c(pob_names, unique(counties))
    } else {
      pob_names <- c(pob_names, NA)
    }
  }
  return(pob_names)
}



# Combine PRS, PCA, Phenotypes
mergeData <- function(PRS, PCA, Phenotypes){
  merged_data <- merge(PRS, PCA, by="vkood")
  merged_data <- merge(merged_data, Phenotypes, by="vkood")
  return(merged_data)
}


# Define POB-POR groups on the county level
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


# Define POB-POR groups on the city level
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



option_list = list(
  make_option(c("-e", "--pca_est"), type="character", default=NULL,
              help="file with 100 PCs of estonians", metavar="character"),
  make_option(c("-r", "--pca_rus"), type="character", default=NULL,
              help="file with 100 PCs of russians", metavar="character"),
  make_option(c("-f", "--pheno"), type="character", default=NULL,
              help="file with phenotypes", metavar="character"),
  make_option(c("-s", "--score"), type="character", default=NULL,
              help="file with vkoods and PRS", metavar="character"),
  make_option(c("-k", "--kin"), type="character", default=NULL,
              help="file with kinship table", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="results.tsv",
              help="output directory", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



# Add PRS from EA4 excluding 23andMe and EstBB
PRS_tmp <- fread(opt$score)
colnames(PRS_tmp) <- c("vkood", "EA4")
PRS <- PRS_tmp
prs_codes <- "EA4"

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



# Select the pairs of siblings from KING kinship table
k <- fread(opt$kin)
k <- k[InfType=="FS" & Kinship > 0.177 & Kinship < 0.354, .(ID1, ID2)]
kvariables <- c("vkood", "YoB", "PoB", "PoR", "Nat", prs_codes, paste0("PC", 1:100))
k <- merge(k, ebb[Nat=="Eestlane", ..kvariables], by.x="ID1", by.y="vkood")
k <- merge(k, ebb[Nat=="Eestlane", ..kvariables], by.x="ID2", by.y="vkood")
# Keep only siblings born in the same county
k <- k[PoB.x==PoB.y,]

# find individuals with more than one sibling
mult_sib <- names(table(c(k$ID1, k$ID2)))[table(c(k$ID1, k$ID2))>1]

# create families of several siblings (make the same FID for all the pairs of a sibship)
k$FID <- 1:nrow(k)
fid_excl <- c()
for(i in 2:nrow(k)){
  # if somebody from the pair has more than one sibling
  if(k[i, ID1] %in% mult_sib | k[i, ID2] %in% mult_sib){
    # keep all the rows with them (and their sibling) from the rows above to a table
    k_tmp <- k[1:(i-1),][ID1 %in% k$ID1[i] | ID1 %in% k$ID2[i] | ID2 %in% k$ID1[i] | ID2 %in% k$ID2[i], ]
    # if there are some
    if(nrow(k_tmp)>0){
      # if some of that families have different FIDs
      if(length(unique(k_tmp$FID))>1){
        # this is useful when there are "linked sibships"
        # for example data.table(ID1=c(1,2,2), ID2=c(4,3,4), FID=1:3)
        # this could wrong for normal sibships in the samples would haven't been sorted
        # for example data.table(ID1=c(1,3,1,3,1,2), ID2=c(2,4,4,2,3,4), FID=1:6)
        # such family would be deleted completely
        # but KING orders the relatives in a good way and we don't have such cases
        print(k_tmp)
        fid_excl <- c(fid_excl, unique(k_tmp$FID))
      }
      k[i, FID := k_tmp[1, FID]]
    }
  }
}
k <- k[!(FID %in% fid_excl),]


# Keep sibships with the numbers of siblings which makes sense
kN <- k[, .(.N), by=FID]
kN <- kN[N %in% c(1, 3, 6, 10, 15, 21), ]
k <- k[FID %in% kN$FID,]


# Exclude sibships where not all individuals are siblings
tr <- k[FID %in% names(table(k[, FID]))[as.vector(table(k[, FID]))>1], ]
famsize <- data.table(n=c(3,4,5,6,7), edgeN=c(3,6,10,15,21))
for(f in unique(tr$FID)){
  tr_tmp <- tr[FID==f, ]
  expected_N <- famsize[edgeN == nrow(tr_tmp), n]
  fam_ids <- unique(c(tr_tmp$ID1, tr_tmp$ID2))
  real_N <- length(fam_ids)
  if(real_N != expected_N){
    print(fam_ids)
    k <- k[FID != f,]
  }
}


# table with future results
results <- data.frame()

prs <- "EA4"

names(k)[names(k) == paste0(prs, ".x")] <- "prs.x"
names(k)[names(k) == paste0(prs, ".y")] <- "prs.y"

# Get within-sibship PGSs
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


# merge sibling data with other phenotype data and annotate groups by POB and POR
ksib <- annotatteMigrGroups(ksib)

# Add geographic info on the parish level (to distinguish the cities)
ebb <- merge(ksib, ebb[, .(vkood, ParishBirth, ParishRes)])
ebb <- annotatteMigrGroups2(ebb)


# Calculate mean and sd values of PGS for migration groups
ebb_counties_summary <- ebb[, .(mean(prs), sd(prs)), by=migr_group_counties]
ebb_counties_N <- ebb[, .N, by=migr_group_counties]
ebb_counties_summary <- merge(ebb_counties_summary, ebb_counties_N)

colnames(ebb_counties_summary) <- c("migr_group", "mean_PRS", "sd_PRS", "N")

ebb_counties_summary$migr_group <- factor(ebb_counties_summary$migr_group, levels = c("far2far", "far2Tartu", "far2Harju", 
                                                                                      "Tartu2far", "Tartu2Tartu", "Tartu2Harju", 
                                                                                      "Harju2far", "Harju2Tartu", "Harju2Harju"))

# Make the plot with bars for sibling design
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
  scale_fill_manual(name="POB -> POR", #"Direction of\nmigration",
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


# Save as RDS file
saveRDS(pl_prs2_counties, paste0(opt$out, "pl_sib_bars_counties.rds"))

