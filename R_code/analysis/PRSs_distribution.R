################################################################################
# Script calculating the mean and sd for PGSs by region for all the subcohorts #
# And for PGSs adjusted for PGS for EA #########################################
################################################################################


Sys.setlocale("LC_CTYPE", "estonian")

my_packages <- c("Rcpp", "data.table", "pROC", "optparse")
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if(length(not_installed)) install.packages(not_installed)

library(data.table)
library(pROC)
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


# Adjustment PRS for PCA and normalization
adjustPRS <- function(tab, prs_codes){
  
  for(prs in prs_codes){
    
    colnames(tab)[which(colnames(tab)==prs)] <- "PRS"
    
    lm_prs2 <- lm(paste0("PRS ~ Sex + YoB + I(YoB^2) + I(Sex*YoB) + ", paste("PC", 1:100, sep = "", collapse = " + ")),
                  data = tab)
    tab$PRS_adj2 <- tab$PRS - predict(object = lm_prs2, newdata = tab)
    
    tab[, PRS := scale(PRS)]
    tab[, PRS_adj2 := scale(PRS_adj2)]

    colnames(tab)[which(colnames(tab)=="PRS")] <- prs
    colnames(tab)[which(colnames(tab)=="PRS_adj2")] <- paste0(prs, "_adj2")

  }

    return(tab)
}


# Calculating mean value ans SD for all the counties for all the subsets
calculateBasicStatistics <- function(ebb, prs_codes, only.all=FALSE){
  
  # table with future results
  results <- data.frame()
  
  for(nat in list("Eestlane", "Venelane")){
    
    ebb_tmp1 <- ebb[Nat %in% nat, ]
    
    if(nat =="Eestlane"){
      sex_list <- list("1", "2", c("1", "2")) 
    } else{
      sex_list <- list(c("1", "2"))
    }
    
    for(sex in sex_list){
      
      ebb_tmp2 <- ebb_tmp1[Sex %in% sex, ]
      if(length(sex) == 2){
        sex <- "All"
      }
      
      if(nat =="Eestlane" & sex == "All"){
        age_list <- list(c(18, 24), c(25, 48), c(49, 64), c(65, 107), c(18, 107)) 
      } else{
        age_list <- list(c(18, 107))
      }
      
      for(age in age_list){
        
        ebb_tmp3 <- ebb_tmp2[Age >= age[1] & Age <= age[2], ]
        print(age)
        if(age[1] == 18 & age[2] == 107){
          age <- "All"
        } else{
          age <- paste(age[1], age[2], sep="_")
        }
        
        if(nat =="Eestlane" & sex == "All" & age == "All"){
          yoa_list <- list(2001:2016, 2017:2021, 2001:2021) 
        } else{
          yoa_list <- list(2001:2021)
        }
        
        for(yoa in yoa_list){
          
          ebb_tmp4 <- ebb_tmp3[YoA %in% yoa, ]
          if(identical(yoa, 2001:2021)){
            yoa <- "All"
          } else{
            yoa <- paste(yoa[1], yoa[length(yoa)], sep="_")
          }

          if(only.all & !(nat=="Eestlane" & sex=="All" & age=="All" & yoa=="All")){
            next
          }
          
          for(prs in c(prs_codes, paste0(prs_codes, "_adj2"))){

            ebb_tmp5 <- ebb_tmp4
            names(ebb_tmp5)[names(ebb_tmp5) == prs] <- "prs"

            for(place in c("PoB", "PoR")){
              
              ebb_tmp6 <- ebb_tmp5
              names(ebb_tmp6)[names(ebb_tmp6) == place] <- "place"

              for(county in c("All", unique(ebb_tmp6[!is.na(place), place]))){
                
                if(county == "All"){
                  ebb_tmp7 <- ebb_tmp6[!is.na(place), ]
                } else{
                  ebb_tmp7 <- ebb_tmp6[place == county, ] 
                }
                
                m <- ebb_tmp7[, mean(prs, na.rm=TRUE)]
                sd <- ebb_tmp7[, sd(prs, na.rm=TRUE)]
                N <- nrow(ebb_tmp7[!is.na(prs),])
                
                output <- c(nat, sex, age, yoa, prs, place, county,
                            round(m, 8), round(sd, 8), N)
                results <- rbind(results, output)
                
              }
            }
          }
        }
      }
    }
  }
  
  colnames(results) <- c("Nat", "Sex", "Age", "YoA", "Variable", "Birth/Residence", "County", 
                         "mean", "sd", "N")
  
  return(results)
  
}


# Test if places of birth and residence can be explain more PC variance than separately
testDifferenceBirthRes <- function(ebb, prs_codes, only.all=FALSE){
  
  # table with future results
  results <- data.frame()
  
  for(nat in list("Eestlane", "Venelane")){
    
    ebb_tmp1 <- ebb[Nat %in% nat, ]
    
    if(nat =="Eestlane"){
      sex_list <- list("1", "2", c("1", "2")) 
    } else{
      sex_list <- list(c("1", "2"))
    }
    
    for(sex in sex_list){
      
      ebb_tmp2 <- ebb_tmp1[Sex %in% sex, ]
      if(length(sex) == 2){
        sex <- "All"
      }
      
      if(nat =="Eestlane" & sex == "All"){
        # age_list <- list(c(18, 24), c(25, 34), c(35, 64), c(65, 105), c(18, 105)) 
        age_list <- list(c(18, 24), c(25, 48), c(49, 64), c(65, 107), c(18, 107)) 
      } else{
        age_list <- list(c(18, 107))
      }
      
      for(age in age_list){
        
        ebb_tmp3 <- ebb_tmp2[Age >= age[1] & Age <= age[2], ]
        print(age)
        if(age[1] == 18 & age[2] == 107){
          age <- "All"
        } else{
          age <- paste(age[1], age[2], sep="_")
        }
        
        if(nat =="Eestlane" & sex == "All" & age == "All"){
          yoa_list <- list(2001:2016, 2017:2021, 2001:2021) 
        } else{
          yoa_list <- list(2001:2021)
        }
        
        for(yoa in yoa_list){
          
          ebb_tmp4 <- ebb_tmp3[YoA %in% yoa, ]
          if(identical(yoa, 2001:2021)){
            yoa <- "All"
          } else{
            yoa <- paste(yoa[1], yoa[length(yoa)], sep="_")
          }
          
          if(only.all & !(nat=="Eestlane" & sex=="All" & age=="All" & yoa=="All")){
            next
          }
          
          for(prs in c(prs_codes, paste0(prs_codes, "_adj2"))){
            # for(prs in c(paste0(prs_codes, "_adj3"))){
            
            ebb_tmp5 <- ebb_tmp4
            names(ebb_tmp5)[names(ebb_tmp5) == prs] <- "prs"
            # ebb_tmp4 <- ebb_tmp4[!is.na(prs), ]
            
            lm1 <- lm(paste0("prs ~ PoB"), data=ebb_tmp5)
            lm2 <- lm(paste0("prs ~ PoR"), data=ebb_tmp5)
            lm3 <- lm(paste0("prs ~ PoB + PoR"), data=ebb_tmp5)
            if(summary(lm1)$r.squared > summary(lm2)$r.squared){
              p <- anova(lm3, lm2, test="F")[["Pr(>F)"]][2]
            } else{
              p <- anova(lm3, lm1, test="F")[["Pr(>F)"]][2]
            }
            output <- c(nat, sex, age, yoa, prs, p)
            results <- rbind(results, output)
          }
        }
      }
    }
  }
  
  colnames(results) <- c("Nat", "Sex", "Age", "YoA", "Variable", "P")
  
  return(results)
  
}

# Adjustment PRS for another PRS and normalization
adjustPRS_forEAprs <- function(tab, prs_codes, adj = "6138_1_adj2"){
  
  tab[, adjEA := get(adj)]
  i=1
  for(prs in prs_codes){
    
    print(i)
    i <- i + 1
    
    if(prs == adj){
      
      tab[, PRS_adj2 := 0]
      
    }else{
      
      colnames(tab)[which(colnames(tab)==prs)] <- "PRS"
      
      lm_prs <- lm("PRS ~ adjEA", data = tab)
      tab$PRS_adj2 <- tab$PRS - predict(object = lm_prs, newdata = tab)
      tab[, PRS_adj2 := scale(PRS_adj2)]
      
      colnames(tab)[which(colnames(tab)=="PRS")] <- prs
      
    }
    
    colnames(tab)[which(colnames(tab)=="PRS_adj2")] <- paste0(prs, "_adjPGSea")
    
  }
  
  tab[, adjEA := NULL]
  
  return(tab)
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
  make_option(c("-u", "--nonrel"), type="character", default="non_relatives.tsv", 
              help="list of non-related individuals [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="~/tmp/", 
              help="output directory name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# Making file with PRS
PRS <- NA
# Trait codes from Pan-UKB (Supplementary Table 1)
prs_codes <- c(1920, 1930, 1940, 1950, 1960, 1980, 1990, 2000, 2010, 2020, 2040,
            2090, 2100, 4631, "6138_1", "6138_100", "6138_2", "6138_3", "6138_6",
            6145, 20016, 20023, 20127, 2050, 2070, 738, 845, 30610, 30700, 30720,
            30750, 30760, 30770, "20002_1226", "20003_1140874420", "20004_1479", 1697, 30040, 30070, 30080,
            30100, 3147, 4125, 5096, "eGFR",
            30600, 30620, 30630, 30640, 30650, 30680, 30690, 30710, 30730, 30740, 30780, 30810, 30830,
            30860, 30870, 30880, 30890, 1210, "1428_0", "1428_1", 1468, 1538, 1747, "20002_1138",
            "20003_1140861998", "20003_1140874744", "20003_1140881856", "20003_1140923346",
            "20003_1141194794", "20004_1320", "20004_1478", "20004_1480", "20004_1483", "20107",
            "20110", "20111", "20116", "20544", "2129", "2257", "2463", "2654", "3606", "41200_K633",
            "41200_V544", "41246", "5556", 6139, 6144, 6146, 6147, "6149_100", "6149_3", "6149_6",
            6152, 6159, "6160_1", "6160_100", "6160_2", "6160_4", 6162, 1050, 1070, 1110, 1160,
            1180, 1220, 1239, 1249, 1329, 135, 137, 1408, 1458, 1478, 1488, 1498, 1518, 1528, 1558,
            1717, 1757, 20512, 20519, 21001, 2139, 2178, 2217, 2237, 2267, 23104, 30060, 30130, 30180,
            30200, 30210, 399, 4194, 4195, 50, 924, "AG", "LDLC", "PP", "K22", "K40", "K57", "M16",
            "N20", "N81", 300, 317, 327, 471, "550.4", 574, "714.1", 726, "735.3", "743.1", 835,
            960, "alpha1-antagonist_alphablocker", "dexamethasone_aceticacid_neomycin",
            "6138_1_grmAdj", "EA4") # additional PRSs

############################################################
# Use these codes for within-sibship GWAS polygenic scores #
############################################################
# prs_codes <- 4813:4860                                   #
############################################################

# Read and aggregare PRS files
for(code in prs_codes){
 PRS_tmp <- fread(paste0(opt$score, "/", code, "/", code, ".chrALL.sscore"))
 colnames(PRS_tmp) <- c("vkood", code)
 if(is.na(PRS)){
   PRS <- PRS_tmp
 } else{
   PRS <- merge(PRS, PRS_tmp, by = "vkood")
 }
}

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

# Adjustment PRS for PCA and normalization
ebb_est <- adjustPRS(ebb[Nat=="Eestlane", ], prs_codes=prs_codes)
ebb_rus <- adjustPRS(ebb[Nat=="Venelane", ], prs_codes=prs_codes)
ebb <- rbind(ebb_est, ebb_rus)


# Save ebb data file with PGSs, PCs amd phenotypes
write.table(ebb, paste0(opt$out, "PRSs_adj.tsv"), row.names = F, quote = F, sep = "\t")

# calculate and save basic result table
results <- calculateBasicStatistics(ebb = ebb, prs_codes=prs_codes, only.all=FALSE)
write.table(results, paste0(opt$out,"/results_all_PRS_ages2_yoa.tsv"), row.names = F, quote = F, sep = "\t")

# calculate and save the table with statistics on PoB and PoR differences
results_dif <- testDifferenceBirthRes(ebb = ebb, prs_codes=prs_codes, only.all=FALSE)
write.table(results_dif, paste0(opt$out, "results_all_PRS_modelDif_ages2_yoa.tsv"), row.names = F, quote = F, sep = "\t")



# Unrelated individuals
non_rel <- fread(opt$nonrel)

# Calculate and save basic result table
results <- calculateBasicStatistics(ebb = ebb[vkood %in% non_rel$vkood,], prs_codes=prs_codes, only.all=TRUE)
write.table(results, paste0(opt$out,"/results_all_PRS_ages2_yoa_nonrel.tsv"), row.names = F, quote = F, sep = "\t")

# Calculate and save the table with statistics on PoB and PoR differences
results_dif <- testDifferenceBirthRes(ebb = ebb[vkood %in% non_rel$vkood,], prs_codes=prs_codes, only.all=TRUE)
write.table(results_dif, paste0(opt$out, "results_all_PRS_modelDif_ages2_yoa_nonrel.tsv"), row.names = F, quote = F, sep = "\t")



######################################
# Correlations between adjusted PGSs #
######################################

cor_cols <- which( colnames(ebb) %in% paste0(prs_codes, "_adj2") )
cor_matrix <- cor(ebb[vkood %in% non_rel$vkood, ..cor_cols])
write.table(cor_matrix, paste0(opt$out, "PRS_adj2_nonrel_cor_matrix.tsv"), col.names = T, row.names = T, quote = F, sep = "\t") 



################################################
# Var(county) for PGSs adjusted for PGS for EA #
################################################

# Adjustment PRS for PCA and normalization
ebb_est <- adjustPRS(ebb[Nat=="Eestlane", ], prs_codes="EA4")
ebb_rus <- adjustPRS(ebb[Nat=="Venelane", ], prs_codes="EA4")
ebb <- rbind(ebb_est, ebb_rus)

# Adjustment PRS for PGS for EA4
ebb_est <- adjustPRS_forEAprs(ebb[Nat=="Eestlane", ], prs_codes=paste0(prs_codes, "_adj2"), adj = "EA4_adj2")
ebb_rus <- adjustPRS_forEAprs(ebb[Nat=="Venelane", ], prs_codes=paste0(prs_codes, "_adj2"), adj = "EA4_adj2")
ebb <- rbind(ebb_est, ebb_rus)

# calculate and save basic result table
results <- calculateBasicStatistics(ebb = ebb, prs_codes=paste0(prs_codes, "_adj2_adjPGSea"), only.all=TRUE)
write.table(results, paste0(opt$out,"/results_all_PRS_pgsEA4adj_ages2_yoa.tsv"), row.names = F, quote = F, sep = "\t")

results_dif <- testDifferenceBirthRes(ebb = ebb, prs_codes=paste0(prs_codes, "_adj2_adjPGSea"), only.all=TRUE)
write.table(results_dif, paste0(opt$out, "/results_all_PRS_pgsEA4adj_modelDif_ages2_yoa.tsv"), row.names = F, quote = F, sep = "\t")


# Remove the PGSs adjusted for PGS for EA
cols_to_remove <- paste0(prs_codes, "_adj2_adjPGSea")
dt[, (cols_to_remove) := NULL]

# Adjustment PRS for PGSea
ebb_est <- adjustPRS_forEAprs(ebb[Nat=="Eestlane", ], prs_codes=paste0(prs_codes, "_adj2"), adj = "6138_1_adj2")
ebb_rus <- adjustPRS_forEAprs(ebb[Nat=="Venelane", ], prs_codes=paste0(prs_codes, "_adj2"), adj = "6138_1_adj2")
ebb <- rbind(ebb_est, ebb_rus)

# calculate and save basic result table
results <- calculateBasicStatistics(ebb = ebb, prs_codes=paste0(prs_codes, "_adj2_adjPGSea"), only.all=FALSE)
write.table(results, paste0(opt$out,"/results_all_PRS_pgsEAadj_ages2_yoa.tsv"), row.names = F, quote = F, sep = "\t")

results_dif <- testDifferenceBirthRes(ebb = ebb, prs_codes=paste0(prs_codes, "_adj2_adjPGSea"), only.all=FALSE)
write.table(results_dif, paste0(opt$out, "results_all_PRS_pgsEAadj_modelDif_ages2_yoa.tsv"), row.names = F, quote = F, sep = "\t")





#################################
# Cumulative adjustment for PCs #
#################################

# Keep only normalized UKBB EA PRS
ebb <- ebb[, c("vkood", "6138_1")]
colnames(ebb)[2] <- "PRS"

# Combine PRS, PCA, Phenotypes
ebb <- mergeData(PRS=PRS, PCA=PCA, Phenotypes=ebb)

# Keep unrelated individuals
ebb = ebb[vkood %in% non_rel$vkood, ]

# Adjust PRS for demographic covariates (0 PCs)
lm_prs <- lm("PRS ~ Sex + YoB + I(YoB^2) + I(Sex*YoB)", data = ebb)
ebb$PRS_adj <- ebb$PRS - predict(object = lm_prs, newdata = ebb)
ebb[, PRS_adj := scale(PRS_adj)]
colnames(ebb)[which(colnames(ebb)=="PRS_adj")] <- "PRS_0"

# Adjust PRS for domographic covariates and first N PCs,
# Where N is from 1 to 100
for(i in 1:100){
  lm_prs <- lm(paste0("PRS ~ Sex + YoB + I(YoB^2) + I(Sex*YoB) + ", paste("PC", 1:i, sep = "", collapse = " + ")),
               data = ebb)
  ebb$PRS_adj <- ebb$PRS - predict(object = lm_prs, newdata = ebb)
  ebb[, PRS_adj := scale(PRS_adj)]
  colnames(ebb)[which(colnames(ebb)=="PRS_adj")] <- paste0("PRS_", i)
}

# calculate and save basic result table
results <- calculateBasicStatistics(ebb = ebb, prs_codes=c("PRS", paste0("PRS_", 0:100)), only.all=TRUE)
write.table(results, paste0(opt$out,"/results_6138_1_nonrel_cumulative_adjustment.tsv"), row.names = F, quote = F, sep = "\t")

# calculate and save the table with statistics on PoB and PoR differences
results_dif <- testDifferenceBirthRes(ebb = ebb, prs_codes=c("PRS", paste0("PRS_", 0:100)), only.all=TRUE)
write.table(results_dif, paste0(opt$out, "results_6138_1_modelDif_nonrel_cumulative_adjustment.tsv"), row.names = F, quote = F, sep = "\t")
