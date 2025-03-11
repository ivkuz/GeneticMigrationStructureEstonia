###################################################################################
# Script estimating the effects of PGSs and PCs on migration from ORE in siblings #
###################################################################################


Sys.setlocale("LC_CTYPE", "estonian")

my_packages <- c("Rcpp", "data.table", "pROC", "optparse")
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if(length(not_installed)) install.packages(not_installed)

library(data.table)
library(pROC)
library(Rcpp)
library(optparse)
library(lme4)


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


# Adjustment PRS for PCA and normalization
adjustPRS <- function(tab, prs_codes){
  
  for(prs in prs_codes){
    
    colnames(tab)[which(colnames(tab)==prs)] <- "PRS"
    
    lm_prs2 <- lm(paste0("PRS ~ Sex + YoB + I(YoB^2) + I(Sex*YoB) + ", paste("PC", 1:100, sep = "", collapse = " + ")),
                  data = tab)
    
    tab$PRS_adj2 <- tab$PRS - predict(object = lm_prs2, newdata = tab)
    tab[, PRS_adj2 := scale(PRS_adj2)]
    
    colnames(tab)[which(colnames(tab)=="PRS")] <- prs
    colnames(tab)[which(colnames(tab)=="PRS_adj2")] <- paste0(prs, "_adj2")
    
  }
  
  return(tab)
  
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



# function to estimate the effects of PGSs or other variables on migration within sibships
# k - data table with all the required information on the individuals
# prs - the name of the analysed variable
# npc - number of PCs as covariates
# FamilyRandom - whether FID must be treated as random effects variable
# pop - if TRUE, one random individual is chosen from each family, effects of the variable in "unrelated" individuals
# pop - if FALSE, effects in sibships, with within- and between-sibship variable components
sibshipEffects <- function(k, prs, npc, FamilyRandom, pop = F){
  
  # unify the name of the variable
  names(k)[names(k) == paste0(prs, ".x")] <- "prs.x"
  names(k)[names(k) == paste0(prs, ".y")] <- "prs.y"
  
  # calculate mean-sibship value
  k_m <- k[, .(mean(prs.x), mean(prs.y)), by=FID]
  k_m[, prs := rowMeans(k_m[, c("V1", "V2")])]
  k <- merge(k, k_m[, .(FID, prs)], by="FID")
  
  # choose the number of PCs required to be included in the model as covariates
  pcs <- paste0("PC", 1:max(10, npc))
  pcs <- pcs[pcs != prs]
  
  # calculate a within-sibship deviation of the variable and transform the table to a long form
  k[, prs1 := prs.x - prs]
  k[, prs2 := prs.y - prs]
  ksib1 <- k[, c("ID1", "FID", "PoB.x", "PoR.x", "ParishRes.x", "prs.x", "prs1", "prs", "Age.x", "Sex.x", "YOB1", "old.x", paste0(pcs, ".x")), with = FALSE]
  ksib2 <- k[, c("ID2", "FID", "PoB.y", "PoR.y", "ParishRes.y", "prs.y", "prs2", "prs", "Age.y", "Sex.y", "YOB2", "old.y", paste0(pcs, ".y")), with = FALSE]
  colnames(ksib1) <- c("vkood", "FID", "PoB", "PoR", "ParishRes", "prs_pop", "prs", "prs_fam", "Age", "Sex", "yob", "old", pcs)
  colnames(ksib2) <- c("vkood", "FID", "PoB", "PoR", "ParishRes", "prs_pop", "prs", "prs_fam", "Age", "Sex", "yob", "old", pcs)
  ksib <- rbind(ksib1, ksib2)
  ksib <- unique(ksib)
  
  # Scale variables for better LMM fitting
  ksib[, Age := scale(Age)]

  # Filter individuals born in ORE and obtain the migration phenotype
  ksib <- ksib[!(PoB %in% c("Harju", "Tartu")) & (ParishRes %in% c("Tartu linn", "Tallinn") | !(PoR %in% c("Harju", "Tartu"))),]
  ksib[!(ParishRes %in% c("Tartu linn", "Tallinn")), PoR_macro := 0]
  ksib[ParishRes %in% c("Tartu linn", "Tallinn"), PoR_macro := 1]
  
  # aggregate string with PC covariates
  if( npc == 0){
    pc_cov <- ""
  } else{
    pc_cov <- paste0(" + ", paste(pcs, sep = "", collapse = " + ") )
  }
  
  
  # Fitting the models
  if( FamilyRandom ){
    
    if(pop){
      lm_sib <- glmer(paste0("PoR_macro ~ prs_pop + (1|FID) + Sex + Age", pc_cov), #   + Age2 + SxA + SxA2 + old
                      data = ksib, family = "binomial")
    } else{
      lm_sib <- glmer(paste0("PoR_macro ~ prs_fam + (1|FID) + prs + Sex + Age", pc_cov), #   + Age2 + SxA + SxA2 + old
                      data = ksib, family = "binomial")
    }
    
  } else{
    
    if(pop){
      lm_sib <- glm(paste0("PoR_macro ~ prs_pop + Sex + Age", pc_cov), #  + Age2 + SxA + SxA2 + old
                      data = ksib[, .SD[sample(.N, 1)], by = FID], family = "binomial")
    } else{
      lm_sib <- glm(paste0("PoR_macro ~ prs_fam + prs + Sex + Age", pc_cov), #  + Age2 + SxA + SxA2 + old
                    data = ksib, family = "binomial")
    }
    
  }
  
  return(summary(lm_sib)$coefficients)
  
}




# function to estimate the effects of PGSs or other variables on migration within sibships
# in presence of EA phenotype covariates
# k - data table with all the required information on the individuals
# prs - the name of the analysed variable
sibshipEffectsEAadjusted <- function(k, prs, npc){
  
  # unify the name of the variable
  names(k)[names(k) == paste0(prs, ".x")] <- "prs.x"
  names(k)[names(k) == paste0(prs, ".y")] <- "prs.y"
  
  # calculate mean-sibship value
  k_m <- k[, .(mean(prs.x), mean(prs.y)), by=FID]
  k_m[, prs := rowMeans(k_m[, c("V1", "V2")])]
  k <- merge(k, k_m[, .(FID, prs)], by="FID")
  
  
  # calculate a within-sibship deviation of the variable and transform the table to a long form
  k[, prs1 := prs.x - prs]
  k[, prs2 := prs.y - prs]
  ksib1 <- k[, c("ID1", "FID", "PoB.x", "PoR.x", "ParishRes.x", "prs1", "prs", "Age.x", "Sex.x", "YOB1", "old.x", "EduYears.x", "EA.x"), with = FALSE]
  ksib2 <- k[, c("ID2", "FID", "PoB.y", "PoR.y", "ParishRes.y", "prs2", "prs", "Age.y", "Sex.y", "YOB2", "old.y", "EduYears.y", "EA.y"), with = FALSE]
  colnames(ksib1) <- c("vkood", "FID", "PoB", "PoR", "ParishRes", "prs", "prs_fam", "Age", "Sex", "yob", "old", "EduYears", "EA", pcs)
  colnames(ksib2) <- c("vkood", "FID", "PoB", "PoR", "ParishRes", "prs", "prs_fam", "Age", "Sex", "yob", "old", "EduYears", "EA", pcs)
  ksib <- rbind(ksib1, ksib2)
  ksib <- unique(ksib)
  
  # Scale variables for better LMM fitting
  ksib[, Age := scale(Age)]
  
  # Filter individuals born in ORE and obtain the migration phenotype
  ksib <- ksib[!(PoB %in% c("Harju", "Tartu")) & (ParishRes %in% c("Tartu linn", "Tallinn") | !(PoR %in% c("Harju", "Tartu"))),]
  ksib[!(ParishRes %in% c("Tartu linn", "Tallinn")), PoR_macro := 0]
  ksib[ParishRes %in% c("Tartu linn", "Tallinn"), PoR_macro := 1]
  ksib[, EduYears := scale(EduYears)]
  ksib[, Age := scale(Age)]
  ksib[, SxA := scale(SxA)]
  

  # Fitting the models
  lm_sib <- glm("PoR_macro ~ prs_fam + prs + Sex + Age",
                      data = ksib, family = "binomial")
  lm_sib_EduYears <- glm("PoR_macro ~ prs_fam + prs + EduYears + Sex + Age",
                               data = ksib, family = "binomial")
  lm_sib_EA <- glm("PoR_macro ~ prs_fam + prs + factor(EA) + Sex + Age",
                         data = ksib, family = "binomial") 

  
  return(list(mixed = summary(lm_mix_sib)$coefficients,
              EduYears = summary(lm_mix_sib_EduYears)$coefficients,
              EA = summary(lm_mix_sib_EA)$coefficients))

}




option_list = list(
  make_option(c("-e", "--ebb"), type="character", default="PRSs_adj.tsv", 
              help="file with adjusted PGSs, PCs amd phenotypes", metavar="character"),
  make_option(c("-k", "--kin"), type="character", default=NULL,
              help="file with kinship table", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="~/tmp/",
              help="output directory name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

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
               "EA4") # Add PRS from EA4 excluding 23andMe and EstBB



# Upload ebb data file with adjusted PGSs, PCs amd phenotypes
ebb <- fread(opt$ebb)


# Keep Estonian subsample, scale variables
ebb <- ebb[Nat=="Eestlane",]
ebb[, (prs_codes) := lapply(.SD, scale), .SDcols = prs_codes]
ebb <- adjustPRS(tab = ebb, prs_codes = "EA4")
pc_columns <- paste0("PC", 1:100)
ebb[, (pc_columns) := lapply(.SD, scale), .SDcols = pc_columns]

# delete unnecessary columns
ebb[, c(prs_codes, paste0(prs_codes, "_adj"), paste0(prs_codes, "_adj3")) := NULL]

# Adjustment PRS for PGSea
ebb <- adjustPRS_forEAprs(ebb, prs_codes=paste0(prs_codes, "_adj2"), adj = "EA4_adj2")


# Select the pairs of siblings from KING kinship table
k <- fread(opt$kin)
k <- k[InfType=="FS" & Kinship > 0.177 & Kinship < 0.354, .(ID1, ID2)]
# kvariables <- c("vkood", "Age", "YoB", "Sex", "PoB", "PoR", "ParishRes", "Nat", "EA", "EduYears", prs_codes, paste0(prs_codes, "_adj2"), paste0(prs_codes, "_adj2_adjPGSea"), paste0("PC", 1:100))
kvariables <- c("vkood", "Age", "YoB", "Sex", "PoB", "PoR", "ParishRes", "Nat", "EA", "EduYears", paste0(prs_codes, "_adj2"), paste0(prs_codes, "_adj2_adjPGSea"), paste0("PC", 1:100))
k <- merge(k, ebb[Nat=="Eestlane", ..kvariables], by.x="ID1", by.y="vkood")
k <- merge(k, ebb[Nat=="Eestlane", ..kvariables], by.x="ID2", by.y="vkood")
# Keep only siblings born in the same county
k <- k[PoB.x==PoB.y,]


# find individuals with more than one sibling
mult_sib <- names(table(c(k$ID1, k$ID2)))[table(c(k$ID1, k$ID2))>1]

# create families of several siblings
k$FID <- 1:nrow(k)
fid_excl <- c()
# if somebody from the pair has more than one sibling
for(i in 2:nrow(k)){
  if(k[i, ID1] %in% mult_sib | k[i, ID2] %in% mult_sib){
    # keep all the rows with them (and their siblings) from the rows above to a table
    k_tmp <- k[1:(i-1),][ID1 %in% k$ID1[i] | ID1 %in% k$ID2[i] | ID2 %in% k$ID1[i] | ID2 %in% k$ID2[i], ]
    # if there are some
    if(nrow(k_tmp)>0){
      # if some of that families have different FIDs
      if(length(unique(k_tmp$FID))>1){
        # this should be useful when first two parts of the family appears independently
        # for example data.table(ID1=c(1,3,1,3,1,2), ID2=c(2,4,4,2,3,4), FID=1:6)
        # but the procedure should be changed
        # currently, this family would be deleted completely
        # luckly, king orders the relatives in a good way and we don't have such cases
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


# define if an individual is the first born sib among present in the data
k[, old1 := YoB.x >= YoB.y]
k[, old2 := YoB.y >= YoB.x]
older <- rbind(k[, .(ID1, old1)], k[, .(ID2, old2)], use.names=F)
older <- unique(older)
colnames(older) <- c("ID", "old")
dupl <- unique(older$ID[duplicated(older$ID)])
older[ID %in% dupl, old := FALSE]
older[, old := as.integer(old)]
older <- unique(older)

k[, old1 := NULL]
k[, old2 := NULL]
k <- merge(k, older, by.x = "ID1", by.y = "ID")
k <- merge(k, older, by.x = "ID2", by.y = "ID")


# Calculate deviation of YoB from family mean
k_m <- k[, .(mean(YoB.x), mean(YoB.y)), by=FID]
k_m[, mean_yob := rowMeans(k_m[, c("V1", "V2")])]
k <- merge(k, k_m[, .(FID, mean_yob)], by="FID")

k[, YOB1 := YoB.x-mean_yob]
k[, YOB2 := YoB.y-mean_yob]





# ORE-cities migration within- and between-sibship PGS coefficients
prs_codes_list <- list(paste0("PC", 1:100), paste0(prs_codes, "_adj2"), paste0(prs_codes[which(prs_codes != "EA4")], "_adj2_adjPGSea"))
suffixes <- c("PCs", "PRSs", "PRSs_adjEAadj")
prefixes <- c("fam_effects_0pc_LMM_cities_",
              "fam_effects_0pc_LM_cities_",
              "pop_effects_0pc_LMM_cities_",
              "pop_effects_0pc_LM_cities_")

for(j in 1:3){
  prs_codes_tmp <- prs_codes_list[[j]]
  suff <- suffixes[j]

  lm_sib_list_rand <- list()
  lm_sib_list <- list()
  lm_sib_list_pop_rand <- list()
  lm_sib_list_pop <- list()
  for (i in seq_along(prs_codes)){
    
    lm_sib_list_rand[[prs_codes[i]]] <- sibshipEffects(k=k, prs=prs_codes[i], npc=0, FamilyRandom = T, pop = F)
    lm_sib_list[[prs_codes[i]]] <- sibshipEffects(k=k, prs=prs_codes[i], npc=0, FamilyRandom = F, pop = F)
    lm_sib_list_pop_rand[[prs_codes[i]]] <- sibshipEffects(k=k, prs=prs_codes[i], npc=0, FamilyRandom = T, pop = T)
    lm_sib_list_pop[[prs_codes[i]]] <- sibshipEffects(k=k, prs=prs_codes[i], npc=0, FamilyRandom = F, pop = T)
    
  }

  saveRDS(lm_sib_list_rand, paste0(opt$out, "sibship_PRS2.0_", prefixes[1], suff, ".rds"))
  saveRDS(lm_sib_list, paste0(opt$out, "sibship_PRS2.0_", prefixes[2], suff, ".rds"))
  saveRDS(lm_sib_list_pop_rand, paste0(opt$out, "sibship_PRS2.0_", prefixes[3], suff, ".rds"))
  saveRDS(lm_sib_list_pop, paste0(opt$out, "sibship_PRS2.0_", prefixes[4], suff, ".rds"))
  
}


# PGS adjusted for EA phenotype (for Table 1)
prs_codes <- "6138_1_adj2"
lm_sib_list <- list()
for (i in seq_along(prs_codes)){

  lm_sib_list[[prs_codes[i]]] <- sibshipEffectsEAadjusted(k=k, prs=prs_codes[i])

}
 
saveRDS(lm_sib_list, paste0(opt$out, "sibship_PRS_effects_EAadj_0pc_cities_final.rds"))


