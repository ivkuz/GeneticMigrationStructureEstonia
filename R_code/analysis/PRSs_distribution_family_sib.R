################################################################################
# Script calculating the mean and sd for PGSs by region for all the subcohorts #
# in sibling design (deviation of individual's PGS from the sibship mean PGS) ##
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


# Calculating mean value ans SD for all the counties for all the subsets
# Special function for the sibling design
calculateBasicStatistics <- function(k, prs_codes, only.all=FALSE){
  
  k_orig <- k
  # table with future results
  results <- data.frame()

  for(prs in prs_codes){
    
    k <- k_orig
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
    ksib1 <- k[, .(ID1, PoB.x, PoR.x, prs1)]
    ksib2 <- k[, .(ID2, PoB.y, PoR.y, prs2)]
    colnames(ksib1) <- c("vkood", "PoB", "PoR", "prs")
    colnames(ksib2) <- c("vkood", "PoB", "PoR", "prs")
    ksib <- rbind(ksib1, ksib2)
    ksib <- unique(ksib)
    
    
    ebb_tmp <- ksib
    

    # Calculate Fstat and Var(county)    
    for(place in c("PoB", "PoR")){
      
      ebb_tmp2 <- ebb_tmp
      names(ebb_tmp2)[names(ebb_tmp2) == place] <- "place"

      for(county in c("All", unique(ebb_tmp2[, place]))){
        
        if(county == "All"){
          ebb_tmp3 <- ebb_tmp2
        } else{
          ebb_tmp3 <- ebb_tmp2[place == county, ] 
        }
        
        m <- ebb_tmp3[, mean(prs, na.rm=TRUE)]
        sd <- ebb_tmp3[, sd(prs, na.rm=TRUE)]
        N <- nrow(ebb_tmp3)

        output <- c(prs, place, county,
                    round(m, 8), round(sd, 8), N)
        results <- rbind(results, output)
        
      }
    }
  }
  
  colnames(results) <- c("Variable", "Birth/Residence", "County", 
                         "mean", "sd", "N")
  
  return(as.data.table(results))
  
}


# Check if significantly different from zero
anovaFstat <- function(means, sds, ns){
  SSE = sum((ns-1)*sds^2)
  n <- sum(ns)
  k <- length(ns)
  global_mean <- sum((means*ns))/sum(ns)
  SSA = sum(ns*(means-global_mean)^2)
  # SST = SSE + SSA
  # sd_global <- sqrt(SST/(sum(ns)-1))
  MSA <- SSA/(k-1)
  MSE <- SSE/(n-k)
  Fstat <- MSA/MSE
  df1 <- k-1
  df2 <- n-k
  pval <- pf(q = Fstat, df1 = df1, df2 = df2, lower.tail = FALSE)
  prop <- SSA/(SSA+SSE)
  
  out <- list(Fstat = Fstat, df1 = df1, df2 = df2, pval = pval, prop=prop)
  
  return(out)
  
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

# Adjustment PRS for PGSea
ebb <- adjustPRS_forEAprs(ebb, prs_codes=paste0(prs_codes, "_adj2"), adj = "EA4_adj2")

# Select the pairs of siblings from KING kinship table
k <- fread(opt$kin)
k <- k[InfType=="FS" & Kinship > 0.177 & Kinship < 0.354, .(ID1, ID2)]
kvariables <- c("vkood", "YoB", "PoB", "PoR", "Nat", prs_codes, paste0(prs_codes, "_adj2_adjPGSea"), paste0("PC", 1:100))
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

# Save sibling table
write.table(k, paste0(opt$out, "siblings.tsv"), row.names = F, quote = F, sep = "\t")

# calculate and save basic result table
results <- calculateBasicStatistics(k = k, prs_codes=c(prs_codes, paste0("PC", 1:100)), only.all=TRUE)

results[, p := pnorm(-abs(as.numeric(mean)), 0, as.numeric(sd)/sqrt(as.numeric(N)))]
write.table(results, paste0(opt$out, "results_family_sib_est.tsv"), row.names = F, quote = F, sep = "\t")



# Calculate Variance explained by regions and the respective F-statistic
results <- fread(paste0(opt$out, "results_family_sib_est.tsv"))
fstat_tab <- data.frame()
for(trait in unique(results$Variable)){
  for(BorR in c("PoB", "PoR")){
    results_tmp <- results[Variable==trait &
                             `Birth/Residence`==BorR & 
                             County != "All", ]
    if(nrow(results_tmp)>0){
      fstat <- unlist(anovaFstat(means = results_tmp$mean, results_tmp$sd, results_tmp$N))
      fstat_tab <- rbind(fstat_tab, c(trait, BorR, fstat))
    }
  }
}
colnames(fstat_tab) <- c("Trait", "BorR", "Fstat", "df1", "df2", "pval", "prop_var")

write.table(fstat_tab, paste0(opt$out, "Fstat_family_sib_est.tsv"), quote = F, sep = "\t", row.names = F)



###########################################
# Repeat for PGSs adjusted for PGS for EA #
###########################################
results <- calculateBasicStatistics(k = k, prs_codes=paste0(prs_codes, "_adj2_adjPGSea"), only.all=TRUE) # [which(prs_codes != "EA4")]

results[, p := pnorm(-abs(as.numeric(mean)), 0, as.numeric(sd)/sqrt(as.numeric(N)))]
write.table(results, paste0(opt$out, "results_family_sib_est_adjEA4prs.tsv"), row.names = F, quote = F, sep = "\t")

results <- fread(paste0(opt$out, "results_family_sib_est_adjEA4prs.tsv"))
fstat_tab <- data.frame()
for(trait in unique(results$Variable)){
  for(BorR in c("PoB", "PoR")){
    results_tmp <- results[Variable==trait &
                             `Birth/Residence`==BorR & 
                             County != "All", ]
    if(nrow(results_tmp)>0){
      fstat <- unlist(anovaFstat(means = results_tmp$mean, results_tmp$sd, results_tmp$N))
      fstat_tab <- rbind(fstat_tab, c(trait, BorR, fstat))
    }
  }
}
colnames(fstat_tab) <- c("Trait", "BorR", "Fstat", "df1", "df2", "pval", "prop_var")

write.table(fstat_tab, paste0(opt$out, "Fstat_family_sib_est_adjEA4prs.tsv"), quote = F, sep = "\t", row.names = F)

