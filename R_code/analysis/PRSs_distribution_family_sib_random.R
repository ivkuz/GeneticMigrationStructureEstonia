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



# Calculating mean value ans SD for all the counties for all the subsets
# Special function for the sibling design
calculateBasicStatistics <- function(k, prs_codes, only.all=FALSE){
  
  k_orig <- k
  # table with future results
  results <- data.frame()
  results_dif <- data.frame()
  
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
    ksib1 <- k[, .(ID1, PoB.x, PoR.x, prs1, Sex.x, YOB1, old.x)]
    ksib2 <- k[, .(ID2, PoB.y, PoR.y, prs2, Sex.y, YOB2, old.y)]
    colnames(ksib1) <- c("vkood", "PoB", "PoR", "prs", "Sex", "yob", "old")
    colnames(ksib2) <- c("vkood", "PoB", "PoR", "prs", "Sex", "yob", "old")
    ksib <- rbind(ksib1, ksib2)
    ksib <- unique(ksib)
    
    # Adjust the random variable for the covariates
    lm_prs <- lm("prs ~ Sex + yob + old + I(Sex*yob) + I(Sex*old)", data = ksib)
    ksib$prs_adj <- ksib$prs - predict(object = lm_prs, newdata = ksib)
    ksib[, prs_adj := scale(prs_adj)]
    
    ebb_tmp <- ksib
    
    
    # Calculate Fstat and Var(county)
    # Var(county) for POB is always zero
    for(place in c("PoR")){
      
      ebb_tmp2 <- ebb_tmp
      names(ebb_tmp2)[names(ebb_tmp2) == place] <- "place"
      
      for(county in c("All", unique(ebb_tmp2[, place]))){
        
        if(county == "All"){
          ebb_tmp3 <- ebb_tmp2
        } else{
          ebb_tmp3 <- ebb_tmp2[place == county, ] 
        }
        
        m <- ebb_tmp3[, mean(prs_adj, na.rm=TRUE)]
        sd <- ebb_tmp3[, sd(prs_adj, na.rm=TRUE)]
        N <- nrow(ebb_tmp3)
        
        output <- c(prs, place, county,
                    round(m, 8), round(sd, 8), N)
        results <- rbind(results, output)
        
      }
    }
    
  }
  
  colnames(results) <- c("Variable", "BorR", "County", 
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
  make_option(c("-f", "--pheno"), type="character", default=NULL, 
              help="file with phenotypes", metavar="character"),
  make_option(c("-k", "--kin"), type="character", default=NULL, 
              help="file with kinship table", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="~/tmp/",
              help="output directory name [default= %default]", metavar="character")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



# File with phenotypes
ebb <- fread(opt$pheno)



# Select the pairs of siblings from KING kinship table
k <- fread(opt$kin)
k <- k[InfType=="FS" & Kinship > 0.177 & Kinship < 0.354, .(ID1, ID2)]
n <- nrow(ebb[vkood %in% c(k$ID1, k$ID2), ])
# Choose the number of the random variables (pseudo-PGS)
r = 10140

# genetate the random variables
set.seed(1)
PRS <- matrix(data = rnorm(n*r), nrow = n)
PRS <- cbind(ebb[vkood %in% c(k$ID1, k$ID2), vkood], 
             as.data.table(PRS))
colnames(PRS)[1] <- "vkood"
ebb <- merge(ebb, PRS)

kvariables <- c("vkood", "YoB", "Sex", "PoB", "PoR", "Nat", paste0("V", 1:r))
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



# calculate and save basic result table
basicStatistics <- calculateBasicStatistics(k = k, prs_codes=paste0("V", 1:r), only.all=TRUE)

results <- as.data.table(basicStatistics)
results[, p := pnorm(-abs(as.numeric(mean)), 0, as.numeric(sd)/sqrt(as.numeric(N)))]
write.table(results, paste0(opt$out, "results_family_sib_est_adj_random.tsv"), row.names = F, quote = F, sep = "\t")


# Calculate Variance explained by regions and the respective F-statistic
results <- fread(paste0(opt$out, "results_family_sib_est_adj_random.tsv"))
fstat_tab <- data.frame()
for(trait in unique(results$Variable)){
  for(BR in c("PoB", "PoR")){
    results_tmp <- results[Variable==trait &
                             BorR==BR &
                             County != "All", ]
    if(nrow(results_tmp)>0){
      print(trait)
      fstat <- unlist(anovaFstat(means = results_tmp$mean, results_tmp$sd, results_tmp$N))
      fstat_tab <- rbind(fstat_tab, c(trait, BorR, fstat))
    }
  }
}
colnames(fstat_tab) <- c("Trait", "BorR", "Fstat", "df1", "df2", "pval", "prop_var")

write.table(Fstat, paste0(opt$out, "Fstat_family_sib_est_adj_random.tsv"), quote = F, sep = "\t", row.names = F)


