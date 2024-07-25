###############################################################################################
# Script calculating the correlations between polygenic scores from population-bsed GWAS data #
###############################################################################################


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
    
    lm_prs <- lm(paste0("PRS ~ ", paste("PC", 1:100, sep = "", collapse = " + ")),
                 data = tab)
    tab$PRS_adj <- tab$PRS - predict(object = lm_prs, newdata = tab)
    
    lm_prs2 <- lm(paste0("PRS ~ Sex + YoB + I(YoB^2) + I(Sex*YoB) + ", paste("PC", 1:100, sep = "", collapse = " + ")),
                  data = tab)
    tab$PRS_adj2 <- tab$PRS - predict(object = lm_prs2, newdata = tab)
    
    tab[, PRS := (PRS-mean(PRS))/sd(PRS)]
    tab[, PRS_adj := (PRS_adj-mean(PRS_adj))/sd(PRS_adj)]
    tab[, PRS_adj2 := (PRS_adj2-mean(PRS_adj2))/sd(PRS_adj2)]
    
    colnames(tab)[which(colnames(tab)=="PRS")] <- prs
    colnames(tab)[which(colnames(tab)=="PRS_adj")] <- paste0(prs, "_adj")
    colnames(tab)[which(colnames(tab)=="PRS_adj2")] <- paste0(prs, "_adj2")
    
  }
  

    return(tab)
}


option_list = list(
  make_option(c("-s", "--scores"), type="character", default=NULL, 
              help="file with vkoods and all the adjusted PRSs", metavar="character"),
  make_option(c("-u", "--nonrel"), type="character", default="non_relatives.tsv", 
              help="list of non-related individuals [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="results.tsv", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

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


# File with PRS
ebb <- fread(opt$scores)

# unrelated individuals
non_rel <- fread(opt$nonrel)
ebb<- ebb[vkood %in% non_rel$vkood,]

# Correlations with EA4
correlations <- data.frame()
for(prs_name in prs_codes){
  prs <- ebb[, .SD,.SDcols=which(colnames(ebb)==paste0(prs_name, "_adj2"))]
  cor_val_adj2 <- cor(prs, ebb$`EA4_adj2`)
  correlations <- rbind(correlations, c(prs_name, "EA4", cor_val_adj2))
}

colnames(correlations) <- c("Variable1", "Variable2", "cor_adj2")
write.table(correlations, opt$out, row.names = F, quote = F, sep = "\t")
