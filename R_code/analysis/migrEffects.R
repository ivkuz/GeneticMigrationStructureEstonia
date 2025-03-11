########################################################################################
# Script estimating the effects of PGSs on migration from ORE in unrelated individuals #
########################################################################################


my_packages <- c("Rcpp", "data.table", "pROC", "optparse")
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if(length(not_installed)) install.packages(not_installed)

library(data.table)
library(pROC)
library(Rcpp)
library(optparse)


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
  make_option(c("-e", "--ebb"), type="character", default="PRSs_adj.tsv", 
              help="file with adjusted PGSs, PCs amd phenotypes", metavar="character"),
  make_option(c("-u", "--nonrel"), type="character", default="non_relatives.tsv", 
              help="list of non-related individuals [default= %default]", metavar="character"),
  make_option(c("-m", "--migr"), type="character", default="non_relatives.tsv", 
              help="unrelated individuals born in ORE with migration phenotype", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="~/tmp/", 
              help="output directory name [default= %default]", metavar="character")
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
               960, "alpha1-antagonist_alphablocker", "dexamethasone_aceticacid_neomycin", "EA4")



# Upload ebb data file with adjusted PGSs, PCs amd phenotypes
ebb <- fread(opt$ebb)


# Adjustment PRS for PGSea4
ebb <- adjustPRS_forEAprs(ebb, prs_codes=paste0(prs_codes, "_adj2"), adj = "EA4_adj2")


# Migration phenotypes for individuals born in ORE
ore_born_unrel <- fread(opt$migr, header = F)
ore_born_unrel <-  ore_born_unrel[, c(1,3)]
colnames(ore_born_unrel) <- c("vkood", "PoR_macro")
ebb <- merge(ebb, ore_born_unrel, by = "vkood")


# PGS (adjusted) or PC effects on migration from ORE
prs_codes_tmp_list <- list(paste0(prs_codes, "_adj2"), paste0("PC", 1:100), paste0(prs_codes[prs_codes != "EA4"], "_adj2_adjPGSea"))
save_file_list <- c("population2.0_PRS_effects_0pc_cities_PRSs.rds", "population2.0_PRS_effects_0pc_cities_PCs.rds",
                    "population2.0_PRS_effects_0pc_cities_PRSs_adjEAadj.rds ")

for(t in 1:3){
  
  lm_pop_list <- list()
  prs_codes_tmp <- prs_codes_tmp_list[[t]]
  
  for (i in seq_along(prs_codes_tmp)){
    
    print(c(i, prs_codes_tmp[i]))
    
    names(ebb)[names(ebb) == prs_codes_tmp[i]] <- "prs"
    lm_pop <- glm("PoR_macro ~ prs + Sex + Age",
                  data = ebb, family = "binomial")
    lm_pop_list[[prs_codes_tmp[i]]] <- summary(lm_pop)$coefficients
    names(ebb)[names(ebb) == "prs"] <- prs_codes_tmp[i]
    
  }
  
  saveRDS(lm_pop_list, paste0(opt$out, save_file_list[t]))
  
}
