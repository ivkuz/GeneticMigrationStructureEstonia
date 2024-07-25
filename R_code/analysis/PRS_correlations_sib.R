#############################################################################################
# Script calculating the correlations between polygenic scores from sibship-based GWAS data #
#############################################################################################


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


prs_codes <- 4813:4860

# File with PRS
ebb <- fread(opt$scores)

# unrelated individuals
non_rel <- fread(opt$nonrel)
ebb<- ebb[vkood %in% non_rel$vkood,]

# Deal with population and sibship GWAS results separately
correlations <- data.frame()
for(prs_name in prs_codes){
  if(prs_name %% 2 == 0){
    prs_ref_name <- 4836
  } else{
    prs_ref_name <- 4835
  }
  prs_ref <- ebb[, .SD,.SDcols=which(colnames(ebb)==prs_ref_name)]
  prs <- ebb[, .SD,.SDcols=which(colnames(ebb)==prs_name)]
  cor_val <- cor(prs, prs_ref)
  prs_ref <- ebb[, .SD,.SDcols=which(colnames(ebb)==paste0(prs_ref_name, "_adj"))]
  prs <- ebb[, .SD,.SDcols=which(colnames(ebb)==paste0(prs_name, "_adj"))] 
  cor_val_adj <- cor(prs, prs_ref)
  prs_ref <- ebb[, .SD,.SDcols=which(colnames(ebb)==paste0(prs_ref_name, "_adj2"))]
  prs <- ebb[, .SD,.SDcols=which(colnames(ebb)==paste0(prs_name, "_adj2"))]
  cor_val_adj2 <- cor(prs, prs_ref)
  correlations <- rbind(correlations, c(prs_name, prs_ref_name, cor_val, cor_val_adj, cor_val_adj2))
}

colnames(correlations) <- c("Variable1", "Variable2", "cor", "cor_adj", "cor_adj2")
write.table(correlations, opt$out, row.names = F, quote = F, sep = "\t")