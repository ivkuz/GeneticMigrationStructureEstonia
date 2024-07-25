###############################################################################################################
# Script calculating the proportion of variance explained by region and related statistics for PGSs, PCs etc. #
# It uses output file from "PCs/PRSs_distribution.R" as its input #############################################
###############################################################################################################


Sys.setlocale("LC_CTYPE", "estonian")

my_packages <- c("Rcpp", "data.table", "pROC", "optparse")
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if(length(not_installed)) install.packages(not_installed)

library(data.table)
library(reshape2)
library(optparse)


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
  make_option(c("-r", "--results"), type="character", default=NULL,
              help="basic result table from the 'distribution' script", metavar="character"),
  make_option(c("-d", "--diff"), type="character", default=NULL,
              help="result table from the 'distribution' script with statistics on PoB and PoR differences", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="results.tsv",
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


results <- fread(opt$results)

# shorten county names
results[, County := unlist(strsplit(County, " maakond"))]

# Calculate Variance explained by regions and the respective F-statistic
fstat_tab <- data.frame()
for(nat in unique(results$Nat)){
  for(sex in unique(results$Sex)){
    for(age in unique(results$Age)){
      for(yoa in unique(results$YoA)){
        for(trait in unique(results$Variable)){
          for(BorR in c("PoB", "PoR")){
            results_tmp <- results[Nat==nat & Sex==sex & Age == age &
                                     Variable==trait & YoA == yoa &
                                     `Birth/Residence`==BorR & 
                                     County != "All", ]
            if(nrow(results_tmp)>0){
              fstat <- unlist(anovaFstat(means = results_tmp$mean, results_tmp$sd, results_tmp$N))
              fstat_tab <- rbind(fstat_tab, c(nat, sex, age, yoa, trait, BorR, fstat))
            }
          }
        }
      }
    }
  }
}
colnames(fstat_tab) <- c("Nat", "Sex", "Age", "YoA", "Trait", "BorR", "Fstat", "df1", "df2", "pval", "prop_var")

# Table with statistics on PoB and PoR differences
results_dif <- fread(opt$diff)
colnames(results_dif)[5] <- "Trait"

fstat_tab <- merge(fstat_tab, results_dif, by = c("Nat", "Sex", "Age", "YoA", "Trait"))

# Change notations for Nationality and Sex
fstat_tab[Nat=="Eestlane", Nat := "E"]
fstat_tab[Nat=="Venelane", Nat := "R"]
fstat_tab[Sex==1, Sex := "M"]
fstat_tab[Sex==2, Sex := "F"]

write.table(fstat_tab, out$out, quote = F, sep = "\t", row.names = F)
