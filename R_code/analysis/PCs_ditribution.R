##################################################################################@
# Script calculating the mean and sd for 100 PCs by region for all the subcohorts #
##################################################################################@


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



# Combine PCA, Phenotypes
mergeData <- function(PCA, Phenotypes){
  merged_data <- merge(PCA, Phenotypes, by="vkood")
  return(merged_data)
}


# Calculating mean value ans SD for all the counties for all the subsets
calculateBasicStatistics <- function(ebb, only.all=FALSE){
  
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
          
          for(pc in paste0("PC", 1:100)){
            
            ebb_tmp5 <- ebb_tmp4
            names(ebb_tmp5)[names(ebb_tmp5) == pc] <- "pc"
            ebb_tmp5 <- ebb_tmp5[!is.na(pc), ]
            
            for(place in c("PoB", "PoR")){
              
              ebb_tmp6 <- ebb_tmp5
              names(ebb_tmp6)[names(ebb_tmp6) == place] <- "place"
              ebb_tmp6 <- ebb_tmp6[!is.na(place), ]
              
              for(county in unique(ebb_tmp6[!is.na(place), place])){
                
                ebb_tmp7 <- ebb_tmp6[place == county, ] 
                
                m <- ebb_tmp7[, mean(pc, na.rm=TRUE)]
                sd <- ebb_tmp7[, sd(pc, na.rm=TRUE)]
                N <- nrow(ebb_tmp7[!is.na(pc),])
                
                output <- c(nat, sex, age, yoa, pc, place, county, 
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
testDifferenceBirthRes <- function(ebb, only.all = FALSE){
  
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
          
          for(pc in paste0("PC", 1:100)){
            
            ebb_tmp5 <- ebb_tmp4
            names(ebb_tmp5)[names(ebb_tmp5) == pc] <- "pc"
            ebb_tmp5 <- ebb_tmp5[!is.na(pc), ]
            
            lm1 <- lm(paste0("pc ~ PoB"), data=ebb_tmp5)
            lm2 <- lm(paste0("pc ~ PoR"), data=ebb_tmp5)
            lm3 <- lm(paste0("pc ~ PoB + PoR"), data=ebb_tmp5)
            if(summary(lm1)$r.squared > summary(lm2)$r.squared){
              p <- anova(lm3, lm2, test="F")[["Pr(>F)"]][2]
            } else{
              p <- anova(lm3, lm1, test="F")[["Pr(>F)"]][2]
            }
            output <- c(nat, sex, age, yoa, pc, p)
            results <- rbind(results, output)
          }
        }
      }
    }
  }
  
  colnames(results) <- c("Nat", "Sex", "Age", "YoA", "Variable", "P")
  
  return(results)
  
}


option_list = list(
  make_option(c("-e", "--pca_est"), type="character", default=NULL, 
              help="file with 100 PCs of estonians", metavar="character"),
  make_option(c("-r", "--pca_rus"), type="character", default=NULL, 
              help="file with 100 PCs of russians", metavar="character"),
  make_option(c("-f", "--pheno"), type="character", default=NULL, 
              help="file with phenotypes", metavar="character"),
  make_option(c("-u", "--nonrel"), type="character", default="non_relatives.tsv", 
              help="list of non-related individuals [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="results.tsv", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# File with PCs
PCA_est <- fread(opt$pca_est)
PCA_rus <- fread(opt$pca_rus)
PCA <- rbind(PCA_est, PCA_rus)
PCA <- PCA[, FID:=NULL]
colnames(PCA)[1] <- "vkood"

# File with phenotypes
ebb <- fread(opt$pheno)

# unrelated individuals
non_rel <- fread(opt$nonrel)

# Combine PCA, Phenotypes
ebb <- mergeData(PCA=PCA, Phenotypes=ebb)

# calculate and save basic result table
results <- calculateBasicStatistics(ebb = ebb[vkood %in% non_rel$vkood,], only.all=TRUE)
write.table(results, paste0(opt$out, "results_PCA_ages2_yoa_nonrel.tsv"), row.names = F, quote = F, sep = "\t")

# calculate and save the table with statistics on PoB and PoR differences
results_dif <- testDifferenceBirthRes(ebb = ebb[vkood %in% non_rel$vkood,], only.all=TRUE)
write.table(results_dif, paste0(opt$out, "results_PCA_modelDif_ages2_yoa_nonrel.tsv"), row.names = F, quote = F, sep = "\t")

