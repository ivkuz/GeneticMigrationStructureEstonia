#########################################################################
# Supplementary Figure 29 - non-rabdom migration and assortative mating #
#########################################################################



library(data.table)
library(ggplot2)


parent = function(x,output){
  
  ID1 = x[1]
  ID2 = x[2]
  YoB1= as.numeric(x[3])
  YoB2= as.numeric(x[4])
  
  if(YoB2-YoB1>16){
    return(ID1)
  }else if(YoB1-YoB2>16){
    return(ID2)
  }else{
    return(NA)
  }
  
}

offspring = function(x,output){
  
  ID1 = x[1]
  ID2 = x[2]
  YoB1= as.numeric(x[3])
  YoB2= as.numeric(x[4])
  
  if(YoB2-YoB1>16){
    return(ID2)
  }else if(YoB1-YoB2>16){
    return(ID1)
  }else{
    return(NA)
  }
  
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
  make_option(c("-u", "--unrel"), type="character", default=NULL, 
              help="file with a subset of unrelated individuals", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output directory", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)




PRS <- NULL
prs_codes <- c("6138_1", "50", "21001", "EA4")

for(code in prs_codes){
  print(code)
  PRS_tmp <- fread(paste0(opt$score, "/", code, "/", code, ".chrALL.sscore"))
  colnames(PRS_tmp) <- c("vkood", code)
  if(is.null(PRS)){
    PRS <- PRS_tmp
  } else{
    PRS <- merge(PRS, PRS_tmp, by = "vkood")
  }
}

# File with phenotypes
ebb <- fread(opt$pheno)


# File with PCs
PCA_est <- fread(opt$pca_est)
PCA_rus <- fread(opt$pca_rus)
PCA <- rbind(PCA_est, PCA_rus)
PCA <- PCA[, FID:=NULL]
colnames(PCA)[1] <- "vkood"



################################
# Here we start work with data #
################################

# Combine PRS, PCA, Phenotypes
# ebb <- merge(PRS, ebb, by="vkood")
ebb <- mergeData(PRS=PRS, PCA=PCA, Phenotypes=ebb)


# Adjustment PRS for PCA and normalization
ebb <- ebb[Nat=="Eestlane", ]

ebb_yob <- ebb[, .(vkood, YoB)]

colnames(ebb)[2] <- "PRS"

ebb <- adjustPRS(tab = ebb, prs_codes = prs_codes)



# selecting independendent pairs having a child
k <- fread(opt$kin)

# exclude MZ - take a random twin as PGS is the same
k_mz <- k[InfType=="Dup/MZ",]
k <- k[!(ID1 %in% k_mz$ID1 | ID2 %in% k_mz$ID1)]

# select Parent-Offspring pairs
k <- k[InfType=="PO" & Kinship > 0.177 & Kinship < 0.354, .(ID1, ID2)]
k <- merge(k, ebb[, .(vkood, YoB)], by.x="ID1", by.y="vkood")
k <- merge(k, ebb[, .(vkood, YoB)], by.x="ID2", by.y="vkood")
k <- k[, .(ID1, ID2, YoB.x, YoB.y)]

# determine who is parent in each pair
k[, P := unlist(apply(as.data.frame(k), 1 , parent))]
k[, O := unlist(apply(as.data.frame(k), 1 , offspring))]
k <- k[!is.na(P) & !is.na(O), ]

# keep offsprings with two parents
withParents <- unique(k$O[which(duplicated(k$O))])
k <- k[O %in% withParents,]

# merge parents of the same offspring
k_parents1 <-k[duplicated(O),]
k_parents2 <-k[!duplicated(O),]
k_parents <- merge(k_parents1, k_parents2, by="O")

# keep pairs to get rid of pairs including the same individuals
parents <- c(k_parents$P.x, k_parents$P.y)
parents_selected <- !duplicated(parents)
parents_selected <- parents_selected[1:(length(parents_selected)/2)] * 
  parents_selected[(length(parents_selected)/2+1):length(parents_selected)]
parents_selected <- as.logical(parents_selected)
k_parents <- k_parents[parents_selected, .(P.x, P.y)]
colnames(k_parents) <- c("ID1", "ID2")

# correlations between spouses
kvariables <- c("vkood", "PoB", "PoR", paste0(prs_codes, "_adj2"))
k_parents <- merge(k_parents, ebb[, ..kvariables], by.x="ID1", by.y="vkood")
k_parents <- merge(k_parents, ebb[, ..kvariables], by.x="ID2", by.y="vkood")

cor_spouses <- k_parents[, .(cor(`6138_1_adj2.x`, `6138_1_adj2.y`), cor(`EA4_adj2.x`, `EA4_adj2.y`), 
                             cor(`50_adj2.x`, `50_adj2.y`), cor(`21001_adj2.x`, `21001_adj2.y`))]  
colnames(cor_spouses) <- c("6138_1", "EA4", "50", "21001")
cor_spouses_bootstrap <- c()
set.seed(100)
for(i in 1:1000){
  s <- sample(1:nrow(k_parents), nrow(k_parents), replace = T)
  cor_spouses_bootstrap <- rbind(cor_spouses_bootstrap,
                             k_parents[s, .(cor(`6138_1_adj2.x`, `6138_1_adj2.y`), 
                                            cor(`EA4_adj2.x`, `EA4_adj2.y`), 
                                            cor(`50_adj2.x`, `50_adj2.y`),
                                            cor(`21001_adj2.x`, `21001_adj2.y`))])
  
}
colnames(cor_spouses_bootstrap) <- c("6138_1", "EA4", "50", "21001")


# independent individuals
non_rel <- fread(opt$unrel)
ebb <- ebb[vkood %in% non_rel$vkood, ]

# completely random pairs
if((nrow(ebb) %% 2) != 0){
  ebb <- ebb[-nrow(ebb),]
}

cor_random_distr <- c()
set.seed(100)
for(i in 1:1000){
  R1 <- sample(x = 1:nrow(ebb), size = as.integer(nrow(ebb)/2), replace = F)
  R2 <- (1:nrow(ebb))[which(!(1:nrow(ebb) %in% R1))]
  pairs_rand <- as.data.table(
    cbind(ebb[R1, `6138_1_adj2`], ebb[R2, `6138_1_adj2`], 
          ebb[R1, EA4_adj2], ebb[R2, EA4_adj2], 
          ebb[R1, `50_adj2`], ebb[R2, `50_adj2`],
          ebb[R1, `21001_adj2`], ebb[R2, `21001_adj2`])
  )

  cor_random_distr <- rbind(cor_random_distr,
                            pairs_rand[, .(cor(V1, V2), cor(V3, V4), cor(V5, V6), cor(V7, V8))])
}
colnames(cor_random_distr) <- c("6138_1", "EA4", "50", "21001")


# random pairs wit the same PoB or PoR
counties <- unique(ebb$PoB)

cor_pob_distr <- c()
cor_por_distr <- c()
set.seed(100)
for(i in 1:1000){
  print(i)
  pairs_pob <- data.frame()
  pairs_por <- data.frame()
  for(region in counties){
    ebb_pob <- ebb[PoB==region,]
    if((nrow(ebb_pob) %% 2) != 0){ebb_pob <- ebb_pob[-nrow(ebb_pob),]}
    R1 <- sample(x = 1:nrow(ebb_pob), size = as.integer(nrow(ebb_pob)/2), replace = F)
    R2 <- (1:nrow(ebb_pob))[which(!(1:nrow(ebb_pob) %in% R1))]
    pairs_pob <- rbind(pairs_pob,
                       as.data.table(
                         cbind(ebb_pob[R1, .(`6138_1_adj2`, EA4_adj2, `50_adj2`, `21001_adj2`)],
                               ebb_pob[R2, .(`6138_1_adj2`, EA4_adj2, `50_adj2`, `21001_adj2`)])
                       )
    )
    
    ebb_por <- ebb[PoR==region,]
    if((nrow(ebb_por) %% 2) != 0){ebb_por <- ebb_por[-nrow(ebb_por),]}
    R1 <- sample(x = 1:nrow(ebb_por), size = as.integer(nrow(ebb_por)/2), replace = F)
    R2 <- (1:nrow(ebb_por))[which(!(1:nrow(ebb_por) %in% R1))]
    pairs_por <- rbind(pairs_por,
                       as.data.table(
                         cbind(ebb_por[R1, .(`6138_1_adj2`, EA4_adj2, `50_adj2`, `21001_adj2`)],
                               ebb_por[R2, .(`6138_1_adj2`, EA4_adj2, `50_adj2`, `21001_adj2`)])
                       )
    )
    
  }
  
  colnames(pairs_pob) <- paste0("V", 1:8)
  colnames(pairs_por) <- paste0("V", 1:8)
  cor_pob_distr <- rbind(cor_pob_distr,
                            pairs_pob[, .(cor(V1, V5), cor(V2, V6), cor(V3, V7), cor(V4, V8))])
  cor_por_distr <- rbind(cor_por_distr,
                         pairs_por[, .(cor(V1, V5), cor(V2, V6), cor(V3, V7), cor(V4, V8))])
  
  
}

colnames(cor_pob_distr) <- c("6138_1", "EA4", "50", "21001")
colnames(cor_por_distr) <- c("6138_1", "EA4", "50", "21001")



cor_spouses_bootstrap[, names(cor_spouses_bootstrap) := lapply(.SD, sort)]
cor_random_distr[, names(cor_random_distr) := lapply(.SD, sort)]
cor_pob_distr[, names(cor_pob_distr) := lapply(.SD, sort)]
cor_por_distr[, names(cor_por_distr) := lapply(.SD, sort)]


cor_spouses_bootstrap_sel <- cor_spouses_bootstrap[c(26, 500, 975),]
cor_spouses_bootstrap_sel[2, ] <- cor_spouses
cor_random_distr_sel <- cor_random_distr[c(26, 500, 975),]
cor_pob_distr_sel <- cor_pob_distr[c(26, 500, 975),]
cor_por_distr_sel <- cor_por_distr[c(26, 500, 975),]
cor_results <- data.table()
for(cor_sel in list(cor_spouses_bootstrap_sel, cor_random_distr_sel, cor_pob_distr_sel, cor_por_distr_sel)){
  colnames(cor_sel) <- c("EA", "EA4", "Height", "BMI")
  cor_sel[, Est := c("low", "mean", "high")]
  cor_sel <- as.data.table(melt(cor_sel, id.vars = "Est"))
  cor_results <- rbind(cor_results, cor_sel)
}
cor_results[, type := rep(c("Spouses", "Random", "POB", "POR"), each = 12)]
cor_results[, type := factor(type, levels = c("Random", "POB", "POR", "Spouses"))]
colnames(cor_results)[2] <- "PGS"
cor_results <- dcast(cor_results, PGS + type ~ Est, value.var = "value")

cor_list <- list(cor_random_distr, cor_pob_distr,
                 cor_por_distr, cor_spouses_bootstrap,
                 cor_results)
save(cor_list, file = paste0(opt$out, "am_4_PGS_adj2_list.Rdata"))


# Make figure
pdf(paste0(opt$out, "am_4_PGS_adj2.pdf"), width=5, height = 4)
ggplot(cor_results, aes(x=type, y = mean, color = PGS)) + 
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.9)) + 
  theme_bw() + theme(text = element_text(size=10),
                     panel.grid.major.x = element_blank()) +
  geom_errorbar(aes(ymin =  low, ymax = high), width = 0.3,
                position = position_dodge(width = 0.9)) +
  ylab("Correlation") + xlab("Pairing")
dev.off()

