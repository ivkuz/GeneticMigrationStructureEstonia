#############################################
# Supplementary Figure 5 - migration matrix #
#############################################



Sys.setlocale("LC_CTYPE", "estonian")
library(data.table)
library(reshape2)
library(ggplot2)
library(stringr)
library(tidyr)


plotMigrationMatrix <- function(migr_ebb, census, title){
  residence <- migr_ebb[, .(sum(N)), by = .(PoR)]
  colnames(residence) <- c("PoR", "N_ebb")
  pop_sample <- merge(residence, census)
  sampling_stat <- c(title, sum(pop_sample$N_ebb), round(sum(pop_sample$N_ebb)/sum(pop_sample$N), 2))
  pop_sample[, sampling := N_ebb/N]
  migr <- merge(migr_ebb, pop_sample[, .(PoR, sampling)], 
                by = "PoR")
  migr[, N := round(N/sampling)]
  migr[, N_birth := N/sum(N), by = .(PoB)]
  migr[, N_birth_diag:= N_birth]
  migr <- as.data.table(complete(migr, PoR, PoB, 
                                 fill = list(N_birth = 0, N_birth_diag = 0)))
  migr[PoB==PoR, N_birth_diag:= NA]
  migr$PoB <- factor(x = migr$PoB,
                     levels = unique(migr$PoB[rev(order(migr$PoB))]), 
                     ordered = TRUE)
  
  
  print(ggplot(migr, aes(y = PoB, x = PoR, fill = N_birth)) +
          geom_tile() + theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
                axis.text.y = element_text(size=12),
                axis.title = element_text(size=15),
                title = element_text(size=20),
                legend.text = element_text(size=12),
                legend.title = element_text(size=15)) +
          labs(title=title, fill = "",
               y ="County of birth", x = "County of residence") +
          geom_text(aes(label=round(N_birth, 2))) +
          scale_fill_gradientn(limits = c(0, 1), colours = terrain.colors(10))) # ggtitle(i) + 
  
  return(sampling_stat)
  
}



ebb <- fread("phenotypes/EstBB_filtered.tsv", encoding = 'UTF-8')

# numbers of migrants
nat_bothsexes_18andolder <- ebb[, .N, by = c("Nat", "PoR", "PoB")]
est_sex_18andolder <- ebb[Nat=="Eestlane", .N, by = c("Sex", "PoR", "PoB")]

est_bothsexes <- ebb[Nat=="Eestlane" & Age<=24, .N, by = c("PoR", "PoB")]
est_bothsexes[, Age:="18-24"]
est_bothsexes_ages <- data.frame(est_bothsexes)

est_bothsexes <- ebb[Nat=="Eestlane" & Age>=25 & Age<=34, .N, by = c("PoR", "PoB")]
est_bothsexes[, Age:="25-34"]
est_bothsexes_ages <- rbind(est_bothsexes_ages, data.frame(est_bothsexes))

est_bothsexes <- ebb[Nat=="Eestlane" & Age>=35 & Age<=64, .N, by = c("PoR", "PoB")]
est_bothsexes[, Age:="35-64"]
est_bothsexes_ages <- rbind(est_bothsexes_ages, data.frame(est_bothsexes))

est_bothsexes <- ebb[Nat=="Eestlane" & Age>=65, .N, by = c("PoR", "PoB")]
est_bothsexes[, Age:="65+"]
est_bothsexes_ages <- rbind(est_bothsexes_ages, data.frame(est_bothsexes))

est_bothsexes_ages <- data.table(est_bothsexes_ages)



est_bothsexes <- ebb[Nat=="Eestlane" & YoA<=2015, .N, by = c("PoR", "PoB")]
est_bothsexes[, YoA:="phase1"]
est_bothsexes_yoa <- data.frame(est_bothsexes)

est_bothsexes <- ebb[Nat=="Eestlane" & YoA>2015, .N, by = c("PoR", "PoB")]
est_bothsexes[, YoA:="phase2"]
est_bothsexes_yoa <- rbind(est_bothsexes_yoa, data.frame(est_bothsexes))

est_bothsexes_yoa <- data.table(est_bothsexes_yoa)




pdf("figures/migration_matrix.pdf", width = 9, height = 8)

census <- fread("census2021_nat_sex_reg_N.csv", encoding = 'UTF-8')
census <- census[, c(2, 3, 5, 6)]
census <- census[, .(sum(`Males and females Estonian`), sum(`Males and females Russian`)), by = "Place of residence"]
census <- census[-1, ]
census[, `Place of residence` := unlist(strsplit(`Place of residence`, " COUNTY"))]
census$`Place of residence` <- str_to_title(census$`Place of residence`)
colnames(census) <- c("PoR", "Eestlane", "Venelane")
census_m <- melt(census, id.vars = "PoR")
census_m <- census_m[, c(2, 1, 3)]
colnames(census_m) <- c("Nat", "PoR", "N")

sampling_table <- data.frame()
sampling_stat <- plotMigrationMatrix(migr_ebb = nat_bothsexes_18andolder[Nat == "Eestlane", ], 
                                     census = as.data.table(census_m)[Nat == "Eestlane", ],
                                     title = "Estonians")
sampling_table <- rbind(sampling_table, sampling_stat)
colnames(sampling_table) <- c("Group", "N", "Proportion of the population")

sampling_stat <- plotMigrationMatrix(migr_ebb = nat_bothsexes_18andolder[Nat == "Venelane", ], 
                                     census = as.data.table(census_m)[Nat == "Venelane", ],
                                     title = "Russians")
sampling_table <- rbind(sampling_table, sampling_stat)

sampling_stat <- plotMigrationMatrix(migr_ebb = est_bothsexes_yoa[YoA=="phase1", ], 
                                     census = as.data.table(census_m)[Nat == "Eestlane", ],
                                     title = "Estonians until 2015")
sampling_table <- rbind(sampling_table, sampling_stat)
sampling_stat <- plotMigrationMatrix(migr_ebb = est_bothsexes_yoa[YoA=="phase2", ], 
                                     census = as.data.table(census_m)[Nat == "Eestlane", ],
                                     title = "Estonians since 2016")
sampling_table <- rbind(sampling_table, sampling_stat)




census <- fread("data/census2021_nat_sex_reg_N.csv", encoding = 'UTF-8')
census <- census[, c(2, 3, 8, 11)]
census <- census[, .(sum(`Males Estonian`), sum(`Females Estonian`)), by = "Place of residence"]
census <- census[-1, ]
census[, `Place of residence` := unlist(strsplit(`Place of residence`, " COUNTY"))]
census$`Place of residence` <- str_to_title(census$`Place of residence`)
colnames(census) <- c("PoR", "Male", "Female")
census_m <- melt(census, id.vars = "PoR")
census_m <- census_m[, c(2, 1, 3)]
colnames(census_m) <- c("Sex", "PoR", "N")

sampling_stat <- plotMigrationMatrix(migr_ebb = est_sex_18andolder[Sex == "1", ], 
                                      census = as.data.table(census_m)[Sex == "Male", ],
                                      title = "Estonian men")
sampling_table <- rbind(sampling_table, sampling_stat)
sampling_stat <- plotMigrationMatrix(migr_ebb = est_sex_18andolder[Sex == "2", ], 
                                      census = as.data.table(census_m)[Sex == "Female", ],
                                      title = "Estonian women")
sampling_table <- rbind(sampling_table, sampling_stat)



census <- fread("data/census2021_est_bothSexes_age_N.csv", encoding = 'UTF-8')
census <- census[, c(2, 3, 4)]
colnames(census) <- c("Age", "PoR", "N")
census_cast <- as.data.table(dcast(census, PoR~Age))
census_cast[, `18-24` := `0-4` + `5-9` + `10-14` + `15-19` + `20-24` - `0-17`]
census_cast[, `25-34` := `25-24` + `30-34`]
census_cast[, `35-64` := `35-39` + `40-44` + `45-49` + `50-54` + `55-59` + `60-64`]
census_cast[, `65+` := `65-69` + `70-74` + `75-79` + `80-84` + `85 and older`]
census_cast <- census_cast[, .(PoR, `18-24`, `25-34`, `35-64`, `65+`)]
census_cast[, `PoR` := unlist(strsplit(`PoR`, " COUNTY"))]
census_cast$`PoR` <- str_to_title(census_cast$`PoR`)
census_m <- melt(census_cast, id.vars = "PoR")
census_m <- census_m[, c(2, 1, 3)]
colnames(census_m) <- c("Age", "PoR", "N")

sampling_stat <- plotMigrationMatrix(migr_ebb = est_bothsexes_ages[Age == "18-24", ], 
                                      census = as.data.table(census_m)[Age == "18-24", ],
                                      title = "Estonian 18-24")
sampling_table <- rbind(sampling_table, sampling_stat)
sampling_stat <- plotMigrationMatrix(migr_ebb = est_bothsexes_ages[Age == "25-34", ], 
                                      census = as.data.table(census_m)[Age == "25-34", ],
                                      title = "Estonian 25-34")
sampling_table <- rbind(sampling_table, sampling_stat)
sampling_stat <- plotMigrationMatrix(migr_ebb = est_bothsexes_ages[Age == "35-64", ], 
                                      census = as.data.table(census_m)[Age == "35-64", ],
                                      title = "Estonian 35-64")
sampling_table <- rbind(sampling_table, sampling_stat)
sampling_stat <- plotMigrationMatrix(migr_ebb = est_bothsexes_ages[Age == "65+", ], 
                                      census = as.data.table(census_m)[Age == "65+", ],
                                      title = "Estonian 65+")
sampling_table <- rbind(sampling_table, sampling_stat)


dev.off()



pdf("figures/sample.pdf", width = 9, height = 8)
sampling_table$`Proportion of the population` <- as.numeric(sampling_table$`Proportion of the population`)
sampling_table$N <- as.numeric(sampling_table$N)
sampling_table_melt <- as.data.table(melt(sampling_table, id.vars = "Group"))
sampling_table_melt$Group <- factor(sampling_table_melt$Group, 
                                    levels = c("Estonians", "Russians", "Estonian men", "Estonian women",
                                               "Estonian 18-24", "Estonian 25-34", "Estonian 35-64", 
                                               "Estonian 65+", "Estonians until 2015",
                                               "Estonians since 2016"), 
                                    ordered = TRUE)
sampling_table_melt[, Type := ""]
sampling_table_melt[Group %in% c("Estonians", "Russians"), Type := "Ethnicity"]
sampling_table_melt[Group %in% c("Estonian men", "Estonian women"), Type := "Sex"]
sampling_table_melt[Group %in% c("Estonian 18-24", "Estonian 25-34", "Estonian 35-64", "Estonian 65+"), Type := "Age"]
sampling_table_melt[Group %in% c("Estonians until 2015", "Estonians since 2016"), Type := "Cohort"]


sampling_table_melt[variable=="N", variable:="Number"]
ggplot(sampling_table_melt, aes(x = Group, y = value, fill=Type)) +
  geom_bar(stat="identity", position="identity") + # , color = "#00BFC4"
  theme_minimal() + facet_grid(variable~., switch="y", scales='free') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.text.y = element_text(size=12),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=12),
        legend.title = element_text(size=15)) +
  labs(fill = "Divided by", y ="", x = "")

ggplot(ebb, aes(x = Age)) + geom_histogram(color="black", binwidth = 2) + theme_minimal() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15)) +
  labs(y = "Count")
ggplot(ebb, aes(x = YoA)) + geom_histogram(color="black", binwidth = 1) + theme_minimal() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15)) +
  labs(x = "Year of agreement", y = "Count")


dev.off()

