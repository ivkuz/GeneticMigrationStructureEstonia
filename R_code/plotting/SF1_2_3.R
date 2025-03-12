###########################################################################
# Script making the PDFs for Supplementary Note 1, describing the dataset #
# Supplementary Figures 1-3 ###############################################
###########################################################################


library(data.table)
library(stringr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(sf)
library(geos)
library(scales)


# Function that calculates representativeness based on census data
calcRepresentativeness <- function(migr_ebb, census, title){
  
  residence <- migr_ebb[, .(sum(N)), by = .(PoR)]
  colnames(residence) <- c("PoR", "N_ebb")
  pop_sample <- merge(residence, census)
  sampling_stat <- c(title, sum(pop_sample$N_ebb), round(sum(pop_sample$N_ebb)/sum(pop_sample$N), 2))
  
  return(sampling_stat)
  
}


option_list = list(
  make_option(c("-p", "--pheno"), type="character", default=NULL, 
              help="phenotype directory", metavar="character"),
  make_option(c("-u", "--nonrel"), type="character", default="non_relatives.tsv", 
              help="list of non-related individuals [default= %default]", metavar="character"),
  make_option(c("-m", "--map"), type="character", default="maakond_shp/maakond_20230201.shp", 
              help="map with counties' borders (.shp object) [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="results.tsv", 
              help="output directory name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)




###########################################################################
# Maps showing the representativeness of the main nationalities by region #
###########################################################################


# Read filtered data with phenotypes
ebb <- fread(paste0(opt$pheno, "/EstBB_filtered.tsv"), encoding = 'UTF-8')

# Download and process census data
census_reg_dat <- fread(paste0(opt$pheno, "/census2021_nat_sex_reg_N.csv"), encoding = 'UTF-8')
census_reg_dat <- census_reg_dat[, c(2, 3, 5, 6)]
census_reg_dat <- census_reg_dat[, .(sum(`Males and females Estonian`), sum(`Males and females Russian`)), by = "Place of residence"]
census_reg_dat <- census_reg_dat[-1, ]
census_reg_dat[, `Place of residence` := unlist(strsplit(`Place of residence`, " COUNTY"))]
census_reg_dat$`Place of residence` <- str_to_title(census_reg_dat$`Place of residence`)
colnames(census_reg_dat) <- c("PoR", "Eestlane", "Venelane")
census_reg_dat_m <- melt(census_reg_dat, id.vars = "PoR")
census_reg_dat_m <- census_reg_dat_m[, c(2, 1, 3)]
colnames(census_reg_dat_m) <- c("Nat", "PoR", "N_census")

# Number of individuals by Nationality and place of residence
ebb_reg_dat <- ebb[, .N, by = c("Nat", "PoR")]
colnames(ebb_reg_dat) <- c("Nat", "PoR", "N_ebb")

# Merge EstBB and census data
reg_dat <- merge(ebb_reg_dat, census_reg_dat_m, by = c("Nat", "PoR"))
reg_dat[, prop := round(N_ebb/N_census, 2)]
colnames(reg_dat)[c(3, 5)] <- c("N", "P")

# Read spatial data on Estonian regions
maakonds <-  st_read(opt$map)
maakonds$MNIMI <- gsub(x=maakonds$MNIMI, pattern = " maakond", replacement = "")
  
# Maps for Estonian representativeness
sub_res <- reg_dat[Nat=="Eestlane", ]
data_s <- merge(maakonds, sub_res, by.x="MNIMI", by.y="PoR")

pl_1 <- ggplot(data = data_s)+
  geom_sf(aes(geometry = geometry, fill = N))# +

pl_2 <- ggplot(data = data_s)+
  geom_sf(aes(geometry = geometry, fill = P)) #+

# Maps for Russian representativeness
sub_res <- reg_dat[Nat=="Venelane", ]
data_s <- merge(maakonds, sub_res, by.x="MNIMI", by.y="PoR")

pl_3 <- ggplot(data = data_s)+
  geom_sf(aes(geometry = geometry, fill = N)) #+

pl_4 <- ggplot(data = data_s)+
  geom_sf(aes(geometry = geometry, fill = P)) #+

plots <- list(pl_1, pl_2, pl_3, pl_4)

# Add the coordinates of the cities
city <- data.frame(
  name = c("Tartu", "Tallinn"),
  x_p = c(660578, 541078), 
  y_p = c(6471700, 6588700),
  x_t = c(689578, 500078), 
  y_t = c(6490700, 6601700)
)

# Add the universal parameters to the plots
for(i in 1:length(plots)){
  plots[[i]] <- plots[[i]] + 
    theme_classic() + scale_fill_viridis_c(option = "A") + 
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=8, hjust = 0.2, vjust = -0.8), 
          legend.key.width = unit(1, 'cm'),
          legend.key.height = unit(0.1, 'cm'),
          legend.position = "bottom",
          text = element_text(size = 7),
          legend.title = element_text(margin = margin(t = -10))) +
    geom_point(data=city, aes(x=x_p, y=y_p), shape=21, stroke=1.2) +
    geom_label(data=city, aes(x=x_t, y=y_t, label=name), 
               size=2.5, label.padding=unit(0.5, "mm"))
}

# Save a PDF with the plots
pdf(paste0(opt$out, "/SN1_1.pdf"), width=7/3*2, height=5)
print(grid.arrange(grobs = plots, layout_matrix = rbind(c(1, 2),
                                                        c(3, 4))))
grid.text("Estonian", x = 0.098, y = 0.975, gp = gpar(fontsize=14, fontface = "italic" ))
grid.text("A", x = 0.048, y = 0.925, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.548, y = 0.925, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("Russian", x = 0.098, y = 0.475, gp = gpar(fontsize=14, fontface = "italic"))
grid.text("C", x = 0.048, y = 0.425, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("D", x = 0.548, y = 0.425, gp = gpar(fontsize=14, fontface = "bold"))
dev.off()




############################################
# Histograms for Date of agreement and Age #
############################################

# Read filtered data with phenotypes
ebb <- fread(paste0(opt$pheno, "/EstBB_filtered.tsv"), encoding = 'UTF-8')

pl_1 <- ggplot(data = ebb, aes(x = YoA)) + geom_histogram(binwidth=1, colour='black',size=0.2) + theme_bw() +
  theme(text = element_text(size=10)) + xlab("Date of agreement") + ylab("Number")
pl_2 <- ggplot(data = ebb, aes(x = Age)) + geom_histogram(binwidth=2, colour='black',size=0.1) + theme_bw() +
  theme(text = element_text(size=10)) + xlab("Age") + ylab("Number")

plots <- list(pl_1, pl_2)
pdf(paste0(opt$out, "/SN1_2.pdf"), width=7, height=2.5)
print(grid.arrange(grobs = plots, layout_matrix = rbind(c(1, 2))))
grid.text("A", x = 0.018, y = 0.965, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("B", x = 0.518, y = 0.965, gp = gpar(fontsize=14, fontface = "bold"))
dev.off()




################################
# Representativeness by groups #
################################


# Read filtered data with phenotypes
ebb <- fread(paste0(opt$pheno, "/EstBB_filtered.tsv"), encoding = 'UTF-8')
# Unrelated Estonian subsample
non_rel <- fread(opt$nonrel)
ebb_unrel <- ebb[vkood %in% non_rel$vkood,]

# Calculate numbers of the individuals in each subcohort

# By sex
nat_bothsexes_18andolder <- ebb[, .N, by = c("Nat", "PoR", "PoB")]
est_sex_18andolder <- ebb[Nat=="Eestlane", .N, by = c("Sex", "PoR", "PoB")]

# Unrelated
nat_bothsexes_18andolder_unrel <- ebb_unrel[, .N, by = c("Nat", "PoR", "PoB")]

# By Age group
est_bothsexes <- ebb[Nat=="Eestlane" & Age<=24, .N, by = c("PoR", "PoB")]
est_bothsexes[, Age:="18-24"]
est_bothsexes_ages <- data.frame(est_bothsexes)

est_bothsexes <- ebb[Nat=="Eestlane" & Age>=25 & Age<=48, .N, by = c("PoR", "PoB")]
est_bothsexes[, Age:="25-48"]
est_bothsexes_ages <- rbind(est_bothsexes_ages, data.frame(est_bothsexes))

est_bothsexes <- ebb[Nat=="Eestlane" & Age>=49 & Age<=64, .N, by = c("PoR", "PoB")]
est_bothsexes[, Age:="49-64"]
est_bothsexes_ages <- rbind(est_bothsexes_ages, data.frame(est_bothsexes))

est_bothsexes <- ebb[Nat=="Eestlane" & Age>=65, .N, by = c("PoR", "PoB")]
est_bothsexes[, Age:="65+"]
est_bothsexes_ages <- rbind(est_bothsexes_ages, data.frame(est_bothsexes))

est_bothsexes_ages <- data.table(est_bothsexes_ages)

# By recruitment period (phase)
est_bothsexes <- ebb[Nat=="Eestlane" & YoA<=2016, .N, by = c("PoR", "PoB")]
est_bothsexes[, YoA:="phase1"]
est_bothsexes_yoa <- data.frame(est_bothsexes)

est_bothsexes <- ebb[Nat=="Eestlane" & YoA>2016, .N, by = c("PoR", "PoB")]
est_bothsexes[, YoA:="phase2"]
est_bothsexes_yoa <- rbind(est_bothsexes_yoa, data.frame(est_bothsexes))

est_bothsexes_yoa <- data.table(est_bothsexes_yoa)


# Use census data on sex and nationality
census <- fread(paste0(opt$pheno, "/census2021_nat_sex_reg_N.csv"), encoding = 'UTF-8')
census <- census[, c(2, 3, 5, 6)]
census <- census[, .(sum(`Males and females Estonian`), sum(`Males and females Russian`)), by = "Place of residence"]
census <- census[-1, ]
census[, `Place of residence` := unlist(strsplit(`Place of residence`, " COUNTY"))]
census$`Place of residence` <- str_to_title(census$`Place of residence`)
colnames(census) <- c("PoR", "Eestlane", "Venelane")
census_m <- melt(census, id.vars = "PoR")
census_m <- census_m[, c(2, 1, 3)]
colnames(census_m) <- c("Nat", "PoR", "N")

# Calculate representativeness for nationalities and phases
sampling_table <- data.frame()
sampling_stat <- calcRepresentativeness(migr_ebb = nat_bothsexes_18andolder[Nat == "Eestlane", ], 
                                     census = as.data.table(census_m)[Nat == "Eestlane", ],
                                     title = "Estonian")
sampling_table <- rbind(sampling_table, sampling_stat)
colnames(sampling_table) <- c("Group", "N", "Proportion of the population")

sampling_stat <- calcRepresentativeness(migr_ebb = nat_bothsexes_18andolder_unrel[Nat == "Eestlane", ], 
                                     census = as.data.table(census_m)[Nat == "Eestlane", ],
                                     title = "Unrelated Estonian")
sampling_table <- rbind(sampling_table, sampling_stat)
colnames(sampling_table) <- c("Group", "N", "Proportion of the population")

sampling_stat <- calcRepresentativeness(migr_ebb = nat_bothsexes_18andolder[Nat == "Venelane", ], 
                                     census = as.data.table(census_m)[Nat == "Venelane", ],
                                     title = "Russian")
sampling_table <- rbind(sampling_table, sampling_stat)

sampling_stat <- calcRepresentativeness(migr_ebb = est_bothsexes_yoa[YoA=="phase1", ], 
                                     census = as.data.table(census_m)[Nat == "Eestlane", ],
                                     title = "2001-2016")
sampling_table <- rbind(sampling_table, sampling_stat)
sampling_stat <- calcRepresentativeness(migr_ebb = est_bothsexes_yoa[YoA=="phase2", ], 
                                     census = as.data.table(census_m)[Nat == "Eestlane", ],
                                     title = "2017-2021")
sampling_table <- rbind(sampling_table, sampling_stat)


# Use census data on sex and nationality
census <- fread(paste0(opt$pheno, "/census2021_nat_sex_reg_N.csv"), encoding = 'UTF-8')
census <- census[, c(2, 3, 8, 11)]
census <- census[, .(sum(`Males Estonian`), sum(`Females Estonian`)), by = "Place of residence"]
census <- census[-1, ]
census[, `Place of residence` := unlist(strsplit(`Place of residence`, " COUNTY"))]
census$`Place of residence` <- str_to_title(census$`Place of residence`)
colnames(census) <- c("PoR", "Male", "Female")
census_m <- melt(census, id.vars = "PoR")
census_m <- census_m[, c(2, 1, 3)]
colnames(census_m) <- c("Sex", "PoR", "N")

# Calculate representativeness for sexes
sampling_stat <- calcRepresentativeness(migr_ebb = est_sex_18andolder[Sex == "1", ], 
                                     census = as.data.table(census_m)[Sex == "Male", ],
                                     title = "Male")
sampling_table <- rbind(sampling_table, sampling_stat)
sampling_stat <- calcRepresentativeness(migr_ebb = est_sex_18andolder[Sex == "2", ], 
                                     census = as.data.table(census_m)[Sex == "Female", ],
                                     title = "Female")
sampling_table <- rbind(sampling_table, sampling_stat)


# Use census data on age in Estonians
census <- fread(paste0(opt$pheno, "/census2021_est_bothSexes_age_N.csv"), encoding = 'UTF-8')
census <- census[, c(2, 3, 4)]
colnames(census) <- c("Age", "PoR", "N")
census_cast <- as.data.table(dcast(census, PoR~Age))
census_cast[, `18-24` := `0-4` + `5-9` + `10-14` + `15-19` + `20-24` - `0-17`]
census_cast[, `25-48` := `25-24` + `30-34` + `35-39`  + `40-44` + `45-49`]
census_cast[, `49-64` := `50-54` + `55-59` + `60-64`]
census_cast[, `65+` := `65-69` + `70-74` + `75-79` + `80-84` + `85 and older`]
census_cast <- census_cast[, .(PoR, `18-24`, `25-48`, `49-64`, `65+`)]
census_cast[, `PoR` := unlist(strsplit(`PoR`, " COUNTY"))]
census_cast$`PoR` <- str_to_title(census_cast$`PoR`)
census_m <- melt(census_cast, id.vars = "PoR")
census_m <- census_m[, c(2, 1, 3)]
colnames(census_m) <- c("Age", "PoR", "N")

# Calculate representativeness for age groups
sampling_stat <- calcRepresentativeness(migr_ebb = est_bothsexes_ages[Age == "18-24", ], 
                                     census = as.data.table(census_m)[Age == "18-24", ],
                                     title = "18-24")
sampling_table <- rbind(sampling_table, sampling_stat)
sampling_stat <- calcRepresentativeness(migr_ebb = est_bothsexes_ages[Age == "25-48", ], 
                                     census = as.data.table(census_m)[Age == "25-48", ],
                                     title = "25-48")
sampling_table <- rbind(sampling_table, sampling_stat)
sampling_stat <- calcRepresentativeness(migr_ebb = est_bothsexes_ages[Age == "49-64", ], 
                                     census = as.data.table(census_m)[Age == "49-64", ],
                                     title = "49-64")
sampling_table <- rbind(sampling_table, sampling_stat)
sampling_stat <- calcRepresentativeness(migr_ebb = est_bothsexes_ages[Age == "65+", ], 
                                     census = as.data.table(census_m)[Age == "65+", ],
                                     title = "65+")
sampling_table <- rbind(sampling_table, sampling_stat)


# Prepare data for visualisation
sampling_table$`Proportion of the population` <- as.numeric(sampling_table$`Proportion of the population`)
sampling_table$N <- as.numeric(sampling_table$N)
sampling_table_melt <- as.data.table(melt(sampling_table, id.vars = "Group"))
sampling_table_melt$Group <- factor(sampling_table_melt$Group, 
                                    levels = c("Estonian", "Unrelated Estonian", "Russian", "Male", "Female",
                                               "18-24", "25-48", "49-64", 
                                               "65+", "2001-2016", "2017-2021"), 
                                    ordered = TRUE)
sampling_table_melt[, Type := ""]
sampling_table_melt[Group %in% c("Estonian", "Unrelated Estonian", "Russian"), Type := "Ethnicity"]
sampling_table_melt[Group %in% c("Male", "Female"), Type := "Sex"]
sampling_table_melt[Group %in% c("18-24", "25-48", "49-64", "65+"), Type := "Age"]
sampling_table_melt[Group %in% c("2001-2016", "2017-2021"), Type := "Agreement"]
sampling_table_melt$Type <- factor(sampling_table_melt$Type, 
                                    levels = c("Ethnicity", "Sex", "Age", "Agreement"), 
                                    ordered = TRUE)

# Make the plots
sampling_table_melt[variable=="N", variable:="Number"]
sampling_table_melt[, variable := factor(variable, levels = c("Number", "Proportion of the population group"))]
pdf("figures/SN1_3.pdf", width=7, height=5)
ggplot(sampling_table_melt, aes(x = Group, y = value, fill=Type)) +
  geom_bar(stat="identity", position="identity", colour="black", linewidth=0.2) + # , color = "#00BFC4"
  theme_minimal() + facet_grid(variable~., switch="y", scales='free') +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(fill = "Divided by", y ="", x = "")

dev.off()

