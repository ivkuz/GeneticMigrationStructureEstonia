##############################################################################
### Script for removing duplicated genetic data and filtering by nationality #
##############################################################################


Sys.setlocale("LC_CTYPE", "estonian")

library(data.table)
library(reshape2)
library(stringr)
library(optparse)


# Transfor nationality field to Eestlane (Estonian), Venelane (Russian) and Other
transformNat <- function(Nat_vector){
  nat <- strsplit(as.character(Nat_vector), "(;| )+")
  nat_names <- c()
  for(i in nat){
    nationalities <- c()
    for(ncode in i){
      nationalities <- c(nationalities, ncode)
    }
    if(length(unique(nationalities)) == 1){
      if(unique(nationalities) %in% c("Eestlane", "Venelane")){
        nat_names <- c(nat_names, unique(nationalities))
      } else {
        nat_names <- c(nat_names, "Other")
      }
    } else if(length(unique(nationalities)) == 0){
      nat_names <- c(nat_names, "")
    } else{
      nat_names <- c(nat_names, "Other")
    }
  }
  return(nat_names)
}


option_list = list(
  make_option(c("-c", "--codes"), type="character", default="svcodes.tsv", 
              help="table with scodes and vcodes [default= %default]", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default="EstBB_chr21.psam", 
              help=".psam file with the vcodes with imputed genomes [default= %default]", metavar="character"),
  make_option(c("-a", "--ancestry"), type="character", default="ancestry.txt", 
              help="file with inferred ancestry prifiles [default= %default]", metavar="character"),
  make_option(c("-p", "--pheno"), type="character", default="query1.tsv", 
              help="phenotype table [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="EstBB_filtered.tsv", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# load table with skoods and vkoods
codes <- fread(opt$codes)

# load information on samples with imputed genomes (vkoods)
geno <- fread(opt$sample)

# Exclude non-Europeans (restricted to Nothern Europe)
anc <- fread(opt$ancestry)
keep_anc <- c("Europe_(East)", "Europe_(North_West)", "Finland")
anc <- anc[Assigned_group %in% keep_anc, ]

# Extract list of vkoods with genomes, one vkood per skood
codes <- codes[`Sample vkood` %in% geno$`#IID` &
                 `Sample vkood` %in% anc$sample.ID, ]
codes <- unique(codes, by = "Person skood")

# load file with phenotypes and keep only skoods with genomes
ebb <- fread(opt$pheno)
ebb <- ebb[`Person skood` %in% codes$`Person skood`]
ebb <- merge(codes, ebb, by = "Person skood")

# Extract Estonian and Russian subpopulations
ebb <- ebb[`PersonPortrait nationality name` %in% c("Eestlane", "Venelane")]
ebb[, NatQuest := transformNat(`CONCATSTR(Nationality nationality name)`)]
ebb <- ebb[`PersonPortrait nationality name` == NatQuest | NatQuest == "",]

# Make sure there are no age errors (Should not have participants younger than 18)
ebb <- ebb[`Person ageAtAgreement` >= 18,]

# Calculate age of death for dead individuals
ebb[, lastAge := 2022 - `Person birthYear`]
ebb[!is.na(`Person deathYear`), lastAge := `Person deathYear` - `Person birthYear`]

# Select phenotypes to keep further
selected_columns <- c("Sample vkood", "Person skood",
                      "Person birthYear", "Person gender code",
                      "PersonLocation birthCounty name",
                      "PersonLocation residencyCounty name",
                      "PersonPortrait nationality name",
                      "lastAge", "Person ageAtAgreement")

ebb <- ebb[, ..selected_columns]
colnames(ebb) <- c("vkood", "skood", "YoB", "Sex", "PoB", "PoR", "Nat", "Age", "AgeAtAgr")
ebb[, YoA := YoB+AgeAtAgr]

# Delete individuals with no data on place of birth or residence
ebb <- ebb[!(PoB=="" | PoR==""),]

# Make counties' names shorter
ebb[, PoR := unlist(strsplit(PoR, " maakond"))]
ebb[, PoB := unlist(strsplit(PoB, " maakond"))]

# Saved filtered data
write.table(ebb, opt$out, row.names = F, quote = F, sep = "\t")

