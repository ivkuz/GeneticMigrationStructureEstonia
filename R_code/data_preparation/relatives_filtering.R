library(data.table)

####################################################################
### Script for creating the list of unrelated Estonian individuals #
####################################################################


option_list = list(
  make_option(c("-f", "--pheno"), type="character", default="EstBB_filtered.tsv", 
              help="file with phenotypes for filtered individuals [default= %default]", metavar="character"),
  make_option(c("-k", "--kin"), type="character", default="king.kin0", 
              help="kinship table from KING [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="~/tmp/", 
              help="output directory name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Read kinship and phenotype data
king <- opt$kin
estbb_filtered <- fread(opt$pheno)
rel <- fread(king)
individuals <- estbb_filtered[Nat=="Eestlane", vkood]
# Keep individuals present in the phenotype table
rel <- rel[ID1 %in% individuals & ID2 %in% individuals, .(ID1, ID2)]

# Exclude one by one individuals with the largest number of relatives
exclude <- c()
n=2
while(n>1){
  # Calculate the number of relatives for each individual having at least one relative
  rel_tab <- as.data.table(table(c(rel$ID1, rel$ID2)))
  # Find the individual with the largest number of relatives
  rel_tab <- rel_tab[order(N, decreasing=T),]
  n <- rel_tab[1, N]
  hub <- rel_tab[N == n, V1][1]
  # Record this individual to the exclusion list
  exclude <- c(exclude, hub)
  # Exclude from the table of relatives
  rel <- rel[!(ID1 %in% hub) & !(ID2 %in% hub),]
  
}
exclude <- c(exclude, rel$ID1)

# Save the list of the individials excleded
exclude_table <- data.table(vkood=exclude)
write.table(exclude_table, paste0(opt$out, "/excluded_relatives_est_perfect.tsv"), quote = F, row.names = F, sep = "\t")

# Save the list of the individials kept
non_relatives <- data.table(vkood = estbb_filtered[Nat=="Eestlane" & !(vkood %in% exclude), vkood])
write.table(non_relatives, paste0(opt$out, "/non_relatives_est_perfect.tsv"), quote = F, row.names = F, sep = "\t")