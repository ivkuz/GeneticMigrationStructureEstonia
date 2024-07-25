# GeneticMigrationStructureEstonia

These custom R code has been used in the study 
"Assessing the impact of 20th century internal migrations on the genetic structure of Estonia"
Ivan A. Kuznetsov, Estonian Biobank Research Team, Mait Metspalu, Uku Vainik, Luca Pagani, Francesco Montinaro, Vasili Pankratov
bioRxiv 2023.10.25.564036; doi: https://doi.org/10.1101/2023.10.25.564036

There are three folders with the scripts. You can find a short description of what each of the scripts does below:


./data_preparation/:

	sampleFiltering.R	- Script for removing duplicated genetic data and filtering by nationality;

	relatives_filtering.R	- Script for creating the list of unrelated Estonian individuals;
	

./analysis/:

	PCs_ditribution.R	- Script calculating the mean and sd for 100 PCs by region for all the subcohorts;

	PRSs_distribution.R	- Script calculating the mean and sd for PGSs by region for all the subcohorts;

	PRS_correlations.R	- Script calculating the correlations between polygenic scores from population-based GWAS data;

	PRS_correlations_sib.R	- Script calculating the correlations between polygenic scores from sibship-based GWAS data;

	PRSs_distribution_and_cor_family_sib.R	- Script calculating the mean and sd for PGSs by region for all the subcohorts
						in sibling design (deviation of individual's PGS from the sibship mean PGS),
						and calculating the correlations between polygenic scores in sibling design;

	resultsAnalysis_calculate_fstat.R	- Script calculating the proportion of variance explained by region 
						and related statistics for PGSs, PCs etc.,
						It uses output file from "PCs/PRSs_distribution.R" as its input;


./plotting/:

	plotSN1.R		- Script making the PDFs for Supplementary Note 1, describing the dataset
				(Supplementary Figures 1-3)

	plotMaps.R		- Script making the PDFs with maps - mean PGS and EA in the regions
				(Figure 3; Supplementary Figures 14-25);

	plotMapTartu-Tallinn.R	- Script making the PDFs with maps - differences in mean PGS and EA between migrants 
				to Tallinn or Tartu in the regions (Supplementary Figures 32-37);

	plotMigrationBars_final.R	- Script making the PDFs with barplots - mean PGS and EA in the migration groups
					(Figure 4; Supplementary Figures 26-31, 38-63),
					It also makes PDF with plot with mean PRS in the groups of 
					born in or migrated to the cities,
					It also makes PDF with plot with PRS distributions 
					in different migration groups for Q&A;

	plotMigrationBars_family_sib.R	- Script making the PDF with barplots for sibling design PGS;
