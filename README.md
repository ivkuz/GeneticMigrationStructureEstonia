# GeneticMigrationStructureEstonia
Archived version: https://doi.org/10.5281/zenodo.17497533

This custom R code has been used in the study 
"Assessing the impact of 20th century internal migrations on the genetic structure of Estonia"
Ivan A. Kuznetsov, Estonian Biobank Research Team, Mait Metspalu, Uku Vainik, Luca Pagani, Francesco Montinaro, Vasili Pankratov
bioRxiv 2023.10.25.564036; doi: https://doi.org/10.1101/2023.10.25.564036

There are three folders with the scripts in the "R_code" directory. 
Find a short description of what each of the scripts does below:


./data_preparation/:

	sampleFiltering.R	- Script for removing duplicated genetic data and filtering by nationality;

	relatives_filtering.R	- Script for creating the list of unrelated Estonian individuals;
	

./analysis/:

	PCs_ditribution.R	- Script calculating the mean and sd for 100 PCs by region for all the subcohorts;

	PRSs_distribution.R	- Script calculating the mean and sd for PGSs by region for all the subcohorts
 				and for PGSs adjusted for PGS for EA;

	PRSs_distribution_family_sib.R		- Script calculating the mean and sd for PGSs by region for all the subcohorts
						in sibling design (deviation of individual's PGS from the sibship mean PGS);

	PRSs_distribution_family_sib_random.R		- Script estimating the effects of PGSs and PCs on migration from ORE in siblings;

 	migrEffects.R	- Script estimating the effects of PGSs on migration from ORE in unrelated individuals;

  	migrEffects_sib.R	- Script estimating the effects of PGSs and PCs on migration from ORE in siblings;

	resultsAnalysis_calculate_fstat.R	- Script calculating the proportion of variance explained by region 
						and related statistics for PGSs, PCs etc.,
						It uses output file from "PCs/PRSs_distribution.R" as its input;


./plotting/:

	figure1_2.R	- Figure 1. Fraction of the inter-individual variance of (A) PCs, (B) PGSs and 
 			(С) PGSs additionally adjusted for PGSEA, explained by county of birth (POB) and county of residence (POR);
    			Figure 2. Fraction of the inter-individual variance of the deviation of individual’s value from the sibship’s mean for
       			(A) PCs, (B) PGSs, and (C) PGSs additionally adjusted for PGSEA4, explained by county of birth (POB) and county of residence (POR);
	  		Supplementary Figures with EA4 PGS and with unrelated subsample;

   	figure3.R	- Figure 3. PGSEA (PGSEA4) landscape in Estonia
    			(maps with mean PGS per county and forest plots for the difference POR-POB in population and in siblings);

	figure4_sib_panel.R	- Prepares .rds file with Figure 4C, input for figure4.R;

 	figure4_7.R	- Figure 4. PGSEA (PGSEA4) in migration groups defined by combination of place of birth (POB) and residence (POR); 
  			And Supplementary Figures with bar plots for migration groups;
     			Figure 7. Difference in average PGSEA between cities (Tallinn and Tartu combined) and ORE across birth year bins;
			Also makes PDF with plot with PRS distributions in different migration groups for SM;

   	figure5.R	- Figure 5. The difference in mean PGSEA and EA (years of education) between residents of Tallinn and Tartu City by county of birth;
     
	figure6.R	- Figure 6. PGSs as predictors of migration from ORE to the major cities (Tallinn or Tartu);
 			And the corresponding Supplementary Figures;
  
 	SF1_2_3.R	- Supplementary Figure 1. The distribution of EstBB participants after filtering by (A) date of agreement (year of recruitment) and (B) age;
  			Supplementary Figure 2. Geographic distribution of EstBB participants by county of residence;
     			Supplementary Figure 3. Absolute number of participants in each group (top) and the same, 
			normalized by the size of the corresponding group in the general population (bottom);

	SF5.R	- Supplementary Figure 5. Internal migration in Estonia;

 	SF29_AM.R	- Supplementary Figure 29. Mating by proximity and assortative mating in Estonia;
