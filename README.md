# tomato_N_timecourse


## Overview

This repository contains the R code and data associated with the following manuscript:

Dixon M M, Rohrbaugh C R, Delgado J A, Manter D K, Vivanco J M. Rhizosphere microbiome stabilization following initial fluctuations is delayed with nitrogen additions in tomato seedlings.

In this article, we elucide changes occurring tomato-root-associated microbial assembly throughout the duration of the vegetative development stage and how this development is disrupted with additions of nitrogen fertilizers.


## File descriptions

1. "ReadInData.Rmd" includes code to read in data for downstream microbiome analysis. It also includes codes to convert EMU files for indivial sequencing runs and compiles them into one phyloseq object. 
  - Function inputs:
    - "DManter_CustomFunctions.R": Includes custom functions that allow for conversion from EMU to phyloseq ojects.
	- Data inputs:  
		- "Carley_F_clean.xls": Includes sequence reads for plate A.  
		- "Carley_G_clean.xlsx": Includes sequence reads for plate B.  
		- "EMU_rel.abund_A..csv": Taxonomic relative abundance information for plate A.
		- "EMU_rel.abund_B.csv.csv": Taxonomic relative abundance information for plate B.
	- Outputs:  
		- "P.rel.RDS" merged phyloseq object that combines information for plates A and B.
		- "P.count.RDS" A phyloseq object that includes count data for microbial abundance, combined for all plates.
    
2. "WeeklyDiffAbund.Rmd" includes code to determine which bacteria change in abundance week-to-week.
  - Data inputs:  
		- "P.count.RDS": phyloseq object showing bacterial abundance count data (culled from code in "ReadInData.Rmd")
    - "EMU_database.GIBBs.KO.PICRUST2.xlsx" Shows information from PICRUST2-mapped EMU database that shows functional abundance of the mapped bacteria.
	- Outputs:  
		- "DA_CK_1v2.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between unfertilized weeks 1-2. 
		- "DA_CK_2v3.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between unfertilized weeks 2-3. 
		- "DA_CK_3v4.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between unfertilized weeks 3-4. 
		- "DA_CK_4v5.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between unfertilized weeks 4-5. 
		- "DA_CK_5v6.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between unfertilized weeks 5-6. 
		- "DA_CK_7v8.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between unfertilized weeks 7-8. 
		- "DA_LU_1v2.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between low urea weeks 1-2. 
		- "DA_LU_2v3.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between low urea weeks 2-3. 
		- "DA_LU_3v4.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between low urea weeks 3-4. 
		- "DA_LU_4v5.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between low urea weeks 4-5. 
		- "DA_LU_5v6.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between low urea weeks 5-6. 
		- "DA_LU_6v7.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between low urea weeks 6-7. 
		- "DA_LU_7v8.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between low urea weeks 7-8. 
		- "DA_HU_1v2.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between high urea weeks 1-2. 
		- "DA_HU_2v3.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between high urea weeks 2-3. 
		- "DA_HU_3v4.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between high urea weeks 3-4. 
		- "DA_HU_4v5.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between high urea weeks 4-5. 
		- "DA_HU_5v6.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between high urea weeks 5-6. 
		- "DA_HU_7v8.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between high urea weeks 7-8. 
   	- "DA_ESN_1v2.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between ESN weeks 1-2. 
		- "DA_ESN_2v3.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between ESN weeks 2-3. 
		- "DA_ESN_3v4.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between ESN weeks 3-4. 
		- "DA_ESN_4v5.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between ESN weeks 4-5. 
		- "DA_ESN_5v6.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between ESN weeks 5-6. 
		- "DA_ESN_6v7.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between ESN weeks 6-7. 
		- "DA_ESN_7v8.csv" Intermediate data file containing results of DA analysis showing the log fold change in bacteria abundance between ESN weeks 7-8.
    - "DA_all.csv" combined data file containing results from all DA weekly comparisons.
    - Results corresponding to Supplemental Table 1
    - Results corresponding to Supplemental Table 2
    - Results corresponding to Supplemental Table 3
    - Results corresponding to Supplemental Table 4
		- Results corresponding to Figure 1
    - "DA-fert-byweek.png" Figure 1, differential abundance of bacteria week-to-week for each fertilizer treatment.
    - "DA_early-late.csv" Results of the DA between early (weeks 1-4) and late (weeks 5-8).
    - Results corresponding to Figure 3b
    - "DA-Volcano-early-late.png" Figure 3b, a volcano plot showing the amount of bacteria that are differentially abundant between early (weeks 1-4) and late (weeks 5-8).
    - Results corresponding to Figure 2
    -  "DA-functions-byweek.png" Figure 2, a volcono plot showing the variation in functions of bacteria that are differentially abundance week-to-week.

3. "DA_fertilization.Rmd" includes code to analyze the differential abundance of bacteria between the fertilizer treatments.
  - Data inputs:  
		- "P.count.RDS": phyloseq object showing bacterial abundance count data (culled from code in "ReadInData.Rmd")
    - "EMU_database.GIBBs.KO.PICRUST2.xlsx" Shows information from PICRUST2-mapped EMU database that shows functional abundance of the mapped bacteria.
  - Outputs:
    - "_fertilizationDA.csv" results from DA analysis showing the bacteria that vary between each fertilizer treatment.
    - Results corresponding to Table 2
      
4. "dbRDA.R" includes code to assess rhizosphere community structure for samples. 
	- Data inputs:  
		- "P.rel.RDS": merged phyloseq object showing relative abundance files (generated from "ReadInData.Rmd")
	- Outputs:  
		- "rdabiplot.png" Figure 3a, Distance-based redundancy analysis (db-RDA) showing clustering based on Bray-Curtis dissimilarity of the bacterial community structure in the tomato rhizosphere. 
		- Results corresponding to  figure 3a
    - Results corresponding to Table 1
		- "pcoa.png" Supplemental Figure 2, Principal coordinate analysis showing unconstrained clustering of the bacterial community
		- Results corresponding to Supplemental Figure 2
   
5. "biomass.Rmd" includes code to analyze differences in relative root and shoot biomass accumulation. It also includes the analysis for the P stress factor. 
	- Data inputs:  
		- "Dry_Biomass.csv": Raw dry shoot and root biomass values.       
	- Data outputs:  
		- "BiomassShoot.png": Supplemental Figure 1a, Boxplot of shoot biomass over time for different fertilizer treatments a 
		- Results corresponding to Supplemental Figure 1a
		- "BiomassRoot.png": Supplemental Figure 1b, Boxplot of root biomass over time for different fertilizer treatments a 
		- Results corresponding to Supplemental Figure 1b
   
6. "betadisper.Rmd" includes code to determine differences beta dispersion as a function of weekly development in tomato.  
	- Data inputs:  
		- "P.rel.RDS": merged phyloseq object showing relative abundance files (generated from "ReadInData.Rmd")
	- Outputs: 
		- "betadisper.csv" Dataframe showing distances to centroid for each sample
    - "betadisper.png" Supplemental Figure 3, Boxplot of beta dispersion as a function of week and fertilizer
		- Results corresponding to Supplemental Figure 3
