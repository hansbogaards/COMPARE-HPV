# HPV-COMPARE
R-code to evaluate the cost-effectiveness profile of 9v vs 2v HPV vaccination in the Netherlands. \
Used in 'Health and economic effects of nonavalent versus bivalent HPV vaccination in the Netherlands: a data-driven analysis' (link to paper).

# .R files
HPVcompare_functions.R \
HPVcompare_RiskReduction_functions.R \
HPVcompare_ParametersData.R \
HPVcompare_RiskReductions.R \
HPVcompare_Computation.R 

# Data
The Inc_ files contain cancer incidence data from the Netherlands Cancer Registry for the various cancer types \
The Surv_ files contain cancer survival data from the Netherlands Cancer Registry for the various cancer types \
The TablePop_ files contain overall survival data of the general population in 2020 from the Statistics Netherlands database \
The Pop files contain data on the population size stratified by age group from the Statistics Netherlands database \
The BirthsPer1000_ file contains data on the average number of births in the population (per age of the mother) from the Statistics Netherlands database

# Example
HPVcompare_Example.R contains code with which an example scenario can be run:
- It first loads all the required functions by loading HPVcompare_functions.R and HPVcompare_RiskReduction_functions.R
- It inputs all the data and simulates the parameters by loading HPVcompare_ParametersData.R
- It defines the example scenario
- It computes the corresponding risk reductions by loading HPVcompare_RiskReductions.R
- It runs the analysis by loading HPVcompare_Computation.R
