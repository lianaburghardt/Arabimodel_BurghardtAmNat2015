This file explains the files submitted to Dryad to support the paper "Modeling the influence of genetic and environmental variation on the expression of plant life cycles across landscapes"

##################Environmental data folder ######################
Files:
environmental_data_Norwich_60yrs.csv
environmental_data_Valencia_60yrs.csv
environmental_data_Halle_60yrs.csv
environmental_data_Oulu_60yrs.csv

These are environmental data files used to run the integrated life cycle model for Arabidopsis thaliana. They provides simulated moisture, temperature, day length, and daylight fraction for each hour for 60 years in each of four European locations (Valencia, Spain; Norwich, England; Halle, Germany; and Oulu, Finland (525600 x 5 matrix for each location). These data files are randomly assembled from 20 possible simulated years of temperature data paired with stochastically generated of moisture data and a single photoperiod profile as described in Appendix B of the associated paper. These environmental series were used to generate Figures 3, 4, and 5. 

Columns-
year - which of the 20 years of temperature data the temperatures came from
day_length - the hours of sunlight that particular day
day_hours - the fraction of that particular hour that it is light outside (values are between 0 and 1)
temperature - simulated average temperatures each hour in degrees celsius 
moisture - simulated soul moisture level each hour in mPa

Rows- each reflects a different sequential hour of weather data

#################Source Folder###################################
File 1: Arabi_Base_final.R

This file contains all the functions necessary to cycle through the entire A. thaliana life cycle for multiple years and to create the summarized output of results provided below for each location and genotype combination. The code is annotated so details on the specific function can be found there.

File 2: Simulation_Code.R
This script contains code to run the 24 main simulations highlighted in the paper and create the summarized output available below for each genotype x environment combination. You will need to source in the Arabi_base_final.R code to run the actual models. 


##################Results folder######################
All of the summarized results span from year 15 of the simulation to year 60 of the simulation. X in the file names below refers to whatever parameter level was used for Fi (.598 or .737) and psi_mean (0, 1.25, or 2.5). The first part of the file name tells you what environmental data file above was used to run the simulation 

Files:
Norwich60year_FiX_psimnX.RData  (6 files, one for each genotype)
Valencia60year_FiX_psimnX.RData  (6 files, one for each genotype)
Halle60year_FiX_psimnX.RData  (6 files, one for each genotype)
Oulu60year_FiX_psimnX.RData  (6 files, one for each genotype)

When you load the data into R you will find that each list is named after the location where the simulation occurred and the genotype that was simulated. This is so you can load multiple lists in at the same time. Below are the list names for each file. You can also view these through the ls() command after you load the Rdata

List name for each file once loaded into R:
Low Floral Rep (Fi=.598)/ Low Dormancy (psi_mean=0):	Valencia.life	, Halle.life, Oulu.life, Norwich.life	
Low Floral Rep (Fi=.598)/ Med Dormancy (psi_mean=1.25): Valencia.life.halfDor, Halle.life.halfDor, Oulu.life.halfDor, Norwich.life.halfDor
Low Floral Rep (Fi=.598)/ High Dormancy (psi_mean=2.5): Valencia.life.Dor	, Halle.life.Dor, Oulu.life.Dor, Norwich.life.Dor
Low Floral Rep (Fi=.737)/ LowDormancy (psi_mean=0):Valencia.life.FRI	, Halle.life.FRI, Oulu.life.FRI, Norwich.life.FRI	
Low Floral Rep (Fi=.737)/ MedDormancy (psi_mean=1.25):Valencia.life.FRI.halfDor, Halle.life.FRI.halfDor, Oulu.life.FRI.halfDor, Norwich.life.FRI.halfDor
Low Floral Rep (Fi=.737)/ High Dormancy (psi_mean=2.5): Valencia.life.FRI.Dor	, Halle.life.FRI.Dor, Oulu.life.FRI.Dor, Norwich.life.FRI.Dor

General data frames included in each list:
disp.mat - matrix of the dispersal day for each seed class (rows) x seed cohort (columns)
dispersalday -vector the length of the number of cohorts giving the dispersal day for each cohort
totperiod - matrix of the generation length of for each seed class (rows) x seed cohort (columns)
seedperiod - matrix of the length of time spent as a seed for each seed class (rows) x seed cohort (columns) 
veggieperiod - matrix of the length of time spent as a rosette for each seed class (rows) x seed cohort (columns) 
floweringperiod - matrix of the length of time spent flowering for each seed class (rows) x seed cohort (columns) 
nonseedperiod -matrix of the length of time spent as a rosette or flowering for each seed class (rows) x seed cohort (columns) 
cohort.size - vector of length the # of seed cohorts indicating the # of individuals in each cohort.
seed.cohort.classes - # of seed classes used in the run	
n.seedcohort - matrix of # of individuals in each seed class (rows) x seed cohort (columns)	
months - indication of which month each dispersal event happened in 

Specialized data frames included in each list:
Note: the X's below stands in for various life stages (flower, seed, vegetative) or transition (germination, flowering, seed dispersal). Tot is short for total and refers to a value measured across the entire life cycle of a plant such as total life-cycle length

- General results for stage length and cohort sizes
avg.X - average (generation, seed stage, vegetative stage, or flowering stage) length for all individuals across the whole simulation
period.all - data frame holding all the life-stage lengths for each individual
Xperiod.all- vector holding the life-stage lengths for each individual in the simulation
Xperiod.cohort - same as above but averaged across cohorts so there is an average life stage length for each cohort not for each individual.
month.cohortsize.X - vector of cohort sizes averaged per month
sumbymonth.X - vector of average stage length (or generation length) for each seed class summarized for each month of the year

-Output that has to do with when particular life stage transitions occur
X.all - number of germination events, flowering events, or dispersal events on every day of the simulation (vector the length of the # of days of your simulation)
X.days - average # of germination, flowering, or dispersal events on each day of the year across the entire simulation
X.weeks - average # of germination, flowering, or dispersal events on each week of the year across the entire simulation
X.months - average # of germination, flowering, or dispersal events on each month of the year across the entire simulation

-Output that has to do with the # of plants in a given life stage a given time
nX.all - number of seeds, rosettes, or flowering plants on every day of the simulation (vector the length of the # of days of your simulation) 
nX.daily- number of seeds, rosettes, or flowering plants one each day of the year (vector of 365 days of the year) 
.nX.weekly - number of seeds, rosettes, or flowering plants one each day of the year (vector of 52 weeks of the year) 
.nX.monthly - number of seeds, rosettes, or flowering plants one each day of the year (vector of 12 months of the year) 

-Outputs that are designed for making histograms of stage lengths
X.hist - all of the summarized output from the hist function in R for the given life stage length
.hist.density - density output from the hist function in R for the given life stage length (paired with the mid values below to create a histogram)
X.hist.mid - mid point values to pair with the density values above from the hist function in R for the given life stage length 
