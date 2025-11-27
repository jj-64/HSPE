# HSPE
Household Survey Parameter Estimation

# 0 - HOw to install
devtools::install_github("jj-64/HSPE")), the user can load:
library(HSPE)

# 1 - Load data
data("SumData")
head(SumData)

#2 - Create summarized data
In the "Execution"" folder, the code "Create Summarized data.R" will read each LIS.dat file from the "DATA" folder.
If "summarize_data()" is run, all files will be summarized. 
The code will save the generated dataframe in a database in "DataProcessed" folder to be used later on easily.
we can call this data using load("DataProcessed/SumData.rda") command.

If only certain files are needed, we can specify an index (as the file number) or a country code.
For more details, go to the function documentation.
