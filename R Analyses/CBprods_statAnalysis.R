# Load libraries for analysis
library(car)
library(stats)

# Set working directory to be directory of current R script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load CSV file for feature that I would like to analyze
featName = "HR"
dataList = read.csv(paste("CBprods_csvReady_", featName, ".csv", sep = ""), header = FALSE)

# Convert the list that read.csv outputs into a matrix
dataMatrix = matrix(unlist(dataList), nrow = 50, ncol = 3)



############ Friedman ################

# Apply the Friedman's rank sum test to the data
friedman.test(dataMatrix)



########## 1-Way Repeated Measures ANOVA #########
# Create levels and then convert to factor and then convert to dataframe
inputLevels = c(1, 2, 3)
inputFactor = as.factor(inputLevels)
inputFrame = data.frame(inputFactor)

# Fit a linear model to the data
dataModel = lm(dataMatrix ~ 1)

# Employ the ANOVA function
analysis = Anova(dataModel, idata = inputFrame, idesign = ~inputFactor)

# Output the results
summary(analysis)

