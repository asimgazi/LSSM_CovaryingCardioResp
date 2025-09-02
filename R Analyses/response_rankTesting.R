# Load libraries for analysis
library(stats)
library(rmcorr)

# To save time, I will hard code the numbers for binomial tests since I just
# proportions of "successes" which is a really easy number to compute from
# the sign data that I output in my MATLAB script
binomResult_real_traum = binom.test(7, 50)
print(binomResult_real_traum)

binomResult_real_stim = binom.test(39, 50) 
print(binomResult_real_stim)

binomResult_sim_traum = binom.test(12, 50)
print(binomResult_sim_traum)

binomResult_sim_stim = binom.test(34, 50)
print(binomResult_sim_stim)



#### Repeated measures correlation of ranks

# Set working directory to be directory of current R script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load CSV file
rankData_real = as.matrix(read.csv("rankData_real.csv", header = FALSE))
rankData_sim = as.matrix(read.csv("rankData_sim.csv", header = FALSE))

# Store number of subjects and repeated measures
num_subjects = nrow(rankData_real)
num_measures = ncol(rankData_real)

# Create a data frame with subject ID, repeated measure index, and the two variables
df = data.frame(
  Subject = rep(1:num_subjects, each = num_measures),
  Measure = rep(1:num_measures, times = num_subjects),
  Real = as.vector(t(rankData_real)),  # Flatten by row
  Sim = as.vector(t(rankData_sim))
)

# Run repeated measures correlation
rmcorr_result = rmcorr(participant = Subject, measure1 = Real, measure2 = Sim, dataset = df)

# Print results
print(rmcorr_result)

# Plot the repeated measures correlation
plot(rmcorr_result)
