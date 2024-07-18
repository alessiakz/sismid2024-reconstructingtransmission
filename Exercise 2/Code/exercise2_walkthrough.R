# Exercise 2: outbreaker 2 for outbreak reconstruction

#setting file path
here::i_am("Exercise 2/Code/exercise2_walkthrough.R")

#libraries
library(ape) #reading, manipulating sequences
library(outbreaker2) #main reconstruction package
library(Hmisc) #general data science package - use for dealing with dates
library(lubridate) #also for dealing with dates
library(EpiEstim) #package for estimation of reproduction numbers

#set seed to allow for reproducibility
set.seed(50)

#outbreaker2 has fake_outbreak data, small simulated outbreak
data(fake_outbreak)

str(fake_outbreak)
#list with 6 variables

counts <- table(fake_outbreak$onset)
#makes table of how many instances of each onset time there are

#make a bar chart to display
barplot(counts, xlab="Onset time", ylab = "Number of cases", main = "Fake outbreak symptom onset times", col="tomato1")

# Running outbreaker2 on the fake outbreak
  #main function in outbreaker2 package is outbreaker -> runs transmission reconstruction

args(outbreaker) #prints possible inputs for the outbreaker function

#this tutorial will work with the defaults for the outbreaker model
  #therefore, the arguments we only need to provide the data, everything else can be left as is

#run outbreaker2 on the fake outbreak data with the default settings

fake_data <- outbreaker_data(dna = fake_outbreak$dna, dates = fake_outbreak$sample,
                             ctd = fake_outbreak$ctd, w_dens = fake_outbreak$w)
#this formats our data for what we need for outbreaker function
  #dates = collection dates
  #dna = DNA sequences in DNAbin format
  #ctd = contact tracing data, provided as a matrix/dataframe of two columns, indicating a reported contact between the two individuals
  #w_dens = indicates generation time distribution
  #f_dens = distribution of the colonization time

result <- outbreaker(data = fake_data)

class(result) #results are stored in a tyep of data frame called outbreaker_chains

#each row of result contains a single step from the MCMC chain
  #a single sample of each parameter
  # at each step, result returns the value of:
    # the log posterior distribution (post), 
    # the log likelihood (like), 
    # the log prior (prior), 
    # and all parameters/augmented data
    # inferred sampled ancestor (alpha)
    # date of infection (t_inf)
    # number of generations between the case and their sampled ancestor (kappa)

str(result) #shows you all elements of result

# analysis and visualization of results from the fake outbreak

#start by looking at the mixing of the MCMC algorithm, ensure our results are trusthworthy

plot(result) #shows a trace plot of the log-posterior values
#we used 10,000 iterations here, but in a full analysis, we would use more to minimize the effect of seed choice on the results

plot(result, "prior") #get trace plots for any other column
plot(result, "t_inf_29")
  #this one doesn't look well mixed
  #discrete jumps between different states
    #suggests to run more iterations, before trusting results
plot(result, "post")
  #transient stage is not ideal -> represents a chain which has not yet converged
  #discard this stage by defining a 'burn-in'
    #we want to keep as many samples as we can whilst removing that un-coverged stage
    #we want the trace plot to look like a "fuzzy caterpillar"

plot(result, burnin = 1000)
#now it looks like our MCMC mixed fairly well, so we can start getting preliminary results
  # real analysis would want you to run more iterations, but this is time consuming

#plot distribution of the parameters we estimated
  # mutation rate: mu
  # reporting probability: pi
  # contact reporting coverages: eps
  # non-infectious contact rate: lambda

#visualizing mutation rate
plot(result, "mu", type = "hist", burnin = 1000)
plot(result, "mu", type = "density", burnin = 1000)

#visualize inferred ancestries and transmission network
#who-infected-who: the size of the dot pertains to the inferred probability that that person (row) was the infector of each case (column)
plot(result, type = "alpha", burnin = 1000)

#distribution of the inferred infection time of each case
plot(result, type = "t_inf", burnin = 1000)

#how many generations were there between cases?
  # size of dot = probability
  # e.g. 1 generation = no unsampled cases, sampled ancestor is the sampled descendant

plot(result, type = "kappa", burnin = 1000)

#inferred transmission network!
  #thickness of line corresponds to probability of that transmission link
    #proportion of MCMC steps in which that pair were seen together
plot(result, type="network", burnin=1000, min_support=0.01)
  #this opens an interactive plot
  #min_support -> removes all edges with a probability less than its value
    # in this case 0.01, so you don't have a mess of lines, and can see the strongest inferred transmission pairs

#print summary of the outbreak inference
summary(result, burnin = 1000)
  #summary provides:
    #$step -> recap of which steps in the MCMC were used to create the summary
    #$post, $like, $prior, $mu, $pi -> distributional stats for these quantities
    #$tree -> consensus tree, tree made up of the most frequent infector/infectee pairs among MCMC samples
    #$support -> tells you the percentage of sampled trees which that pair occurred in
    #$generations -> most frequent value of kappa among those trees