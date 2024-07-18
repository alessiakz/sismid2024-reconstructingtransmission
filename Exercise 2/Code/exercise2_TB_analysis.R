# Exercise 2: outbreaker 2 for outbreak reconstruction
  #self-analysis of TB data

#setting file path
here::i_am("Exercise 2/Code/exercise2_walkthrough.R")

#libraries
library(ape) #reading, manipulating sequences
library(outbreaker2) #main reconstruction package
library(Hmisc) #general data science package - use for dealing with dates
library(lubridate) #also for dealing with dates
library(EpiEstim) #package for estimation of reproduction numbers

library(tidyverse)

#about the data
  #86 genomes from an outbreak of TB in Hamburg, Germany during late 1990s/early 2000s
  #hard to recreate transmission using epidemiological information alone since they often occur over long time frames
  #86 genomes belong to largest strain cluster

#.fasta -> aligned sequences
#.txt -> dates

#goal: use outbreaker2 to infer transmission networks for the TB outbreaks
  #repeat the analysis we did for the fake outbreak

#loading the data
myseq <- read.FASTA(file = here::here("Exercise 2/Data/roetzer2013.fasta"))
dates <- read.table(file = here::here("Exercise 2/Data/roetzer_dates.txt"))

#set dates to be in data format
dates_formatted <- dates %>% 
  mutate(date = as.Date(V1)) %>% 
  select(-V1)

#set everything to be in months
# For a vector of dates called 'dates'
dates_clean$months <- (year(dates_formatted) - 1997)*12 + month(dates_formatted) + day(dates_formatted)/monthDays(dates_formatted)

#dna is in the .fasta file, no contact tracing data, need dates, w_dens, f_dens
data <- fake_data <- outbreaker_data(dna = myseq, dates = dates,
                                      w_dens = fake_outbreak$w)


