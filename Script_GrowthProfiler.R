#Growth profiler script
#Developed by Simone Zaghen, 18/june/2020
#Updated by Simone Zaghen 07/april/2022

#load packages
library(reshape2)
library(plyr, include.only = "revalue")
library(tidyverse)
library(growthrates)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set wd
rawdata <- read.csv2("Data.csv", header = F) #import data from growth profiler
number_of_replicates <- 3 #set number of replicates used here

#generate well number, sample number, and replicate number
col_name <- c("time", #first column header needs to be the time
              paste0( #from here generate WellNumber_SampleNumber_Replicate_Number
                rep(LETTERS[1:8], each=12), rep(1:12, 8), #generate well number
                "_",
                "s", #s for sample
                rep(1:(96/number_of_replicates), each = number_of_replicates), #generate sample number
                "_",
                rep(1:number_of_replicates))) #generate replicate number
colnames(rawdata) = col_name[1:ncol(rawdata)] #give names to columns

df <- melt(rawdata, 
           id=c("time"), 
           variable.name="sample", 
           value.name="OD", 
           na.rm = TRUE) %>%
  mutate(sample = as.character(sample)) %>%
  separate(sample, c("well", "sample", "replicate"), "_")

ggplot(data = df) +
  geom_line(aes(x = time, 
                y = OD, 
                group = replicate, 
                color = replicate)) +
  facet_wrap(~ sample,
             scales = "free") +
  theme_bw()

#Use this section of code to remove outliers
df <- 
  df %>%
  filter(well != "A1")#if you need to remove more than one sample, remove ) and place , 
         #well != "A2") #and uncomment this line
#you can re-run lines 34-41 to check if you removed all outliers

#create new data frame to store mean and sd of the replicates
average_replicates <- df %>% 
  select(time, sample, OD) %>%
  group_by(time, sample) %>%
  summarise(OD_mean = mean(OD, na.rm = T),
            OD_sd = sd(OD, na.rm= T))
average_replicates[1] <- average_replicates[1]/60 #convert minutes in hours

#rename the samples
average_replicates$sample <- revalue(average_replicates$sample, c("s1"= "Strain 1_Condition 1",
                                                                  "s2"= "Strain 2_Condition 1",
                                                                  "s3"= "Strain 3_Condition 2",
                                                                  "s4"= "Strain 3_Condition 1"))

#"_" is used to split the names into more columns
average_replicates <- average_replicates %>% separate(sample, c("strain", "condition"), "_")

#create plot
plot <- ggplot(average_replicates,
               mapping = aes(x = time,
                             y = OD_mean,
                             group = strain,
                             color = strain)) +
  geom_line() +
  geom_ribbon(aes(x=time, ymin = OD_mean - OD_sd, ymax = OD_mean + OD_sd, fill = strain), 
              alpha = 0.5, colour = NA) +
  facet_wrap(~ condition,
             scales = "free")
plot

#calculate lag phase (as time when OD is over certain threshold [approximation!])
lag_phase_all_samples <- average_replicates %>%
  group_by(strain, condition) %>% 
  filter(OD_mean > 0.25) %>% #set the threshold here 
  group_by(strain, condition) %>% 
  summarise(end_lag_phase = time[1])

#calculate maximum OD
OD_all_samples <- average_replicates %>%
  group_by(strain, condition) %>% 
  filter(OD_mean > 0.25) %>% #to filter out strains that did not grow
  filter(OD_mean == max(OD_mean)) %>%
  rename(max_OD = OD_mean) %>%
  select(strain, condition, max_OD)

#calculate the growth rate with the "growthrate" package - cite them if you use it! 
all_easylinear_fits <- all_easylinear(OD_mean ~ time | strain + condition, 
                                      data = average_replicates, spar = 0.5)
par(mfrow=c(2,2)) #set the parameters in order to display all the samples!
plot(all_easylinear_fits) #check if the growth rate was calculated in the correct part of the curve!

growth_rate <- as.data.frame(coef(all_easylinear_fits)) %>% 
  select(mumax) %>% 
  rename(growth_rate = mumax) %>% 
  rownames_to_column("sample") %>%
  separate(sample, c("strain", "condition"), ":")  

#merge and export the data
parameters <- merge(growth_rate, lag_phase_all_samples) %>% 
  merge(OD_all_samples)
write.csv2(parameters, "parameters.csv", row.names = FALSE)
