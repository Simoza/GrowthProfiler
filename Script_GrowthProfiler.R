#Growth profiler script
#Developed by Simone Zaghen, 18/june/2020

#load packages
library(readxl)
library(tidyr)
library(reshape2)
library(Rmisc)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
library(deSolve)
library(dplyr)
library(tibble)
library(growthrates)


setwd("/Users/zaghen/OneDrive - Chalmers/z.LAB/Templates/GrowthProfiler/")
rawdata <- read_excel("Data.xlsx")

#create the df and plot all the replicates 
sampleindex = 1
for (i in 2:ncol(rawdata)) {
  replicateindex = (i-1) %% 3 #3 for triplicates, 4 for quadruplicates and so on
  
  if(replicateindex == 0){
    replicateindex = 3} #3 for triplicates, 4 for quadruplicates and so on
  if(replicateindex == 1 && i!=2){
    sampleindex = sampleindex + 1}
  samplename = paste0("s",sampleindex)
  colnames(rawdata)[i] = paste0(colnames(rawdata)[i],"_",
                                samplename,
                                "_",
                                replicateindex)}

df_to_visualise_triplicates <- melt(rawdata, 
                                    id=c("time"), 
                                    variable.name="sample", 
                                    value.name="OD", 
                                    na.rm = TRUE) %>%
  mutate(sample = as.character(sample)) %>%
  separate(sample, c("well", "sample", "replicate"), "_")

ggplot(data = df_to_visualise_triplicates) +
  geom_line(aes(x = time, 
                y = OD, 
                group = replicate, 
                color = replicate)) +
  facet_wrap(~ sample,
             scales = "free") +
  theme_bw()

#create new data frame to store mean and sd of the replicates
#you will receive a warning because summarise_each has been deprecated.
#I will fix this sooner or later - but it still works. 
average_replicates <- df_to_visualise_triplicates %>% 
  select(time, sample, OD) %>%
  group_by(time, sample) %>%
  summarise_each(funs(mean, sd)) %>%
  rename(OD = mean)
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
                             y = OD,
                             group = strain,
                             color = strain)) +
  geom_line() +
  geom_ribbon(aes(x=time, ymin = OD - sd, ymax = OD + sd, fill = strain), alpha = 0.5, colour = NA) +
  facet_wrap(~ condition,
             scales = "free")
plot

#calculate lag phase (as time when OD is over certain threshold [approximation!])
lag_phase_all_samples <- average_replicates %>%
  group_by(strain, condition) %>% 
  filter(OD > 0.25) %>% #set the threshold here 
  group_by(strain, condition) %>% 
  summarise(end_lag_phase = time[1])

#calculate maximum OD
OD_all_samples <- average_replicates %>%
  group_by(strain, condition) %>% 
  filter(OD > 0.25) %>% #to filter out strains that did not grow
  filter(OD == max(OD)) %>%
  rename(max_OD = OD) %>%
  select(strain, condition, max_OD)

#calculate the growth rate with the "growthrate" package - cite them if you use it! 
all_easylinear_fits <- all_easylinear(OD ~ time | strain + condition, data = average_replicates, spar = 0.5)
par(mfrow=c(2,2)) #set the parameters in order to display all the samples!
plot(all_easylinear_fits) #check if the growth rate was calculated in the correct part of the curve!

growth_rate <- as.data.frame(coef(all_easylinear_fits)) %>% 
  select(mumax) %>% 
  rename(growth_rate = mumax) %>% 
  rownames_to_column("sample") %>%
  separate(sample, c("strain", "condition"), ":")  

#merge and export the dataframes
parameters <- merge(growth_rate, lag_phase_all_samples) %>% 
  merge(OD_all_samples)
write.csv2(parameters, "parameters.csv", row.names = FALSE)
