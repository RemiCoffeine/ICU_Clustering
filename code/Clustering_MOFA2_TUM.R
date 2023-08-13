

#Author: Johannes Krell, Klinikum Rechts der Isar

#implementation of MOFA2 as robust clustering tool which can handle missing at random  data without interpolation
#source and citation: https://biofam.github.io/MOFA2/

######################################################################################################################

library(openxlsx)
library(MOFA2)
library(ggplot2)
library(tidyverse)
library(factoextra)
library(psych)
library(reticulate)


#your working dir
setwd("")

#folder with input data
data_location <- "input/"

#name of data file
data_file <- "12_08_23_SEPIS_ICU_ANONYM.csv"


#  check FAQ!
#  https://biofam.github.io/MOFA2/faq.html


#read data
sample_data <- read.csv( paste0(data_location,paste0('/',data_file)))

#short explanation: view= type of data, feature = name of feature, value= numeric measurement, sample = ID of patient
#temporal covariates are possible but not included in this script
#feature-sample duplicates are not possible!
#data is real ICU data, but anomynised and scrambled + normalized

###################################################################################################

#setup model, check homepage for different input formats
#currently use a tall dataframe 

model_setup<- create_mofa(data = sample_data, extract_metadata = T)


#plots the data structure
plot_data_overview(model_setup,
                   show_covariate = F,
                   show_dimensions = T)



#set model options
data_opts <- get_default_data_options(model_setup)
data_opts$center_groups <- FALSE

model_opts <- get_default_model_options(model_setup)

#set number of factors
model_opts$num_factors <- 5

#set training options
train_opts <- get_default_training_options(model_setup)
#train_opts$drop_factor_threshold <- 0.01
train_opts$convergence_mode <- "slow"
train_opts$verbose <- T
train_opts$gpu_mode <- F
train_opts$stochastic <- F

train_opts$seed <- 2021


#bind it all together
model_setup <- prepare_mofa(
  object = model_setup,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
) 


#name trained model file
outfile <- "MOFA_SEPSIS_sampledata.hdf5"

#run training
MOFAobject <- run_mofa(model_setup, outfile = outfile, use_basilisk = T)


##############################################################
#If it CRASHES: otherwise ignore

#it connects via the "reticulate" package to your python library of coice

# "pip install mofapy2" before execution

#insert path to your python environment here
use_python("/Users/***/anaconda3/bin/python3")

MOFAobject <- run_mofa(model_setup, outfile = outfile, use_basilisk = F)
###############################################################

#load pretrained model if needed
MOFAobject <- load_model(outfile, load_interpol_Z = F)

#check if factor corr is low
plot_factor_cor(MOFAobject)

#look at variance in model
p <- plot_variance_explained(MOFAobject, 
                             #x = "group",
                             #y = "view",
                             # split_by = NA,
                             plot_total = T) 

#variance in detail
p[[1]]+ theme(axis.text.x = element_text(angle = 90))

#overview
p[[2]]+ theme(axis.text.x = element_text(angle = 90))


######################################
#clustering of factors

#produce elbow plot
#Fetch factors in matrix format (a list, one matrix per group)
factors <- get_factors(MOFAobject)

# Concatenate groups
factors <- do.call("rbind",factors)

#elbow
fviz_nbclust(factors, kmeans, method = "wss")

#run kmeans
clusters <- cluster_samples(MOFAobject, k=5, factors=1:5)

#exporting clusters as ID-Cluster table
cluster_export <- data.frame(clusters$cluster)
write.csv(cluster_export, "yourfilename.csv")


#adds clusters to metadata##########
samples_metadata(MOFAobject) <- left_join(samples_metadata(MOFAobject),
                                          data.frame(clusters$cluster, sample = names(clusters$cluster)),
                                          by = c("sample")) 

#plot factors 
plot_factors(MOFAobject, 
             factors = c(1,2), 
             dot_size = 2.5,
             color_by = "clusters.cluster"
)


#Analysis:
#join clusters to original rawdata to figure out the deeper meaning and explore the factors with MOFA2



#optional:
#join metadata by patient ID, must have identical row number to "samples_metadata(MOFAobject)"
nrow(samples_metadata(MOFAobject))

samples_metadata(MOFAobject) <- left_join(samples_metadata(MOFAobject),
                                          your_metadata_file,
                                          by = c("sample")) 





