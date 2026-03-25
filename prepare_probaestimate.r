rm(list=ls(all=TRUE)) 

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(brms)
library(kableExtra)


######################################
## prepare the data for fits etc.
#####################################

exp <- 15; # specify which experiment to run
data <- read.csv(paste0("./data/DATA_EXPChoiceReliab_Exp", as.character(exp),".csv"))

## Generate the variables proba_1,.., proba_6 and color_1,...,color_6
data <- data %>%
  separate(Sample_Reliability, into = paste0("proba_", 1:6), sep = "\\s+") %>%
  mutate(across(starts_with("proba_"), as.numeric)) %>%     
  mutate(color = str_extract_all(Sample_Color, "blue|red"))  %>%
  unnest_wider(color, names_sep = "_") 

write.csv(data, paste0("./data/data_reliability_exp", as.character(exp),".csv"))
save(data, file = paste0("./data/data_reliability_exp", as.character(exp),".rdata"))



##################################################
##        PREPARE THE DATA       
##################################################

load(file = paste0("./data/data_reliability_exp", as.character(exp),".rdata"))

##  if ResponseButtonOrder= 1:  A scale is from red to blue
##  if ResponseButtonOrder= 0:  A scale is from blue to red

data <- data %>%
    mutate(PosteriorBlue = PosteriorBlue/100) %>%  
    mutate(est = SliderResponse/100) %>%
    mutate(estBlue = case_when(
      (ResponseButtonOrder == 1) ~ est,
      (ResponseButtonOrder == 0) ~ 1-est
    )) %>%
    mutate(estBlue_safe = case_when(
      (estBlue == 0) ~ 0.001,
      (estBlue == 1) ~ 0.999,
      TRUE ~ estBlue
    )) %>%   
    mutate(lo_proba = log(estBlue_safe / (1 - estBlue_safe))) %>%
    mutate(sim_choice = case_when(
      ## we simulate a binary choice by the subjective estimate 
      (estBlue_safe < 0.5) ~ 2, # 2= red, if the subj estimate was less than 50%, that means that ppt believed that red was more likely
      (estBlue_safe > 0.5) ~ 1, # 1= blue
      (estBlue_safe == 0.5) ~ sample(1:2, size = n(), replace=TRUE)
    )) %>%
    mutate_at(vars(starts_with("color")), ~ ifelse(. == "blue", 1, 2)) %>%
    rowwise() %>%
    mutate(sample_number = sum(!is.na(c_across(starts_with("proba_"))))) %>%  
    ungroup() %>%
    group_by(ParticipantPrivateID) %>% 
    mutate(normestBlue = (  (estBlue - min(estBlue)) / (max(estBlue) - min(estBlue))  ) 
           ) %>%
    mutate(normestBlue_safe = case_when(
      (normestBlue == 0) ~ 0.001,
      (normestBlue == 1) ~ 0.999,
      TRUE ~ normestBlue
    )) %>%   
    mutate(lo_normproba = log(normestBlue_safe / (1 - normestBlue_safe)))

N = length(unique(data$ParticipantPrivateID))
I_max <- max(data$sample_number) ## max number of samples/trial
## compute trials by subject
d <- data %>%
    group_by(ParticipantPrivateID) %>%
    summarise(t_subjs = n())
t_subjs <- d$t_subjs
subjs <- unique(data$ParticipantPrivateID)
T_max = t_subjs[1]
l_obs = log( max(data$proba_1/100) / (1 - max(data$proba_1/100))); # the maxima of the observed reliability

    
## Initialize data arrays
lo_proba  <- array(-1, c(N, T_max))
lo_normproba <- array(-1, c(N, T_max))
simchoice  <- array(-1, c(N, T_max))
posteriorBlue  <- array(-1, c(N, T_max))
color <- array( -1, c(N, T_max, I_max))
proba <- array(-1, c(N, T_max, I_max))
sample <- array(-1, c(N, T_max))
## fill the  arrays
for (n in 1:N) { ## loop through subjects
  t <- t_subjs[n] ## number of trials for subj i
  data_subj <- data %>% filter(ParticipantPrivateID == subjs[n])
  lo_proba[n, 1:t] <- data_subj$lo_proba
  lo_normproba[n, 1:t] <- data_subj$lo_normproba
  simchoice[n, 1:t] <- data_subj$sim_choice
  posteriorBlue[n, 1:t] <- data_subj$PosteriorBlue
  
  for (k in 1:t) { ## loop through trials
      data_subj_t <- data_subj[k,]
      sample[n,k] <- data_subj_t$sample_number
      for (i in 1:data_subj_t$sample_number) {
          color_var <- paste0("color_", i)
          proba_var <- paste0("proba_", i)
          color[n, k, i] <- data_subj[[color_var]][k]
          proba[n, k, i] <- data_subj[[proba_var]][k]/100
      }
  }
}


data_list_proba <- list(
    N = N,
    T_max = T_max,
    I_max = I_max,
    Tsubj = t_subjs,
    color = color,
    proba = proba,
    lo_proba = lo_proba, ##
    posteriorBlue = posteriorBlue,
    sample = sample,
    l_obs = l_obs
)

data_list_normproba <- list(
  N = N,
  T_max = T_max,
  I_max = I_max,
  Tsubj = t_subjs,
  color = color,
  proba = proba,
  lo_proba = lo_normproba, ##
  posteriorBlue = posteriorBlue,
  sample = sample,
  l_obs = l_obs
)

data_list_simchoice <- list(
  N = N,
  T_max = T_max,
  I_max = I_max,
  Tsubj = t_subjs,
  color = color,
  proba = proba,
  choice = simchoice, 
  posteriorBlue = posteriorBlue,
  sample = sample,
  l_obs = l_obs
)

save(data_list_proba, file = paste0('./data/data_list_proba_exp', as.character(exp), '.rdata'))
save(data_list_normproba, file = paste0('./data/data_list_normproba_exp', as.character(exp), '.rdata'))
save(data_list_simchoice, file = paste0('./data/data_list_simchoice_exp', as.character(exp), '.rdata'))

