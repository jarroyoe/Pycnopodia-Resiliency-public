library(dplyr)
library(ggplot2)
library(randomForest)
library(readr)
library(forcats)
set.seed(11002093)

# read data
sim_data_stab <- read_csv("I:/My Drive/Kelp Internship/Project Code/Kelp-Project/sim_data_stab.csv")

# assign column names
colnames(sim_data_stab) <- c('r','K','kD','aA','kS','dA','eD','eU','aD','aU','gU','kA','dU','B','eS','dS','dD','propForest')

##############################################################################
################## with kelp growth parameters ###############################

# create random forest
PRF <- randomForest(propForest~.,data = sim_data_stab)

# create dataframe with parameters and importance
PRF_GSA <-
  data.frame(Parameter = rownames(PRF$importance),
             Importance = PRF$importance[, 1])

# set parameters as factor, order by importance
PRF_GSA$Parameter <-
  factor(PRF_GSA$Parameter, levels = PRF_GSA$Parameter[order(PRF_GSA$Importance,decreasing = TRUE)])

# plot results
PRF_GSA %>% ggplot(aes(Parameter, Importance)) + geom_col() + theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size=16),
    axis.text.x = element_text(angle = 90)
  )

##################
## fancy bar chart

# parameter descriptions
descriptions_PRF <- c(
  '(r) Kelp Growth Rate', #r
  '(K) Kelp Carrying Capacity', #K
  '(\u03BA_D) Urchin Preference for Drift Kelp', #kD
  '(\u03B1_A) Urchin Grazing Rate', #aA
  '(\u03BA_S) Urchin Fear Response to Predators', #kS
  '(\u03B4_A) Kelp to Drift Kelp Conversion Rate', #dA
  '(\u03B5_D) Drift Kelp Retention', #eD
  '(\u03B5_U) Kelp to Urchin Conversion Efficiency', #eU
  '(\u03B1_D) Drift Kelp Consumption Rate', #aD
  '(\u03B1_U) Seastar Predation Rate', #aU
  '(\u03B3_U) Seastar Predation Saturation Constant', #gU
  '(\u03BA_A) Urchin Starvation Effect on Seastar Predation Rate', #kA
  '(\u03B4_U) Urchin Death Rate', #dU
  '(\u03b2) Kelp Density Impact to Nutritional Value of Urchins', #B
  '(\u03B5_S) Urchin to Seastar Conversion Efficiency', #eS
  '(\u03B4_S) Seastar Death Rate', #dS
  '(\u03B4_D) Drift Kelp Escape Rate' #dD
)

PRF_GSA[,1] <- descriptions_PRF

colors_PRF <- rev(c(
  'seagreen', #K
  'darkseagreen', #dA
  'seagreen', #r
  'darkseagreen', #aD
  'coral2', #eS
  'slateblue3', #eU
  'darkseagreen', #eD
  'slateblue3', #aA
  'coral2', #aU
  'coral2', #gU
  'slateblue3', #dU
  'darkseagreen', #dD
  'coral2', #kS
  'seagreen', #B
  'coral2',#dS
  'darkseagreen', #kD
  'slateblue3' #kA
))

PRF_GSA %>% ggplot(aes(
  fct_reorder(Parameter, Importance), 
  Importance,
  fill = fct_reorder(Parameter, Importance))) + 
  geom_bar(stat = 'identity') + 
  coord_flip() +
  scale_fill_manual(values = colors_PRF) +
  labs(
#    title = "Importance of Ecological Factors to Kelp Forest Stable State",
    x = "Ecological Factor",
    y = "Importance"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 24),
    legend.position = "none"
  )

##############################################################################
################## without kelp growth parameters ############################

# kelp growth parameters
kelpGrowth <- c('r','K','dA')

# create dataframe that excludes kelp growth parameters
ex_kelpGrowth <- sim_data_stab[,!(names(sim_data_stab) %in% kelpGrowth)]

# create random forest
TRF <- randomForest(propForest~.,data = ex_kelpGrowth)

# global sensitivity analysis
TRF_GSA <-
  data.frame(Parameter = rownames(TRF$importance),
             Importance = TRF$importance[, 1])
TRF_GSA$Parameter <-
  factor(TRF_GSA$Parameter, levels = TRF_GSA$Parameter[order(TRF_GSA$Importance,decreasing = TRUE)])

# plot results
TRF_GSA %>% ggplot(aes(Parameter, Importance)) + geom_col() + theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size=16),
    axis.text.x = element_text(angle = 90)
  )

##################
## fancy bar chart

# parameter descriptions
descriptions_TRF <- c(
  '(\u03BA_D) Urchin Preference for Drift Kelp', #kD
  '(\u03B1_A) Urchin Grazing Rate', #aA
  '(\u03BA_S) Urchin Fear Response to Predators', #kS
  '(\u03B5_D) Drift Kelp Retention', #eD
  '(\u03B5_U) Kelp to Urchin Conversion Efficiency', #eU
  '(\u03B1_D) Drift Kelp Consumption Rate', #aD
  '(\u03B1_U) Seastar Predation Rate', #aU
  '(\u03B3_U) Seastar Predation Saturation Constant', #gU
  '(\u03BA_A) Urchin Starvation Effect on Seastar Predation Rate', #kA
  '(\u03B4_U) Urchin Death Rate', #dU
  '(\u03b2) Kelp Density Impact to Nutritional Value of Urchins', #B
  '(\u03B5_S) Urchin to Seastar Conversion Efficiency', #eS
  '(\u03B4_S) Seastar Death Rate', #dS
  '(\u03B4_D) Drift Kelp Escape Rate' #dD
)

TRF_GSA[,1] <- descriptions_TRF

colors_TRF <- rev(c(
  'coral2', #gU
  'darkseagreen', #eD
  'darkseagreen', #aD
  'coral2', #eS
  'slateblue3', #aA
  'slateblue3', #eU
  'coral2',#dS
  'slateblue3', #dU
  'coral2', #aU
  'darkseagreen', #dD
  'seagreen', #B
  'coral2', #kS
  'slateblue3', #kA
  'darkseagreen' #kD
))

TRF_GSA %>% ggplot(aes(
  fct_reorder(Parameter, Importance), 
  Importance,
  fill = fct_reorder(Parameter, Importance))) + 
  geom_bar(stat = 'identity') + 
  coord_flip() +
  scale_fill_manual(values = colors_TRF) +
  labs(
#    title = "Importance of Ecological Factors to Kelp Forest Stable State\n(Trophic Chain)",
    x = "Ecological Factor",
    y = "Importance"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 24),
    legend.position = "none"
  )


























