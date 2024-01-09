library(dplyr)
library(ggplot2)
library(randomForest)
library(readr)
library(forcats)
set.seed(11002093)

# import data
sim_data <- read_table("I:/My Drive/Kelp Internship/Project Code/Kelp-Project/sim_data.csv", 
                       col_names = FALSE)
# create vector for column names
columns <-  c('S0', 'U0', 'A0', 'r', 'K', 'k_D', 'a_A', 'k_S', 'd_A', 'e_D', 'e_U', 'a_D', 
              'a_U', 'g_U', 'k_A', 'd_U', 'B', 'e_S', 'd_S', 'd_D', 'Send', 'Uend', 'Aend')

# assign column names
colnames(sim_data) <- columns

# create vector for parameter descriptions
descriptions <- c(
  '(S_0) Initial Seastar Population Density', #S0
  '(U_0) Initial Urchin Population Density', #U0
  '(A_0) Initial Kelp Population Density', #A0
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

############################################################################################
##################### SEASTARS #############################################################
# create random forest
SRF <- randomForest(Send~.,data = sim_data %>% select(-Aend,-Uend))

# create dataframe with parameters and importance
SRF_GSA <-
  data.frame(Parameter = rownames(SRF$importance),
             Importance = SRF$importance[, 1])

# set parameters as factor, order by importance
SRF_GSA$Parameter <-
  factor(SRF_GSA$Parameter, levels = SRF_GSA$Parameter[order(SRF_GSA$Importance,decreasing = TRUE)])

colors_SRF <- rev(c(
  'darkseagreen', #kD
  'seagreen', #B
  'seagreen', #K
  'slateblue3', #dU
  'slateblue3', #kA
  'coral2', #gU
  'slateblue3', #eU
  'darkseagreen', #eD
  'slateblue3', #aA
  'slateblue3', #U0
  'coral2', #S0
  'seagreen', #A0
  'darkseagreen', #dD
  'coral2', #eS
  'seagreen', #r
  'coral2', #kS
  'coral2', #dS
  'coral2', #aU
  'darkseagreen', #aD
  'darkseagreen' #dA
))

# plot results
SRF_GSA %>% ggplot(aes(Parameter, Importance)) + geom_col() + theme_bw() +
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

SRF_GSA[,1] <- descriptions

SRF_GSA %>% ggplot(aes(
  fct_reorder(Parameter, Importance), 
  Importance,
  fill = fct_reorder(Parameter, Importance))) + 
  geom_bar(stat = 'identity') + 
  coord_flip() +
  scale_fill_manual(values = colors_SRF) +
  labs(
#    title = "Importance of Ecological Factors to Seastar Population Density",
    x = "Ecological Factor",
    y = "Importance"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 24),
    legend.position = 'none'
  )

############################################################################################
##################### KELP #################################################################
# create random forest
ARF <- randomForest(Aend~.,data = sim_data %>% select(-Send,-Uend))

# create dataframe with parameters and importance
ARF_GSA <-
  data.frame(Parameter = rownames(ARF$importance),
             Importance = ARF$importance[, 1])

# set parameters as factor, order by importance
ARF_GSA$Parameter <-
  factor(ARF_GSA$Parameter, levels = ARF_GSA$Parameter[order(ARF_GSA$Importance,decreasing = TRUE)])

# plot results
ARF_GSA %>% ggplot(aes(Parameter, Importance)) + geom_col() + theme_bw() +
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

ARF_GSA[,1] <- descriptions

colors_ARF <- rev(c(
  'seagreen', #K
  'coral2', #dS
  'coral2', #gU
  'seagreen', #A0
  'darkseagreen', #dA
  'darkseagreen', #aD
  'darkseagreen', #dD
  'slateblue3', #eU
  'coral2', #eS
  'slateblue3', #dU
  'slateblue3', #U0
  'coral2', #S0
  'slateblue3', #aA
  'darkseagreen', #eD
  'coral2', #aU
  'seagreen', #r
  'darkseagreen', #kD
  'seagreen', #B
  'coral2', #kS
  'slateblue3' #kA
))

ARF_GSA %>% ggplot(aes(
  fct_reorder(Parameter, Importance), 
  Importance,
  fill = fct_reorder(Parameter, Importance))) + 
  geom_bar(stat = 'identity') + 
  coord_flip() +
  scale_fill_manual(values = colors_ARF) +
  labs(
#    title = "Importance of Ecological Factors to Kelp Population Density",
    x = "Ecological Factor",
    y = "Importance"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 24),
    legend.position = 'none'
  )

############################################################################################
##################### URCHINS ##############################################################
# create random forest
URF <- randomForest(Uend~.,data = sim_data %>% select(-Send,-Aend))

# create dataframe with parameters and importance
URF_GSA <-
  data.frame(Parameter = rownames(URF$importance),
             Importance = URF$importance[, 1])

# set parameters as factor, order by importance
URF_GSA$Parameter <-
  factor(URF_GSA$Parameter, levels = URF_GSA$Parameter[order(URF_GSA$Importance,decreasing = TRUE)])

# plot results
URF_GSA %>% ggplot(aes(Parameter, Importance)) + geom_col() + theme_bw() +
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

URF_GSA[,1] <- descriptions

colors_URF <- rev(c(
  'slateblue3', #eU
  'coral2', #aU
  'darkseagreen', #dA
  'seagreen', #A0
  'coral2', #kS
  'slateblue3', #kA
  'darkseagreen', #kD
  'seagreen', #r
  'slateblue3', #U0
  'seagreen', #K
  'darkseagreen', #dD
  'coral2', #dS
  'darkseagreen', #eD
  'slateblue3', #dU
  'darkseagreen', #aD
  'coral2', #gU
  'seagreen', #B
  'coral2', #S0
  'coral2', #eS
  'slateblue3' #aA
))

URF_GSA %>% ggplot(aes(
  fct_reorder(Parameter, Importance), 
  Importance,
  fill = fct_reorder(Parameter, Importance))) + 
  geom_bar(stat = 'identity') + 
  coord_flip() +
  scale_fill_manual(values = colors_URF) +
  labs(
#    title = "Importance of Ecological Factors to Urchin Population Density",
    x = "Ecological Factor",
    y = "Importance"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 24),
    legend.position = 'none'
  )
