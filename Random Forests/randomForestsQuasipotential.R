library(dplyr)
library(ggplot2)
library(randomForest)
library(readr)
set.seed(11002093)

quasi_data <- read.delim("quasi_data.csv",header = FALSE)
colnames(quasi_data) <- c('r','K','kD','aA','kS','dA','eD','eU','aD','aU','gU','kA','dU','B','eS','dS','dD','qp')
quasi_data$qp <- exp(quasi_data$qp)

PRF <- randomForest(qp~.,data = quasi_data)
PRF_GSA <-
  data.frame(Parameter = rownames(PRF$importance),
             Importance = PRF$importance[, 1])
PRF_GSA$Parameter <-
  factor(PRF_GSA$Parameter, levels = PRF_GSA$Parameter[order(PRF_GSA$Importance,decreasing = FALSE)])
PRF_GSA %>% ggplot(aes(Importance,Parameter)) + geom_col() + theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size=16),
    axis.text.x = element_text(angle = 90)
  ) + scale_y_discrete(labels = rev(c("Predator Predation Rate","Urchin Grazing Rate", "Urchin Conversion Efficiency",
                                  "Kelp Carrying Capacity","Drift Kelp Escape Rate","Fear Response of Urchins to Predators",
                                  "Kelp Growth Rate","Drift Kelp Retention","Predation Saturation Constant","Urchins Death Rate",
                                  "Effect of Urchins Starvation on Predation","Drift Kelp Consumption Rate",
                                  "Nutritional Value of Urchins","Predator Death Rate","Urchins Preference of Drift Kelp",
                                  "Drift Kelp Formation Rate","Predator Conversion Efficiency")))
