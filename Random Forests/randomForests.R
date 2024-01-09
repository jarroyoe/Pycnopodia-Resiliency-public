library(dplyr)
library(ggplot2)
library(randomForest)
set.seed(11002093)

sim_data <- read.delim("~/Documentos/Github/kelp-climate-change/DS Analysis/sim_data.csv")

SRF <- randomForest(Send~.,data = sim_data %>% select(-Aend,-Uend))
SRF_GSA <-
  data.frame(Parameter = rownames(SRF$importance),
             Importance = SRF$importance[, 1])
SRF_GSA$Parameter <-
  factor(SRF_GSA$Parameter, levels = SRF_GSA$Parameter[order(SRF_GSA$Importance,decreasing = TRUE)])
SRF_GSA %>% ggplot(aes(Parameter, Importance)) + geom_col() + theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size=16),
    axis.text.x = element_text(angle = 90)
  )

ARF <- randomForest(Aend~.,data = sim_data %>% select(-Send,-Uend))
ARF_GSA <-
  data.frame(Parameter = rownames(ARF$importance),
             Importance = ARF$importance[, 1])
ARF_GSA$Parameter <-
  factor(ARF_GSA$Parameter, levels = ARF_GSA$Parameter[order(ARF_GSA$Importance,decreasing = TRUE)])
ARF_GSA %>% ggplot(aes(Parameter, Importance)) + geom_col() + theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size=16),
    axis.text.x = element_text(angle = 90)
  )

URF <- randomForest(Uend~.,data = sim_data %>% select(-Send,-Aend))
URF_GSU <-
  data.frame(Parameter = rownames(URF$importance),
             Importance = URF$importance[, 1])
URF_GSU$Parameter <-
  factor(URF_GSU$Parameter, levels = URF_GSU$Parameter[order(URF_GSU$Importance,decreasing = TRUE)])
URF_GSU %>% ggplot(aes(Parameter, Importance)) + geom_col() + theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size=16),
    axis.text.x = element_text(angle = 90)
  )