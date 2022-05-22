## TODO: add header

# Packages
library(tidyverse)
library(FSA)
library(psych)
library(car)
library(rcompanion)
library(MASS)

# Import utils script
source(".\\analysis\\utils.R")

# Constants
DATA_PATH = '.\\data\\pbrt-samples2k.csv'

# Load data
Data <- read.csv(DATA_PATH, header = TRUE)
Data$Environment <- factor(Data$Environment,
                           levels=unique(Data$Environment))

Data <- Data %>% 
  select(everything()) %>%
  mutate(across(c(Environment, Scene, Parameters), as.factor))


# Verify integrity
headTail(Data)
str(Data)
summary(Data)


# Check number of samples per specific scenario
scenario_counts <- Data %>% 
  group_by(Environment, Scene, Parameters) %>%
  tally()

# Generate interaction plot for each factor
# Environment and Scene
interaction.plot(
  x.factor = Data$Scene,
  trace.factor = Data$Environment,
  response = Data$Time,
  fun = mean,
  type = "b",
  col = c("black", "red", "green", "blue", "purple"),
  pch = c(19, 17, 15),
  fixed = TRUE,
  leg.bty = "o"
)

# Environment and Parameters
interaction.plot(
  x.factor = Data$Parameters,
  trace.factor = Data$Environment,
  response = Data$Time,
  fun = mean,
  type = "b",
  col = c("black", "red", "green", "blue", "purple"),
  pch = c(19, 17, 15),
  fixed = TRUE,
  leg.bty = "o"
)


# Test ANOVA assumptions
model <- lm(Time ~ Environment*Scene*Parameters,
           data = Data)
test_assumptions(model, "Original", hist_xlim = c(-15,15))


# Neither of the assumptions (i.e. normality, homoscedasticity) are met
# We'll start using transformations to correct this

# Square root

T_sqrt <- sqrt(Data$Time)
sqrt_model <- lm(T_sqrt ~ Environment*Scene*Parameters,
                data = Data)
test_assumptions(sqrt_model, "Square Root", hist_xlim = c(-0.4, 0.7))

# Assumptions still not met


# Cubic root

T_cubic <- sign(Data$Time) * abs(Data$Time)^(1/3)
cubic_model <- lm(T_cubic ~ Environment*Scene*Parameters,
                data = Data)
test_assumptions(cubic_model, "Cubic Root", hist_xlim = c(-0.07, 0.15))


# Logarithmic

T_log <- log(Data$Time)
log_model <- lm(T_log ~ Environment*Scene*Parameters,
                  data = Data)
test_assumptions(log_model, "Logarithmic", hist_xlim = c(-0.03, 0.06))

Anova(cubic_model, type="II")


# Tukey's ladder of powers
T_tuk <- transformTukey(Data$Time,
                        plotit = FALSE)
tuk_model <- lm(T_tuk ~ Environment*Scene*Parameters,
                data = Data)
test_assumptions(tuk_model, "Tukey", hist_xlim = c(-0.0018, 0.003))


# Box-Cox
Box <- boxcox(Data$Time ~ 1,              # Transform Turbidity as a single vector
             lambda = seq(-6,6,0.1)      # Try values -6 to 6 by 0.1
)

Cox <- data.frame(Box$x, Box$y)

Cox2 <- Cox[with(Cox, order(-Cox$Box.y)),]
Cox2[1,]

lambda <- Cox2[1, "Box.x"]

T_box <- (Data$Time ^ lambda - 1) / lambda
box_model <- lm(T_box ~ Environment*Scene*Parameters,
                data = Data)
test_assumptions(box_model, "Box-Cox", hist_xlim = c(-0.01, 0.02))

leveneTest(Time ~ Environment*Scene*Parameters, data = Data)
leveneTest(T_sqrt ~ Environment*Scene*Parameters, data = Data)
leveneTest(T_cubic ~ Environment*Scene*Parameters, data = Data)
leveneTest(T_log ~ Environment*Scene*Parameters, data = Data)
leveneTest(T_tuk ~ Environment*Scene*Parameters, data = Data)
leveneTest(T_box ~ Environment*Scene*Parameters, data = Data)
