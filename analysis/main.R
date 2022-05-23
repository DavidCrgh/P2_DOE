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

#Data <- Data %>% 
 # select(everything()) %>%
  #mutate(across(c(Environment, Scene, Parameters), as.factor))

Data$Environment <- as.factor(Data$Environment)
Data$Scene <- as.factor(Data$Scene)
Data$Parameters <- as.factor(Data$Parameters)


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

Anova(model, type="II")


# Neither of the assumptions (i.e. normality, homoscedasticity) are met
# We'll start using transformations to correct this

# Square root

T_sqrt <- sqrt(Data$Time)
sqrt_model <- lm(T_sqrt ~ Environment*Scene*Parameters,
                data = Data)
test_assumptions(sqrt_model, "Square Root", hist_xlim = c(-0.4, 0.7))


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

Anova(log_model, type="II")

# Visual analysis show slight but acceptable deviations from normality
# Homoscedasticity looks OK with a few outliers, no observable pattern


# # Tukey's ladder of powers
# T_tuk <- transformTukey(Data$Time,
#                         plotit = FALSE)
# tuk_model <- lm(T_tuk ~ Environment*Scene*Parameters,
#                 data = Data)
# test_assumptions(tuk_model, "Tukey", hist_xlim = c(-0.0018, 0.003))
# 
# 
# # Box-Cox
# Box <- boxcox(Data$Time ~ 1,              
#              lambda = seq(-6,6,0.1)     
# )
# 
# Cox <- data.frame(Box$x, Box$y)
# 
# Cox2 <- Cox[with(Cox, order(-Cox$Box.y)),]
# Cox2[1,]
# 
# lambda <- Cox2[1, "Box.x"]
# 
# T_box <- (Data$Time ^ lambda - 1) / lambda
# box_model <- lm(T_box ~ Environment*Scene*Parameters,
#                 data = Data)
# test_assumptions(box_model, "Box-Cox", hist_xlim = c(-0.01, 0.02))


# Perform Levene Tests to formally test homogeneity of variance

leveneTest(Time ~ Environment*Scene*Parameters, data = Data)
leveneTest(T_sqrt ~ Environment*Scene*Parameters, data = Data)
leveneTest(T_cubic ~ Environment*Scene*Parameters, data = Data)
#leveneTest(T_log ~ Environment*Scene*Parameters, data = Data) 
#leveneTest(T_tuk ~ Environment*Scene*Parameters, data = Data)
#leveneTest(T_box ~ Environment*Scene*Parameters, data = Data)


# Test assumptions without outliers
#Data_no_outliers <- Data %>%
 # filter(!row_number() %in% c(172,170,258))

#model2 <- lm(Time ~ Environment*Scene*Parameters,
 #           data = Data_no_outliers)

#test_assumptions(model2, "Original No outliers", hist_xlim = c(-15,15))

#T_log2 <- log(Data_no_outliers$Time)
#log_model2 <- lm(T_log2 ~ Environment*Scene*Parameters,
 #               data = Data_no_outliers)

#test_assumptions(log_model2, "Logarithmic No Outliers", hist_xlim = c(-0.03, 0.06))

# leveneTest(T_log2 ~ Environment*Scene*Parameters, data = Data_no_outliers) 


# Generate boxplots for Time in response to Environment and other factors
# Environment ==========================================================
Sum = Summarize(T_log ~ Environment,
                data = Data,
                digits = 3)

# Compute the standard error (se)
Sum$se = Sum$sd / sqrt(Sum$n)
Sum$se = signif(Sum$se, digits = 3)

pd = position_dodge(.2)
ggplot(Sum, aes(x = Environment,
                y = mean,
                color = Environment)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se),
                width = .2, size = 0.7, position = pd) +
  geom_point(aes(shape = Environment), size = 5, position = pd) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 20,
                                  hjust = 0.5),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(face = "bold", size = 15),
        plot.caption = element_text(hjust = 0),
        legend.text = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 20),
        legend.justification = c(1,0),
        legend.position = "none") +
  ylab("Mean of Time (log-transformed)") +
  ggtitle("Time vs Environment")

# De-transform the data
Sum$mean = exp(Sum$mean)
Sum$sd = exp(Sum$sd)

# Compute the standard error (se)
Sum$se = Sum$sd / sqrt(Sum$n)
Sum$se = signif(Sum$se, digits = 3)

ggplot(Sum, aes(x = Environment,
                y = mean,
                color = Environment)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se),
                width = .2, size = 0.7, position = pd) +
  geom_point(aes(shape = Environment), size = 5, position = pd) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 20,
                                  hjust = 0.5),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(face = "bold", size = 15),
        plot.caption = element_text(hjust = 0),
        legend.text = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 20),
        legend.justification = c(1,0),
        legend.position = "none") +
  ylab("Mean of time (s)") +
  ggtitle("Time vs Environment")

# Environment and scenes =========================================

Sum = Summarize(T_log ~ Environment + Scene,
                data = Data,
                digits = 3)

# Compute the standard error (se)
Sum$se = Sum$sd / sqrt(Sum$n)
Sum$se = signif(Sum$se, digits = 3)

ggplot(Sum, aes(x = Scene,
                y = mean,
                color = Environment)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se),
                width = .2, size = 0.7, position = pd) +
  geom_point(aes(shape = Environment), size = 5, position = pd) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 15,
                                  hjust = 0.5),
        axis.title = element_text(face = "bold", size = 15),
        axis.text = element_text(face = "bold", size = 10),
        plot.caption = element_text(hjust = 0),
        legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 12),
        legend.justification = c(0,1),
        legend.position = c(0.01, 0.99)) +
  xlab("Scene") +
  ylab("Mean of time (s) (log-transformed)") +
  ggtitle("Time vs Enviroment and Scene")
  
# De-transform the data

Sum = Summarize(T_log ~ Environment + Scene,
                data = Data,
                digits = 3)

Sum$mean = exp(Sum$mean)
Sum$sd = exp(Sum$sd)

# Compute the standard error (se)
Sum$se = Sum$sd / sqrt(Sum$n)
Sum$se = signif(Sum$se, digits = 3)

ggplot(Sum, aes(x = Scene,
                y = mean,
                color = Environment)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se),
                width = .2, size = 0.7, position = pd) +
  geom_point(aes(shape = Environment), size = 5, position = pd) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 15,
                                  hjust = 0.5),
        axis.title = element_text(face = "bold", size = 15),
        axis.text = element_text(face = "bold", size = 10),
        plot.caption = element_text(hjust = 0),
        legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 12),
        legend.justification = c(0,1),
        legend.position = c(0.01, 0.99)) +
  xlab("Scene") +
  ylab("Mean of time (s)") +
  ggtitle("Time vs Enviroment and Scene")


# Environment and parameters =========================================
Sum = Summarize(T_log ~ Environment + Parameters,
                data = Data,
                digits = 3)

# Compute the standard error (se)
Sum$se = Sum$sd / sqrt(Sum$n)
Sum$se = signif(Sum$se, digits = 3)

ggplot(Sum, aes(x = Parameters,
                y = mean,
                color = Environment)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se),
                width = .2, size = 0.7, position = pd) +
  geom_point(aes(shape = Environment), size = 5, position = pd) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 15,
                                  hjust = 0.5),
        axis.title = element_text(face = "bold", size = 15),
        axis.text = element_text(face = "bold", size = 10),
        plot.caption = element_text(hjust = 0),
        legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 12),
        legend.justification = c(0,1),
        legend.position = c(0.01, 0.99)) +
  xlab("Parameters") +
  ylab("Mean of time (s) (log-transformed)") +
  ggtitle("Time vs Enviroment and Parameters")

# De-transform the data

Sum = Summarize(T_log ~ Environment + Parameters,
                data = Data,
                digits = 3)

Sum$mean = exp(Sum$mean)
Sum$sd = exp(Sum$sd)

# Compute the standard error (se)
Sum$se = Sum$sd / sqrt(Sum$n)
Sum$se = signif(Sum$se, digits = 3)

ggplot(Sum, aes(x = Parameters,
                y = mean,
                color = Environment)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se),
                width = .2, size = 0.7, position = pd) +
  geom_point(aes(shape = Environment), size = 5, position = pd) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 15,
                                  hjust = 0.5),
        axis.title = element_text(face = "bold", size = 15),
        axis.text = element_text(face = "bold", size = 10),
        plot.caption = element_text(hjust = 0),
        legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 12),
        legend.justification = c(0,1),
        legend.position = c(0.01, 0.99)) +
  xlab("Parameters") +
  ylab("Mean of time (s)") +
  ggtitle("Time vs Enviroment and Parameters")


# Conduct pairwise t-tests to view p-values between groups
# Just environment
pairwise.t.test(T_log, Data$Environment,
                p.adjust.method = "BH")


# Environment and Scene
pairwise.t.test(T_log, Data$Environment:Data$Scene,
                p.adjust.method = "BH")

# Environment and Parameters
pairwise.t.test(T_log, Data$Environment:Data$Parameters,
                p.adjust.method = "BH")