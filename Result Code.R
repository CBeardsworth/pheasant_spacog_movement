# R Code for Analysis on Paper : Spatial cognitive ability is associated with speed of movement in the early exploration of an environment

# Packages required:

library(lme4)
library(ggplot2)
library(sjPlot)

# Data required:

Efficiency <- read.csv("TransitMSData_Efficiency.csv", stringsAsFactors = F)

Efficiency <- Efficiency[Efficiency$Duration < 21600,] # Remove transit paths that were over 6 hours long. 

##### Does Spatial ability, sex or motivation affect speed or straightness of transit? ####

#### Speed ####

#Full Model
m1 <- glmer(Speed ~ ScaledDay * Spatial_Score + Sex + ScaledTestOrder + (1|ATLAS), data=Efficiency, family = Gamma(link="log"))
summary(m1)
drop1(m1, test="Chi")

#drop Sex
m1 <- glmer(Speed ~ ScaledDay * Spatial_Score + ScaledTestOrder + (1|ATLAS), data=Efficiency, family = Gamma(link="log"))
summary(m1)
drop1(m1, test="Chi")

#drop TO
m1 <- glmer(Speed ~ ScaledDay * Spatial_Score + (1|ATLAS), data=Efficiency, family = Gamma(link="log"))
summary(m1)
drop1(m1, test="Chi")

plot_model(m1, "int")+
  scale_x_continuous(expand = c(0, 0), breaks = c(-1.0172767,0.2919676,	1.6012120), labels=c(10,20,30), limits = c(-2, 1.6012120)) + # Cannot change plot to have non scaled days within the model as confidence intervals change slightly, however can match days to scaled days manually on the x axis so that model stays the same.
  scale_y_continuous(expand = c(0, 0), limits=c(0,0.11)) +
  scale_colour_manual(values=c("red", "deepskyblue3"))+
  scale_fill_manual(values=c("red","deepskyblue3"))+
  ylab("Speed (m/s)") +
  xlab("Date (September)") +
  expand_limits(x = 0, y = -2) +
  theme_classic()+
  ggtitle("")+
  theme(legend.position="none")# This automatically chooses ind with 9 (highest score) as comparison but has no conf intervals near end. Instead, plot using groups for low, med, high performance.

#ggsave("SpeedGLMM_int.jpeg", units="cm",width=15, height = 10, dpi= 600)


#### straightness ####

m1 <- glmer(Straightness ~ ScaledDay * Spatial_Score + Sex + ScaledTestOrder + scale(N) + (1|ATLAS), data=Efficiency, family = Gamma(link="log"))
#dredge(m1)
summary(m1)
drop1(m1, test="Chi")

m1 <- glmer(Straightness ~ ScaledDay * Spatial_Score + Sex + scale(N) + (1|ATLAS), data=Efficiency, family = Gamma(link="log"))
summary(m1)
drop1(m1, test="Chi")

m1 <- glmer(Straightness ~ ScaledDay * Spatial_Score + scale(N) + (1|ATLAS), data=Efficiency, family = Gamma(link="log"))
summary(m1)
drop1(m1, test="Chi")

m1 <- glmer(Straightness ~ ScaledDay + Spatial_Score + scale(N)+ (1|ATLAS), data=Efficiency, family = Gamma(link="log"))
summary(m1)
drop1(m1, test="Chi")

m1 <- glmer(Straightness ~ ScaledDay + scale(N) + (1|ATLAS), data=Efficiency, family = Gamma(link="log"))
summary(m1)
drop1(m1, test="Chi")

plot_model(m1, "pred")$ScaledDay +
  theme_classic()+
  scale_x_continuous(expand = c(0, 0),breaks = c(-1.0172767,0.2919676,	1.6012120), labels=c(10,20,30), limits = c(-2, 1.6012120)) +
  xlab("Date (September)") +
  scale_y_continuous(expand = c(0, 0), limits=c(0,1)) +
  ggtitle("")

#ggsave("StraightGLMM_eff.jpeg", units="cm",width=15, height = 10, dpi= 600)

