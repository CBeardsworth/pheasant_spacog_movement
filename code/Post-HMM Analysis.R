# R Code for Analysis on Paper : Spatial cognitive ability is associated with speed of movement in the early exploration of an environment

# Packages required:

library(lme4)
library(ggplot2)
library(sjPlot)
library(cowplot)

setwd("D:/OneDrive - NIOZ/99_PhD/1_Manuscripts/Movement and Cognition/Submitted Version_7RSOS/Data/data_for_analysis")

# Did the birds learn?
realVsim <- read.csv("obs_exp.csv")
sum((realVsim$Real- realVsim$expected)^2/realVsim$expected) #chi2
result <- chisq.test(realVsim$Real, p = realVsim$prop_expected, simulate.p.value = T) #simulate p values
result

fig_data <- rbind(data.frame(type= "real", score=realVsim$Score, freq= realVsim$Real), data.frame(type= "expected", score=realVsim$Score, freq= realVsim$expected)) #make df for figure

#make and save figure
pdf("D:/OneDrive - NIOZ/99_PhD/1_Manuscripts/Movement and Cognition/Submitted Version_7RSOS/Figs/Fig4_SimulatedVsRealSpaScore.pdf", width= 9, height=5)

hist <- ggplot(fig_data, aes(x=score, y=freq, fill=type, group = type))+
  geom_histogram(position="dodge", stat="identity")+
  scale_fill_manual(values = c("turquoise4", "grey60"))+
  scale_y_continuous(breaks = seq(0,35,5), limits = c(0,35), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0,12,1), limits = c(-0.5,12.5)) +
  labs(x="Score \n(maximum number of consecutive zero-error trails)", y = "Number of Birds")+
  theme(panel.background = element_rect(fill = "white"),
        #panel.grid.major = element_line(linetype="dashed", colour = "grey90"), 
        panel.border= element_rect(colour="black", fill="NA"))
 hist
 
 dev.off()

#Import straightness, speed + cognitive data
Efficiency <- read.csv("TransitMSData_Efficiency.csv", stringsAsFactors = F)
Efficiency <- Efficiency[Efficiency$Duration < 21600,] # Remove transit paths that were over 6 hours long. 

##### Does Spatial ability, sex or motivation affect speed or straightness of transit? ####

#### Speed ####

#Full Model
m1 <- glmer(Speed ~ ScaledDay * Spatial_Score + Sex + ScaledTestOrder + (1|ATLAS), data=Efficiency, family = Gamma(link="log"), na.action="na.fail")
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

speed <- plot_model(m1, "int")+
  scale_x_continuous(expand = c(0, 0), breaks = c(-1.0172767,0.2919676,	1.6012120), labels=c(10,20,30), limits = c(-2, 1.6012120)) + # Cannot change plot to have non scaled days within the model as confidence intervals change slightly, however can match days to scaled days manually on the x axis so that model stays the same.
  scale_y_continuous(expand = c(0, 0), limits=c(0,0.11)) +
  scale_colour_manual(values=c("red", "deepskyblue3"))+
  scale_fill_manual(values=c("red","deepskyblue3"))+
  ylab("Speed (m/s)") +
  xlab("Date (September)") +
  expand_limits(x = 0, y = -2) +
  theme_classic()+
  ggtitle("")#+
  theme(legend.position="none")# 

#### straightness ####

m1 <- glmer(Straightness ~ ScaledDay * Spatial_Score + Sex + ScaledTestOrder + scale(N) + (1|ATLAS), data=Efficiency, family = Gamma(link="log"), na.action="na.fail")
plot(m1)
dredge(m1)
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

straight <- plot_model(m1, "pred")$ScaledDay +
  theme_classic()+
  scale_x_continuous(expand = c(0, 0),breaks = c(-1.0172767,0.2919676,	1.6012120), labels=c(10,20,30), limits = c(-2, 1.6012120)) +
  xlab("Date (September)") +
  scale_y_continuous(expand = c(0, 0), limits=c(0,1)) +
  ggtitle("")

pdf("D:/OneDrive - NIOZ/99_PhD/1_Manuscripts/Movement and Cognition/Submitted Version_7RSOS/Figs/RSOS_Fig4.pdf", width= 9, height=5)
mve <- plot_grid(speed,straight,labels="auto")
mve
dev.off()
