# R Code for Analysis on Paper : spatial cognitive ability is associated with the development of movement strategies

# Packages required:

library(lme4)
library(ggplot2)
library(sjPlot)
library(Rmisc)
library(tidyr)
library(cowplot)

# Data required:

Efficiency <- read.csv("data/TransitMSData_Efficiency.csv", stringsAsFactors = F)
Efficiency <- Efficiency[Efficiency$Duration < 21600,] # Remove transit paths that were over 6 hours long.
Efficiency$Date <- as.Date(Efficiency$Date, origin= "1970-01-01", format= "%d/%m/%Y")

# ESM_S3: Tag failures throughout september, show that there wasn't a difference between tags of different spatial ratings.----- 

h <- data.frame(table(Efficiency$SpatialRating,Efficiency$Date))
dates <- data.frame(Var2=rep(seq.Date(as.Date("2017-09-01"),as.Date("2017-09-30"), by=1), each=3), Var1 = rep(c("Good","Ok","Bad"), by=30), Accuracy = rep(c("High","Medium","Low"), by=30))
h$Var1 <- as.character(h$Var1)
h$Var2 <- as.Date(h$Var2)
h2 <- merge(h,dates, by=c("Var2", "Var1"), all.y=T)
h2[is.na(h2$Freq),]$Freq <-0

h2$Accuracy <- factor(h2$Accuracy, levels = c("High","Medium","Low") )


s3 <- ggplot(data=h2, aes(x=as.Date(Var2), y=Freq, col=Accuracy))+
  geom_line(lwd=1.1)+
  xlab("Date")+
  ylab("Number of Birds")+
  scale_color_manual(values=c("royalblue","orange", "red3"))+
  theme(legend.position="none",
        panel.background = element_rect(fill = "white"),
        panel.border= element_rect(colour="black", fill="NA"))

s3

#ggsave("figs/esm_s3_tagSignalsSeptember_bySpatial.jpeg", units="cm",width=16.5,height=10, dpi=600)

##### SI Result 4: Did birds disperse within the study period #### ------------------------------------------------------------

######

data <- read.csv("data/2017_ATLASdata_fullyfiltered.csv", stringsAsFactors = F) #this data is not provided, since it is very large.
rp <- read_sf("data/ReleasePen2.shp")

data <- st_as_sf(data, coords = c("medianEast","medianNort"), crs=27700)
data$rp_distance <- as.numeric(st_distance(data$geometry, rp[1]))
data$rp_within <- as.numeric(st_within(data$geometry, rp[1]))
data[is.na(data$rp_within),]$rp_within <- 0
data[data$rp_within==1,]$rp_distance <- 0
pData <- aggregate(data$rp_distance, list(data$ATLAS, as.Date(data$Date)), mean)
pData2 <- aggregate(data$rp_distance, list(as.Date(data$Date)), mean)
names(pData) <- c("Tag", "Date" , "rp_distance")
names(pData2) <- c("Date" , "rp_distance")


s4 <- ggplot()+
  geom_line(data= pData, aes(y=rp_distance,x = Date, group=Tag),alpha=0.2)+
  geom_point(data= pData, aes(y = rp_distance, x = Date, group=Tag),alpha=0.2)+
  geom_line(data= pData2, aes(y=rp_distance,x = Date), col="purple", size=1, alpha=0.95)+
  scale_x_date(date_breaks = "1 month", date_labels = "%b %d", expand=c(0,0), limits = c(as.Date("2017-08-15"), as.Date("2017-12-01")))+
  scale_y_continuous(expand=c(0,0))+
  # coord_cartesian(ylim=c(0, 500))+
  labs(y="Distance from Release Pen (m)")+
  theme(panel.background = element_rect(fill = "white"),
        #panel.grid.major = element_line(linetype="dashed", colour = "grey90"), 
        panel.border= element_rect(colour="black", fill="NA"))
s4

#ggsave("figs/esm_s4_moveDistancefromRP.jpeg", width=16.5,height=10,units="cm", dpi=600,device="jpeg")


##### SI S5: unique detections on each day to demonstrate sample sizes throughout study period. -------------------------------
# unique birds per day

data <- read.csv("data/2017_ATLASdata_fullyfiltered.csv", stringsAsFactors = F) #this data is not provided, since it is very large.
data$Date <- factor(data$Date, as.character(seq.Date(as.Date("2017-07-22"),as.Date("2018-02-28"), by=1)))

#calculate number of unique tags localised.
all <- as.data.frame(table(data$ATLAS, data$Date))%>%
    mutate(Freq = ifelse(Freq > 0, 1,0), Date=Var2)%>%
    group_by(Date)%>%
    summarise(unique_id = sum(Freq))
all$Date <- as.Date(all$Date)

s5 <- ggplot(all,(aes(y= unique_id, x=Date)))  +
  geom_rect(mapping=aes(xmin=as.Date("2017-09-01"), xmax=as.Date("2017-09-30"), ymin=0,ymax=70), fill="thistle1", col="orchid1", alpha=0.1)+
  geom_line(lwd=1, col="gray4")+
  scale_x_date(date_labels = '%b', date_breaks="1 month")+
  scale_y_continuous(limits=c(0,70), expand=c(0,0), breaks=seq(0,70,10))+
  ylab("Number of birds \nLocalised")+
  xlab("Date (2017-2018)")+
  theme(panel.background = element_rect(fill = "white"),
        #panel.grid.major = element_line(linetype="dashed", colour = "grey90"), 
        panel.border= element_rect(colour="black", fill="NA"))

s5

#ggsave("figs/esm_s5_unique_localisations.jpeg", width=16.5,height=10,units="cm", dpi=600,device="jpeg")


##### SI S7: Density plots --------------------------------------------------------------------------

data <- read.csv("data/bird_data.csv")

d1 <- readRDS("data/pheas2017.hmmOutput1_190510.RData")
d2 <- readRDS("data/pheas2017.hmmOutput2_190510.RData")
d3 <- readRDS("data/pheas2017.hmmOutput3_190510.RData")
d4 <- readRDS("data/pheas2017.hmmOutput4_190510.RData")
d5 <- readRDS("data/pheas2017.hmmOutput5_190510.RData")

hmm <- unique(c(d1,d2,d3,d4,d5))

# # # Vector of maximum log-likelihoods 
# # # (negative sign because moveHMM works with the negative log-likelihood)
llk <- unlist(lapply(hmm, function(m) return(-m$mod$minimum)))
aicM <- unlist(lapply(hmm, AIC))

# # Index of model with largest maximum log-likelihood
whichBest <- which.max(llk)

# # Best model (from fitted models)
bestm <- hmm[[whichBest]]
#plot(bestm)
states <- viterbi(bestm) 

# Colours for states
nstate <- 3

# Estimated step length parameters
stepMean <- bestm$mle$stepPar["mean",]
stepSD <- bestm$mle$stepPar["sd",]

# Estimated turning angle parameters
angleMean <- bestm$mle$anglePar["mean",]
angleCon <- bestm$mle$anglePar["concentration",]

stepShape <- stepMean^2/stepSD^2
stepRate <- stepMean/stepSD^2

df <- data.frame(State = c("State 1","State 2","State 3"), StepShape=stepShape, StepRate = stepRate)
# Grid of step length values, to plot densities
stepgrid <- seq(min(data$step, na.rm = TRUE),
                max(data$step, na.rm = TRUE),
                length = 1000)
mycols <- c(rgb(1,0.8,0.2,0.6),rgb(0.4,0.8,1,0.7),rgb(0.8,0,0.7,0.7))
mycols2 <- c(rgb(0.8,0,0.7,1), rgb(0.4,0.8,1,1),rgb(1,0.8,0.2,1))

# Loop over states

png("figs/esm_s7_distribution_StepandAnglesHMM.png", units="in", width=9,height=7, res=600)
par(mfrow=c(1, 2))

# S7a:
for(s in c(1,3,2)){
        # Indices of observations in state s (excluding steps of length 0)
        ind <- which(states == s & data$step != 0)
        # Histogram of step lengths in state s
        if(s==1){
          hist(data$step[ind], col = mycols[s], border = 1,main="",xaxs="i",yaxs="i",
               xlab = "Step Length (m)", probability = TRUE,
               xlim = c(0,150),
               ylim = c(0,0.2),
               cex.lab=1,
               mgp=c(1.5,0.5,0),
               breaks = seq(0, max(data$step[ind]),length = max(data$step[ind])/10), 
               add=F)}
        #title(xlab="Step Length (m)", line=2, cex.lab=1.2, family="Calibri Light")}
        
        else{
          hist(data$step[ind], col = mycols[s], border = 1,main="",xaxs="i",yaxs="i",
               xlab = "Step length (m)", probability = TRUE,
               xlim = c(0,150),
               ylim = c(0,0.2),
               breaks = seq(0, max(data$step[ind]),length = max(data$step[ind])/10), 
               add=T)}}

      for(s in c(1,3,2)){
        # Indices of observations in state s (excluding steps of length 0)
        ind <- which(states == s & data$step != 0)
        points(stepgrid,
               dgamma(stepgrid, shape = stepShape[s], rate = stepRate[s]),
               col = mycols2[s], type = "l", lwd=3)
      }
mtext(paste("a"),side=3,adj=-0.15, line=0.5, cex=1.5)


# Estimated gamma density for state s

# Grid of turning angle values
anglegrid <- seq(-pi, pi, length = 1000)
#png("AngleGridDensity.png", units="in", width=3.5,height=5, res=300)
for(s in c(1,3,2)) { # Loop over states
  # Indices of observations in state s
  ind <- which(states == s)
  # Histogram of turning angles in state s
  if(s==1){
    hist(data$angle[ind], col =  mycols[s], border = 1, main = "",xaxs="i",yaxs="i",
         xlab = "Turning Angle (Radians)", probability = TRUE,
         mgp=c(1.5,0.5,0),
         ylim=c(0,0.3),
         xlim = c(-pi, pi), breaks = seq(-pi, pi, length = 30))
  }
  else{
    hist(data$angle[ind], col =  mycols[s], border = 1, main = "",xaxs="i",yaxs="i",
         xlab = "Turning Angle (Radians)", probability = TRUE,
         xlim = c(-pi, pi), breaks = seq(-pi, pi, length = 30), add=T)}}
# Estimated von Mises density for state s
for(s in c(1,3,2)){
  points(anglegrid,
         dvm(anglegrid, mu = angleMean[s], kappa = angleCon[s]),
         col = mycols2[s], type = "l",lwd=3)
}
mtext(paste("b"),side=3,adj=-0.15, line= 0.5, cex=1.5)
dev.off()
#ggsave("figs/esm_s7_distribution_StepandAnglesHMM.jpeg", width=16.5,height=10,units="cm", dpi=600,device="jpeg")


##### SI S8: All birds cog scores-------------------------------------------------------------------------------------------

cog <- read.csv("data/cognition_scores.csv")
cog <- cog[!is.na(cog$score),] #remove NAs

s8 <- ggplot(cog, aes(x=score))+
  geom_histogram(fill = "turquoise4", binwidth=1, col="white")+
  scale_y_continuous(breaks = seq(0,60,10), limits = c(0,60), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0,12,1), limits = c(-0.5,12.5)) +
  labs(x="Score \n(maximum number of consecutive zero-error trails)", y = "Number of Birds")+
  theme(panel.background = element_rect(fill = "white"),
        #panel.grid.major = element_line(linetype="dashed", colour = "grey90"), 
        panel.border= element_rect(colour="black", fill="NA"))

#ggsave("figs/esm_s8_distribution_allbirdCog.jpeg", width=16.5,height=10,units="cm", dpi=600,device="jpeg")


