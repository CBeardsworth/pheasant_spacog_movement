# R Code (1) for Analysis on Paper : Spatial cognitive ability is associated with speed of movement in the early exploration of an environment

# This script takes filtered data, subsets it and runs HMMs. 

library(moveHMM)

dat<- read.csv("data/2017_ATLASdata_fullyfiltered.csv", stringsAsFactors = F) #this data is not provided, since it is very large but can be requested. 
d <- dat[dat$ATLAS!=44001004033,] # very questionable track, bird may have been buried. 
d <- d[d$Day==T,]
d <- d[d$Date>=as.Date("2017-09-01")& d$Date<as.Date("2017-10-01"),]

# check how many days have more than six hours data to create a list of birds
t <- as.data.frame(table(d$ATLAS,as.Date(d$RealTime)))
t <- t[t$Freq!=0,]
t<- t[t$Freq >=72,] 
d2 <- as.data.frame(table(t$Var1))
birds2 <- as.character(d2[d2$Freq >=7,]$Var1) #list of birds used for analysis.
d <- d[d$ATLAS%in%birds2,] # subset appropriate birds
#check all birds still have 7 days of data

d <- d[d$ATLAS!=44001004088,] # 4088 has only 6 days once the short duration tracks are removed, easier to explain n decrease in flowchart if this is removed before this step. 
d <- subset(d,medianNort > 98000 | medianEast > 266020) #remove large outliers
# create regular time sets (cannot have many minutes of missing data) - translate into subtracks per individual
d$REF <- d$ATLAS
#plot(d$medianEast,d$medianNort)

# make sure birds have continuous data (no breaks more than 1 hours)
all <- NULL

for(i in unique(d$ATLAS)){
    b <- subset(d,ATLAS==i)
    b <- b[order(b$End),]
    ref <- 1
    start <- b$End[1]
    b$REF[1] <- paste0(b$REF[1],"-",ref)
    for(j in 2:nrow(b)){
        if((b$End[j]-start)<=3600){ # must be less than 1 hours break otherwise it is a new reference ID
            b$REF[j] <- paste0(b$REF[j],"-",ref)
            start <- b$End[j] }
        else{
            ref <- ref+1
            b$REF[j] <- paste0(b$REF[j],"-",ref)
            start <- b$End[j] 
        }}
    all <- rbind(all, b)
}
all <- all[,c(63,1:62)]
d <- all


d<- d[d$REF %in% names(table(d$REF))[table(d$REF)>=72],] # only use tracks that are at least 6 hours long
#remove obvious errors (see figures in transit folder, clear jumps to around F35 then leap back)
length(unique(d$ATLAS))

# Prep and Run HMM ---------------------------------------------------------------------------------------------------------

names(d)[1] <- "ID"
d <- d[order(d$ID,d$TimeGroup), 1:21]
hmm_data <- prepData(d, type="UTM",coordNames=c("medianEast","medianNort")) # 0.5 filter
names(hmm_data)[15] <- "step_to"

allm3 <- list()#list to keep models in

# Number of runs - 25 is standard for HMM papers but it is first necessary to try a number of different initial values to get a 'feel' for the data. Then randomise around these numbers. We ran 5 of these on a few R sessions simultaneously to save time, rather than running 25 on one R session. 
ntries <- 1

for(try in 1:ntries) {
    cat("Iteration", try, "...\n")
    print(paste("start time is:",Sys.time()))
    # Generate starting values randomly for 3 state model
    # initial parameters for gamma and von Mises distributions
    mu0 <- runif(3, min = c(5, 10, 20), max = c(10, 20, 30)) # step mean (two parameters: one for each state)
    sigma0 <- runif(3, min = c(2, 2, 4), max = c(7,10,10))
    zeromass0 <- c(0.01,0.01,0.01) # step zero-mass
    stepPar0 <- c(mu0,sigma0,zeromass0)
    angleMean0 <- runif(3, min = c(0, 0,0), max =c(0.9, 0.9,0.5)) # angle mean
    kappa0 <- runif(3, min = 0, max =1) # angle mean # angle concentration
    anglePar0 <- c(angleMean0,kappa0)
    
    # Fit model
    m <- fitHMM(data=hmm_data,
                nbStates=3,
                stepPar0=stepPar0,
                anglePar0=anglePar0,
                formula=~1)
    
    
    # Save fitted model
    allm3[[try]] <- m
    print(paste("end time is:",Sys.time()))
    #saveRDS(allm3,file="pheas2017.hmmOutput3_190510.RData") #save outputs to compare later.
    
}

#write.csv(hmm_data, "data/bird_data.csv", row.names = F)
