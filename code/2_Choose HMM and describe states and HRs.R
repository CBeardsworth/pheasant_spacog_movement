# R Code (2) for Analysis on Paper : Spatial cognitive ability is associated with speed of movement in the early exploration of an environment

#This script combines all 25 HMM models to find the 'best' model (see Michelot et al. 2016). Also describes step length and trajectory durations for MS.

#Packages required
library(moveHMM)
library(ks)
library(OpenStreetMap)
library(ggplot2)
library(cowplot)
library(adehabitatHR)
library(sp)
library(sf)
library(ggspatial)
library(ggsn)
library(amt)

cog <- read.csv("data/cognition_scores_50HMMbirds.csv") 
tag_ids <- unique(cog$ATLAS) #get IDs of birds that are in the analysis (not all birds)

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
table(round(llk, 1))
which.min(aicM)

# # Best model (from fitted models)
bestm <- hmm[[whichBest]]
bestm

#######
bird_data <- read.csv("data/bird_data.csv") #import filtered movement data (that was used in creating the HMMs)
bird_data$State <- viterbi(bestm) # add state from hmm model to bird data
bird_data <- merge(bird_data,cog, by="ATLAS") # add cognition data 
bird_data <- bird_data[bird_data$ATLAS %in% tag_ids,] # remove birds that are not used in the analysis
bird_data$time_period <- "beginning" # make time periods for HR analysis (methods), this is  to ensure that the birds are still using the same core areas and can therefore use their experience to learn. 
bird_data[as.integer(bird_data$Date)>15, ]$time_period <- "end"

# Calculate Home Ranges and make figs ---------------------------------------------------------------------------------------

kd_info <- NULL

for(i in tag_ids){
    if(nrow(bird_data[bird_data$ATLAS==i & bird_data$time_period=="beginning", ]) > 0){
        if(nrow(bird_data[bird_data$ATLAS==i & bird_data$time_period=="end", ]) > 0){
            
            # beginning and end of september KDE
            df1<- bird_data[bird_data$ATLAS==i, 
                    c("time_period", "x", "y")]
            coordinates(df1) <- ~x+y
            kd <- kernelUD(df1, h = "href",same4all=T, grid=600)
            
            temp <- data.frame(ATLAS = i, 
                               prob_overlap_50 = kerneloverlaphr(kd, meth="VI", percent= 50)[2,1],
                               prob_overlap_95 = kerneloverlaphr(kd, meth="VI", percent= 95)[2,1],
                               beg_area_50 = kernel.area(kd, percent =  50, unin="m", unout="ha")[1,1], 
                               beg_area_95 = kernel.area(kd, percent =  95, unin="m", unout="ha")[1,1], 
                               end_area_50 = kernel.area(kd, percent =  50,unin="m", unout="ha")[1,2], 
                               end_area_95 = kernel.area(kd, percent =  95,unin="m", unout="ha")[1,2])
            
            kd_info <- rbind(kd_info,temp)
            
        }
        if(nrow(bird_data[bird_data$ATLAS==i & bird_data$time_period=="end", ]) == 0){
            
            # beginning of september KDE only
            df1<- bird_data[bird_data$ATLAS==i, 
                            c("time_period", "x", "y")]
            coordinates(df1) <- ~x+y
            kd <- kernelUD(df1, h = "href",same4all=T, grid=600)
            
            temp <- data.frame(ATLAS = i, 
                               prob_overlap_50 = NA,
                               prob_overlap_95 = NA,
                               beg_area_50 = kernel.area(kd, percent =  50, unin="m", unout="ha")[1,1], 
                               beg_area_95 = kernel.area(kd, percent =  95, unin="m", unout="ha")[1,1], 
                               end_area_50 = NA, 
                               end_area_95 = NA)
            
            kd_info <- rbind(kd_info,temp)
            
            
        }}}

kd_info <- merge(kd_info, cog, by= "ATLAS")

# start total home range size (Ha)
mean(kd_info$beg_area_95)
sd(kd_info$beg_area_95)

# end total home range size (Ha)
mean(kd_info$end_area_95, na.rm= T)
sd(kd_info$end_area_95, na.rm=T)

# 
mean(kd_info$prob_overlap_50, na.rm= T)
sd(kd_info$prob_overlap_50, na.rm= T)


# Maps -----------------------------------------------------------------------------------------------------
LAT1 =  50.7660 ; LAT2 = 50.7785
LON1 = -3.908 ; LON2 = -3.890

map <- openmap(c(LAT2,LON1), c(LAT1,LON2), zoom = NULL,
               type = c("osm", "stamen-toner", "stamen-terrain","stamen-watercolor", "esri","esri-topo")[6],
               mergeTiles = TRUE)

map <- openproj(map, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
release_pen <- read_sf("data/ReleasePen2.shp")%>%
    st_transform(4326)%>%
    st_coordinates()%>%
    as.data.frame()

#plotting function

plot_conts_foraging <- function(data, i, map){

        if(nrow(data[data$ATLAS==i & data$time_period=="beginning", ]) > 0){
                if(nrow(data[data$ATLAS==i & data$time_period=="end", ]) > 0){
                    
                    #beginning of month contour
                    df1<- data[data$ATLAS==i, c("time_period", "Long", "Lat")]
                    coordinates(df1) <- ~Long+Lat
                    kd <- kernelUD(df1, h = "href", same4all=T, grid=600)
            
                    # Beginning
                    contour_50_beg <- getverticeshr(kd$beginning, percent = 50) %>%
                        st_as_sf()%>%
                        st_set_crs(4326)%>%
                        st_coordinates()%>%
                        as.data.frame()
                    
                    contour_95_beg <- getverticeshr(kd$beginning, percent = 95) %>%
                        st_as_sf()%>%
                        st_set_crs(4326)%>%
                        st_coordinates()%>%
                        as.data.frame()
                    
                    #End
                    contour_50_end <- getverticeshr(kd$end, percent = 50) %>%
                        st_as_sf()%>%
                        st_set_crs(4326)%>%
                        st_coordinates()%>%
                        as.data.frame()
                    
                    contour_95_end <- getverticeshr(kd$end, percent = 95) %>%
                        st_as_sf()%>%
                        st_set_crs(4326)%>%
                        st_coordinates()%>%
                        as.data.frame()
            
                    # Plot
                    title <- as.character(i)
                    autoplot(map)+
                        ggtitle(title)+
                        xlab("Longitude")+
                        ylab("Latitude")+
                        geom_polygon(aes(X,Y, group = L2, col="1-15 Sept"), data=contour_50_beg, fill="blue") +
                        geom_polygon(aes(X,Y, group = L2, col="1-15 Sept"), data=contour_95_beg, fill="blue", alpha=0.2)+
                        geom_polygon(aes(X,Y, group = L2,col = "16-30 Sept"), data=contour_50_end, fill="purple", alpha=0.7) +
                        geom_polygon(aes(X,Y, group = L2, col = "16-30 Sept"), data=contour_95_end, fill="purple", alpha=0.2)+
                        geom_polygon(aes(X,Y, col="Release pen"), data=release_pen, fill=NA,  size=1)+
                        annotation_north_arrow(location="tl",height= unit(1,"cm"), width= unit(0.9, "cm"), 
                                               style=north_arrow_orienteering(text_size=6))+
                        scalebar(x.min = LON1, x.max= LON2, y.min = LAT1, y.max= LAT2, location= "bottomright",  dist = 100,   dist_unit = "m", transform = T,  model = "WGS84", st.bottom=F, border.size = 0.25, st.size=3)+
                        scale_color_manual(values = c("Release pen" = "gold", "1-15 Sept" = "blue", "16-30 Sept"= "purple"),
                                           name = NULL) +
                        theme(panel.background = element_rect(fill = "white"),
                              panel.grid.major = element_line(linetype="dashed", colour = "grey90"), 
                              panel.border= element_rect(colour="black", fill="NA"),
                              legend.key=element_blank(),
                              legend.spacing.y = unit(-0.12, "cm"), 
                              legend.background = element_rect(fill="white",colour="black"), 
                              legend.key.height = unit(0.7, "line"),
                              legend.position = c(0.9,0.9))+
                        guides(color=guide_legend(override.aes=list(fill=NA, size = 0.7)))
                               
                }
            else{
                #beginning of month contour
                df1<- data[data$ATLAS==i, c("time_period", "Long", "Lat")]
                coordinates(df1) <- ~Long+Lat
                kd <- kernelUD(df1, h = "href", same4all=T, grid=600)
                
                # Beginning
                contour_50_beg <- getverticeshr(kd$beginning, percent = 50) %>%
                    st_as_sf()%>%
                    st_set_crs(4326)%>%
                    st_coordinates()%>%
                    as.data.frame()
                
                contour_95_beg <- getverticeshr(kd$beginning, percent = 95) %>%
                    st_as_sf()%>%
                    st_set_crs(4326)%>%
                    st_coordinates()%>%
                    as.data.frame()
                
                # Plot
                title <- as.character(i)
                
                autoplot(map)+
                    ggtitle(title)+
                    xlab("Longitude")+
                    ylab("Latitude")+
                    geom_polygon(aes(X,Y, group = L2, col="1-15 Sept"), data=contour_50_beg, fill="blue") +
                    geom_polygon(aes(X,Y, group = L2, col="1-15 Sept"), data=contour_95_beg, fill="blue", alpha=0.2)+
                    geom_polygon(aes(X,Y, col="Release pen"), data=release_pen, fill=NA,  size=1)+
                    annotation_north_arrow(location="tl",height= unit(1,"cm"), width= unit(0.9, "cm"), 
                                           style=north_arrow_orienteering(text_size=6))+
                    scalebar(x.min = LON1, x.max= LON2, y.min = LAT1, y.max= LAT2, location= "bottomright",  dist = 100,   dist_unit = "m", transform = T,  model = "WGS84", st.bottom=F, border.size = 0.25, st.size=3)+
                    scale_color_manual(values = c("Release pen" = "gold", "1-15 Sept" = "blue"),
                                       name = NULL) +
                    theme(panel.background = element_rect(fill = "white"),
                          panel.grid.major = element_line(linetype="dashed", colour = "grey90"), 
                          panel.border= element_rect(colour="black", fill="NA"),
                          legend.key=element_blank(),
                          legend.spacing.y = unit(-0.12, "cm"), 
                          legend.background = element_rect(fill="white",colour="black"), 
                          legend.position = c(0.9,0.9), legend.key.height = unit(0.7, "line"))+
                    guides(color=guide_legend(override.aes=list(fill=NA)))

            }}}

my_plots <- lapply(tag_ids[order(tag_ids)], FUN=plot_conts_foraging, data=bird_data, map = map)


pdf(file = "figs/ForReviewer_KDE_all.pdf", #for reviewer
    width = 18, height = 200)
p <- cowplot::plot_grid(plotlist = my_plots, ncol=2, nrow=25)
p
dev.off()

#combine some examples 
pdf(file = "figs/Fig2_KDEexamples.pdf",
    width = 15, height = 10)
p <- cowplot::plot_grid(plotlist = my_plots[c(9,21,22,39)], ncol=2, nrow=2)
p
dev.off()

# Calculate length of time in each state ----------------------------------------------------------------------------------------

# Get summary information on different states:Steps
mean(bird_data[bird_data$State==1,]$step, na.rm=T)
sd(bird_data[bird_data$State==1,]$step, na.rm=T)
nrow(bird_data[bird_data$State==1,])

mean(bird_data[bird_data$State==2,]$step, na.rm=T)
sd(bird_data[bird_data$State==2,]$step, na.rm=T)
nrow(bird_data[bird_data$State==2,])

mean(bird_data[bird_data$State==3,]$step, na.rm=T)
sd(bird_data[bird_data$State==3,]$step, na.rm=T)
nrow(bird_data[bird_data$State==3,])

bird_data$RealTime <- as.POSIXct(bird_data$RealTime, tz="GMT")
data <- bird_data
new.data2 <- NULL

for(i in unique(data$ID)){
    b<- subset(data,ID==i) #subset each section of activity 
    b <- b[order(b$RealTime),]
    b$seq <- sequence(rle(as.character(b$State))$lengths) # count from 1 :  transition to another behaviour
    b$seqTotal <- rep(rle(as.character(b$State))$lengths,rle(as.character(b$State))$lengths) #end of behaviour
    Start <- b[b$seq==1,] # first rows
    Finish <- b[b$seq==b$seqTotal,]$RealTime + 299# last row of behaviours' time
    new.data <- cbind(Start,Finish) #start info plus end time of behaviour
    # print(Start$RealTime)
    # print(Finish)
    new.data$StateTime <- as.numeric(as.POSIXct(new.data$Finish))-as.numeric(as.POSIXct(new.data$RealTime)) # total time in this state, ignoring small gaps as they were accounted for before the HMM occurred.
    new.data$SectionLength <- sum(new.data$StateTime) # total amount of time in detected states (cannot take away final from first time as it seems some gaps are present)
    new.data2 <- rbind(new.data2, new.data)
    
}

##### Foraging #####

#foraging only 
foraging <- subset(new.data2, State==3)
foraging$Date <- as.Date(foraging$Date)

# get only foraging paths from original data

all <- NULL
count <- 0
foraging <- foraging[foraging$seqTotal >=3, ] # make sure the traj has 3 or more points

for(i in 1:nrow(foraging)){
    count <- count + 1
    x <- subset(data, ATLAS==foraging$ATLAS[i] & RealTime > foraging$RealTime[i] & RealTime <= foraging$Finish[i])
    x$trajRef <- count
    all <- rbind(all,x)
}

sin <- NULL
for(i in unique(all$trajRef)){
    x <- subset(all, trajRef==i)
    x <- x[order(x$RealTime),]
    if(nrow(x)>=3){
        track1 <- track(x$x,x$y,x$RealTime, crs=("+init=espg:27700"))
        y <- data.frame(ATLAS = x$ATLAS[1], N = nrow(track1), Date = unique(as.Date(track1$t_)), Duration = (as.numeric(track1[nrow(track1),3])-as.numeric(track1[1,3])), Distance = tot_dist(track1), mean.TurnAngle = tac(track1), Sinuosity = sinuosity(track1), Straightness=straightness(track1))
        #print(y)
        sin <- rbind(sin,y)
    }}

dat_Foraging <- sin

##### Transit #####

#transit only 
transit <- subset(new.data2, State==2)
transit$Date <- as.Date(transit$Date)

# get only transit paths from original data

all <- NULL
count <- 0
transit <- transit[transit$seqTotal >=3, ] # make sure the traj has 3 or more points

for(i in 1:nrow(transit)){
    count <- count + 1
    x <- subset(data, ATLAS==transit$ATLAS[i] & RealTime > transit$RealTime[i] & RealTime <= transit$Finish[i])
    x$trajRef <- count
    all <- rbind(all,x)
}

sin <- NULL
for(i in unique(all$trajRef)){
    x <- subset(all, trajRef==i)
    x <- x[order(x$RealTime),]
    if(nrow(x)>=3){
        track1 <- track(x$x,x$y,x$RealTime, crs=("+init=espg:27700"))
        y <- data.frame(ATLAS = x$ATLAS[1], N = nrow(track1), Date = unique(as.Date(track1$t_)), Duration = (as.numeric(track1[nrow(track1),3]) - as.numeric(track1[1,3])), Distance = tot_dist(track1), mean.TurnAngle = tac(track1), Sinuosity = sinuosity(track1), Straightness = straightness(track1))
        #print(y)
        sin <- rbind(sin,y)
    }}

dat_Transit <- sin

#write.csv(dat_Transit, "data/TransitMSData_Efficiency.csv") # This file was later modified to merge individual information on cognition score and mass to this dataset. 

#### Resting #####

#resting only 
resting <- subset(new.data2, State==1)
resting$Date <- as.Date(resting$Date)

# get only resting paths from original data

all <- NULL
count <- 0
resting <- resting[resting$seqTotal >=3, ] # make sure the traj has 3 or more points

for(i in 1:nrow(resting)){
    count <- count + 1
    x <- subset(data, ATLAS==resting$ATLAS[i] & RealTime > resting$RealTime[i] & RealTime <= resting$Finish[i])
    x$trajRef <- count
    all <- rbind(all,x)
}

sin <- NULL

for(i in unique(all$trajRef)){
    x <- subset(all, trajRef==i)
    x <- x[order(x$RealTime),]
    if(nrow(x)>=3){
        track1 <- track(x$x,x$y,x$RealTime, crs=("+init=espg:27700"))
        y <- data.frame(ATLAS = x$ATLAS[1], N = nrow(track1), Date = unique(as.Date(track1$t_)), Duration = (as.numeric(track1[nrow(track1),3])-as.numeric(track1[1,3])), Distance = tot_dist(track1), mean.TurnAngle = tac(track1), Sinuosity = sinuosity(track1), Straightness=straightness(track1))
        #print(y)
        sin <- rbind(sin,y)
    }}

dat_Resting <- sin



# Get summary information on different states:Paths
mean(dat_Resting$Distance)
sd(dat_Resting$Distance)
nrow(dat_Resting)

mean(dat_Foraging$Distance)
sd(dat_Foraging$Distance)
nrow(dat_Foraging)

mean(dat_Transit$Distance)
sd(dat_Transit$Distance)
nrow(dat_Transit)

# Outliers? 
quantile(dat_Transit$Duration/3600, 0.95)
boxplot(dat_Transit$Duration/3600)
summary(dat_Transit$Duration/3600)
sd(dat_Transit$Duration/3600)
nrow(dat_Transit[dat_Transit$Duration/3600 >=6,])
length(unique(dat_Transit$ATLAS))


