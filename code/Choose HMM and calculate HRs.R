library(moveHMM)
library(ks)
library(OpenStreetMap)
library(ggplot2)
library(cowplot)
library(adehabitatHR)
library(sp)
library(sf)
library(ggspatial)

setwd("D:/OneDrive - NIOZ/99_PhD/1_Manuscripts/Movement and Cognition/Code and Data")
cog <- read.csv("cognition_scores_50HMMbirds.csv") 
tag_ids <- unique(cog$ATLAS) 

d1 <- readRDS("pheas2017.hmmOutput1_190510.RData")
d2 <- readRDS("pheas2017.hmmOutput2_190510.RData")
d3 <- readRDS("pheas2017.hmmOutput3_190510.RData")
d4 <- readRDS("pheas2017.hmmOutput4_190510.RData")
d5 <- readRDS("pheas2017.hmmOutput5_190510.RData")

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
bird_data <- read.csv("birdDat_190510HMM.csv")[,c(1:22)] #import filtered movement data (that was used in creating the HMMs)
bird_data$State <- viterbi(bestm) # add state from hmm model to bird data
bird_data <- merge(bird_data,cog, by="ATLAS") # add cognition data 
bird_data <- bird_data[bird_data$ATLAS %in% tag_ids,] # remove birds that are not used in the analysis
bird_data$time_period <- "beginning" # make time periods for HR analysis (methods), this is  to ensure that the birds are still using the same core areas and can therefore use their experience to learn. 
bird_data[as.integer(bird_data$Date)>15, ]$time_period <- "end"

# Calculate Home Ranges -------------------------------------------------------------------------------------------

kd_info <- NULL

for(i in tag_ids){
    if(nrow(bird_data[bird_data$ATLAS==i & bird_data$time_period=="beginning", ]) > 0){
        if(nrow(bird_data[bird_data$ATLAS==i & bird_data$time_period=="end", ]) > 0){
            
            # beginning and end of september KDE
            df1<- bird_data[bird_data$ATLAS==i, 
                    c("time_period", "medianEast", "medianNort")]
            coordinates(df1) <- ~medianEast+medianNort
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
                            c("time_period", "medianEast", "medianNort")]
            coordinates(df1) <- ~medianEast+medianNort
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
release_pen <- read_sf("ReleasePen2.shp")%>%
    st_transform(4326)%>%
    st_coordinates()%>%
    as.data.frame()

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

setwd("D:/OneDrive - NIOZ/99_PhD/1_Manuscripts/Movement and Cognition/Submitted Version_7RSOS/Figs")

pdf(file = "ESM_KDE_all.pdf",
    width = 18, height = 200)
p <- plot_grid(plotlist = my_plots, ncol=2, nrow=25)
p
dev.off()

#combine some examples 
pdf(file = "Fig2_KDEexamples.pdf",
    width = 15, height = 10)
p <- plot_grid(plotlist = my_plots[c(9,21,22,39)], ncol=2, nrow=2)
p
dev.off()


