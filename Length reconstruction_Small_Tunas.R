####################################
# Small Tunas Length Reconstruction
# by: Matheus Lourenco
###################################

#cleanning workspace and devices
rm(list = ls())
#dev.off()
#graphics.off()

# R packages
#please check if you have all packages below..

######@> Package list...
#install.packages("dplyr")
library(dplyr)

#install.packages('ggplot2')
library(ggplot2)

#install.packages("gghalves")
library(gghalves)  #violin plot

#install.packages('fitdistrplus')
library(fitdistrplus) #fit distributions

#install.packages("lmtest")
library(lmtest)

#install.packages('actuar')
library(actuar) #density distributions

#install.packages('mgcv')
library(mgcv) #GAM fit

# the development version from github LBSPR
#install.packages("devtools")
#devtools::install_github("AdrianHordyk/LBSPR",force = TRUE)
library(LBSPR)


#directory-Change the directory to the folder where the script and the spreadsheets are located
setwd("C:/Matheus/Universidade/Doutorado/Length reconstruction_Small_Tunas_GAM_GLM")

#dataset

#ICCAT size small tunas
sz1 <- read.csv("t2sz_smttunas_1950-20_sub.csv", header = TRUE, sep = ",")

#Size updated small tunas #REVIZEE, BNDA..
sz2<-read.csv("sizedata_updated(9-12) - s_iccat.csv", header = TRUE, sep = ",")

# Mean Lengths 
sz3<-read.csv("HistoricalLengths_SmallTunas.csv", header = TRUE, sep = ",",encoding = 'latin1')

#Biological parameters
pars=read.csv('tabela_protuna_v5.csv', sep = ',',dec = '.')
#------------------------------------------------------------------------------


#First ICCAT database
#Subsetting for small tunas for Brazilian flag
sz1sub<- sz1 %>%
  dplyr::filter(SpeciesCode %in% c('BLF','BLT','BON','BRS','FRI',
                                   'KGM','LTA','WAH','DOL'), 
                SampAreaCode=="BIL96")


#Latitude and longitude are not in absolute values
#We have to correct them by quadrant area# Quadrants 3 and 4 are Brazil
lat.q <- dplyr::recode(sz1sub$QuadID, "1" = 1, "2" = -1 , "3" = -1 , "4" = 1,"0"= 0)
table(lat.q)
lon.q <- dplyr::recode(sz1sub$QuadID, "1" = 1, "2" = 1 , "3" = -1 , "4" = -1, "0"= 0)
table(lon.q)

sz1sub= sz1sub %>% 
  dplyr::mutate(lat = Lat*lat.q, lon = Lon*lon.q) %>% 
  dplyr::select(-c(QuadID, Lon,Lat)) 

#get rid of lat and lon zeros
for (i in 1:length(sz1sub$lat)) {
  
  if (sz1sub$lat[i]==0 & sz1sub$lon[i]==0) {
    
    sz1sub$lat[i]=NA; sz1sub$lon[i]=NA
  }
}

#creating region column
sz1sub$region= "ALL"
sz1sub$region[sz1sub$lat>=-1 & sz1sub$lat<=12 & sz1sub$lon<= -45 & sz1sub$lon>= -51]= 'N'
sz1sub$region[sz1sub$lat>=-18 & sz1sub$lat<=8 & sz1sub$lon<=-20 & sz1sub$lon>=-46]= 'NE'
sz1sub$region[sz1sub$lat>=-25 & sz1sub$lat<=-18 & sz1sub$lon>=-48 & sz1sub$lon<=-20]= 'SE'
sz1sub$region[sz1sub$lat< -25]= 'S'


#data frame with frequencies
sz1sub= data.frame(yr=sz1sub$YearC,region=sz1sub$region,codsp=sz1sub$SpeciesCode,
                   codgr=sz1sub$GearGrpCode,class=sz1sub$ClassFrq,n=sz1sub$Nr) 

#expand class frequency to all "Nr" individuals
sz1sub= sz1sub %>% 
  tidyr::uncount(n)

#final data frame with Specie cod, cod gear, region, length and Source
sz1sub= data.frame(yr=sz1sub$yr,region=sz1sub$region,codsp=sz1sub$codsp,codgr=sz1sub$codgr,
                   fl=sz1sub$class,source="ICCAT")
#---------------------------------------------------

#Second BNDA, REVIZEE, GUSTAVO, TATIANA Sources
#Size updated subset with Specie cod, cod gear, region,length and Source
sz2sub=data.frame(codsp=sz2$Cod_ICCAT,fl=sz2$FL.cm.,
                        yr=sz2$Year,source=sz2$Source)


#Adding gear for each project
sz2sub$codgr='UN'
sz2sub$codgr[sz2sub$source=='Revizee N']= 'GN'
sz2sub$codgr[sz2sub$source=='Revizee NE']= 'GN'
sz2sub$codgr[sz2sub$source=='BNDA']= 'LL'
sz2sub$codgr[sz2sub$source=='ODB_Protuna']= 'LL'
sz2sub$codgr[sz2sub$source=='Tatiana_Noronha']= 'HL'
sz2sub$codgr[sz2sub$source=='RS_Gustavo']= 'HL'


#Adding region for each project
sz2sub$region='ALL'
sz2sub$region[sz2sub$source=='Revizee N']= 'N'
sz2sub$region[sz2sub$source=='Revizee NE']= 'NE'
sz2sub$region[sz2sub$source=='BNDA']= 'ALL'
sz2sub$region[sz2sub$source=='ODB_Protuna']= 'ALL'
sz2sub$region[sz2sub$source=='Tatiana_Noronha']= 'NE'
sz2sub$region[sz2sub$source=='RS_Gustavo']= 'S'

#changing program names #tatiana ,Gustavo...
sz2sub$source[sz2sub$source=='Tatiana_Noronha']= 'Noronha'
sz2sub$source[sz2sub$source=='RS_Gustavo']= 'Tubarão Azul'
sz2sub$source[sz2sub$source=='ODB_Protuna']= 'Protuna-ODB'
sz2sub$source[sz2sub$source=='Revizee N']= 'REVIZEE' #standardizing Revizee source
sz2sub$source[sz2sub$source=='Revizee NE']= 'REVIZEE'


#final data frame with Specie cod, cod gear, region, length and Source
sz2sub=data.frame(yr=sz2sub$yr,region=sz2sub$region,
                  codsp=sz2sub$codsp,codgr=sz2sub$codgr,
                  fl=sz2sub$fl,source=sz2sub$source)
#----------------------------------------------------

# Frequency of lengths pooled together (ICCAT+ Size Update)
szfreqobs=rbind(sz1sub,sz2sub)
szfreqobs= szfreqobs[order(szfreqobs$yr),] #order by year

unique(szfreqobs$codsp) #species list

#-----------------------------------------------------

#table with mean lengths from projects (PROTUNA, Revizee, ICCAT) and from literature

#sumarizing observed mean fork length by year, region, specie and gear
sz3sub = sz3 %>% 
  group_by(year, region, codsp, codgr) %>%
  summarize(
    mean = round(mean(flmean..cm.), 2),
    sd = round(sd(flmean..cm.), 2),
    cv = round(sd(flmean..cm.) / mean(flmean..cm.), 2),
    min = round(min(flmean..cm.), 2),
    max = round(max(flmean..cm.), 2)
  ) %>%
  mutate(source = 'Literature') %>% 
     rename(
       yr = year,
       region = region,
       codsp = codsp,
       codgr = codgr,
       mean = mean,
       sd = sd,
       source = source
      )



#creating the final observed mean length (Lengths from literature,Protuna,BNDA,ICCAT,Gustavo...)
szmeanobs = szfreqobs %>%
  group_by(yr, region, codsp, codgr, source) %>%
  summarize(
    mean = round(mean(fl), 2),
    sd = round(sd(fl), 2),
    cv = round(sd(fl) / mean(fl), 2),
    min = round(min(fl), 2),
    max = round(max(fl), 2)
                          ) %>%
              bind_rows(sz3sub) %>%
              arrange(yr) %>%
              mutate(mean = ifelse(is.nan(mean), NA, mean))


#assuming that inicial (1950) lengths when fishing was incipient
szmeansub=data.frame(yr=1950,region="ALL",
                     codsp= unique(szmeanobs$codsp[szmeanobs$codsp!=""]), 
                     codgr="UN",source='Literature',
                     mean=NA, 
                     sd=NA,cv=NA,min=NA,max=NA)

#Final observed mean size 
szmeanobs= szmeanobs %>%
  dplyr::bind_rows(szmeansub) %>%
  dplyr::arrange(yr)
#---------------------------------------------------------------


#adding ICCAT codes for each specie in 'par'
pars= pars %>%
  dplyr::mutate(codsp=plyr::revalue(Species,
                                    c( "Thunnus atlanticus"= "BLF",
                                       "Auxis thazard" = "FRI",
                                       "Auxis rochei"= "BLT",
                                       "Euthynnus alleteratus" = "LTA",
                                       "Sarda sarda" = "BON",
                                       "Scomberomorus brasiliensis" = "BRS",
                                       "Scomberomorus regalis" = "CER",
                                       "Scomberomorus cavalla" = "KGM",
                                       "Acanthocybium solandri"= "WAH",
                                       "Coryphaena hippurus"= "DOL")))
#-----------------------------------------------------------------------


#-------------------------------
#exploratory analysis of lengths 
#-------------------------------
#distribution type plot

p1<-ggplot2::ggplot(szfreqobs, aes(x=fl, fill=codsp)) + 
            geom_density(alpha=.5,aes(colour=codsp),linewidth=1)+
  labs(x="Fork Length (cm)",y="Density",fill="",color="")+
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))

p1

ggsave("Length_Density_Small_Tunas_1.png", plot = p1, device = "png", units = "cm",
       width = 20, height = 12)

#kernel density plot
sp<- unique(szfreqobs$codsp)
#colours=heat.colors(length(sp),alpha = 0.5,rev = TRUE)
colours=rainbow(length(sp),alpha = 0.65,rev = FALSE)

#colours=c('tomato1','deepskyblue4','springgreen','firebrick',
#          'deepskyblue','black','dimgray', 'gold')

d=list(BLF= density(szfreqobs$fl[szfreqobs$codsp=='BLF']),
     FRI= density(szfreqobs$fl[szfreqobs$codsp=='FRI']),
     WAH= density(szfreqobs$fl[szfreqobs$codsp=='WAH']),
     BRS= density(szfreqobs$fl[szfreqobs$codsp=='BRS']),
     KGM= density(szfreqobs$fl[szfreqobs$codsp=='KGM']),
     CER= density(szfreqobs$fl[szfreqobs$codsp=='CER']),
     LTA= density(szfreqobs$fl[szfreqobs$codsp=='LTA']),
     DOL= density(szfreqobs$fl[szfreqobs$codsp=='DOL']))


{jpeg("Length_Density_Small_Tunas_2.jpg",
     width=17,height=12,
     units="cm",
     res=500,quality=150,
     antialias = "cleartype")

#Graphical parameters
par(mar=c(3.3,3.2,0.1,0.2),oma=c(0.3,0.3,0.5,0.5),xpd=F,
    mgp=c(2.3,1,0), cex.lab=1,cex.axis=0.8, bty="n",las=1,lwd=2.5,
    col.axis='gray30',col.lab='gray10')

for (i in 1:length(sp)) {
plot(d[[i]],
     xlab = 'Fork Length (cm)', ylab = "Density",
     xlim=c(0,195), ylim=c(0,0.13),
     main = '',bty='L',
     col=colours[i])
  polygon(d[[i]], col=colours[i], border=colours[i])
  par(new=TRUE)
}

legend('topright', legend = sp, fill=colours,bty='n',border =colours) 

dev.off()
}
#-----------------------------

#Half violin plot- Densities of lengths for each year

# simple half violin plot (gghalve+ggplot)

p3<- ggplot() +
  gghalves::geom_half_violin(data=szfreqobs, 
                             aes(x=codsp, y=fl,linetype='Density'),col='gray30',fill='gray30',
                   side = "r", 
                   #show.legend = TRUE,
                   #fill=,
                   nudge = 0,
                   scale = 'width',
                   na.rm = TRUE,
                   lwd=1,
                   trim = TRUE)+
                   #draw_quantiles = c(0.5))+
  geom_point(data=dplyr::filter(szmeanobs,codsp!=""),
               aes(x=codsp,y=mean,col=source))+
  #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Specie",y="Fork length (cm)",fill="",color="",linetype='')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p3

ggsave("Length_Density_Small_Tunas_3.png", plot = p3, device = "png", units = "cm",
       width = 23, height = 13)
#------------------------------

#Half violin plot
p4<- ggplot() +
  gghalves::geom_half_violin(data=szfreqobs, 
                             aes(x=yr,y=fl,linetype='Density'),col='gray20',fill='gray20',
                             alpha=0.8,
                             side = "r", 
                             #show.legend = TRUE,
                             position = "identity",
                             #fill=,
                             #nudge = 0,
                             scale = 'count',
                             na.rm = TRUE,
                             lwd=1,
                             trim = TRUE)+
  facet_wrap(~codsp,scales= "free_y")+
  #draw_quantiles = c(0.5))+
  geom_point(data=dplyr::filter(szmeanobs,codsp!=""),
              aes(x=yr,y=mean,col=source),size=2)+
  #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)",fill="",color="",linetype='')+
  scale_y_continuous() +
  scale_x_continuous(breaks = seq(1950, 2022, 30),limits = c(1950, 2022)) +
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p4

ggsave("Length_Density_Small_Tunas_4.png", plot = p4, device = "png", units = "cm",
       width = 24, height = 14)
#-------------------------------


#boxplot type + Mean lengths (dots)
p5<- ggplot() +
  geom_point(data=dplyr::filter(szmeanobs,codsp!="" & codsp!="CER"),
             aes(x=factor(yr),y=mean,col=source),size=1.5)+
  geom_boxplot(data=dplyr::filter(szfreqobs,codsp!="" & codsp!="CER"), 
    aes(x=factor(yr), y=fl,fill=source),
    #show.legend = TRUE)+
    #position = "identity",
    lwd=0.55,
    outlier.shape=8,
    outlier.size = 0.5,
    outlier.alpha = 0.5)+
  geom_point(data=dplyr::filter(szmeanobs,codsp!="" & codsp!="CER"),
    aes(x=factor(yr),y=mean,col=source),size=1.5)+
  facet_wrap(~codsp, scales = "free_y") +
  labs(x="Year",y="Fork length (cm)",color="Mean Length",fill="Frequency")+
  scale_y_continuous() +
  #scale_x_discrete(beaks = seq(1950, 2022, 10),limits = c(1950, 2022)) +
  scale_x_discrete(breaks=c(factor(as.character(seq(1950,2024,10)))),expand = c(0,0))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p5

ggsave("Length_boxplot_Small_Tunas_1.png", plot = p5, device = "png", units = "cm",
       width = 31, height = 16)
#------------------------------

#mean and sd plot (error bar)

p6<-ggplot(dplyr::filter(szmeanobs,codsp!=""), 
       aes(x=yr,y=mean,fill=source,col=source),size=2) + 
  #geom_line(size=1,aes(fill=source)) +
  geom_point(size=1.8)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=2,
                position="identity",size=1.2)+
  facet_wrap(~codsp,scales = "free_y")+
  labs(y="Mean Length (cm)",x="Year",col='',fill="")+
  scale_y_continuous() +
  scale_x_continuous(breaks = seq(1950, 2022, 30),limits = c(1950, 2022)) +
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p6

ggsave("Length_mean_sd_Small_Tunas.png", plot = p6, device = "png", units = "cm",
       width = 32, height = 18)
#-------------------------------


#functional relationship/ colinearity analysis (all species)

#shape (all sp) by source
p7<-ggplot(szfreqobs, 
           aes(x=fl,fill=source),size=2) + 
  #geom_line(size=1,aes(fill=source)) +
  geom_histogram(size=2,bins = 30)+
  facet_wrap(~codsp,scales = "free")+
  labs(y="Frequency",x="Fork Length (cm)",col='',fill="")+
  scale_y_continuous() +
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p7

ggsave("Length_histogram_Small_Tunas_1.png", plot = p7, device = "png", units = "cm",
       width = 24, height = 14)
#------------------------------

#shape (all sp) by gear
p8<-ggplot(szfreqobs, 
           aes(x=fl,fill=codgr),size=2) + 
  #geom_line(size=1,aes(fill=source)) +
  geom_histogram(size=2,bins = 30)+
  facet_wrap(~codsp,scales = "free")+
  labs(y="Frequency",x="Fork Length (cm)",col='',fill="")+
  scale_y_continuous() +
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p8

ggsave("Length_histogram_Small_Tunas_2.png", plot = p8, device = "png", units = "cm",
       width = 24, height = 14)
#------------------------------


#shape (all sp) by year
p9<-ggplot(szfreqobs, 
           aes(x=fl,fill=factor(yr)),size=2) + 
  #geom_line(size=1,aes(fill=source)) +
  geom_histogram(size=2,bins = 30)+
  facet_wrap(~codsp,scales = "free")+
  labs(y="Frequency",x="Fork Length (cm)",col='',fill="")+
  scale_y_continuous() +
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p9

ggsave("Length_histogram_Small_Tunas_3.png", plot = p9, device = "png", units = "cm",
       width = 24, height = 14)
#------------------------------

#shape (all sp) by region
p10<-ggplot(szfreqobs, 
           aes(x=fl,fill=region),size=2) + 
  #geom_line(size=1,aes(fill=source)) +
  geom_histogram(size=2,bins = 30)+
  facet_wrap(~codsp,scales = "free")+
  labs(y="Frequency",x="Fork Length (cm)",col='',fill="")+
  scale_y_continuous() +
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p10

ggsave("Length_histogram_Small_Tunas_4.png", plot = p10, device = "png", units = "cm",
       width = 24, height = 14)
#--------------------------------


# exploratory for mean lengths of all species
#Boxplot (all sp) by source
p11<-ggplot(dplyr::filter(szmeanobs,codsp!=""), 
           aes(x=source,y=mean,col=source),size=2) + 
  geom_boxplot(lwd=1,na.rm = TRUE)+
  geom_jitter(color="black", size=0.5, alpha=0.9)+
  facet_wrap(~codsp,scales = "free")+
  labs(y='Mean Length (cm)',x="Source",col='',fill="")+
  scale_y_continuous() +
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p11

ggsave("Length_boxplot_Small_Tunas_2.png", plot = p11, device = "png", units = "cm",
       width = 32, height = 18)

#boxplot (all sp) by gear
p12<-ggplot(dplyr::filter(szmeanobs,codsp!=""), 
            aes(x=codgr,y=mean,col=codgr),size=2) + 
  geom_boxplot(lwd=1,na.rm = TRUE)+
  geom_jitter(color="black", size=0.5, alpha=0.9)+
  facet_wrap(~codsp,scales = "free")+
  labs(y='Mean Length (cm)',x="Gear",col='',fill="")+
  scale_y_continuous() +
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p12

ggsave("Length_boxplot_Small_Tunas_3.png", plot = p12, device = "png", units = "cm",
       width = 32, height = 18)


#boxplot (all sp) by year
p13<-ggplot(dplyr::filter(szmeanobs,codsp!=""), 
            aes(x=factor(yr),y=mean),size=2) + 
  geom_boxplot(lwd=1,na.rm = TRUE)+
  geom_jitter(color="black", size=0.5, alpha=0.9)+
  facet_wrap(~codsp,scales = "free")+
  labs(y='Mean Length (cm)',x="Year",col='',fill="")+
  scale_y_continuous() +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p13

ggsave("Length_boxplot_Small_Tunas_4.png", plot = p13, device = "png", units = "cm",
       width = 24, height = 16)
#--------------------------------

#Functional relationship between variables

#scatter plot (all sp) by year
p14<-ggplot(dplyr::filter(szmeanobs,codsp!=""), 
            aes(x=yr,y=mean),size=2) + 
  geom_point(size=1,na.rm = TRUE)+
  geom_smooth(method = 'loess',  color="red")+
  facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='',fill="")+
  scale_y_continuous() +
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p14

ggsave("Length_functional_Small_Tunas_1.png", plot = p14,device = "png", units = "cm",
       width = 24, height = 16)
#------------------------------


#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
# models by species
# BLF model - Length Reconstruction 
#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
sp='BLF'

#-------------------------------------
#find the best distribution for specie
#-------------------------------------
dist=data.frame(yr=NULL,codsp=NULL,dist=NULL,codgr=NULL, #empty frame
                mean=NULL,shape=NULL,cv=NULL,aic=NULL) #Shape for gamma is like the mean

for (i in sp) {
  
  years=unique(szfreqobs$yr[szfreqobs$codsp==i])
  
  for (j in years) {
    
    #frequency for each year  
    sz= szfreqobs$fl[szfreqobs$codsp==i & szfreqobs$yr==j]
    gear= unique(szfreqobs$codgr[szfreqobs$codsp==i & szfreqobs$yr==j])[1]
    
    if (length(sz)>5 & (length(unique(sz))==1)==FALSE) {
      
      #graphical analysis
      #fitdistrplus::descdist(sz, discrete = FALSE,boot = 1000)
      
      #fitting distributions
      fit_1 <- fitdistrplus::fitdist(sz, "norm")
      fit_2<- fitdistrplus::fitdist(sz, "invgauss",start=list(mean=mean(sz),shape= 1))
      fit_3 <- fitdistrplus::fitdist(sz, "gamma")
      
      dat=data.frame(yr=j,codsp=i,
                     dist=c(fit_1$distname,fit_2$distname,fit_3$distname),
                     codgr=gear,
                     mean=c(fit_1$estimate[1],fit_2$estimate[1],fit_3$estimate[1]), 
                     shape=c(fit_1$estimate[2],fit_2$estimate[2],fit_3$estimate[2]),
                     cv=c(fit_1$estimate[2]/fit_1$estimate[1],
                          fit_2$estimate[2]/fit_2$estimate[1],
                          fit_3$estimate[2]/fit_3$estimate[1]),
                     aic=c(fit_1$aic,fit_2$aic,fit_3$aic))
      
      dist=rbind(dist,dat)
    }
  }
}

dist=dist %>%
  dplyr::group_by(yr,codsp,codgr) %>% 
  slice_min(aic, n = 1)  %>%
  dplyr::arrange(dist) 

#Table with the best distribution
table(dist$codsp,dist$dist)  
distribution='gamma'

#----------------------------------------------------
#predict 1950 length when there was incipient fishing
#---------------------------------------------------
#Assuming 10% higher than the current mean length level
avg1950=1.1 * mean(szmeanobs$mean[szmeanobs$codsp==sp], na.rm = TRUE)

#put in the observed mean lengths
szmeanobs$mean[szmeanobs$codsp==sp & szmeanobs$yr==1950]= round(avg1950,2)
#----------------------------------

#selecting variables .. subseting
szmeanobssub= szmeanobs[,c("yr",'region','codsp','codgr','mean')]

#removing extreme values
szmeanobssub$mean[szmeanobssub$codsp=='BLF' & szmeanobssub$mean<50] =NA
szmeanobssub$mean[szmeanobssub$codsp=='DOL' & szmeanobssub$mean<75] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean>100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean<60] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean<100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean>180] =NA

#get only non-NA values
szmeanobssub=szmeanobssub[is.na(szmeanobssub$mean)==FALSE & szmeanobssub$codsp==sp,]

#grouping factors most likely (Fredou et al., 2021)
szmeanobssub$region[szmeanobssub$region=='ALL']='SE'
szmeanobssub$region[szmeanobssub$region=='S']='SE'
szmeanobssub$codgr[szmeanobssub$codgr=='BB']='UN'

#transforming columns in to factors
szmeanobssub$region= factor(szmeanobssub$region)
szmeanobssub$codgr= factor(szmeanobssub$codgr)
szmeanobssub$codsp= factor(szmeanobssub$codsp)


# model GAM selection----------------
mod1<- mgcv::gam(mean ~ s(yr),
                   data= szmeanobssub,
                   method = "REML" ,
                   family =Gamma(link='log'))

mod2<- mgcv::gam(mean ~ s(yr)+region,
                 data= szmeanobssub,
                 method = "REML" ,
                 family =Gamma(link='log'))

mod3<- mgcv::gam(mean ~ s(yr)+codgr, 
                 data= szmeanobssub,
                 method = "REML" ,
                 family =Gamma(link='log'))

mod4<- mgcv::gam(mean ~ s(yr)+region+codgr,
                 data= szmeanobssub,
                 method = "REML" ,
                 family =Gamma(link='log'))
#aic check
mod1$aic
mod2$aic
mod3$aic
mod4$aic

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)


#final model
modgam<- mgcv::gam(mean ~ s(yr)+codgr, 
                   sp=0.01,
                   data= szmeanobssub,
                   method = "REML" ,
                   family =Gamma(link='log'))

#assign GAM model to the specie
modgam_blf=modgam

#diagnostic
summary(modgam)

par(mfrow=c(2,2))
gam.check(modgam)# model check

#residual 
#homocedasticity
lmtest::bptest(modgam)

#normality
shapiro.test(resid(modgam))

# Check overall concurvity
round(concurvity(modgam, full = TRUE),2)

plot.gam(modgam, # plot partial effects
     rug=TRUE,       #puts the x variable
     #select = c(1), #select some partial effects
     pages = 1,
     all.terms = TRUE,
     residuals = TRUE,
     #cex=2,
     #lwd=2,
     se=TRUE,
     #shade = TRUE,
     #shade.col = "lightblue", #Shading intervals
     shift = coef(modgam)[1]) #shift the scale by


#----------------------
#predictions "year"
newdataexpand= expand.grid(yr=seq(1950,2022,1),
                          region="ALL",
                          codgr=szmeanobssub$codgr)

pred<- data.frame(predict(modgam, newdata = newdataexpand,
                          type = "response", se.fit=TRUE))

pred$up<-pred$fit+qnorm(0.975)*pred$se.fit #confidence intervals
pred$lw<-pred$fit-qnorm(0.975)*pred$se.fit

gam=cbind(type='GAM',
          codsp=sp,
          family=modgam$family$family,
          link=modgam$family$link,
          formula=paste0(c(paste(modgam$formula)[2],
                           paste(modgam$formula)[1],
                           paste(modgam$formula)[3]),
                         collapse = ' '),
          newdataexpand,
          pred,
          aic=modgam$aic,
          r.sq=summary(modgam)$r.sq,
          dev.ex=summary(modgam)$dev.expl)


#extract mean of all levels
szmeanpred = gam %>%
  group_by(type, codsp, family, link, formula, yr) %>%
  summarize(across(c('fit', 'se.fit', 'up', 'lw', 'aic', 'r.sq', 'dev.ex'),
    ~ round(mean(.), 2)))

#creating data frame
szmeanpred_blf=szmeanpred #size mean predicted for all factors combined
#----------------------------------------------

#partial effects 
newdata= szmeanobssub   #desired partial effect

pareffect<- data.frame(predict(modgam, newdata = newdata,
                               type = "terms", se.fit=TRUE)) #terms gets the effects

#intervals for the partial effect
pareffect$up.codgr<-pareffect$fit.codgr+qnorm(0.975)*pareffect$se.fit.codgr #confidence intervals
pareffect$lw.codgr<-pareffect$fit.codgr-qnorm(0.975)*pareffect$se.fit.codgr
pareffect$codgr= szmeanobssub$codgr

# Extract mean of all levels
pareffect = pareffect[grepl("codgr", names(pareffect))]
pareffect = pareffect %>%  # select only the desired partial effect
  group_by(codgr) %>%
  summarize(across(c("fit.codgr", "se.fit.codgr", "up.codgr", "lw.codgr"),
    ~ round(mean(.), 2))) %>%
  rename(effect = codgr,
    fit = fit.codgr, 
    se.fit = se.fit.codgr,
    up = up.codgr,
    lw = lw.codgr) %>%
  mutate(codsp = sp, var = "Gear")

pareffect_blf = pareffect


#----------------------------------
#simulating frequency distributions
#----------------------------------

#estimating FL distributions from reconstructed mean length
flloop=data.frame(yr=NULL,codsp=NULL,fl=NULL,source=NULL) #empty data frame
sp=sp  #list of species for loop

#coefficient of variation for the corresponding distribution 
cv=mean(dist$cv[dist$codsp==sp & dist$dist==distribution]) #the shape is like mean for gamma, so the cv is corrupted

for (j in 1:length(sp)) {  #loop through all species
    #codes
  years=1950:2022
  n=100 #number of simulate individuals
  
     for (i in 1:length(years)) { #loop through each year of available information
    
    l= data.frame(yr=rep(years[i],n),
                  codsp=rep(sp[j],n),
                  fl=rgamma(n,  #distribution 
                           shape=szmeanpred$fit[szmeanpred$yr==years[i] & szmeanpred$codsp==sp], 
                           scale=szmeanpred$fit[szmeanpred$yr==years[i] & szmeanpred$codsp==sp]*cv),
                  source=rep(as.character('Pred'),n))
    
     flloop=rbind(flloop,l)
    }
 }


# size frequency predicted data frame
szfreqpred= flloop
szfreqpred_blf=szfreqpred
#-------------------------------------


#-----------------------------------------
#Observed SPR vs Reconstructed SPR (LBSPR)
#-----------------------------------------
mypars <- new("LB_pars",verbose=F) #blank parameters object
mypars@Species <- sp
mypars@Linf <- pars$linf[pars$codsp==sp]
mypars@L50 <-  pars$l50[pars$codsp==sp] 
mypars@L95 <- 1.1*pars$l50[pars$codsp==sp] 
mypars@MK <-  0.94/pars$k[pars$codsp==sp] #M=0.94 (Freire et al., 2005)
mypars@L_units <- "cm"


#SPR from observed length frequencies 
#getting size classes and counts
szfreqobssub= szfreqobs[szfreqobs$codsp==sp,]
szfreqobssub= szfreqobssub[szfreqobssub$yr!= 2017,] #removing year with no data

szloop=data.frame(codsp=NULL,yr=NULL,lmids=NULL,n=NULL)

  years=unique(szfreqobssub$yr[szfreqobssub$codsp==sp])
  
  for (i in years) {
    h=hist(szfreqobssub$fl[szfreqobssub$yr==i],
           breaks=seq(from=min(szfreqobssub$fl[szfreqobssub$codsp==sp])-5,
                      to=max(szfreqobssub$fl[szfreqobssub$codsp==sp])+5, by=5), #binwidth
           xlab = 'Fork length (cm)',
           main = paste(j))
    freq=data.frame(codsp=sp,yr=i,lmids=h$mids,n=h$counts)    
    szloop=rbind(szloop,freq)
  }

  #LBSPR lengths format
  szfreqobssub =szloop %>% #removing first row 
        tidyr::spread(., yr, n) %>%  #years in columns
          `rownames<-`(.[,2]) %>%   #years in rownames
        dplyr::select(-codsp,-lmids) 
  
#observed length frequencies object
mylengthsobs <- new("LB_lengths",LB_pars=mypars,dataType="freq") #blank length object

mylengthsobs@LMids<- as.numeric(rownames(szfreqobssub)) #mid points
mylengthsobs@LData<- as.matrix(szfreqobssub)    #frequency table
mylengthsobs@L_units<- "cm"  #units
mylengthsobs@Years<- as.numeric(colnames(szfreqobssub)) #list of years
mylengthsobs@NYears<- length(as.numeric(colnames(szfreqobssub))) #number of years

#plotSize(mylengthsobs)
fitlbspr= LBSPRfit(mypars, mylengthsobs,verbose = FALSE)
#plotSize(fitlbspr)
#plotEsts(fitlbspr)

lbsprobs= data.frame(yr=fitlbspr@Years,  #years 
                   codsp=sp,
                   rawsl50=fitlbspr@SL50, #raw sl50
                   sl50sd=sqrt(fitlbspr@Vars[,1]),  #sl50 sd
                   sl50= fitlbspr@Ests[,1],   #smooth sl50
                   rawsl95=fitlbspr@SL95,
                   sl95sd=sqrt(fitlbspr@Vars[,2]),
                   sl95= fitlbspr@Ests[,2],
                   rawfm=fitlbspr@FM,
                   fmsd=sqrt(fitlbspr@Vars[,3]),
                   fm= fitlbspr@Ests[,3],
                   rawspr=fitlbspr@SPR,
                   sprsd=sqrt(fitlbspr@Vars[,4]),
                   spr= fitlbspr@Ests[,4],
                   source="Obs") #observed length frequency

#----------------------------------------------------------

#SPR from Reconstructed length frequencies 
#getting size classes and counts
szfreqpredsub=szfreqpred[szfreqpred$codsp==sp,]
#szfreqrecsub= szfreqrecsub[szfreqrecsub$yr!= 2017,] #removing year with no data

szloop=data.frame(codsp=NULL,yr=NULL,lmids=NULL,n=NULL)

years=unique(szfreqpredsub$yr[szfreqpredsub$codsp==sp])

for (i in years) {
  h=hist(szfreqpredsub$fl[szfreqpredsub$yr==i],
         breaks=seq(from=min(szfreqpredsub$fl[szfreqpredsub$codsp==sp])-5,
                    to=max(szfreqpredsub$fl[szfreqpredsub$codsp==sp])+5, by=5), #binwidth
         xlab = 'Fork length (cm)',
         main = paste(j))
  freq=data.frame(codsp=sp,yr=i,lmids=round(h$mids,1),n=h$counts)    
  szloop=rbind(szloop,freq)
}

#LBSPR lengths format
szfreqpredsub =szloop %>% #removing first row 
  tidyr::spread(., yr, n) %>%  #years in columns
  `rownames<-`(.[,2]) %>%   #years in rownames
  dplyr::select(-codsp,-lmids) 

#observed length frequencies object
mylengthspred <- new("LB_lengths",LB_pars=mypars,dataType="freq",verbose=F) #blank length object

mylengthspred@LMids<- as.numeric(rownames(szfreqpredsub)) #mid points
mylengthspred@LData<- as.matrix(szfreqpredsub)    #frequency table
mylengthspred@L_units<- "cm"  #units
mylengthspred@Years<- as.numeric(colnames(szfreqpredsub)) #list of years
mylengthspred@NYears<- length(as.numeric(colnames(szfreqpredsub))) #number of years

#plotSize(mylengthsobs)

fitlbspr= LBSPRfit(mypars, mylengthspred,verbose = FALSE)
#plotSize(fitlbspr)
#plotEsts(fitlbspr)


lbsprpred= data.frame(yr=fitlbspr@Years,  #years 
                    codsp=sp,
                    rawsl50=fitlbspr@SL50, #raw sl50
                    sl50sd=sqrt(fitlbspr@Vars[,1]),  #sl50 sd
                    sl50= fitlbspr@Ests[,1],   #smooth sl50
                    rawsl95=fitlbspr@SL95,
                    sl95sd=sqrt(fitlbspr@Vars[,2]),
                    sl95= fitlbspr@Ests[,2],
                    rawfm=fitlbspr@FM,
                    fmsd=sqrt(fitlbspr@Vars[,3]),
                    fm= fitlbspr@Ests[,3],
                    rawspr=fitlbspr@SPR,
                    sprsd=sqrt(fitlbspr@Vars[,4]),
                    spr= fitlbspr@Ests[,4],
                    source="Pred") #Predicted length frequency


#binding spr results from observed and reconstructed lengths (data frame)
lbspr= rbind(lbsprobs,lbsprpred)
lbspr_blf= lbspr #general object for plotting


#------------------------------------------------------
#length frequency reconstructed vs LBSPR Pop Simulation 
#------------------------------------------------------

szfreqsim=data.frame(yr=NULL,codsp=NULL,lmids=NULL,nsim=NULL,npred=NULL,res=NULL,rmse=NULL,
                     smd=NULL,ks=NULL,t=NULL,x2=NULL) #empty frame for simulated frequencies

#getting subset for current specie of predicted sizes 
szfreqpredsub=szfreqpred[szfreqpred$codsp==sp,]


for (j in 1:length(sp)) {  #loop through all species codes
  
  years=1950:2022
  
  for (i in 1:length(years)) { #loop through each year of available information
    
    
    #------------------------------------------------------
    #biological parameters for LBSPR simulation
    mypars <- new("LB_pars",verbose=F) #blank parameters object for simulation
    mypars@Species <- sp
    mypars@Linf <- pars$linf[pars$codsp==sp]
    mypars@L50 <-  pars$l50[pars$codsp==sp] 
    mypars@L95 <- 1.1*pars$l50[pars$codsp==sp] 
    mypars@MK <-  0.94/pars$k[pars$codsp==sp] #M=0.94 (Freire et al., 2005)
    mypars@L_units <- "cm"
    
    #Get selectivity and mortality parameters (SL50,SL95, FM) from LBSPR estimation
    mypars@SL50 <- lbsprpred$rawsl50[lbsprpred$yr==years[i]] 
    mypars@SL95 <- lbsprpred$rawsl95[lbsprpred$yr==years[i]]
    mypars@FM <-   lbsprpred$rawfm[lbsprpred$yr==years[i]]
    mypars@BinMax<- 1.4*pars$linf[pars$codsp==sp]
    mypars@BinMin<-  0
    mypars@BinWidth <- 5
    
    #LBSPR simulation
    mysim <- LBSPRsim(mypars)
    
    #normalize (min-max) function to make frequencies comparable
    normalize <- function(x, na.rm = TRUE) {
      return((x- min(x)) /(max(x)-min(x)))
    }
    
    
    lsim= data.frame(yr=rep(years[i], length(mysim@pLPop[,1])),
                    codsp=rep(sp[j],length(mysim@pLPop[,1])),
                    lmids=mysim@pLPop[,1],  
                    nsim=normalize(mysim@pLPop[,5]), #surviving in terms of probability
                    npred=NA)
    
    #-----------------------------------------------------
    # transforming predicted frequencies in length classes
    #we need to cut intervals to match with the simulated sizes
    #creating intervals between min and max length bins 
    lpred= data.frame(fl=szfreqpredsub$fl[szfreqpredsub$codsp==sp  
                                           & szfreqpredsub$yr==years[i]])
    freqtable <- data.frame(table(cut(lpred$fl,
                                  seq(min(mysim@pLPop[,1]),max(mysim@pLPop[,1])+5,5))))
    
    intervals=as.character(freqtable$Var1) #character vector
    intervals=gsub("[()]", "", intervals)       #remove parenthesis
    intervals=sapply(strsplit(intervals, split=',', fixed=TRUE), `[`, 1) #get the first class
    freqtable$Var1=intervals
    
    #normalized predicted frequency
    lsim$npred= normalize(freqtable$Freq)
    #------------------------------------------------------
    
    #Residual length predicted and simulated frequencies
    lsim$res= lsim$nsim-lsim$npred
    lsim$res[lsim$npred==0]=NA #avoiding Residual length in classes with no data
    
    #Root Mean Squared Error- RMSE
    rmse= function(x, y, na.rm = TRUE) {
      return(sqrt(sum((x- y)^2) /length(x)))
    }
    
    lsim$rmse= rmse(lsim$npred,lsim$nsim)
    
    #standardized mean reserence (SMD)
    smd= function(x, y, na.rm = TRUE) {
    return( abs(mean(x)- mean(y)) / sqrt((sd(x)^2+sd(y)^2)/2) )
    }
    
    lsim$smd= smd(lsim$npred,lsim$nsim)
    
    
    #Kolmogorov-smirnov test (comparison of cumulative distributions)
    #we need to convert probabilities in sample distributions to perform KS test
    lpredfreq=lpred$fl #predicted length distribution
    
    #simulated frequency distribution 
    lsimfreq= data.frame(lmids=lsim$lmids,n=round(mysim@pLPop[,5]*100))
    lsimfreq= lsimfreq %>%
      tidyr::uncount(n)
    lsimfreq=lsimfreq$lmids
    
    # Kolmogorov smirnov test to compare if the samples come from the same probability distributions
    lsim$ks=suppressWarnings(ks.test(lpredfreq,lsimfreq)$p.value) #kolmogorov-smirnov test
    
    #t-student test to compare if the means of distributions are the same
    lsim$t=suppressWarnings( t.test(lpredfreq,lsimfreq)$p.value)
    
    # chi-square goodness-of-fit test
    lsim$x2= suppressWarnings(chisq.test(sample(lpredfreq, 95,replace = FALSE), 
                    sample(lsimfreq, 95,replace = FALSE))$p.value)
    
    #bind into size frequency simulate data frame
    szfreqsim=rbind(szfreqsim,lsim)
  }
}

#get comparison of frequencies data frame 
szfreqsim_blf= szfreqsim
#--------------------------------------


#-------------
#plot section
#-------------
#plot GAM mean size prediction
p15<-ggplot(data=szmeanpred,aes(x=yr,y=fit))+
  geom_ribbon(data=dplyr::filter(szmeanpred,codsp==sp),
                      aes(ymin = lw, ymax = up,linetype='Pred'), 
               col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(szmeanpred,codsp==sp),
                      aes(x=yr,y=fit,linetype='Pred'),col="gray60",size=2)+
  #geom_line(data=dplyr::filter(gam,codsp==sp),aes(x=yr,y=fit,col=codgr,linetype=region),size=1.2)+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
              aes(x=yr,y=mean),size=2,col="gray40")+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='')+
  scale_y_continuous() +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p15


ggsave(paste(c("Length_model_prediction_Small_Tunas_",sp,".png"),collapse=''), plot = p15, device = "png", units = "cm",
       width = 24, height = 14)

# plot GAM partial effects
p16<-ggplot(dplyr::filter(pareffect,codsp==sp), 
           aes(x=effect,y=fit,fill=var,col="Gear"),size=2) + 
  #geom_line(aes(x=effect,y=up, col="Gear",group=1)) +
  geom_point(aes(x=effect,y=fit,col="Gear"),size=3.5)+
  geom_errorbar(aes(ymin=lw, ymax=up,col="Gear"), width=0.5,
                position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Effect",x="Gear",col='',fill="")+
  scale_y_continuous() +
  scale_x_discrete() +
  scale_colour_manual(values="gray60") +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p16

ggsave(paste(c("Effects_model_Small_Tunas_",sp,".png"),collapse=''), plot = p16, device = "png", units = "cm",
       width = 24, height = 14)


#plot GAM boxplot size freq prediction
p17<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqpred,codsp==sp), 
               aes(x=factor(yr), y=fl),
               col="gray40",
               fill="gray60",
               alpha=0.6,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=factor(yr),y=mean),size=2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Fork length (cm)",col="Region",linetype='',fill='Gear')+
  scale_y_continuous() +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p17


ggsave(paste(c("Length_pred_boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p17, device = "png", units = "cm",
       width = 24, height = 15)


#plot half violin for the predicted frequencies
p18<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(szfreqpred,codsp==sp), 
                             aes(x=factor(yr),y=fl),
                             fill="gray20",
                             col="gray20",
                             alpha=0.7,
                             side = "r")+ 
  #show.legend = TRUE,
  #position = "identity",
  #fill=,
  #nudge = 0,
  #scale = 'area',
  #na.rm = TRUE,
  #lwd=1,
  #trim = TRUE)+
  #facet_wrap(~codsp,scales= "free_y")+
  #draw_quantiles = c(0.5))+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
            aes(x=factor(yr),y=mean),size=2)+
  # #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p18

ggsave(paste(c("Length_pred_density_Small_Tunas_",sp,".png"),collapse = ''), plot = p18, device = "png", units = "cm",
       width = 35, height = 12)



#plot SPR trajectories comparison (observed vs Predicted)
p19<-ggplot(dplyr::filter(lbspr,codsp==sp), 
            aes(x=yr,y=rawspr,fill=source,col=source),size=2) + 
  #geom_line(size=1,aes(fill=source)) +
  geom_point(size=2)+
  geom_smooth(size=1.2,method = 'loess' ,formula = 'y ~ x',se =FALSE)+
  geom_errorbar(aes(ymin=ifelse(rawspr-sprsd < 0, 0, rawspr-sprsd),
                    ymax=ifelse(rawspr+sprsd > 1, 1, rawspr+sprsd)), width=1.5,
                position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Spawning Potential Ratio (SPR)",x="Year",col='',fill="")+
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks = seq(1950, 2023, 30),limits = c(1950, 2023)) +
  scale_colour_manual(values = c("gray60","gray40")) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p19

ggsave(paste(c("SPR_Small_Tunas_",sp,".png"),collapse=''), plot = p19, device = "png", units = "cm",
       width = 24, height = 14)



#plot GAM boxplot size freq prediction
p20<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqsim,codsp==sp), 
               aes(group=as.character(yr),x=as.character(yr), y=res),
               col="gray40",
               fill="gray60",
               alpha=0.4,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_smooth(data=dplyr::filter(szfreqsim,codsp==sp),
              method = "loess",span=0.2,col="gray20", se = FALSE,
              aes(group=1,x=as.character(yr), y=res,linetype="LOESS"),size=1.5)+
  geom_text(data=data.frame(szfreqsim), aes(x="2014", y=0.9*max(res,na.rm = TRUE), 
                              label=paste("RMSE=",round(mean(rmse,na.rm = TRUE),2))))+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Length residual",col="",linetype='',fill='')+
  scale_y_continuous() +
  scale_fill_manual(values = c("blue", "red"))+
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p20

ggsave(paste(c("Length_res_boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p20, device = "png", units = "cm",
       width = 24, height = 15)
#x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x#




#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
# models by species
# BRS model - Length Reconstruction 
#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
sp='BRS'

#-------------------------------------
#find the best distribution for specie
#-------------------------------------
dist=data.frame(yr=NULL,codsp=NULL,dist=NULL,codgr=NULL, #empty frame
                mean=NULL,shape=NULL,cv=NULL,aic=NULL) #Shape for gamma is like the mean

for (i in sp) {
  
  years=unique(szfreqobs$yr[szfreqobs$codsp==i])
  
  for (j in years) {
    
    #frequency for each year  
    sz= szfreqobs$fl[szfreqobs$codsp==i & szfreqobs$yr==j]
    gear= unique(szfreqobs$codgr[szfreqobs$codsp==i & szfreqobs$yr==j])[1]
    
    if (length(sz)>5 & (length(unique(sz))==1)==FALSE) {
      
      #graphical analysis
      #fitdistrplus::descdist(sz, discrete = FALSE,boot = 1000)
      
      #fitting distributions
      fit_1 <- fitdistrplus::fitdist(sz, "norm")
      fit_2<- fitdistrplus::fitdist(sz, "invgauss",start=list(mean=mean(sz),shape= 1))
      fit_3 <- fitdistrplus::fitdist(sz, "gamma")
      
      dat=data.frame(yr=j,codsp=i,
                     dist=c(fit_1$distname,fit_2$distname,fit_3$distname),
                     codgr=gear,
                     mean=c(fit_1$estimate[1],fit_2$estimate[1],fit_3$estimate[1]), 
                     shape=c(fit_1$estimate[2],fit_2$estimate[2],fit_3$estimate[2]),
                     cv=c(fit_1$estimate[2]/fit_1$estimate[1],
                          fit_2$estimate[2]/fit_2$estimate[1],
                          fit_3$estimate[2]/fit_3$estimate[1]),
                     aic=c(fit_1$aic,fit_2$aic,fit_3$aic))
      
      dist=rbind(dist,dat)
    }
  }
}

dist=dist %>%
  dplyr::group_by(yr,codsp,codgr) %>% 
  slice_min(aic, n = 1)  %>%
  dplyr::arrange(dist) 

#Table with the best distribution
table(dist$codsp,dist$dist)  
distribution='norm'

#----------------------------------------------------
#predict 1950 length when there was incipient fishing
#---------------------------------------------------
#Assuming 10% higher than the current mean length level
avg1950=1.1 * mean(szmeanobs$mean[szmeanobs$codsp==sp], na.rm = TRUE)

#put in the observed mean lengths
szmeanobs$mean[szmeanobs$codsp==sp & szmeanobs$yr==1950]= round(avg1950,2)
#----------------------------------

#selecting variables .. subseting
szmeanobssub= szmeanobs[,c("yr",'region','codsp','codgr','mean')]

#removing extreme values
szmeanobssub$mean[szmeanobssub$codsp=='BLF' & szmeanobssub$mean<50] =NA
szmeanobssub$mean[szmeanobssub$codsp=='DOL' & szmeanobssub$mean<75] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean>100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean<60] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean<100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean>180] =NA

#get only non NA values
szmeanobssub=szmeanobssub[is.na(szmeanobssub$mean)==FALSE & szmeanobssub$codsp==sp,]

#grouping factors most likely (Fredou et al., 2021)
szmeanobssub$region[szmeanobssub$region=='ALL']='NE'
szmeanobssub$codgr[szmeanobssub$codgr=='UN']='TR'

#transforming columns in to factors
szmeanobssub$region= factor(szmeanobssub$region)
szmeanobssub$codgr= factor(szmeanobssub$codgr)
szmeanobssub$codsp= factor(szmeanobssub$codsp)


# model GAM selection
mod1<- mgcv::gam(mean ~ s(yr)+region,
                 data= szmeanobssub,
                 method = "REML" ,
                 family =gaussian(link='log'))

mod2<- mgcv::gam(mean ~ s(yr)+codgr, 
                 data= szmeanobssub,
                 method = "REML" ,
                 family =gaussian(link='log'))

mod3<- mgcv::gam(mean ~ s(yr)+region+codgr,
                 data= szmeanobssub,
                 method = "REML" ,
                 family =gaussian(link='log'))
#aic check
mod1$aic
mod2$aic
mod3$aic

summary(mod1)
summary(mod2)
summary(mod3)


#final model
modgam<- mgcv::gam(mean ~ s(yr)+region, 
                   data= szmeanobssub,
                   method = "REML" ,
                   family =gaussian(link='log'))

#assign GAM model to the specie
modgam_brs=modgam

#diagnostic
summary(modgam)

par(mfrow=c(2,2))
gam.check(modgam)# model check

#residual 
#homocedasticity
lmtest::bptest(modgam)

#normality
shapiro.test(resid(modgam))

# Check overall concurvity
round(concurvity(modgam, full = TRUE),2)

plot.gam(modgam, # plot partial effects
         rug=TRUE,       #puts the x variable
         #select = c(1), #select some partial effects
         pages = 1,
         all.terms = TRUE,
         residuals = TRUE,
         #cex=2,
         #lwd=2,
         se=TRUE,
         #shade = TRUE,
         #shade.col = "lightblue", #Shading intervals
         shift = coef(modgam)[1]) #shift the scale by


#----------------------
#predictions "year"
newdataexpand= expand.grid(yr=seq(1950,2022,1),
                           region=szmeanobssub$region,
                           codgr="UN")

pred<- data.frame(predict(modgam, newdata = newdataexpand,
                          type = "response", se.fit=TRUE))

pred$up<-pred$fit+qnorm(0.975)*pred$se.fit #confidence intervals
pred$lw<-pred$fit-qnorm(0.975)*pred$se.fit

gam=cbind(type='GAM',
          codsp=sp,
          family=modgam$family$family,
          link=modgam$family$link,
          formula=paste0(c(paste(modgam$formula)[2],
                           paste(modgam$formula)[1],
                           paste(modgam$formula)[3]),
                         collapse = ' '),
          newdataexpand,
          pred,
          aic=modgam$aic,
          r.sq=summary(modgam)$r.sq,
          dev.ex=summary(modgam)$dev.expl)


#extract mean of all levels
szmeanpred <- gam %>%
  dplyr::group_by(type, codsp, family, link, formula, yr) %>%
  dplyr::summarize(across(c(fit, se.fit, up, lw, aic, r.sq, dev.ex), ~ round(mean(.), 2)))

#creating data frame
szmeanpred_brs=szmeanpred #size mean predicted for all factors combined
#----------------------------------------------

#partial effects 
newdata= szmeanobssub   #desired partial effect

pareffect<- data.frame(predict(modgam, newdata = newdata,
                               type = "terms", se.fit=TRUE)) #terms gets the effects

#intervals for the partial effect
pareffect$up.region<-pareffect$fit.region+qnorm(0.975)*pareffect$se.fit.region #confidence intervals
pareffect$lw.region<-pareffect$fit.region-qnorm(0.975)*pareffect$se.fit.region
pareffect$region= szmeanobssub$region

#extract mean of all levels
pareffect= pareffect[grepl("region", names(pareffect))]
pareffect <- pareffect %>%
  dplyr::group_by(region) %>%
  dplyr::summarize(across(c(fit.region, se.fit.region, up.region, lw.region), ~ round(mean(.), 2))) %>%
  dplyr::rename(.,  "effect"="region", "fit"= "fit.region", "se.fit"="se.fit.region", "up"="up.region", "lw"="lw.region") %>%
  dplyr::mutate(., codsp = sp, var = "Region")

pareffect_brs=pareffect


#----------------------------------
#simulating frequency distributions
#----------------------------------

#estimating FL distributions from reconstructed mean length
flloop=data.frame(yr=NULL,codsp=NULL,fl=NULL,source=NULL) #empty data frame
sp=sp  #list of species for loop

#coefficient of variation for the corresponding distribution 
cv=mean(dist$cv[dist$codsp==sp & dist$dist==distribution]) #the shape is like mean for gamma, so the cv is corrupted

for (j in 1:length(sp)) {  #loop through all species
  #codes
  years=1950:2022
  n=100 #number of simulate individuals
  
  for (i in 1:length(years)) { #loop through each year of available information
    
    l= data.frame(yr=rep(years[i],n),
                  codsp=rep(sp[j],n),
                  fl=rnorm(n,  #distribution 
                           mean=szmeanpred$fit[szmeanpred$yr==years[i] & szmeanpred$codsp==sp], 
                           sd=szmeanpred$fit[szmeanpred$yr==years[i] & szmeanpred$codsp==sp]*cv),
                  source=rep(as.character('Pred'),n))
    
    flloop=rbind(flloop,l)
  }
}


# size frequency predicted data frame
szfreqpred= flloop
szfreqpred_brs=szfreqpred
#-------------------------------------


#-----------------------------------------
#Observed SPR vs Reconstructed SPR (LBSPR)
#-----------------------------------------
mypars <- new("LB_pars",verbose=F) #blank parameters object
mypars@Species <- sp
mypars@Linf <- pars$linf[pars$codsp==sp]
mypars@L50 <-  pars$l50[pars$codsp==sp] 
mypars@L95 <- 1.1*pars$l50[pars$codsp==sp] 
mypars@MK <-  0.36/0.16 #m=0.36, ~k=0.16(Nobrega, 2002)
mypars@L_units <- "cm"


#SPR from observed length frequencies 
#getting size classes and counts
szfreqobssub= szfreqobs[szfreqobs$codsp==sp,]
#szfreqobssub= szfreqobssub[szfreqobssub$yr!= 2017,] #removing year with no data

szloop=data.frame(codsp=NULL,yr=NULL,lmids=NULL,n=NULL)

years=unique(szfreqobssub$yr[szfreqobssub$codsp==sp])

for (i in years) {
  h=hist(szfreqobssub$fl[szfreqobssub$yr==i],
         breaks=seq(from=min(szfreqobssub$fl[szfreqobssub$codsp==sp])-5,
                    to=max(szfreqobssub$fl[szfreqobssub$codsp==sp])+5, by=5), #binwidth
         xlab = 'Fork length (cm)',
         main = paste(j))
  freq=data.frame(codsp=sp,yr=i,lmids=h$mids,n=h$counts)    
  szloop=rbind(szloop,freq)
}

#LBSPR lengths format
szfreqobssub =szloop %>% #removing first row 
  tidyr::spread(., yr, n) %>%  #years in columns
  `rownames<-`(.[,2]) %>%   #years in rownames
  dplyr::select(-codsp,-lmids) 

#observed length frequencies object
mylengthsobs <- new("LB_lengths",LB_pars=mypars,dataType="freq",verbose=F) #blank length object

mylengthsobs@LMids<- as.numeric(rownames(szfreqobssub)) #mid points
mylengthsobs@LData<- as.matrix(szfreqobssub)    #frequency table
mylengthsobs@L_units<- "cm"  #units
mylengthsobs@Years<- as.numeric(colnames(szfreqobssub)) #list of years
mylengthsobs@NYears<- length(as.numeric(colnames(szfreqobssub))) #number of years

#plotSize(mylengthsobs)
fitlbspr= LBSPRfit(mypars, mylengthsobs,verbose = FALSE)
#plotSize(fitlbspr)
#plotEsts(fitlbspr)

lbsprobs= data.frame(yr=fitlbspr@Years,  #years 
                     codsp=sp,
                     rawsl50=fitlbspr@SL50, #raw sl50
                     sl50sd=sqrt(fitlbspr@Vars[,1]),  #sl50 sd
                     sl50= fitlbspr@Ests[,1],   #smooth sl50
                     rawsl95=fitlbspr@SL95,
                     sl95sd=sqrt(fitlbspr@Vars[,2]),
                     sl95= fitlbspr@Ests[,2],
                     rawfm=fitlbspr@FM,
                     fmsd=sqrt(fitlbspr@Vars[,3]),
                     fm= fitlbspr@Ests[,3],
                     rawspr=fitlbspr@SPR,
                     sprsd=sqrt(fitlbspr@Vars[,4]),
                     spr= fitlbspr@Ests[,4],
                     source="Obs") #observed length frequency

#----------------------------------------------------------

#SPR from Reconstructed length frequencies 
#getting size classes and counts
szfreqpredsub=szfreqpred[szfreqpred$codsp==sp,]
#szfreqrecsub= szfreqrecsub[szfreqrecsub$yr!= 2017,] #removing year with no data

szloop=data.frame(codsp=NULL,yr=NULL,lmids=NULL,n=NULL)

years=unique(szfreqpredsub$yr[szfreqpredsub$codsp==sp])

for (i in years) {
  h=hist(szfreqpredsub$fl[szfreqpredsub$yr==i],
         breaks=seq(from=min(szfreqpredsub$fl[szfreqpredsub$codsp==sp])-5,
                    to=max(szfreqpredsub$fl[szfreqpredsub$codsp==sp])+25, by=5), #binwidth
         xlab = 'Fork length (cm)',
         main = paste(j))
  freq=data.frame(codsp=sp,yr=i,lmids=round(h$mids,1),n=h$counts)    
  szloop=rbind(szloop,freq)
}

#LBSPR lengths format
szfreqpredsub =szloop %>% #removing first row 
  tidyr::spread(., yr, n) %>%  #years in columns
  `rownames<-`(.[,2]) %>%   #years in rownames
  dplyr::select(-codsp,-lmids) 

#observed length frequencies object
mylengthspred <- new("LB_lengths",LB_pars=mypars,dataType="freq",verbose=F) #blank length object
mylengthspred@LMids<- as.numeric(rownames(szfreqpredsub)) #mid points
mylengthspred@LData<- as.matrix(szfreqpredsub)    #frequency table
mylengthspred@L_units<- "cm"  #units
mylengthspred@Years<- as.numeric(colnames(szfreqpredsub)) #list of years
mylengthspred@NYears<- length(as.numeric(colnames(szfreqpredsub))) #number of years

#plotSize(mylengthspred)

fitlbspr= LBSPRfit(mypars, mylengthspred,verbose = FALSE)
#plotSize(fitlbspr)
#plotEsts(fitlbspr)


lbsprpred= data.frame(yr=fitlbspr@Years,  #years 
                      codsp=sp,
                      rawsl50=fitlbspr@SL50, #raw sl50
                      sl50sd=sqrt(fitlbspr@Vars[,1]),  #sl50 sd
                      sl50= fitlbspr@Ests[,1],   #smooth sl50
                      rawsl95=fitlbspr@SL95,
                      sl95sd=sqrt(fitlbspr@Vars[,2]),
                      sl95= fitlbspr@Ests[,2],
                      rawfm=fitlbspr@FM,
                      fmsd=sqrt(fitlbspr@Vars[,3]),
                      fm= fitlbspr@Ests[,3],
                      rawspr=fitlbspr@SPR,
                      sprsd=sqrt(fitlbspr@Vars[,4]),
                      spr= fitlbspr@Ests[,4],
                      source="Pred") #Predicted length frequency


#binding spr results from observed and reconstructed lengths (data frame)
lbspr= rbind(lbsprobs,lbsprpred)
lbspr_brs= lbspr #general object for plotting


#------------------------------------------------------
#length frequency reconstructed vs LBSPR Pop Simulation 
#------------------------------------------------------

szfreqsim=data.frame(yr=NULL,codsp=NULL,lmids=NULL,nsim=NULL,npred=NULL,res=NULL,rmse=NULL,
                     smd=NULL,ks=NULL,t=NULL,x2=NULL) #empty frame for simulated frequencies

#getting subset for current specie of predicted sizes 
szfreqpredsub=szfreqpred[szfreqpred$codsp==sp,]


for (j in 1:length(sp)) {  #loop through all species codes
  
  years=1950:2022
  
  for (i in 1:length(years)) { #loop through each year of available information
    
    
    #------------------------------------------------------
    #biological parameters for LBSPR simulation
    mypars <- new("LB_pars",verbose=F) #blank parameters object for simulation
    mypars@Species <- sp
    mypars@Linf <- pars$linf[pars$codsp==sp]
    mypars@L50 <-  pars$l50[pars$codsp==sp] 
    mypars@L95 <- 1.1*pars$l50[pars$codsp==sp] 
    mypars@MK <-  0.36/0.16 #m=0.36, ~k=0.16(Nobrega, 2002)
    mypars@L_units <- "cm"
    
    #Get selectivity and mortality parameters (SL50,SL95, FM) from LBSPR estimation
    mypars@SL50 <- lbsprpred$rawsl50[lbsprpred$yr==years[i]] 
    mypars@SL95 <- lbsprpred$rawsl95[lbsprpred$yr==years[i]]
    mypars@FM <-   lbsprpred$rawfm[lbsprpred$yr==years[i]]
    mypars@BinMax<- 1.4*pars$linf[pars$codsp==sp]
    mypars@BinMin<-  0
    mypars@BinWidth <- 5
    
    #LBSPR simulation
    mysim <- LBSPRsim(mypars)
    
    #normalize (min-max) function to make frequencies comparable
    normalize <- function(x, na.rm = TRUE) {
      return((x- min(x)) /(max(x)-min(x)))
    }
    
    
    lsim= data.frame(yr=rep(years[i], length(mysim@pLPop[,1])),
                     codsp=rep(sp[j],length(mysim@pLPop[,1])),
                     lmids=mysim@pLPop[,1],  
                     nsim=normalize(mysim@pLPop[,5]), #surviving in terms of probability
                     npred=NA)
    
    #-----------------------------------------------------
    # transforming predicted frequencies in length classes
    #we need to cut intervals to match with the simulated sizes
    #creating intervals between min and max length bins 
    lpred= data.frame(fl=szfreqpredsub$fl[szfreqpredsub$codsp==sp  
                                          & szfreqpredsub$yr==years[i]])
    freqtable <- data.frame(table(cut(lpred$fl,
                                      seq(min(mysim@pLPop[,1]),max(mysim@pLPop[,1])+5,5))))
    
    intervals=as.character(freqtable$Var1) #character vector
    intervals=gsub("[()]", "", intervals)       #remove parenthesis
    intervals=sapply(strsplit(intervals, split=',', fixed=TRUE), `[`, 1) #get the first class
    freqtable$Var1=intervals
    
    #normalized predicted frequency
    lsim$npred= normalize(freqtable$Freq)
    #------------------------------------------------------
    
    #Residual length predicted and simulated frequencies
    lsim$res= lsim$nsim-lsim$npred
    lsim$res[lsim$npred==0]=NA #avoiding Residual length in classes with no data
    
    #Root Mean Squared Error- RMSE
    rmse= function(x, y, na.rm = TRUE) {
      return(sqrt(sum((x- y)^2) /length(x)))
    }
    
    lsim$rmse= rmse(lsim$npred,lsim$nsim)
    
    #standardized mean reserence (SMD)
    smd= function(x, y, na.rm = TRUE) {
      return( abs(mean(x)- mean(y)) / sqrt((sd(x)^2+sd(y)^2)/2) )
    }
    
    lsim$smd= smd(lsim$npred,lsim$nsim)
    
    
    #Kolmogorov-smirnov test (comparison of cumulative distributions)
    #we need to convert probabilities in sample distributions to perform KS test
    lpredfreq=lpred$fl #predicted length distribution
    
    #simulated frequency distribution 
    lsimfreq= data.frame(lmids=mysim@pLPop[,1],n=round(mysim@pLPop[,5]*100))
    lsimfreq= lsimfreq %>%
      tidyr::uncount(n)
    lsimfreq=lsimfreq$lmids
    
    # Kolmogorov smirnov test to compare if the samples come from the same probability distributions
    lsim$ks=suppressWarnings(ks.test(lpredfreq,lsimfreq)$p.value) #kolmogorov-smirnov test
    
    #t-student test to compare if the means of distributions are the same
    lsim$t= suppressWarnings(t.test(lpredfreq,lsimfreq)$p.value)
    
    # chi-square goodness-of-fit test
    lsim$x2= suppressWarnings(chisq.test(sample(lpredfreq, 95,replace = FALSE), sample(lsimfreq, 95,replace = FALSE))$p.value)
    
    #bind into size frequency simulate data frame
    szfreqsim=rbind(szfreqsim,lsim)
  }
}


#get comparison of frequencies data frame 
szfreqsim_brs= szfreqsim
#--------------------------------------


#-------------
#plot section
#-------------
#plot GAM mean size prediction
p21<-ggplot(data=szmeanpred,aes(x=yr,y=fit))+
  geom_ribbon(data=dplyr::filter(szmeanpred,codsp==sp),
              aes(ymin = lw, ymax = up,linetype='Pred'), 
              col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(szmeanpred,codsp==sp),
            aes(x=yr,y=fit,linetype='Pred'),col="gray60",size=2)+
  #geom_line(data=dplyr::filter(gam,codsp==sp),aes(x=yr,y=fit,col=codgr,linetype=region),size=1.2)+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=yr,y=mean),size=2,col="gray40")+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='')+
  scale_y_continuous() +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p21


ggsave(paste(c("Length_model_prediction_Small_Tunas_",sp,".png"),collapse=''), plot = p21, device = "png", units = "cm",
       width = 24, height = 14)

# plot GAM partial effects
p22<-ggplot(dplyr::filter(pareffect,codsp==sp), 
            aes(x=effect,y=fit,fill=var,col="Region"),size=2) + 
  #geom_line(aes(x=effect,y=up, col="Gear",group=1)) +
  geom_point(aes(x=effect,y=fit,col="Region"),size=3.5)+
  geom_errorbar(aes(ymin=lw, ymax=up,col="Region"), width=0.5,
                position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Effect",x="Region",col='',fill="")+
  scale_y_continuous() +
  scale_x_discrete() +
  scale_colour_manual(values="gray60") +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p22

ggsave(paste(c("Effects_model_Small_Tunas_",sp,".png"),collapse=''), plot = p22, device = "png", units = "cm",
       width = 24, height = 14)


#plot GAM boxplot size freq prediction
p23<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqpred,codsp==sp), 
               aes(x=factor(yr), y=fl),
               col="gray40",
               fill="gray60",
               alpha=0.6,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=factor(yr),y=mean),size=2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Fork length (cm)",col="Region",linetype='',fill='Gear')+
  scale_y_continuous() +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p23

ggsave(paste(c("Length_pred_boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p23, device = "png", units = "cm",
       width = 24, height = 15)


#plot half violin for the predicted frequencies
p24<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(szfreqpred,codsp==sp), 
                             aes(x=factor(yr),y=fl),
                             fill="gray20",
                             col="gray20",
                             alpha=0.7,
                             side = "r")+ 
  #show.legend = TRUE,
  #position = "identity",
  #fill=,
  #nudge = 0,
  #scale = 'area',
  #na.rm = TRUE,
  #lwd=1,
  #trim = TRUE)+
  #facet_wrap(~codsp,scales= "free_y")+
  #draw_quantiles = c(0.5))+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=factor(yr),y=mean),size=2)+
  # #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p24

ggsave(paste(c("Length_pred_density_Small_Tunas_",sp,".png"),collapse = ''),plot = p24, device = "png", units = "cm",
       width = 35, height = 12)


#plot SPR trajectories comparison (observed vs Predicted)
p25<-ggplot(dplyr::filter(lbspr,codsp==sp), 
            aes(x=yr,y=rawspr,fill=source,col=source),size=2) + 
  #geom_line(size=1,aes(fill=source)) +
  geom_point(size=2)+
  geom_smooth(size=1.2,method = 'loess' ,formula = 'y ~ x',se =FALSE,span = 0.4)+
  geom_errorbar(aes(ymin=ifelse(rawspr-sprsd < 0, 0, rawspr-sprsd),
                    ymax=ifelse(rawspr+sprsd > 1, 1, rawspr+sprsd)), width=1.5,
                position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Spawning Potential Ratio (SPR)",x="Year",col='',fill="")+
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks = seq(1950, 2023, 30),limits = c(1950, 2023)) +
  scale_colour_manual(values = c("gray60","gray40")) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p25

ggsave(paste(c("SPR_Small_Tunas_",sp,".png"),collapse=''), plot = p25, device = "png", units = "cm",
       width = 24, height = 14)



#plot GAM boxplot size freq prediction
p26<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqsim,codsp==sp), 
               aes(group=as.character(yr),x=as.character(yr), y=res),
               col="gray40",
               fill="gray60",
               alpha=0.4,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_smooth(data=dplyr::filter(szfreqsim,codsp==sp),
              method = "loess",span=0.1,col="gray20", se = FALSE,
              aes(group=1,x=as.character(yr), y=res,linetype="LOESS"),size=1.5)+
  geom_text(data=dplyr::filter(data.frame(szfreqsim),codsp==sp), 
                               aes(x="2014", y=0.9*max(res,na.rm = TRUE), 
                                    label=paste("RMSE=",round(mean(rmse,na.rm = TRUE),2))))+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Length residual",col="",linetype='',fill='')+
  scale_y_continuous() +
  scale_fill_manual(values = c("blue", "red"))+
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p26

ggsave(paste(c("Length_res_boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p26, device = "png", units = "cm",
       width = 24, height = 15)
#x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x#



#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
# models by species
# DOL model - Length Reconstruction 
#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
sp='DOL'

#-------------------------------------
#find the best distribution for specie
#-------------------------------------
dist=data.frame(yr=NULL,codsp=NULL,dist=NULL,codgr=NULL, #empty frame
                mean=NULL,shape=NULL,cv=NULL,aic=NULL) #Shape for gamma is like the mean

for (i in sp) {
  
  years=unique(szfreqobs$yr[szfreqobs$codsp==i])
  
  for (j in years) {
    
    #frequency for each year  
    sz= szfreqobs$fl[szfreqobs$codsp==i & szfreqobs$yr==j]
    gear= unique(szfreqobs$codgr[szfreqobs$codsp==i & szfreqobs$yr==j])[1]
    
    if (length(sz)>5 & (length(unique(sz))==1)==FALSE) {
      
      #graphical analysis
      #fitdistrplus::descdist(sz, discrete = FALSE,boot = 1000)
      
      #fitting distributions
      fit_1 <- fitdistrplus::fitdist(sz, "norm")
      fit_2<- fitdistrplus::fitdist(sz, "invgauss",start=list(mean=mean(sz),shape= 1))
      fit_3 <- fitdistrplus::fitdist(sz, "gamma")
      
      dat=data.frame(yr=j,codsp=i,
                     dist=c(fit_1$distname,fit_2$distname,fit_3$distname),
                     codgr=gear,
                     mean=c(fit_1$estimate[1],fit_2$estimate[1],fit_3$estimate[1]), 
                     shape=c(fit_1$estimate[2],fit_2$estimate[2],fit_3$estimate[2]),
                     cv=c(fit_1$estimate[2]/fit_1$estimate[1],
                          fit_2$estimate[2]/fit_2$estimate[1],
                          fit_3$estimate[2]/fit_3$estimate[1]),
                     aic=c(fit_1$aic,fit_2$aic,fit_3$aic))
      
      dist=rbind(dist,dat)
    }
  }
}

dist=dist %>%
  dplyr::group_by(yr,codsp,codgr) %>% 
  slice_min(aic, n = 1)  %>%
  dplyr::arrange(dist) 

#Table with the best distribution
table(dist$codsp,dist$dist)  
distribution='invgauss'

#----------------------------------------------------
#predict 1950 length when there was incipient fishing
#---------------------------------------------------
#Assuming 10% higher than the current mean length level
avg1950=1.1 * mean(szmeanobs$mean[szmeanobs$codsp==sp], na.rm = TRUE)

#put in the observed mean lengths
szmeanobs$mean[szmeanobs$codsp==sp & szmeanobs$yr==1950]= round(avg1950,2)
#----------------------------------

#selecting variables .. subseting
szmeanobssub= szmeanobs[,c("yr",'region','codsp','codgr','mean')]

#removing extreme values
szmeanobssub$mean[szmeanobssub$codsp=='BLF' & szmeanobssub$mean<50] =NA
szmeanobssub$mean[szmeanobssub$codsp=='DOL' & szmeanobssub$mean<75] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean>100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean<60] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean<100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean>180] =NA

#get only non NA values
szmeanobssub=szmeanobssub[is.na(szmeanobssub$mean)==FALSE & szmeanobssub$codsp==sp,]

#transforming columns in to factors
szmeanobssub$region= factor(szmeanobssub$region)
szmeanobssub$codgr= factor(szmeanobssub$codgr)
szmeanobssub$codsp= factor(szmeanobssub$codsp)


# model GAM selection
mod1<- mgcv::gam(mean ~ s(yr)+region,
                 data= szmeanobssub,
                 method = "REML" ,
                 family =inverse.gaussian(link='log'))

mod2<- mgcv::gam(mean ~ s(yr)+codgr, 
                 data= szmeanobssub,
                 method = "REML" ,
                 family =inverse.gaussian(link='log'))

mod3<- mgcv::gam(mean ~ s(yr)+region+codgr,
                 data= szmeanobssub,
                 method = "REML" ,
                 family =inverse.gaussian(link='log'))
#aic check
mod1$aic
mod2$aic
mod3$aic

summary(mod1)
summary(mod2)
summary(mod3)


#final model
modgam<- mgcv::gam(mean ~ s(yr)+codgr, 
                   data= szmeanobssub,
                   method = "REML" ,
                   family =inverse.gaussian(link='log'))

#assign GAM model to the specie
modgam_dol=modgam

#diagnostic
summary(modgam)

par(mfrow=c(2,2))
gam.check(modgam)# model check

#residual 
#homocedasticity
lmtest::bptest(modgam)

#normality
shapiro.test(resid(modgam))

# Check overall concurvity
round(concurvity(modgam, full = TRUE),2)

plot.gam(modgam, # plot partial effects
         rug=TRUE,       #puts the x variable
         #select = c(1), #select some partial effects
         pages = 1,
         all.terms = TRUE,
         residuals = TRUE,
         #cex=2,
         #lwd=2,
         se=TRUE,
         #shade = TRUE,
         #shade.col = "lightblue", #Shading intervals
         shift = coef(modgam)[1]) #shift the scale by


#----------------------
#predictions "year"
newdataexpand= expand.grid(yr=seq(1950,2022,1),
                           region="ALL",
                           codgr=szmeanobssub$codgr)

pred<- data.frame(predict(modgam, newdata = newdataexpand,
                          type = "response", se.fit=TRUE))

pred$up<-pred$fit+qnorm(0.975)*pred$se.fit #confidence intervals
pred$lw<-pred$fit-qnorm(0.975)*pred$se.fit

gam=cbind(type='GAM',
          codsp=sp,
          family=modgam$family$family,
          link=modgam$family$link,
          formula=paste0(c(paste(modgam$formula)[2],
                           paste(modgam$formula)[1],
                           paste(modgam$formula)[3]),
                         collapse = ' '),
          newdataexpand,
          pred,
          aic=modgam$aic,
          r.sq=summary(modgam)$r.sq,
          dev.ex=summary(modgam)$dev.expl)


#extract mean of all levels
szmeanpred <- gam %>%
  dplyr::group_by(type, codsp, family, link, formula, yr) %>%
  dplyr::summarize(across(c(fit, se.fit, up, lw, aic, r.sq, dev.ex), ~ round(mean(.), 2)))

#creating data frame
szmeanpred_dol=szmeanpred #size mean predicted for all factors combined
#----------------------------------------------

#partial effects 
newdata= szmeanobssub   #desired partial effect

pareffect<- data.frame(predict(modgam, newdata = newdata,
                               type = "terms", se.fit=TRUE)) #terms gets the effects

#intervals for the partial effect
pareffect$up.codgr<-pareffect$fit.codgr+qnorm(0.975)*pareffect$se.fit.codgr #confidence intervals
pareffect$lw.codgr<-pareffect$fit.codgr-qnorm(0.975)*pareffect$se.fit.codgr
pareffect$codgr= szmeanobssub$codgr

#extract mean of all levels
pareffect= pareffect[grepl("codgr", names(pareffect))]
pareffect <- pareffect %>%
  dplyr::group_by(codgr) %>%
  dplyr::summarize(across(c(fit.codgr, se.fit.codgr, up.codgr, lw.codgr), ~ round(mean(.), 2))) %>%
  dplyr::rename(., "effect" = codgr, "fit" = "fit.codgr", "se.fit" = "se.fit.codgr", "up" = "up.codgr", "lw" = "lw.codgr") %>%
  dplyr::mutate(., codsp = sp, var = "Gear")

pareffect_dol=pareffect


#----------------------------------
#simulating frequency distributions
#----------------------------------

#estimating FL distributions from reconstructed mean length
flloop=data.frame(yr=NULL,codsp=NULL,fl=NULL,source=NULL) #empty data frame
sp=sp  #list of species for loop

#coefficient of variation for the corresponding distribution 
cv=mean(dist$cv[dist$codsp==sp & dist$dist==distribution]) #the shape is like mean for gamma, so the cv is corrupted

for (j in 1:length(sp)) {  #loop through all species
  #codes
  years=1950:2022
  n=100 #number of simulate individuals
  
  for (i in 1:length(years)) { #loop through each year of available information
    
    l= data.frame(yr=rep(years[i],n),
                  codsp=rep(sp[j],n),
                  fl=rinvgauss(n,  #distribution 
                           mean=szmeanpred$fit[szmeanpred$yr==years[i] & szmeanpred$codsp==sp], 
                           shape=szmeanpred$fit[szmeanpred$yr==years[i] & szmeanpred$codsp==sp]*cv),
                  source=rep(as.character('Pred'),n))
    
    flloop=rbind(flloop,l)
  }
}


# size frequency predicted data frame
szfreqpred= flloop
szfreqpred_dol=szfreqpred
#-------------------------------------


#-----------------------------------------
#Observed SPR vs Reconstructed SPR (LBSPR)
#-----------------------------------------
mypars <- new("LB_pars",verbose=F) #blank parameters object
mypars@Species <- sp
mypars@Linf <- pars$linf[pars$codsp==sp]
mypars@L50 <-  pars$l50[pars$codsp==sp] 
mypars@L95 <- 1.1*pars$l50[pars$codsp==sp] 
mypars@MK <- pars$m[pars$codsp==sp]/pars$k[pars$codsp==sp] #m=0.36, ~k=0.16(Nobrega, 2002)
mypars@L_units <- "cm"


#SPR from observed length frequencies 
#getting size classes and counts
szfreqobssub= szfreqobs[szfreqobs$codsp==sp,]
# szfreqobssub= szfreqobssub[szfreqobssub$yr!= 2020,] #removing year with no data


szloop=data.frame(codsp=NULL,yr=NULL,lmids=NULL,n=NULL)

years=unique(szfreqobssub$yr[szfreqobssub$codsp==sp])

for (i in years) {
  h=hist(szfreqobssub$fl[szfreqobssub$yr==i],
         breaks=seq(from=min(szfreqobssub$fl[szfreqobssub$codsp==sp])-5,
                    to=max(szfreqobssub$fl[szfreqobssub$codsp==sp])+5, by=5), #binwidth
         xlab = 'Fork length (cm)',
         main = paste(j))
  freq=data.frame(codsp=sp,yr=i,lmids=h$mids,n=h$counts)    
  szloop=rbind(szloop,freq)
}

#LBSPR lengths format
szfreqobssub =szloop %>% #removing first row 
  tidyr::spread(., yr, n) %>%  #years in columns
  `rownames<-`(.[,2]) %>%   #years in rownames
  dplyr::select(-codsp,-lmids) 

#observed length frequencies object
mylengthsobs <- new("LB_lengths",LB_pars=mypars,dataType="freq",verbose=F) #blank length object

mylengthsobs@LMids<- as.numeric(rownames(szfreqobssub)) #mid points
mylengthsobs@LData<- as.matrix(szfreqobssub)    #frequency table
mylengthsobs@L_units<- "cm"  #units
mylengthsobs@Years<- as.numeric(colnames(szfreqobssub)) #list of years
mylengthsobs@NYears<- length(as.numeric(colnames(szfreqobssub))) #number of years

#plotSize(mylengthsobs)
fitlbspr= LBSPRfit(mypars, mylengthsobs,verbose = FALSE)
#plotSize(fitlbspr)
#plotEsts(fitlbspr)

lbsprobs= data.frame(yr=fitlbspr@Years,  #years 
                     codsp=sp,
                     rawsl50=fitlbspr@SL50, #raw sl50
                     sl50sd=sqrt(fitlbspr@Vars[,1]),  #sl50 sd
                     sl50= fitlbspr@Ests[,1],   #smooth sl50
                     rawsl95=fitlbspr@SL95,
                     sl95sd=sqrt(fitlbspr@Vars[,2]),
                     sl95= fitlbspr@Ests[,2],
                     rawfm=fitlbspr@FM,
                     fmsd=sqrt(fitlbspr@Vars[,3]),
                     fm= fitlbspr@Ests[,3],
                     rawspr=fitlbspr@SPR,
                     sprsd=sqrt(fitlbspr@Vars[,4]),
                     spr= fitlbspr@Ests[,4],
                     source="Obs") #observed length frequency

#----------------------------------------------------------

#SPR from Reconstructed length frequencies 
#getting size classes and counts
szfreqpredsub=szfreqpred[szfreqpred$codsp==sp,]
#szfreqrecsub= szfreqrecsub[szfreqrecsub$yr!= 2017,] #removing year with no data

szloop=data.frame(codsp=NULL,yr=NULL,lmids=NULL,n=NULL)

years=unique(szfreqpredsub$yr[szfreqpredsub$codsp==sp])

for (i in years) {
  h=hist(szfreqpredsub$fl[szfreqpredsub$yr==i],
         breaks=seq(from=min(szfreqpredsub$fl[szfreqpredsub$codsp==sp])-5,
                    to=max(szfreqpredsub$fl[szfreqpredsub$codsp==sp])+5, by=5), #binwidth
         xlab = 'Fork length (cm)',
         main = paste(j))
  freq=data.frame(codsp=sp,yr=i,lmids=round(h$mids,1),n=h$counts)    
  szloop=rbind(szloop,freq)
}

#LBSPR lengths format
szfreqpredsub =szloop %>% #removing first row 
  tidyr::spread(., yr, n) %>%  #years in columns
  `rownames<-`(.[,2]) %>%   #years in rownames
  dplyr::select(-codsp,-lmids) 

#observed length frequencies object
mylengthspred <- new("LB_lengths",LB_pars=mypars,dataType="freq",verbose=F) #blank length object
mylengthspred@LMids<- as.numeric(rownames(szfreqpredsub)) #mid points
mylengthspred@LData<- as.matrix(szfreqpredsub)    #frequency table
mylengthspred@L_units<- "cm"  #units
mylengthspred@Years<- as.numeric(colnames(szfreqpredsub)) #list of years
mylengthspred@NYears<- length(as.numeric(colnames(szfreqpredsub))) #number of years

#plotSize(mylengthspred)

fitlbspr= LBSPRfit(mypars, mylengthspred,verbose = FALSE)
#plotSize(fitlbspr)
#plotEsts(fitlbspr)


lbsprpred= data.frame(yr=fitlbspr@Years,  #years 
                      codsp=sp,
                      rawsl50=fitlbspr@SL50, #raw sl50
                      sl50sd=sqrt(fitlbspr@Vars[,1]),  #sl50 sd
                      sl50= fitlbspr@Ests[,1],   #smooth sl50
                      rawsl95=fitlbspr@SL95,
                      sl95sd=sqrt(fitlbspr@Vars[,2]),
                      sl95= fitlbspr@Ests[,2],
                      rawfm=fitlbspr@FM,
                      fmsd=sqrt(fitlbspr@Vars[,3]),
                      fm= fitlbspr@Ests[,3],
                      rawspr=fitlbspr@SPR,
                      sprsd=sqrt(fitlbspr@Vars[,4]),
                      spr= fitlbspr@Ests[,4],
                      source="Pred") #Predicted length frequency


#binding spr results from observed and reconstructed lengths (data frame)
lbspr= rbind(lbsprobs,lbsprpred)
lbspr_dol= lbspr #general object for plotting


#------------------------------------------------------
#length frequency reconstructed vs LBSPR Pop Simulation 
#------------------------------------------------------
szfreqsim=data.frame(yr=NULL,codsp=NULL,lmids=NULL,nsim=NULL,npred=NULL,res=NULL,rmse=NULL,
                     smd=NULL,ks=NULL,t=NULL,x2=NULL) #empty frame for simulated frequencies

#getting subset for current specie of predicted sizes 
szfreqpredsub=szfreqpred[szfreqpred$codsp==sp,]


for (j in 1:length(sp)) {  #loop through all species codes
  
  years=1950:2022
  
  for (i in 1:length(years)) { #loop through each year of available information
    
    
    #------------------------------------------------------
    #biological parameters for LBSPR simulation
    mypars <- new("LB_pars",verbose=F) #blank parameters object for simulation
    mypars@Species <- sp
    mypars@Linf <- pars$linf[pars$codsp==sp]
    mypars@L50 <-  pars$l50[pars$codsp==sp] 
    mypars@L95 <- 1.1*pars$l50[pars$codsp==sp] 
    mypars@MK <- pars$m[pars$codsp==sp]/pars$k[pars$codsp==sp] #
    mypars@L_units <- "cm"
    
    #Get selectivity and mortality parameters (SL50,SL95, FM) from LBSPR estimation
    mypars@SL50 <- lbsprpred$rawsl50[lbsprpred$yr==years[i]] 
    mypars@SL95 <- lbsprpred$rawsl95[lbsprpred$yr==years[i]]
    mypars@FM <-   lbsprpred$rawfm[lbsprpred$yr==years[i]]
    mypars@BinMax<- 1.4*pars$linf[pars$codsp==sp]
    mypars@BinMin<-  0
    mypars@BinWidth <- 5
    
    #LBSPR simulation
    mysim <- LBSPRsim(mypars)
    
    #normalize (min-max) function to make frequencies comparable
    normalize <- function(x, na.rm = TRUE) {
      return((x- min(x)) /(max(x)-min(x)))
    }
    
    
    lsim= data.frame(yr=rep(years[i], length(mysim@pLPop[,1])),
                     codsp=rep(sp[j],length(mysim@pLPop[,1])),
                     lmids=mysim@pLPop[,1],  
                     nsim=normalize(mysim@pLPop[,5]), #surviving in terms of probability
                     npred=NA)
    
    #-----------------------------------------------------
    # transforming predicted frequencies in length classes
    #we need to cut intervals to match with the simulated sizes
    #creating intervals between min and max length bins 
    lpred= data.frame(fl=szfreqpredsub$fl[szfreqpredsub$codsp==sp  
                                          & szfreqpredsub$yr==years[i]])
    freqtable <- data.frame(table(cut(lpred$fl,
                                      seq(min(mysim@pLPop[,1]),max(mysim@pLPop[,1])+5,5))))
    
    intervals=as.character(freqtable$Var1) #character vector
    intervals=gsub("[()]", "", intervals)       #remove parenthesis
    intervals=sapply(strsplit(intervals, split=',', fixed=TRUE), `[`, 1) #get the first class
    freqtable$Var1=intervals
    
    #normalized predicted frequency
    lsim$npred= normalize(freqtable$Freq)
    #------------------------------------------------------
    
    #Residual length predicted and simulated frequencies
    lsim$res= lsim$nsim-lsim$npred
    lsim$res[lsim$npred==0]=NA #avoiding Residual length in classes with no data
    
    #Root Mean Squared Error- RMSE
    rmse= function(x, y, na.rm = TRUE) {
      return(sqrt(sum((x- y)^2) /length(x)))
    }
    
    lsim$rmse= rmse(lsim$npred,lsim$nsim)
    
    #standardized mean reserence (SMD)
    smd= function(x, y, na.rm = TRUE) {
      return( abs(mean(x)- mean(y)) / sqrt((sd(x)^2+sd(y)^2)/2) )
    }
    
    lsim$smd= smd(lsim$npred,lsim$nsim)
    
    
    #Kolmogorov-smirnov test (comparison of cumulative distributions)
    #we need to convert probabilities in sample distributions to perform KS test
    lpredfreq=lpred$fl #predicted length distribution
    
    #simulated frequency distribution 
    lsimfreq= data.frame(lmids=mysim@pLPop[,1],n=round(mysim@pLPop[,5]*100))
    lsimfreq= lsimfreq %>%
      tidyr::uncount(n)
    lsimfreq=lsimfreq$lmids
    
    # Kolmogorov smirnov test to compare if the samples come from the same probability distributions
    lsim$ks=suppressWarnings(ks.test(lpredfreq,lsimfreq)$p.value) #kolmogorov-smirnov test
    
    #t-student test to compare if the means of distributions are the same
    lsim$t=suppressWarnings(t.test(lpredfreq,lsimfreq)$p.value)
    
    # chi-square goodness-of-fit test
    lsim$x2=suppressWarnings(chisq.test(sample(lpredfreq, 95,replace = FALSE), sample(lsimfreq, 95,replace = FALSE))$p.value)
    
    #bind into size frequency simulate data frame
    szfreqsim=rbind(szfreqsim,lsim)
  }
}

#get comparison of frequencies data frame 
szfreqsim_dol= szfreqsim
#--------------------------------------


#-------------
#plot section
#-------------
#plot GAM mean size prediction
p27<-ggplot(data=szmeanpred,aes(x=yr,y=fit))+
  geom_ribbon(data=dplyr::filter(szmeanpred,codsp==sp),
              aes(ymin = lw, ymax = up,linetype='Pred'), 
              col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(szmeanpred,codsp==sp),
            aes(x=yr,y=fit,linetype='Pred'),col="gray60",size=2)+
  #geom_line(data=dplyr::filter(gam,codsp==sp),aes(x=yr,y=fit,col=codgr,linetype=region),size=1.2)+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=yr,y=mean),size=2,col="gray40")+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='')+
  scale_y_continuous() +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p27


ggsave(paste(c("Length_model_prediction_Small_Tunas_",sp,".png"),collapse=''), plot = p27, device = "png", units = "cm",
       width = 24, height = 14)

# plot GAM partial effects
p28<-ggplot(dplyr::filter(pareffect,codsp==sp), 
            aes(x=effect,y=fit,fill=var,col="Gear"),size=2) + 
  #geom_line(aes(x=effect,y=up, col="Gear",group=1)) +
  geom_point(aes(x=effect,y=fit,col="Gear"),size=3.5)+
  geom_errorbar(aes(ymin=lw, ymax=up,col="Gear"), width=0.5,
                position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Effect",x="Gear",col='',fill="")+
  scale_y_continuous() +
  scale_x_discrete() +
  scale_colour_manual(values="gray60") +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p28

ggsave(paste(c("Effects_model_Small_Tunas_",sp,".png"),collapse=''), plot = p28, device = "png", units = "cm",
       width = 24, height = 14)


#plot GAM boxplot size freq prediction
p29<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqpred,codsp==sp), 
               aes(x=factor(yr), y=fl),
               col="gray40",
               fill="gray60",
               alpha=0.6,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=factor(yr),y=mean),size=2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Fork length (cm)",col="Region",linetype='',fill='Gear')+
  scale_y_continuous() +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p29

ggsave(paste(c("Length_pred_boxplot_Small_Tunas_",sp,".png"),collapse = ''),plot=p29, device = "png", units = "cm",
       width = 24, height = 15)


#plot half violin for the predicted frequencies
p30<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(szfreqpred,codsp==sp), 
                             aes(x=factor(yr),y=fl),
                             fill="gray20",
                             col="gray20",
                             alpha=0.7,
                             side = "r")+ 
  #show.legend = TRUE,
  #position = "identity",
  #fill=,
  #nudge = 0,
  #scale = 'area',
  #na.rm = TRUE,
  #lwd=1,
  #trim = TRUE)+
  #facet_wrap(~codsp,scales= "free_y")+
  #draw_quantiles = c(0.5))+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=factor(yr),y=mean),size=2)+
  # #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p30

ggsave(paste(c("Length_pred_density_Small_Tunas_",sp,".png"),collapse = ''),plot=p30, device = "png", units = "cm",
       width = 35, height = 12)


#plot SPR trajectories comparison (observed vs Predicted)
p31<-ggplot(dplyr::filter(lbspr,codsp==sp), 
            aes(x=yr,y=rawspr,fill=source,col=source),size=2) + 
  #geom_line(size=1,aes(fill=source)) +
  geom_point(size=2)+
  geom_smooth(size=1.2,method = 'loess' ,formula = 'y ~ x',se =FALSE,span = 0.4)+
  geom_errorbar(aes(ymin=ifelse(rawspr-sprsd < 0, 0, rawspr-sprsd),
                    ymax=ifelse(rawspr+sprsd > 1, 1, rawspr+sprsd)), width=1.5,
                position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Spawning Potential Ratio (SPR)",x="Year",col='',fill="")+
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks = seq(1950, 2023, 30),limits = c(1950, 2023)) +
  scale_colour_manual(values = c("gray60","gray40")) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p31

ggsave(paste(c("SPR_Small_Tunas_",sp,".png"),collapse=''),plot=p31, device = "png", units = "cm",
       width = 24, height = 14)



#plot GAM boxplot size freq prediction
p32<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqsim,codsp==sp), 
               aes(group=as.character(yr),x=as.character(yr), y=res),
               col="gray40",
               fill="gray60",
               alpha=0.4,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_smooth(data=dplyr::filter(szfreqsim,codsp==sp),
              method = "loess",span=0.1,col="gray20", se = FALSE,
              aes(group=1,x=as.character(yr), y=res,linetype="LOESS"),size=1.5)+
  geom_text(data=dplyr::filter(data.frame(szfreqsim),codsp==sp), 
            aes(x="2014", y=0.9*max(res,na.rm = TRUE), 
                label=paste("RMSE=",round(mean(rmse,na.rm = TRUE),2))))+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Length residual",col="",linetype='',fill='')+
  scale_y_continuous() +
  scale_fill_manual(values = c("blue", "red"))+
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p32

ggsave(paste(c("Length_res_boxplot_Small_Tunas_",sp,".png"),collapse = ''),plot=p32, device = "png", units = "cm",
       width = 24, height = 15)
#x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x#




#X##X##X##X##X##X##X##X##X##X#X##X##X##X##
# models by species
# FRI model - Length Reconstruction 
#X##X##X##X##X##X##X##X##X##X#X##X##X##X##
sp='FRI'

#-------------------------------------
#find the best distribution for specie
#-------------------------------------
dist=data.frame(yr=NULL,codsp=NULL,dist=NULL,codgr=NULL, #empty frame
                mean=NULL,shape=NULL,cv=NULL,aic=NULL) #Shape for gamma is like the mean

for (i in sp) {
  
  years=unique(szfreqobs$yr[szfreqobs$codsp==i])
  
  for (j in years) {
    
    #frequency for each year  
    sz= szfreqobs$fl[szfreqobs$codsp==i & szfreqobs$yr==j]
    gear= unique(szfreqobs$codgr[szfreqobs$codsp==i & szfreqobs$yr==j])[1]
    
    if (length(sz)>5 & (length(unique(sz))==1)==FALSE) {
      
      #graphical analysis
      #fitdistrplus::descdist(sz, discrete = FALSE,boot = 1000)
      
      #fitting distributions
      fit_1 <- fitdistrplus::fitdist(sz, "norm")
      fit_2<- fitdistrplus::fitdist(sz, "invgauss",start=list(mean=mean(sz),shape= 1))
      fit_3 <- fitdistrplus::fitdist(sz, "gamma")
      
      dat=data.frame(yr=j,codsp=i,
                     dist=c(fit_1$distname,fit_2$distname,fit_3$distname),
                     codgr=gear,
                     mean=c(fit_1$estimate[1],fit_2$estimate[1],fit_3$estimate[1]), 
                     shape=c(fit_1$estimate[2],fit_2$estimate[2],fit_3$estimate[2]),
                     cv=c(fit_1$estimate[2]/fit_1$estimate[1],
                          fit_2$estimate[2]/fit_2$estimate[1],
                          fit_3$estimate[2]/fit_3$estimate[1]),
                     aic=c(fit_1$aic,fit_2$aic,fit_3$aic))
      
      dist=rbind(dist,dat)
    }
  }
}

dist=dist %>%
  dplyr::group_by(yr,codsp,codgr) %>% 
  slice_min(aic, n = 1)  %>%
  dplyr::arrange(dist) 

#Table with the best distribution
table(dist$codsp,dist$dist)  
distribution='norm'

#----------------------------------------------------
#predict 1950 length when there was incipient fishing
#---------------------------------------------------
#Assuming 10% higher than the current mean length level
avg1950=1.1 * mean(szmeanobs$mean[szmeanobs$codsp==sp], na.rm = TRUE)

#put in the observed mean lengths
szmeanobs$mean[szmeanobs$codsp==sp & szmeanobs$yr==1950]= round(avg1950,2)
#----------------------------------

#selecting variables .. subseting
szmeanobssub= szmeanobs[,c("yr",'region','codsp','codgr','mean')]

#removing extreme values
szmeanobssub$mean[szmeanobssub$codsp=='BLF' & szmeanobssub$mean<50] =NA
szmeanobssub$mean[szmeanobssub$codsp=='DOL' & szmeanobssub$mean<75] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean>100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean<60] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean<100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean>180] =NA
#szmeanobssub$mean[szmeanobssub$codsp=='FRI' & szmeanobssub$mean>44] =NA


#get only non NA values
szmeanobssub=szmeanobssub[is.na(szmeanobssub$mean)==FALSE & szmeanobssub$codsp==sp,]

#grouping factors most likely (Fredou et al., 2021)
szmeanobssub$region[szmeanobssub$region=='ALL']='SE'
szmeanobssub$region[szmeanobssub$region=='S']='SE'
szmeanobssub$codgr[szmeanobssub$codgr=='UN']='BB'

#transforming columns in to factors
szmeanobssub$region= factor(szmeanobssub$region)
szmeanobssub$codgr= factor(szmeanobssub$codgr)
szmeanobssub$codsp= factor(szmeanobssub$codsp)

# model GAM selection
mod1<- mgcv::gam(mean ~ s(yr)+region,
                 data= szmeanobssub,
                 method = "REML" ,
                 family =gaussian(link='log'))

mod2<- mgcv::gam(mean ~ s(yr)+codgr, 
                 data= szmeanobssub,
                 method = "REML" ,
                 family =gaussian(link='log'))

mod3<- mgcv::gam(mean ~ s(yr)+region+codgr,
                 data= szmeanobssub,
                 method = "REML" ,
                 family =gaussian(link='log'))
#aic check
mod1$aic
mod2$aic
mod3$aic

summary(mod1)
summary(mod2)
summary(mod3)


#final model
modgam<- mgcv::gam(mean ~ s(yr)+codgr, 
                   #sp=0.1,
                   data= szmeanobssub,
                   method = "REML" ,
                   #select = TRUE,
                   family =gaussian(link='log'))

#assign GAM model to the specie
modgam_fri=modgam

#diagnostic
summary(modgam)

par(mfrow=c(2,2))
gam.check(modgam)# model check

#residual 
#homocedasticity
lmtest::bptest(modgam)

#normality
shapiro.test(resid(modgam))

# Check overall concurvity
round(concurvity(modgam, full = TRUE),2)

plot.gam(modgam, # plot partial effects
         rug=TRUE,       #puts the x variable
         #select = c(1), #select some partial effects
         pages = 1,
         all.terms = TRUE,
         residuals = TRUE,
         #cex=2,
         #lwd=2,
         se=TRUE,
         #shade = TRUE,
         #shade.col = "lightblue", #Shading intervals
         shift = coef(modgam)[1]) #shift the scale by


#----------------------
#predictions "year"
newdataexpand= expand.grid(yr=seq(1950,2022,1),
                           region="ALL",
                           codgr=szmeanobssub$codgr)


pred<- data.frame(predict(modgam, newdata = newdataexpand,
                          type = "response", se.fit=TRUE))

pred$up<-pred$fit+qnorm(0.975)*pred$se.fit #confidence intervals
pred$lw<-pred$fit-qnorm(0.975)*pred$se.fit

gam=cbind(type='GAM',
          codsp=sp,
          family=modgam$family$family,
          link=modgam$family$link,
          formula=paste0(c(paste(modgam$formula)[2],
                           paste(modgam$formula)[1],
                           paste(modgam$formula)[3]),
                         collapse = ' '),
          newdataexpand,
          pred,
          aic=modgam$aic,
          r.sq=summary(modgam)$r.sq,
          dev.ex=summary(modgam)$dev.expl)


# Extract mean of all levels
szmeanpred <- gam %>%
  dplyr::group_by(type, codsp, family, link, formula, yr) %>%
  dplyr::summarize(across(c(fit, se.fit, up, lw, aic, r.sq, dev.ex), ~ round(mean(.), 2) ))

#creating data frame
szmeanpred_fri=szmeanpred #size mean predicted for all factors combined
#----------------------------------------------

#partial effects 
newdata= szmeanobssub   #desired partial effect

pareffect<- data.frame(predict(modgam, newdata = newdata,
                               type = "terms", se.fit=TRUE)) #terms gets the effects

#intervals for the partial effect
pareffect$up.codgr<-pareffect$fit.codgr+qnorm(0.975)*pareffect$se.fit.codgr #confidence intervals
pareffect$lw.codgr<-pareffect$fit.codgr-qnorm(0.975)*pareffect$se.fit.codgr
pareffect$codgr= szmeanobssub$codgr

#extract mean of all levels
pareffect= pareffect[grepl("codgr", names(pareffect))]
pareffect <- pareffect %>%
  dplyr::group_by(codgr) %>%
  dplyr::summarize(across(c(fit.codgr, se.fit.codgr, up.codgr, lw.codgr), ~ round(mean(.), 2))) %>%
  dplyr::rename(effect = codgr, 
    fit = fit.codgr, 
    se.fit = se.fit.codgr, 
    up = up.codgr, 
    lw = lw.codgr) %>%
  dplyr::mutate(codsp = sp, var = "Gear")

pareffect_fri=pareffect


#----------------------------------
#simulating frequency distributions
#----------------------------------

#estimating FL distributions from reconstructed mean length
flloop=data.frame(yr=NULL,codsp=NULL,fl=NULL,source=NULL) #empty data frame
sp=sp  #list of species for loop

#coefficient of variation for the corresponding distribution 
cv=mean(dist$cv[dist$codsp==sp & dist$dist==distribution]) #the shape is like mean for gamma, so the cv is corrupted

for (j in 1:length(sp)) {  #loop through all species
  #codes
  years=1950:2022
  n=100 #number of simulate individuals
  
  for (i in 1:length(years)) { #loop through each year of available information
    
    l= data.frame(yr=rep(years[i],n),
                  codsp=rep(sp[j],n),
                  fl=rnorm(n,  #distribution 
                           mean=szmeanpred$fit[szmeanpred$yr==years[i] & szmeanpred$codsp==sp], 
                           sd=szmeanpred$fit[szmeanpred$yr==years[i] & szmeanpred$codsp==sp]*cv),
                  source=rep(as.character('Pred'),n))
    
    flloop=rbind(flloop,l)
  }
}


# size frequency predicted data frame
szfreqpred= flloop
szfreqpred_fri=szfreqpred
#-------------------------------------


#-----------------------------------------
#Observed SPR vs Reconstructed SPR (LBSPR)
#-----------------------------------------
mypars <- new("LB_pars",verbose=F) #blank parameters object
mypars@Species <- sp
mypars@Linf <- pars$linf[pars$codsp==sp]
mypars@L50 <-  pars$l50[pars$codsp==sp] 
mypars@L95 <- 1.1*pars$l50[pars$codsp==sp] 
mypars@MK <- pars$m[pars$codsp==sp]/pars$k[pars$codsp==sp] #m=0.36, ~k=0.16(Nobrega, 2002)
mypars@L_units <- "cm"


#SPR from observed length frequencies 
#getting size classes and counts
szfreqobssub= szfreqobs[szfreqobs$codsp==sp,]
# szfreqobssub= szfreqobssub[szfreqobssub$yr!= 2002,] #removing year with no data
# szfreqobssub= szfreqobssub[szfreqobssub$yr!= 2003,] #removing year with no data
# szfreqobssub= szfreqobssub[szfreqobssub$yr!= 2004,] #removing year with no data
# szfreqobssub= szfreqobssub[szfreqobssub$yr!= 2013,] #removing year with no data
# szfreqobssub= szfreqobssub[szfreqobssub$yr!= 2014,] #removing year with no data
# szfreqobssub= szfreqobssub[szfreqobssub$yr!= 2015,] #removing year with no data
# szfreqobssub= szfreqobssub[szfreqobssub$yr!= 2020,] #removing year with no data


szloop=data.frame(codsp=NULL,yr=NULL,lmids=NULL,n=NULL)

years=unique(szfreqobssub$yr[szfreqobssub$codsp==sp])

for (i in years) {
  h=hist(szfreqobssub$fl[szfreqobssub$yr==i],
         breaks=seq(from=min(szfreqobssub$fl[szfreqobssub$codsp==sp])-5,
                    to=max(szfreqobssub$fl[szfreqobssub$codsp==sp])+5, by=5), #binwidth
         xlab = 'Fork length (cm)',
         main = paste(j))
  freq=data.frame(codsp=sp,yr=i,lmids=h$mids,n=h$counts)    
  szloop=rbind(szloop,freq)
}

#LBSPR lengths format
szfreqobssub =szloop %>% #removing first row 
  tidyr::spread(., yr, n) %>%  #years in columns
  `rownames<-`(.[,2]) %>%   #years in rownames
  dplyr::select(-codsp,-lmids) 

#observed length frequencies object
mylengthsobs <- new("LB_lengths",LB_pars=mypars,dataType="freq",verbose=F) #blank length object

mylengthsobs@LMids<- as.numeric(rownames(szfreqobssub)) #mid points
mylengthsobs@LData<- as.matrix(szfreqobssub)    #frequency table
mylengthsobs@L_units<- "cm"  #units
mylengthsobs@Years<- as.numeric(colnames(szfreqobssub)) #list of years
mylengthsobs@NYears<- length(as.numeric(colnames(szfreqobssub))) #number of years

#plotSize(mylengthsobs)
fitlbspr= LBSPRfit(mypars, mylengthsobs,verbose = FALSE)
#plotSize(fitlbspr)
#plotEsts(fitlbspr)

lbsprobs= data.frame(yr=fitlbspr@Years,  #years 
                     codsp=sp,
                     rawsl50=fitlbspr@SL50, #raw sl50
                     sl50sd=sqrt(fitlbspr@Vars[,1]),  #sl50 sd
                     sl50= fitlbspr@Ests[,1],   #smooth sl50
                     rawsl95=fitlbspr@SL95,
                     sl95sd=sqrt(fitlbspr@Vars[,2]),
                     sl95= fitlbspr@Ests[,2],
                     rawfm=fitlbspr@FM,
                     fmsd=sqrt(fitlbspr@Vars[,3]),
                     fm= fitlbspr@Ests[,3],
                     rawspr=fitlbspr@SPR,
                     sprsd=sqrt(fitlbspr@Vars[,4]),
                     spr= fitlbspr@Ests[,4],
                     source="Obs") #observed length frequency

#----------------------------------------------------------

#SPR from Reconstructed length frequencies 
#getting size classes and counts
szfreqpredsub=szfreqpred[szfreqpred$codsp==sp,]
#szfreqrecsub= szfreqrecsub[szfreqrecsub$yr!= 2017,] #removing year with no data

szloop=data.frame(codsp=NULL,yr=NULL,lmids=NULL,n=NULL)

years=unique(szfreqpredsub$yr[szfreqpredsub$codsp==sp])

for (i in years) {
  h=hist(szfreqpredsub$fl[szfreqpredsub$yr==i],
         breaks=seq(from=min(szfreqpredsub$fl[szfreqpredsub$codsp==sp])-5,
                    to=max(szfreqpredsub$fl[szfreqpredsub$codsp==sp])+5, by=5), #binwidth
         xlab = 'Fork length (cm)',
         main = paste(j))
  freq=data.frame(codsp=sp,yr=i,lmids=round(h$mids,1),n=h$counts)    
  szloop=rbind(szloop,freq)
}

#LBSPR lengths format
szfreqpredsub =szloop %>% #removing first row 
  tidyr::spread(., yr, n) %>%  #years in columns
  `rownames<-`(.[,2]) %>%   #years in rownames
  dplyr::select(-codsp,-lmids) 

#observed length frequencies object
mylengthspred <- new("LB_lengths",LB_pars=mypars,dataType="freq",verbose=F) #blank length object
mylengthspred@LMids<- as.numeric(rownames(szfreqpredsub)) #mid points
mylengthspred@LData<- as.matrix(szfreqpredsub)    #frequency table
mylengthspred@L_units<- "cm"  #units
mylengthspred@Years<- as.numeric(colnames(szfreqpredsub)) #list of years
mylengthspred@NYears<- length(as.numeric(colnames(szfreqpredsub))) #number of years

#plotSize(mylengthspred)

fitlbspr= LBSPRfit(mypars, mylengthspred,verbose = FALSE)
#plotSize(fitlbspr)
#plotEsts(fitlbspr)


lbsprpred= data.frame(yr=fitlbspr@Years,  #years 
                      codsp=sp,
                      rawsl50=fitlbspr@SL50, #raw sl50
                      sl50sd=sqrt(fitlbspr@Vars[,1]),  #sl50 sd
                      sl50= fitlbspr@Ests[,1],   #smooth sl50
                      rawsl95=fitlbspr@SL95,
                      sl95sd=sqrt(fitlbspr@Vars[,2]),
                      sl95= fitlbspr@Ests[,2],
                      rawfm=fitlbspr@FM,
                      fmsd=sqrt(fitlbspr@Vars[,3]),
                      fm= fitlbspr@Ests[,3],
                      rawspr=fitlbspr@SPR,
                      sprsd=sqrt(fitlbspr@Vars[,4]),
                      spr= fitlbspr@Ests[,4],
                      source="Pred") #Predicted length frequency


#binding spr results from observed and reconstructed lengths (data frame)
lbspr= rbind(lbsprobs,lbsprpred)
lbspr_fri= lbspr #general object for plotting

#------------------------------------------------------
#length frequency reconstructed vs LBSPR Pop Simulation 
#------------------------------------------------------
szfreqsim=data.frame(yr=NULL,codsp=NULL,lmids=NULL,nsim=NULL,npred=NULL,res=NULL,rmse=NULL,
                     smd=NULL,ks=NULL,t=NULL,x2=NULL) #empty frame for simulated frequencies

#getting subset for current specie of predicted sizes 
szfreqpredsub=szfreqpred[szfreqpred$codsp==sp,]


for (j in 1:length(sp)) {  #loop through all species codes
  
  years=1950:2022
  
  for (i in 1:length(years)) { #loop through each year of available information
    
    
    #------------------------------------------------------
    #biological parameters for LBSPR simulation
    mypars <- new("LB_pars",verbose=F) #blank parameters object for simulation
    mypars@Species <- sp
    mypars@Linf <- pars$linf[pars$codsp==sp]
    mypars@L50 <-  pars$l50[pars$codsp==sp] 
    mypars@L95 <- 1.1*pars$l50[pars$codsp==sp] 
    mypars@MK <- pars$m[pars$codsp==sp]/pars$k[pars$codsp==sp] #
    mypars@L_units <- "cm"
    
    #Get selectivity and mortality parameters (SL50,SL95, FM) from LBSPR estimation
    mypars@SL50 <- lbsprpred$rawsl50[lbsprpred$yr==years[i]] 
    mypars@SL95 <- lbsprpred$rawsl95[lbsprpred$yr==years[i]]
    mypars@FM <-   lbsprpred$rawfm[lbsprpred$yr==years[i]]
    mypars@BinMax<- 1.4*pars$linf[pars$codsp==sp]
    mypars@BinMin<-  0
    mypars@BinWidth <- 5
    
    #LBSPR simulation
    mysim <- LBSPRsim(mypars)
    
    #normalize (min-max) function to make frequencies comparable
    normalize <- function(x, na.rm = TRUE) {
      return((x- min(x)) /(max(x)-min(x)))
    }
    
    
    lsim= data.frame(yr=rep(years[i], length(mysim@pLPop[,1])),
                     codsp=rep(sp[j],length(mysim@pLPop[,1])),
                     lmids=mysim@pLPop[,1],  
                     nsim=normalize(mysim@pLPop[,5]), #surviving in terms of probability
                     npred=NA)
    
    #-----------------------------------------------------
    # transforming predicted frequencies in length classes
    #we need to cut intervals to match with the simulated sizes
    #creating intervals between min and max length bins 
    lpred= data.frame(fl=szfreqpredsub$fl[szfreqpredsub$codsp==sp  
                                          & szfreqpredsub$yr==years[i]])
    freqtable <- data.frame(table(cut(lpred$fl,
                                      seq(min(mysim@pLPop[,1]),max(mysim@pLPop[,1])+5,5))))
    
    intervals=as.character(freqtable$Var1) #character vector
    intervals=gsub("[()]", "", intervals)       #remove parenthesis
    intervals=sapply(strsplit(intervals, split=',', fixed=TRUE), `[`, 1) #get the first class
    freqtable$Var1=intervals
    
    #normalized predicted frequency
    lsim$npred= normalize(freqtable$Freq)
    #------------------------------------------------------
    
    #Residual length predicted and simulated frequencies
    lsim$res= lsim$nsim-lsim$npred
    lsim$res[lsim$npred==0]=NA #avoiding Residual length in classes with no data
    
    #Root Mean Squared Error- RMSE
    rmse= function(x, y, na.rm = TRUE) {
      return(sqrt(sum((x- y)^2) /length(x)))
    }
    
    lsim$rmse= rmse(lsim$npred,lsim$nsim)
    
    #standardized mean reserence (SMD)
    smd= function(x, y, na.rm = TRUE) {
      return( abs(mean(x)- mean(y)) / sqrt((sd(x)^2+sd(y)^2)/2) )
    }
    
    lsim$smd= smd(lsim$npred,lsim$nsim)
    
    
    #Kolmogorov-smirnov test (comparison of cumulative distributions)
    #we need to convert probabilities in sample distributions to perform KS test
    lpredfreq=lpred$fl #predicted length distribution
    
    #simulated frequency distribution 
    lsimfreq= data.frame(lmids=mysim@pLPop[,1],n=round(mysim@pLPop[,5]*100))
    lsimfreq= lsimfreq %>%
      tidyr::uncount(n)
    lsimfreq=lsimfreq$lmids
    
    # Kolmogorov smirnov test to compare if the samples come from the same probability distributions
    lsim$ks=suppressWarnings(ks.test(lpredfreq,lsimfreq)$p.value) #kolmogorov-smirnov test
    
    #t-student test to compare if the means of distributions are the same
    lsim$t= suppressWarnings(t.test(lpredfreq,lsimfreq)$p.value)
    
    # chi-square goodness-of-fit test
    lsim$x2= suppressWarnings(chisq.test(sample(lpredfreq, 95,replace = FALSE), sample(lsimfreq, 95,replace = FALSE))$p.value)
    
    #bind into size frequency simulate data frame
    szfreqsim=rbind(szfreqsim,lsim)
  }
}
#get comparison of frequencies data frame 
szfreqsim_fri= szfreqsim
#--------------------------------------


#-------------
#plot section
#-------------
#plot GAM mean size prediction
p33<-ggplot(data=szmeanpred,aes(x=yr,y=fit))+
  geom_ribbon(data=dplyr::filter(szmeanpred,codsp==sp),
              aes(ymin = lw, ymax = up,linetype='Pred'), 
              col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(szmeanpred,codsp==sp),
            aes(x=yr,y=fit,linetype='Pred'),col="gray60",size=2)+
  #geom_line(data=dplyr::filter(gam,codsp==sp),aes(x=yr,y=fit,col=codgr,linetype=region),size=1.2)+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=yr,y=mean),size=2,col="gray40")+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='')+
  scale_y_continuous() +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p33


ggsave(paste(c("Length_model_prediction_Small_Tunas_",sp,".png"),collapse=''), plot = p33, device = "png", units = "cm",
       width = 24, height = 14)

# plot GAM partial effects
p34<-ggplot(dplyr::filter(pareffect,codsp==sp), 
            aes(x=effect,y=fit,fill=var,col="Gear"),size=2) + 
  #geom_line(aes(x=effect,y=up, col="Gear",group=1)) +
  geom_point(aes(x=effect,y=fit,col="Gear"),size=3.5)+
  geom_errorbar(aes(ymin=lw, ymax=up,col="Gear"), width=0.5,
                position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Effect",x="Gear",col='',fill="")+
  scale_y_continuous() +
  scale_x_discrete() +
  scale_colour_manual(values="gray60") +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p34

ggsave(paste(c("Effects_model_Small_Tunas_",sp,".png"),collapse=''), plot = p34, device = "png", units = "cm",
       width = 24, height = 14)


#plot GAM boxplot size freq prediction
p35<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqpred,codsp==sp), 
               aes(x=factor(yr), y=fl),
               col="gray40",
               fill="gray60",
               alpha=0.6,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=factor(yr),y=mean),size=2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Fork length (cm)",col="Region",linetype='',fill='Gear')+
  scale_y_continuous() +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p35

ggsave(paste(c("Length_pred_boxplot_Small_Tunas_",sp,".png"),collapse = ''),plot=p35, device = "png", units = "cm",
       width = 24, height = 15)


#plot half violin for the predicted frequencies
p36<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(szfreqpred,codsp==sp), 
                             aes(x=factor(yr),y=fl),
                             fill="gray20",
                             col="gray20",
                             alpha=0.7,
                             side = "r")+ 
  #show.legend = TRUE,
  #position = "identity",
  #fill=,
  #nudge = 0,
  #scale = 'area',
  #na.rm = TRUE,
  #lwd=1,
  #trim = TRUE)+
  #facet_wrap(~codsp,scales= "free_y")+
  #draw_quantiles = c(0.5))+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=factor(yr),y=mean),size=2)+
  # #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p36

ggsave(paste(c("Length_pred_density_Small_Tunas_",sp,".png"),collapse = ''),plot=p36, device = "png", units = "cm",
       width = 35, height = 12)


#plot SPR trajectories comparison (observed vs Predicted)
p37<-ggplot(dplyr::filter(lbspr,codsp==sp), 
            aes(x=yr,y=rawspr,fill=source,col=source),size=2) + 
  #geom_line(size=1,aes(fill=source)) +
  geom_point(size=2)+
  geom_smooth(size=1.2,method = 'loess' ,formula = 'y ~ x',se =FALSE,span = 0.4)+
  geom_errorbar(aes(ymin=ifelse(rawspr-sprsd < 0, 0, rawspr-sprsd),
                    ymax=ifelse(rawspr+sprsd > 1, 1, rawspr+sprsd)), width=1.5,
                position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Spawning Potential Ratio (SPR)",x="Year",col='',fill="")+
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks = seq(1950, 2023, 30),limits = c(1950, 2023)) +
  scale_colour_manual(values = c("gray60","gray40")) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p37

ggsave(paste(c("SPR_Small_Tunas_",sp,".png"),collapse=''),plot=p37, device = "png", units = "cm",
       width = 24, height = 14)



#plot GAM boxplot size freq prediction
p38<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqsim,codsp==sp), 
               aes(group=as.character(yr),x=as.character(yr), y=res),
               col="gray40",
               fill="gray60",
               alpha=0.4,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_smooth(data=dplyr::filter(szfreqsim,codsp==sp),
              method = "loess",span=0.1,col="gray20", se = FALSE,
              aes(group=1,x=as.character(yr), y=res,linetype="LOESS"),size=1.5)+
  geom_text(data=dplyr::filter(data.frame(szfreqsim),codsp==sp), 
            aes(x="2014", y=0.9*max(res,na.rm = TRUE), 
                label=paste("RMSE=",round(mean(rmse,na.rm = TRUE),2))))+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Length residual",col="",linetype='',fill='')+
  scale_y_continuous() +
  scale_fill_manual(values = c("blue", "red"))+
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p38

ggsave(paste(c("Length_res_boxplot_Small_Tunas_",sp,".png"),collapse = ''),plot=p38, device = "png", units = "cm",
       width = 24, height = 15)
#x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x#




#X##X##X##X##X##X##X##X##X##X#X##X##X##X##
# models by species
# KGM model - Length Reconstruction 
#X##X##X##X##X##X##X##X##X##X#X##X##X##X##
sp='KGM'

#-------------------------------------
#find the best distribution for specie
#-------------------------------------
dist=data.frame(yr=NULL,codsp=NULL,dist=NULL,codgr=NULL, #empty frame
                mean=NULL,shape=NULL,cv=NULL,aic=NULL) #Shape for gamma is like the mean

for (i in sp) {
  
  years=unique(szfreqobs$yr[szfreqobs$codsp==i])
  
  for (j in years) {
    
    #frequency for each year  
    sz= szfreqobs$fl[szfreqobs$codsp==i & szfreqobs$yr==j]
    gear= unique(szfreqobs$codgr[szfreqobs$codsp==i & szfreqobs$yr==j])[1]
    
    if (length(sz)>5 & (length(unique(sz))==1)==FALSE) {
      
      #graphical analysis
      #fitdistrplus::descdist(sz, discrete = FALSE,boot = 1000)
      
      #fitting distributions
      fit_1 <- fitdistrplus::fitdist(sz, "norm")
      fit_2<- fitdistrplus::fitdist(sz, "invgauss",start=list(mean=mean(sz),shape= 1))
      fit_3 <- fitdistrplus::fitdist(sz, "gamma")
      
      dat=data.frame(yr=j,codsp=i,
                     dist=c(fit_1$distname,fit_2$distname,fit_3$distname),
                     codgr=gear,
                     mean=c(fit_1$estimate[1],fit_2$estimate[1],fit_3$estimate[1]), 
                     shape=c(fit_1$estimate[2],fit_2$estimate[2],fit_3$estimate[2]),
                     cv=c(fit_1$estimate[2]/fit_1$estimate[1],
                          fit_2$estimate[2]/fit_2$estimate[1],
                          fit_3$estimate[2]/fit_3$estimate[1]),
                     aic=c(fit_1$aic,fit_2$aic,fit_3$aic))
      
      dist=rbind(dist,dat)
    }
  }
}

dist=dist %>%
  dplyr::group_by(yr,codsp,codgr) %>% 
  slice_min(aic, n = 1)  %>%
  dplyr::arrange(dist) 

#Table with the best distribution
table(dist$codsp,dist$dist)  
distribution='invgauss'

#----------------------------------------------------
#predict 1950 length when there was incipient fishing
#---------------------------------------------------
#Assuming 10% higher than the current mean length level
avg1950=1.1 * mean(szmeanobs$mean[szmeanobs$codsp==sp], na.rm = TRUE)
#put in the observed mean lengths
szmeanobs$mean[szmeanobs$codsp==sp & szmeanobs$yr==1950]= round(avg1950,2)
#----------------------------------

#selecting variables .. subseting
szmeanobssub= szmeanobs[,c("yr",'region','codsp','codgr','mean')]

#removing extreme values
szmeanobssub$mean[szmeanobssub$codsp=='BLF' & szmeanobssub$mean<50] =NA
szmeanobssub$mean[szmeanobssub$codsp=='DOL' & szmeanobssub$mean<75] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean>100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean<60] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean<100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean>180] =NA


#get only non NA values
szmeanobssub=szmeanobssub[is.na(szmeanobssub$mean)==FALSE & szmeanobssub$codsp==sp,]

#grouping factors most likely (Fredou et al., 2021)
szmeanobssub$region[szmeanobssub$region=='ALL']='NE'

#transforming columns in to factors
szmeanobssub$region= factor(szmeanobssub$region)
szmeanobssub$codgr= factor(szmeanobssub$codgr)
szmeanobssub$codsp= factor(szmeanobssub$codsp)


# model GAM selection
mod1<- mgcv::gam(mean ~ s(yr)+region,
                 data= szmeanobssub,
                 method = "REML" ,
                 family =inverse.gaussian(link='log'))

mod2<- mgcv::gam(mean ~ s(yr)+codgr, 
                 data= szmeanobssub,
                 method = "REML" ,
                 family =inverse.gaussian(link='log'))

mod3<- mgcv::gam(mean ~ s(yr)+region+codgr,
                 data= szmeanobssub,
                 method = "REML" ,
                 family =inverse.gaussian(link='log'))
#aic check
mod1$aic
mod2$aic
mod3$aic

summary(mod1)
summary(mod2)
summary(mod3)


#final model
modgam<- mgcv::gam(mean ~ s(yr)+region+codgr,
                   sp=0.001,
                   data= szmeanobssub,
                   method = "REML" ,
                   #select = TRUE,
                   family =inverse.gaussian(link='log'))

#assign GAM model to the specie
modgam_kgm=modgam

#diagnostic
summary(modgam)

par(mfrow=c(2,2))
gam.check(modgam)# model check

#residual 
#homocedasticity
lmtest::bptest(modgam)

#normality
shapiro.test(resid(modgam))

# Check overall concurvity
round(concurvity(modgam, full = TRUE),2)

plot.gam(modgam, # plot partial effects
         rug=TRUE,       #puts the x variable
         #select = c(1), #select some partial effects
         pages = 1,
         all.terms = TRUE,
         residuals = TRUE,
         #cex=2,
         #lwd=2,
         se=TRUE,
         #shade = TRUE,
         #shade.col = "lightblue", #Shading intervals
         shift = coef(modgam)[1]) #shift the scale by


#----------------------
#predictions "year"
newdataexpand= expand.grid(yr=seq(1950,2022,1),
                           region=szmeanobssub$region,
                           codgr=szmeanobssub$codgr)


pred<- data.frame(predict(modgam, newdata = newdataexpand,
                          type = "response", se.fit=TRUE))

pred$up<-pred$fit+qnorm(0.975)*pred$se.fit #confidence intervals
pred$lw<-pred$fit-qnorm(0.975)*pred$se.fit

gam=cbind(type='GAM',
          codsp=sp,
          family=modgam$family$family,
          link=modgam$family$link,
          formula=paste0(c(paste(modgam$formula)[2],
                           paste(modgam$formula)[1],
                           paste(modgam$formula)[3]),
                         collapse = ' '),
          newdataexpand,
          pred,
          aic=modgam$aic,
          r.sq=summary(modgam)$r.sq,
          dev.ex=summary(modgam)$dev.expl)


# Extract mean of all levels
szmeanpred <- gam %>%
  dplyr::group_by(type, codsp, family, link, formula, yr) %>%
  dplyr::summarize(across(c(fit, se.fit, up, lw, aic, r.sq, dev.ex), ~ round(mean(.), 2)))

#creating data frame
szmeanpred_kgm=szmeanpred #size mean predicted for all factors combined
#----------------------------------------------

#partial effects 
newdata= szmeanobssub   #desired partial effect

pareffect<- data.frame(predict(modgam, newdata = newdata,
                               type = "terms", se.fit=TRUE)) #terms gets the effects

#intervals for the partial effect
pareffect$up.codgr<-pareffect$fit.codgr+qnorm(0.975)*pareffect$se.fit.codgr #confidence intervals
pareffect$lw.codgr<-pareffect$fit.codgr-qnorm(0.975)*pareffect$se.fit.codgr
pareffect$codgr= szmeanobssub$codgr

pareffect$up.region<-pareffect$fit.region+qnorm(0.975)*pareffect$se.fit.region #confidence intervals
pareffect$lw.region<-pareffect$fit.region-qnorm(0.975)*pareffect$se.fit.region
pareffect$region= szmeanobssub$region


#extract mean of all levels (gear)
pareffectgear= pareffect[grepl("codgr", names(pareffect))]
pareffectgear <- pareffectgear %>%
  dplyr::group_by(codgr) %>%
  dplyr::summarize(across(c(fit.codgr, se.fit.codgr, up.codgr, lw.codgr), ~ round(mean(.), 2))) %>%
  dplyr::rename(effect = codgr, 
    fit = fit.codgr, 
    se.fit = se.fit.codgr, 
    up = up.codgr, 
    lw = lw.codgr) %>%
  dplyr::mutate(codsp = sp, var = "Gear")

#extract mean of all levels (Region)
pareffectregion= pareffect[grepl("region", names(pareffect))]
pareffectregion <- pareffectregion %>%
  dplyr::group_by(region) %>%
  dplyr::summarize(across(c(fit.region, se.fit.region, up.region, lw.region), ~ round(mean(.), 2))) %>%
  dplyr::rename(effect = region, 
    fit = fit.region, 
    se.fit = se.fit.region, 
    up = up.region, 
    lw = lw.region) %>%
  dplyr::mutate(codsp = sp, var = "Region")

pareffect=rbind(pareffectregion,pareffectgear)

pareffect_kgm=pareffect


#----------------------------------
#simulating frequency distributions
#----------------------------------

#estimating FL distributions from reconstructed mean length
flloop=data.frame(yr=NULL,codsp=NULL,fl=NULL,source=NULL) #empty data frame
sp=sp  #list of species for loop

#coefficient of variation for the corresponding distribution 
cv=mean(dist$cv[dist$codsp==sp & dist$dist==distribution]) #the shape is like mean for gamma, so the cv is corrupted

for (j in 1:length(sp)) {  #loop through all species
  #codes
  years=1950:2022
  n=100 #number of simulate individuals
  
  for (i in 1:length(years)) { #loop through each year of available information
    
    l= data.frame(yr=rep(years[i],n),
                  codsp=rep(sp[j],n),
                  fl=rinvgauss(n,  #distribution 
                           mean=szmeanpred$fit[szmeanpred$yr==years[i] & szmeanpred$codsp==sp], 
                           shape=szmeanpred$fit[szmeanpred$yr==years[i] & szmeanpred$codsp==sp]*cv),
                  source=rep(as.character('Pred'),n))
    
    flloop=rbind(flloop,l)
  }
}


# size frequency predicted data frame
szfreqpred= flloop
szfreqpred_kgm=szfreqpred
#-------------------------------------


#-----------------------------------------
#Observed SPR vs Reconstructed SPR (LBSPR)
#-----------------------------------------
mypars <- new("LB_pars",verbose=F) #blank parameters object
mypars@Species <- sp
mypars@Linf <- pars$linf[pars$codsp==sp]
mypars@L50 <-  pars$l50[pars$codsp==sp] 
mypars@L95 <- 1.1*pars$l50[pars$codsp==sp] 
mypars@MK <- pars$m[pars$codsp==sp]/pars$k[pars$codsp==sp] 
mypars@L_units <- "cm"


#SPR from observed length frequencies 
#getting size classes and counts
szfreqobssub= szfreqobs[szfreqobs$codsp==sp,]

szloop=data.frame(codsp=NULL,yr=NULL,lmids=NULL,n=NULL)

years=unique(szfreqobssub$yr[szfreqobssub$codsp==sp])

for (i in years) {
  h=hist(szfreqobssub$fl[szfreqobssub$yr==i],
         breaks=seq(from=min(szfreqobssub$fl[szfreqobssub$codsp==sp])-5,
                    to=max(szfreqobssub$fl[szfreqobssub$codsp==sp])+5, by=5), #binwidth
         xlab = 'Fork length (cm)',
         main = paste(j))
  freq=data.frame(codsp=sp,yr=i,lmids=h$mids,n=h$counts)    
  szloop=rbind(szloop,freq)
}

#LBSPR lengths format
szfreqobssub =szloop %>% #removing first row 
  tidyr::spread(., yr, n) %>%  #years in columns
  `rownames<-`(.[,2]) %>%   #years in rownames
  dplyr::select(-codsp,-lmids) 

#observed length frequencies object
mylengthsobs <- new("LB_lengths",LB_pars=mypars,dataType="freq",verbose=F) #blank length object

mylengthsobs@LMids<- as.numeric(rownames(szfreqobssub)) #mid points
mylengthsobs@LData<- as.matrix(szfreqobssub)    #frequency table
mylengthsobs@L_units<- "cm"  #units
mylengthsobs@Years<- as.numeric(colnames(szfreqobssub)) #list of years
mylengthsobs@NYears<- length(as.numeric(colnames(szfreqobssub))) #number of years

#plotSize(mylengthsobs)
fitlbspr= LBSPRfit(mypars, mylengthsobs,verbose = FALSE)
#plotSize(fitlbspr)
#plotEsts(fitlbspr)

lbsprobs= data.frame(yr=fitlbspr@Years,  #years 
                     codsp=sp,
                     rawsl50=fitlbspr@SL50, #raw sl50
                     sl50sd=sqrt(fitlbspr@Vars[,1]),  #sl50 sd
                     sl50= fitlbspr@Ests[,1],   #smooth sl50
                     rawsl95=fitlbspr@SL95,
                     sl95sd=sqrt(fitlbspr@Vars[,2]),
                     sl95= fitlbspr@Ests[,2],
                     rawfm=fitlbspr@FM,
                     fmsd=sqrt(fitlbspr@Vars[,3]),
                     fm= fitlbspr@Ests[,3],
                     rawspr=fitlbspr@SPR,
                     sprsd=sqrt(fitlbspr@Vars[,4]),
                     spr= fitlbspr@Ests[,4],
                     source="Obs") #observed length frequency

#----------------------------------------------------------

#SPR from Reconstructed length frequencies 
#getting size classes and counts
szfreqpredsub=szfreqpred[szfreqpred$codsp==sp,]
#szfreqrecsub= szfreqrecsub[szfreqrecsub$yr!= 2017,] #removing year with no data

szloop=data.frame(codsp=NULL,yr=NULL,lmids=NULL,n=NULL)

years=unique(szfreqpredsub$yr[szfreqpredsub$codsp==sp])

for (i in years) {
  h=hist(szfreqpredsub$fl[szfreqpredsub$yr==i],
         breaks=seq(from=min(szfreqpredsub$fl[szfreqpredsub$codsp==sp])-5,
                    to=max(szfreqpredsub$fl[szfreqpredsub$codsp==sp])+5, by=5), #binwidth
         xlab = 'Fork length (cm)',
         main = paste(j))
  freq=data.frame(codsp=sp,yr=i,lmids=round(h$mids,1),n=h$counts)    
  szloop=rbind(szloop,freq)
}

#LBSPR lengths format
szfreqpredsub =szloop %>% #removing first row 
  tidyr::spread(., yr, n) %>%  #years in columns
  `rownames<-`(.[,2]) %>%   #years in rownames
  dplyr::select(-codsp,-lmids) 

#observed length frequencies object
mylengthspred <- new("LB_lengths",LB_pars=mypars,dataType="freq",verbose=F) #blank length object
mylengthspred@LMids<- as.numeric(rownames(szfreqpredsub)) #mid points
mylengthspred@LData<- as.matrix(szfreqpredsub)    #frequency table
mylengthspred@L_units<- "cm"  #units
mylengthspred@Years<- as.numeric(colnames(szfreqpredsub)) #list of years
mylengthspred@NYears<- length(as.numeric(colnames(szfreqpredsub))) #number of years

#plotSize(mylengthspred)

fitlbspr= LBSPRfit(mypars, mylengthspred,verbose = FALSE)
#plotSize(fitlbspr)
#plotEsts(fitlbspr)


lbsprpred= data.frame(yr=fitlbspr@Years,  #years 
                      codsp=sp,
                      rawsl50=fitlbspr@SL50, #raw sl50
                      sl50sd=sqrt(fitlbspr@Vars[,1]),  #sl50 sd
                      sl50= fitlbspr@Ests[,1],   #smooth sl50
                      rawsl95=fitlbspr@SL95,
                      sl95sd=sqrt(fitlbspr@Vars[,2]),
                      sl95= fitlbspr@Ests[,2],
                      rawfm=fitlbspr@FM,
                      fmsd=sqrt(fitlbspr@Vars[,3]),
                      fm= fitlbspr@Ests[,3],
                      rawspr=fitlbspr@SPR,
                      sprsd=sqrt(fitlbspr@Vars[,4]),
                      spr= fitlbspr@Ests[,4],
                      source="Pred") #Predicted length frequency


#binding spr results from observed and reconstructed lengths (data frame)
lbspr= rbind(lbsprobs,lbsprpred)
lbspr_kgm= lbspr #general object for plotting

#------------------------------------------------------
#length frequency reconstructed vs LBSPR Pop Simulation 
#------------------------------------------------------
szfreqsim=data.frame(yr=NULL,codsp=NULL,lmids=NULL,nsim=NULL,npred=NULL,res=NULL,rmse=NULL,
                     smd=NULL,ks=NULL,t=NULL,x2=NULL) #empty frame for simulated frequencies

#getting subset for current specie of predicted sizes 
szfreqpredsub=szfreqpred[szfreqpred$codsp==sp,]


for (j in 1:length(sp)) {  #loop through all species codes
  
  years=1950:2022
  
  for (i in 1:length(years)) { #loop through each year of available information
    
    
    #------------------------------------------------------
    #biological parameters for LBSPR simulation
    mypars <- new("LB_pars",verbose=F) #blank parameters object for simulation
    mypars@Species <- sp
    mypars@Linf <- pars$linf[pars$codsp==sp]
    mypars@L50 <-  pars$l50[pars$codsp==sp] 
    mypars@L95 <- 1.1*pars$l50[pars$codsp==sp] 
    mypars@MK <- pars$m[pars$codsp==sp]/pars$k[pars$codsp==sp] #
    mypars@L_units <- "cm"
    
    #Get selectivity and mortality parameters (SL50,SL95, FM) from LBSPR estimation
    mypars@SL50 <- lbsprpred$rawsl50[lbsprpred$yr==years[i]] 
    mypars@SL95 <- lbsprpred$rawsl95[lbsprpred$yr==years[i]]
    mypars@FM <-   lbsprpred$rawfm[lbsprpred$yr==years[i]]
    mypars@BinMax<- 1.4*pars$linf[pars$codsp==sp]
    mypars@BinMin<-  0
    mypars@BinWidth <- 5
    
    #LBSPR simulation
    mysim <- LBSPRsim(mypars)
    
    #normalize (min-max) function to make frequencies comparable
    normalize <- function(x, na.rm = TRUE) {
      return((x- min(x)) /(max(x)-min(x)))
    }
    
    
    lsim= data.frame(yr=rep(years[i], length(mysim@pLPop[,1])),
                     codsp=rep(sp[j],length(mysim@pLPop[,1])),
                     lmids=mysim@pLPop[,1],  
                     nsim=normalize(mysim@pLPop[,5]), #surviving in terms of probability
                     npred=NA)
    
    #-----------------------------------------------------
    # transforming predicted frequencies in length classes
    #we need to cut intervals to match with the simulated sizes
    #creating intervals between min and max length bins 
    lpred= data.frame(fl=szfreqpredsub$fl[szfreqpredsub$codsp==sp  
                                          & szfreqpredsub$yr==years[i]])
    freqtable <- data.frame(table(cut(lpred$fl,
                                      seq(min(mysim@pLPop[,1]),max(mysim@pLPop[,1])+5,5))))
    
    intervals=as.character(freqtable$Var1) #character vector
    intervals=gsub("[()]", "", intervals)       #remove parenthesis
    intervals=sapply(strsplit(intervals, split=',', fixed=TRUE), `[`, 1) #get the first class
    freqtable$Var1=intervals
    
    #normalized predicted frequency
    lsim$npred= normalize(freqtable$Freq)
    #------------------------------------------------------
    
    #Residual length predicted and simulated frequencies
    lsim$res= lsim$nsim-lsim$npred
    lsim$res[lsim$npred==0]=NA #avoiding Residual length in classes with no data
    
    #Root Mean Squared Error- RMSE
    rmse= function(x, y, na.rm = TRUE) {
      return(sqrt(sum((x- y)^2) /length(x)))
    }
    
    lsim$rmse= rmse(lsim$npred,lsim$nsim)
    
    #standardized mean reserence (SMD)
    smd= function(x, y, na.rm = TRUE) {
      return( abs(mean(x)- mean(y)) / sqrt((sd(x)^2+sd(y)^2)/2) )
    }
    
    lsim$smd= smd(lsim$npred,lsim$nsim)
    
    
    #Kolmogorov-smirnov test (comparison of cumulative distributions)
    #we need to convert probabilities in sample distributions to perform KS test
    lpredfreq=lpred$fl #predicted length distribution
    
    #simulated frequency distribution 
    lsimfreq= data.frame(lmids=mysim@pLPop[,1],n=round(mysim@pLPop[,5]*100))
    lsimfreq= lsimfreq %>%
      tidyr::uncount(n)
    lsimfreq=lsimfreq$lmids
    
    # Kolmogorov smirnov test to compare if the samples come from the same probability distributions
    lsim$ks=suppressWarnings(ks.test(lpredfreq,lsimfreq)$p.value) #kolmogorov-smirnov test
    
    #t-student test to compare if the means of distributions are the same
    lsim$t=suppressWarnings(t.test(lpredfreq,lsimfreq)$p.value)
    
    # chi-square goodness-of-fit test
    lsim$x2=suppressWarnings(chisq.test(sample(lpredfreq, 95,replace = FALSE), sample(lsimfreq, 95,replace = FALSE))$p.value)
    
    #bind into size frequency simulate data frame
    szfreqsim=rbind(szfreqsim,lsim)
  }
}

#get comparison of frequencies data frame 
szfreqsim_kgm= szfreqsim
#--------------------------------------


#-------------
#plot section
#-------------
#plot GAM mean size prediction
p38<-ggplot(data=szmeanpred,aes(x=yr,y=fit))+
  geom_ribbon(data=dplyr::filter(szmeanpred,codsp==sp),
              aes(ymin = lw, ymax = up,linetype='Pred'), 
              col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(szmeanpred,codsp==sp),
            aes(x=yr,y=fit,linetype='Pred'),col="gray60",size=2)+
  #geom_line(data=dplyr::filter(gam,codsp==sp),aes(x=yr,y=fit,col=codgr,linetype=region),size=1.2)+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=yr,y=mean),size=2,col="gray40")+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='')+
  scale_y_continuous() +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p38

ggsave(paste(c("Length_model_prediction_Small_Tunas_",sp,".png"),collapse=''), plot = p38, device = "png", units = "cm",
       width = 24, height = 14)

# plot GAM partial effects
p39<-ggplot(dplyr::filter(pareffect,codsp==sp), 
            aes(x=effect,y=fit,col=var),size=2) + 
  #geom_line(aes(x=effect,y=up, col="Gear",group=1)) +
  geom_point(aes(x=effect,y=fit,col=var),size=3.5)+
  geom_errorbar(aes(ymin=lw, ymax=up,col=var), width=0.5,
                position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Effect",x="Factor",col='',fill="")+
  scale_y_continuous() +
  scale_x_discrete() +
  scale_colour_manual(values=c("gray60","gray30")) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p39

ggsave(paste(c("Effects_model_Small_Tunas_",sp,".png"),collapse=''), plot = p39, device = "png", units = "cm",
       width = 24, height = 14)


#plot GAM boxplot size freq prediction
p40<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqpred,codsp==sp), 
               aes(x=factor(yr), y=fl),
               col="gray40",
               fill="gray60",
               alpha=0.6,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=factor(yr),y=mean),size=2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Fork length (cm)",col="Region",linetype='',fill='Gear')+
  scale_y_continuous() +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p40

ggsave(paste(c("Length_pred_boxplot_Small_Tunas_",sp,".png"),collapse = ''),plot=p40, device = "png", units = "cm",
       width = 24, height = 15)


#plot half violin for the predicted frequencies
p41<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(szfreqpred,codsp==sp), 
                             aes(x=factor(yr),y=fl),
                             fill="gray20",
                             col="gray20",
                             alpha=0.7,
                             side = "r")+ 
  #show.legend = TRUE,
  #position = "identity",
  #fill=,
  #nudge = 0,
  #scale = 'area',
  #na.rm = TRUE,
  #lwd=1,
  #trim = TRUE)+
  #facet_wrap(~codsp,scales= "free_y")+
  #draw_quantiles = c(0.5))+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=factor(yr),y=mean),size=2)+
  # #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p41

ggsave(paste(c("Length_pred_density_Small_Tunas_",sp,".png"),collapse = ''),plot=p41, device = "png", units = "cm",
       width = 35, height = 12)


#plot SPR trajectories comparison (observed vs Predicted)
p42<-ggplot(dplyr::filter(lbspr,codsp==sp), 
            aes(x=yr,y=rawspr,fill=source,col=source),size=2) + 
  #geom_line(size=1,aes(fill=source)) +
  geom_point(size=2)+
  geom_smooth(size=1.2,method = 'loess' ,formula = 'y ~ x',se =FALSE)+
  geom_errorbar(aes(ymin=ifelse(rawspr-sprsd < 0, 0, rawspr-sprsd),
                    ymax=ifelse(rawspr+sprsd > 1, 1, rawspr+sprsd)), width=1.5,
                position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Spawning Potential Ratio (SPR)",x="Year",col='',fill="")+
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks = seq(1950, 2023, 30),limits = c(1950, 2023)) +
  scale_colour_manual(values = c("gray60","gray40")) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p42

ggsave(paste(c("SPR_Small_Tunas_",sp,".png"),collapse=''),plot=p42, device = "png", units = "cm",
       width = 24, height = 14)



#plot GAM boxplot size freq prediction
p43<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqsim,codsp==sp), 
               aes(group=as.character(yr),x=as.character(yr), y=res),
               col="gray40",
               fill="gray60",
               alpha=0.4,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_smooth(data=dplyr::filter(szfreqsim,codsp==sp),
              method = "loess",span=0.1,col="gray20", se = FALSE,
              aes(group=1,x=as.character(yr), y=res,linetype="LOESS"),size=1.5)+
  geom_text(data=dplyr::filter(data.frame(szfreqsim),codsp==sp), 
            aes(x="2014", y=0.9*max(res,na.rm = TRUE), 
                label=paste("RMSE=",round(mean(rmse,na.rm = TRUE),2))))+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Length residual",col="",linetype='',fill='')+
  scale_y_continuous() +
  scale_fill_manual(values = c("blue", "red"))+
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p43

ggsave(paste(c("Length_res_boxplot_Small_Tunas_",sp,".png"),collapse = ''),plot=p43, device = "png", units = "cm",
       width = 24, height = 15)
#x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x#





#X##X##X##X##X##X##X##X##X##X#X##X##X##X##
# models by species
# LTA model - Length Reconstruction 
#X##X##X##X##X##X##X##X##X##X#X##X##X##X##
sp='LTA'

#-------------------------------------
#find the best distribution for specie
#-------------------------------------
dist=data.frame(yr=NULL,codsp=NULL,dist=NULL,codgr=NULL, #empty frame
                mean=NULL,shape=NULL,cv=NULL,aic=NULL) #Shape for gamma is like the mean

for (i in sp) {
  
  years=unique(szfreqobs$yr[szfreqobs$codsp==i])
  
  for (j in years) {
    
    #frequency for each year  
    sz= szfreqobs$fl[szfreqobs$codsp==i & szfreqobs$yr==j]
    gear= unique(szfreqobs$codgr[szfreqobs$codsp==i & szfreqobs$yr==j])[1]
    
    if (length(sz)>5 & (length(unique(sz))==1)==FALSE) {
      
      #graphical analysis
      #fitdistrplus::descdist(sz, discrete = FALSE,boot = 1000)
      
      #fitting distributions
      fit_1 <- fitdistrplus::fitdist(sz, "norm")
      fit_2<- fitdistrplus::fitdist(sz, "invgauss",start=list(mean=mean(sz),shape= 1))
      fit_3 <- fitdistrplus::fitdist(sz, "gamma")
      
      dat=data.frame(yr=j,codsp=i,
                     dist=c(fit_1$distname,fit_2$distname,fit_3$distname),
                     codgr=gear,
                     mean=c(fit_1$estimate[1],fit_2$estimate[1],fit_3$estimate[1]), 
                     shape=c(fit_1$estimate[2],fit_2$estimate[2],fit_3$estimate[2]),
                     cv=c(fit_1$estimate[2]/fit_1$estimate[1],
                          fit_2$estimate[2]/fit_2$estimate[1],
                          fit_3$estimate[2]/fit_3$estimate[1]),
                     aic=c(fit_1$aic,fit_2$aic,fit_3$aic))
      
      dist=rbind(dist,dat)
    }
  }
}

dist=dist %>%
  dplyr::group_by(yr,codsp,codgr) %>% 
  slice_min(aic, n = 1)  %>%
  dplyr::arrange(dist) 

#Table with the best distribution
table(dist$codsp,dist$dist)  
distribution='invgauss'

#----------------------------------------------------
#predict 1950 length when there was incipient fishing
#---------------------------------------------------
#Assuming 10% higher than the current mean length level
avg1950=1.1 * mean(szmeanobs$mean[szmeanobs$codsp==sp], na.rm = TRUE)

#put in the observed mean lengths
szmeanobs$mean[szmeanobs$codsp==sp & szmeanobs$yr==1950]= round(avg1950,2)
#----------------------------------

#selecting variables .. subseting
szmeanobssub= szmeanobs[,c("yr",'region','codsp','codgr','mean')]

#removing extreme values
szmeanobssub$mean[szmeanobssub$codsp=='BLF' & szmeanobssub$mean<50] =NA
szmeanobssub$mean[szmeanobssub$codsp=='DOL' & szmeanobssub$mean<75] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean>100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean<60] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean<100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean>180] =NA


#get only non NA values
szmeanobssub=szmeanobssub[is.na(szmeanobssub$mean)==FALSE & szmeanobssub$codsp==sp,]

#grouping factors most likely (Fredou et al., 2021)
szmeanobssub$region[szmeanobssub$region=='ALL']='SE'

#transforming columns in to factors
szmeanobssub$region= factor(szmeanobssub$region)
szmeanobssub$codgr= factor(szmeanobssub$codgr)
szmeanobssub$codsp= factor(szmeanobssub$codsp)

# model GAM selection
mod1<- mgcv::gam(mean ~ s(yr),
                 data= szmeanobssub,
                 method = "REML" ,
                 family =inverse.gaussian(link='log'))

# mod2<- mgcv::gam(mean ~ s(yr)+region,
#                  data= szmeanobssub,
#                  method = "REML" ,
#                  family =inverse.gaussian(link='log'))
# 
# mod3<- mgcv::gam(mean ~ s(yr)+codgr, 
#                  data= szmeanobssub,
#                  method = "REML" ,
#                  family =inverse.gaussian(link='log'))
# 
# mod4<- mgcv::gam(mean ~ s(yr)+region+codgr,
#                  data= szmeanobssub,
#                  method = "REML" ,
#                  family =inverse.gaussian(link='log'))

#aic check
mod1$aic
mod2$aic
mod3$aic

summary(mod1)
summary(mod2)
summary(mod3)


#final model
modgam<- mgcv::gam(mean ~ s(yr),
                   sp=0.001,
                   data= szmeanobssub,
                   method = "REML" ,
                   family =inverse.gaussian(link='log'))

#assign GAM model to the specie
modgam_lta=modgam

#diagnostic
summary(modgam)

par(mfrow=c(2,2))
gam.check(modgam)# model check

#residual 
#homocedasticity
lmtest::bptest(modgam)

#normality
shapiro.test(resid(modgam))

# Check overall concurvity
round(concurvity(modgam, full = TRUE),2)

plot.gam(modgam, # plot partial effects
         rug=TRUE,       #puts the x variable
         #select = c(1), #select some partial effects
         pages = 1,
         all.terms = TRUE,
         residuals = TRUE,
         #cex=2,
         #lwd=2,
         se=TRUE,
         #shade = TRUE,
         #shade.col = "lightblue", #Shading intervals
         shift = coef(modgam)[1]) #shift the scale by


#----------------------
#predictions "year"
newdataexpand= expand.grid(yr=seq(1950,2022,1),
                           region="ALL",
                           codgr="UN")


pred<- data.frame(predict(modgam, newdata = newdataexpand,
                          type = "response", se.fit=TRUE))

pred$up<-pred$fit+qnorm(0.975)*pred$se.fit #confidence intervals
pred$lw<-pred$fit-qnorm(0.975)*pred$se.fit

gam=cbind(type='GAM',
          codsp=sp,
          family=modgam$family$family,
          link=modgam$family$link,
          formula=paste0(c(paste(modgam$formula)[2],
                           paste(modgam$formula)[1],
                           paste(modgam$formula)[3]),
                         collapse = ' '),
          newdataexpand,
          pred,
          aic=modgam$aic,
          r.sq=summary(modgam)$r.sq,
          dev.ex=summary(modgam)$dev.expl)


# Extract mean of all levels
szmeanpred <- gam %>%
  dplyr::group_by(type, codsp, family, link, formula, yr) %>%
  dplyr::summarize(across(c(fit, se.fit, up, lw, aic, r.sq, dev.ex), ~ round(mean(.), 2)))

#creating data frame
szmeanpred_lta=szmeanpred #size mean predicted for all factors combined
#----------------------------------------------

#partial effects 
# newdata= szmeanobssub   #desired partial effect
# 
# pareffect<- data.frame(predict(modgam, newdata = newdata,
#                                type = "terms", se.fit=TRUE)) #terms gets the effects
# 
# #intervals for the partial effect
# pareffect$up.codgr<-pareffect$fit.codgr+qnorm(0.975)*pareffect$se.fit.codgr #confidence intervals
# pareffect$lw.codgr<-pareffect$fit.codgr-qnorm(0.975)*pareffect$se.fit.codgr
# pareffect$codgr= szmeanobssub$codgr
# 
# pareffect$up.region<-pareffect$fit.region+qnorm(0.975)*pareffect$se.fit.region #confidence intervals
# pareffect$lw.region<-pareffect$fit.region-qnorm(0.975)*pareffect$se.fit.region
# pareffect$region= szmeanobssub$region
# 
# 
# #extract mean of all levels (gear)
# pareffectgear= pareffect[grepl("codgr", names(pareffect))]
# pareffectgear= pareffectgear %>% #select only the desired partial effect
#   dplyr::group_by(codgr) %>%
#   dplyr::summarize_at(c("fit.codgr", "se.fit.codgr","up.codgr","lw.codgr"),
#                       funs(round(mean(.), 2))) %>%
#   dplyr::rename(.,  "effect"="codgr",
#                 "fit"= "fit.codgr", 
#                 "se.fit"="se.fit.codgr",
#                 "up"="up.codgr",
#                 "lw"="lw.codgr")  %>%
#   dplyr::mutate(., codsp=sp, var="Gear")
# 
# 
# #extract mean of all levels (Region)
# pareffectregion= pareffect[grepl("region", names(pareffect))]
# pareffectregion= pareffectregion %>% #select only the desired partial effect
#   dplyr::group_by(region) %>%
#   dplyr::summarize_at(c("fit.region", "se.fit.region","up.region","lw.region"),
#                       funs(round(mean(.), 2))) %>%
#   dplyr::rename(.,  "effect"="region",
#                 "fit"= "fit.region", 
#                 "se.fit"="se.fit.region",
#                 "up"="up.region",
#                 "lw"="lw.region")  %>%
#   dplyr::mutate(., codsp=sp, var="Region")
# 
# pareffect=rbind(pareffectregion,pareffectgear)
# 
# pareffect_lta=pareffect


#----------------------------------
#simulating frequency distributions
#----------------------------------

#estimating FL distributions from reconstructed mean length
flloop=data.frame(yr=NULL,codsp=NULL,fl=NULL,source=NULL) #empty data frame
sp=sp  #list of species for loop

#coefficient of variation for the corresponding distribution 
cv=mean(dist$cv[dist$codsp==sp & dist$dist==distribution]) #the shape is like mean for gamma, so the cv is corrupted

for (j in 1:length(sp)) {  #loop through all species
  #codes
  years=1950:2022
  n=100 #number of simulate individuals
  
  for (i in 1:length(years)) { #loop through each year of available information
    
    l= data.frame(yr=rep(years[i],n),
                  codsp=rep(sp[j],n),
                  fl=rinvgauss(n,  #distribution 
                               mean=szmeanpred$fit[szmeanpred$yr==years[i] & szmeanpred$codsp==sp], 
                               shape=szmeanpred$fit[szmeanpred$yr==years[i] & szmeanpred$codsp==sp]*cv),
                  source=rep(as.character('Pred'),n))
    
    flloop=rbind(flloop,l)
  }
}


# size frequency predicted data frame
szfreqpred= flloop
szfreqpred_lta=szfreqpred
#-------------------------------------


#-----------------------------------------
#Observed SPR vs Reconstructed SPR (LBSPR)
#-----------------------------------------
mypars <- new("LB_pars",verbose=F) #blank parameters object
mypars@Species <- sp
mypars@Linf <- pars$linf[pars$codsp==sp]
mypars@L50 <-  pars$l50[pars$codsp==sp] 
mypars@L95 <- 1.1*pars$l50[pars$codsp==sp] 
mypars@MK <- pars$m[pars$codsp==sp]/pars$k[pars$codsp==sp] 
mypars@L_units <- "cm"


#SPR from observed length frequencies 
#getting size classes and counts
szfreqobssub= szfreqobs[szfreqobs$codsp==sp,]

szloop=data.frame(codsp=NULL,yr=NULL,lmids=NULL,n=NULL)

years=unique(szfreqobssub$yr[szfreqobssub$codsp==sp])

for (i in years) {
  h=hist(szfreqobssub$fl[szfreqobssub$yr==i],
         breaks=seq(from=min(szfreqobssub$fl[szfreqobssub$codsp==sp])-5,
                    to=max(szfreqobssub$fl[szfreqobssub$codsp==sp])+60, by=5), #binwidth
         xlab = 'Fork length (cm)',
         main = paste(j))
  freq=data.frame(codsp=sp,yr=i,lmids=h$mids,n=h$counts)    
  szloop=rbind(szloop,freq)
}

#LBSPR lengths format
szfreqobssub =szloop %>% #removing first row 
  tidyr::spread(., yr, n) %>%  #years in columns
  `rownames<-`(.[,2]) %>%   #years in rownames
  dplyr::select(-codsp,-lmids) 

#observed length frequencies object
mylengthsobs <- new("LB_lengths",LB_pars=mypars,dataType="freq",verbose=F) #blank length object

mylengthsobs@LMids<- as.numeric(rownames(szfreqobssub)) #mid points
mylengthsobs@LData<- as.matrix(szfreqobssub)    #frequency table
mylengthsobs@L_units<- "cm"  #units
mylengthsobs@Years<- as.numeric(colnames(szfreqobssub)) #list of years
mylengthsobs@NYears<- length(as.numeric(colnames(szfreqobssub))) #number of years

#plotSize(mylengthsobs)
fitlbspr= LBSPRfit(mypars, mylengthsobs,verbose = FALSE)
#plotSize(fitlbspr)
#plotEsts(fitlbspr)

lbsprobs= data.frame(yr=fitlbspr@Years,  #years 
                     codsp=sp,
                     rawsl50=fitlbspr@SL50, #raw sl50
                     sl50sd=sqrt(fitlbspr@Vars[,1]),  #sl50 sd
                     sl50= fitlbspr@Ests[,1],   #smooth sl50
                     rawsl95=fitlbspr@SL95,
                     sl95sd=sqrt(fitlbspr@Vars[,2]),
                     sl95= fitlbspr@Ests[,2],
                     rawfm=fitlbspr@FM,
                     fmsd=sqrt(fitlbspr@Vars[,3]),
                     fm= fitlbspr@Ests[,3],
                     rawspr=fitlbspr@SPR,
                     sprsd=sqrt(fitlbspr@Vars[,4]),
                     spr= fitlbspr@Ests[,4],
                     source="Obs") #observed length frequency

#----------------------------------------------------------

#SPR from Reconstructed length frequencies 
#getting size classes and counts
szfreqpredsub=szfreqpred[szfreqpred$codsp==sp,]
#szfreqrecsub= szfreqrecsub[szfreqrecsub$yr!= 2017,] #removing year with no data

szloop=data.frame(codsp=NULL,yr=NULL,lmids=NULL,n=NULL)

years=unique(szfreqpredsub$yr[szfreqpredsub$codsp==sp])

for (i in years) {
  h=hist(szfreqpredsub$fl[szfreqpredsub$yr==i],
         breaks=seq(from=min(szfreqpredsub$fl[szfreqpredsub$codsp==sp])-5,
                    to=max(szfreqpredsub$fl[szfreqpredsub$codsp==sp])+60, by=5), #binwidth
         xlab = 'Fork length (cm)',
         main = paste(j))
  freq=data.frame(codsp=sp,yr=i,lmids=round(h$mids,1),n=h$counts)    
  szloop=rbind(szloop,freq)
}

#LBSPR lengths format
szfreqpredsub =szloop %>% #removing first row 
  tidyr::spread(., yr, n) %>%  #years in columns
  `rownames<-`(.[,2]) %>%   #years in rownames
  dplyr::select(-codsp,-lmids) 

#observed length frequencies object
mylengthspred <- new("LB_lengths",LB_pars=mypars,dataType="freq",verbose=F) #blank length object
mylengthspred@LMids<- as.numeric(rownames(szfreqpredsub)) #mid points
mylengthspred@LData<- as.matrix(szfreqpredsub)    #frequency table
mylengthspred@L_units<- "cm"  #units
mylengthspred@Years<- as.numeric(colnames(szfreqpredsub)) #list of years
mylengthspred@NYears<- length(as.numeric(colnames(szfreqpredsub))) #number of years

#plotSize(mylengthspred)

fitlbspr= LBSPRfit(mypars, mylengthspred,verbose = FALSE)
#plotSize(fitlbspr)
#plotEsts(fitlbspr)


lbsprpred= data.frame(yr=fitlbspr@Years,  #years 
                      codsp=sp,
                      rawsl50=fitlbspr@SL50, #raw sl50
                      sl50sd=sqrt(fitlbspr@Vars[,1]),  #sl50 sd
                      sl50= fitlbspr@Ests[,1],   #smooth sl50
                      rawsl95=fitlbspr@SL95,
                      sl95sd=sqrt(fitlbspr@Vars[,2]),
                      sl95= fitlbspr@Ests[,2],
                      rawfm=fitlbspr@FM,
                      fmsd=sqrt(fitlbspr@Vars[,3]),
                      fm= fitlbspr@Ests[,3],
                      rawspr=fitlbspr@SPR,
                      sprsd=sqrt(fitlbspr@Vars[,4]),
                      spr= fitlbspr@Ests[,4],
                      source="Pred") #Predicted length frequency


#binding spr results from observed and reconstructed lengths (data frame)
lbspr= rbind(lbsprobs,lbsprpred)
lbspr_lta= lbspr #general object for plotting


#------------------------------------------------------
#length frequency reconstructed vs LBSPR Pop Simulation 
#------------------------------------------------------
szfreqsim=data.frame(yr=NULL,codsp=NULL,lmids=NULL,nsim=NULL,npred=NULL,res=NULL,rmse=NULL,
                     smd=NULL,ks=NULL,t=NULL,x2=NULL) #empty frame for simulated frequencies

#getting subset for current specie of predicted sizes 
szfreqpredsub=szfreqpred[szfreqpred$codsp==sp,]


for (j in 1:length(sp)) {  #loop through all species codes
  
  years=1950:2022
  
  for (i in 1:length(years)) { #loop through each year of available information
    
    
    #------------------------------------------------------
    #biological parameters for LBSPR simulation
    mypars <- new("LB_pars",verbose=F) #blank parameters object for simulation
    mypars@Species <- sp
    mypars@Linf <- pars$linf[pars$codsp==sp]
    mypars@L50 <-  pars$l50[pars$codsp==sp] 
    mypars@L95 <- 1.1*pars$l50[pars$codsp==sp] 
    mypars@MK <- pars$m[pars$codsp==sp]/pars$k[pars$codsp==sp] #
    mypars@L_units <- "cm"
    
    #Get selectivity and mortality parameters (SL50,SL95, FM) from LBSPR estimation
    mypars@SL50 <- lbsprpred$rawsl50[lbsprpred$yr==years[i]] 
    mypars@SL95 <- lbsprpred$rawsl95[lbsprpred$yr==years[i]]
    mypars@FM <-   lbsprpred$rawfm[lbsprpred$yr==years[i]]
    mypars@BinMax<- 1.4*pars$linf[pars$codsp==sp]
    mypars@BinMin<-  0
    mypars@BinWidth <- 5
    
    #LBSPR simulation
    mysim <- LBSPRsim(mypars)
    
    #normalize (min-max) function to make frequencies comparable
    normalize <- function(x, na.rm = TRUE) {
      return((x- min(x)) /(max(x)-min(x)))
    }
    
    
    lsim= data.frame(yr=rep(years[i], length(mysim@pLPop[,1])),
                     codsp=rep(sp[j],length(mysim@pLPop[,1])),
                     lmids=mysim@pLPop[,1],  
                     nsim=normalize(mysim@pLPop[,5]), #surviving in terms of probability
                     npred=NA)
    
    #-----------------------------------------------------
    # transforming predicted frequencies in length classes
    #we need to cut intervals to match with the simulated sizes
    #creating intervals between min and max length bins 
    lpred= data.frame(fl=szfreqpredsub$fl[szfreqpredsub$codsp==sp  
                                          & szfreqpredsub$yr==years[i]])
    freqtable <- data.frame(table(cut(lpred$fl,
                                      seq(min(mysim@pLPop[,1]),max(mysim@pLPop[,1])+5,5))))
    
    intervals=as.character(freqtable$Var1) #character vector
    intervals=gsub("[()]", "", intervals)       #remove parenthesis
    intervals=sapply(strsplit(intervals, split=',', fixed=TRUE), `[`, 1) #get the first class
    freqtable$Var1=intervals
    
    #normalized predicted frequency
    lsim$npred= normalize(freqtable$Freq)
    #------------------------------------------------------
    
    #Residual length predicted and simulated frequencies
    lsim$res= lsim$nsim-lsim$npred
    lsim$res[lsim$npred==0]=NA #avoiding Residual length in classes with no data
    
    #Root Mean Squared Error- RMSE
    rmse= function(x, y, na.rm = TRUE) {
      return(sqrt(sum((x- y)^2) /length(x)))
    }
    
    lsim$rmse= rmse(lsim$npred,lsim$nsim)
    
    #standardized mean reserence (SMD)
    smd= function(x, y, na.rm = TRUE) {
      return( abs(mean(x)- mean(y)) / sqrt((sd(x)^2+sd(y)^2)/2) )
    }
    
    lsim$smd= smd(lsim$npred,lsim$nsim)
    
    
    #Kolmogorov-smirnov test (comparison of cumulative distributions)
    #we need to convert probabilities in sample distributions to perform KS test
    lpredfreq=lpred$fl #predicted length distribution
    
    #simulated frequency distribution 
    lsimfreq= data.frame(lmids=mysim@pLPop[,1],n=round(mysim@pLPop[,5]*100))
    lsimfreq= lsimfreq %>%
      tidyr::uncount(n)
    lsimfreq=lsimfreq$lmids
    
    # Kolmogorov smirnov test to compare if the samples come from the same probability distributions
    lsim$ks=suppressWarnings(ks.test(lpredfreq,lsimfreq)$p.value) #kolmogorov-smirnov test
    
    #t-student test to compare if the means of distributions are the same
    lsim$t= suppressWarnings(t.test(lpredfreq,lsimfreq)$p.value)
    
    # chi-square goodness-of-fit test
    lsim$x2=suppressWarnings(chisq.test(sample(lpredfreq, 95,replace = FALSE), sample(lsimfreq, 95,replace = FALSE))$p.value)
    
    #bind into size frequency simulate data frame
    szfreqsim=rbind(szfreqsim,lsim)
  }
}

#get comparison of frequencies data frame 
szfreqsim_lta= szfreqsim
#--------------------------------------


#-------------
#plot section
#-------------
#plot GAM mean size prediction
p44<-ggplot(data=szmeanpred,aes(x=yr,y=fit))+
  geom_ribbon(data=dplyr::filter(szmeanpred,codsp==sp),
              aes(ymin = lw, ymax = up,linetype='Pred'), 
              col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(szmeanpred,codsp==sp),
            aes(x=yr,y=fit,linetype='Pred'),col="gray60",size=2)+
  #geom_line(data=dplyr::filter(gam,codsp==sp),aes(x=yr,y=fit,col=codgr,linetype=region),size=1.2)+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=yr,y=mean),size=2,col="gray40")+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='')+
  scale_y_continuous() +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p44

ggsave(paste(c("Length_model_prediction_Small_Tunas_",sp,".png"),collapse=''),plot=p44, device = "png", units = "cm",
       width = 24, height = 14)

# plot GAM partial effects
# p45<-ggplot(dplyr::filter(pareffect,codsp==sp), 
#             aes(x=effect,y=fit,col=var),size=2) + 
#   #geom_line(aes(x=effect,y=up, col="Gear",group=1)) +
#   geom_point(aes(x=effect,y=fit,col=var),size=3.5)+
#   geom_errorbar(aes(ymin=lw, ymax=up,col=var), width=0.5,
#                 position="identity",size=1.2)+
#   #facet_wrap(~codsp,scales = "free_y")+
#   labs(y="Effect",x="Factor",col='',fill="")+
#   scale_y_continuous() +
#   scale_x_discrete() +
#   scale_colour_manual(values=c("gray60","gray30")) +
#   theme_classic(base_size = 14) %+replace%
#   theme(strip.background=element_blank(),
#         plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
# p45
# 
# ggsave(paste(c("Effects_model_Small_Tunas_",sp,".png"),collapse=''),plot= p45, device = "png", units = "cm",
#        width = 24, height = 14)


#plot GAM boxplot size freq prediction
p46<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqpred,codsp==sp), 
               aes(x=factor(yr), y=fl),
               col="gray40",
               fill="gray60",
               alpha=0.6,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=factor(yr),y=mean),size=2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Fork length (cm)",col="Region",linetype='',fill='Gear')+
  scale_y_continuous() +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p46

ggsave(paste(c("Length_pred_boxplot_Small_Tunas_",sp,".png"),collapse = ''),plot=p46, device = "png", units = "cm",
       width = 24, height = 15)


#plot half violin for the predicted frequencies
p47<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(szfreqpred,codsp==sp), 
                             aes(x=factor(yr),y=fl),
                             fill="gray20",
                             col="gray20",
                             alpha=0.7,
                             side = "r")+ 
  #show.legend = TRUE,
  #position = "identity",
  #fill=,
  #nudge = 0,
  #scale = 'area',
  #na.rm = TRUE,
  #lwd=1,
  #trim = TRUE)+
  #facet_wrap(~codsp,scales= "free_y")+
  #draw_quantiles = c(0.5))+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=factor(yr),y=mean),size=2)+
  # #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p47

ggsave(paste(c("Length_pred_density_Small_Tunas_",sp,".png"),collapse = ''),plot=p47, device = "png", units = "cm",
       width = 35, height = 12)


#plot SPR trajectories comparison (observed vs Predicted)
p48<-ggplot(dplyr::filter(lbspr,codsp==sp), 
            aes(x=yr,y=rawspr,fill=source,col=source),size=2) + 
  #geom_line(size=1,aes(fill=source)) +
  geom_point(size=2)+
  geom_smooth(size=1.2,method = 'loess' ,formula = 'y ~ x',se =FALSE)+
  geom_errorbar(aes(ymin=ifelse(rawspr-sprsd < 0, 0, rawspr-sprsd),
                    ymax=ifelse(rawspr+sprsd > 1, 1, rawspr+sprsd)), width=1.5,
                position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Spawning Potential Ratio (SPR)",x="Year",col='',fill="")+
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks = seq(1950, 2023, 30),limits = c(1950, 2023)) +
  scale_colour_manual(values = c("gray60","gray40")) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p48

ggsave(paste(c("SPR_Small_Tunas_",sp,".png"),collapse=''),plot=p48, device = "png", units = "cm",
       width = 24, height = 14)



#plot GAM boxplot size freq prediction
p49<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqsim,codsp==sp), 
               aes(group=as.character(yr),x=as.character(yr), y=res),
               col="gray40",
               fill="gray60",
               alpha=0.4,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_smooth(data=dplyr::filter(szfreqsim,codsp==sp),
              method = "loess",span=0.1,col="gray20", se = FALSE,
              aes(group=1,x=as.character(yr), y=res,linetype="LOESS"),size=1.5)+
  geom_text(data=dplyr::filter(data.frame(szfreqsim),codsp==sp), 
            aes(x="2014", y=0.9*max(res,na.rm = TRUE), 
                label=paste("RMSE=",round(mean(rmse,na.rm = TRUE),2))))+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Length residual",col="",linetype='',fill='')+
  scale_y_continuous() +
  scale_fill_manual(values = c("blue", "red"))+
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p49

ggsave(paste(c("Length_res_boxplot_Small_Tunas_",sp,".png"),collapse = ''),plot=p49, device = "png", units = "cm",
       width = 24, height = 15)
#x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x#




#X##X##X##X##X##X##X##X##X##X#X##X##X##X
# models by species
# WAH model - Length Reconstruction 
#X##X##X##X##X##X##X##X##X##X#X##X##X##X
sp='WAH'

#-------------------------------------
#find the best distribution for specie
#-------------------------------------
dist=data.frame(yr=NULL,codsp=NULL,dist=NULL,codgr=NULL, #empty frame
                mean=NULL,shape=NULL,cv=NULL,aic=NULL) #Shape for gamma is like the mean

for (i in sp) {
  
  years=unique(szfreqobs$yr[szfreqobs$codsp==i])
  
  for (j in years) {
    
    #frequency for each year  
    sz= szfreqobs$fl[szfreqobs$codsp==i & szfreqobs$yr==j]
    gear= unique(szfreqobs$codgr[szfreqobs$codsp==i & szfreqobs$yr==j])[1]
    
    if (length(sz)>5 & (length(unique(sz))==1)==FALSE) {
      
      #graphical analysis
      #fitdistrplus::descdist(sz, discrete = FALSE,boot = 1000)
      
      #fitting distributions
      fit_1 <- fitdistrplus::fitdist(sz, "norm")
      fit_2<- fitdistrplus::fitdist(sz, "invgauss",start=list(mean=mean(sz),shape= 1))
      fit_3 <- fitdistrplus::fitdist(sz, "gamma")
      
      dat=data.frame(yr=j,codsp=i,
                     dist=c(fit_1$distname,fit_2$distname,fit_3$distname),
                     codgr=gear,
                     mean=c(fit_1$estimate[1],fit_2$estimate[1],fit_3$estimate[1]), 
                     shape=c(fit_1$estimate[2],fit_2$estimate[2],fit_3$estimate[2]),
                     cv=c(fit_1$estimate[2]/fit_1$estimate[1],
                          fit_2$estimate[2]/fit_2$estimate[1],
                          fit_3$estimate[2]/fit_3$estimate[1]),
                     aic=c(fit_1$aic,fit_2$aic,fit_3$aic))
      
      dist=rbind(dist,dat)
    }
  }
}

dist=dist %>%
  dplyr::group_by(yr,codsp,codgr) %>% 
  slice_min(aic, n = 1)  %>%
  dplyr::arrange(dist) 

#Table with the best distribution
table(dist$codsp,dist$dist)  
distribution='norm'

#----------------------------------------------------
#predict 1950 length when there was incipient fishing
#---------------------------------------------------
#Assuming 10% higher than the current mean length level
avg1950=1.1 * mean(szmeanobs$mean[szmeanobs$codsp==sp], na.rm = TRUE)
#put in the observed mean lengths
szmeanobs$mean[szmeanobs$codsp==sp & szmeanobs$yr==1950]= round(avg1950,2)
#----------------------------------

#selecting variables .. subseting
szmeanobssub= szmeanobs[,c("yr",'region','codsp','codgr','mean')]

#removing extreme values
szmeanobssub$mean[szmeanobssub$codsp=='BLF' & szmeanobssub$mean<50] =NA
szmeanobssub$mean[szmeanobssub$codsp=='DOL' & szmeanobssub$mean<75] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean>100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean<60] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean<100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean>180] =NA


#get only non NA values
szmeanobssub=szmeanobssub[is.na(szmeanobssub$mean)==FALSE & szmeanobssub$codsp==sp,]

#transforming columns in to factors
szmeanobssub$region= factor(szmeanobssub$region)
szmeanobssub$codgr= factor(szmeanobssub$codgr)
szmeanobssub$codsp= factor(szmeanobssub$codsp)


# model GAM selection
mod1<- mgcv::gam(mean ~ s(yr),
                 data= szmeanobssub,
                 method = "REML" ,
                 family =gaussian(link='log'))

mod2<- mgcv::gam(mean ~ s(yr)+region,
                 data= szmeanobssub,
                 method = "REML" ,
                 family =gaussian(link='log'))

mod3<- mgcv::gam(mean ~ s(yr)+codgr,
                 data= szmeanobssub,
                 method = "REML" ,
                 family =gaussian(link='log'))

mod4<- mgcv::gam(mean ~ s(yr)+region+codgr,
                 data= szmeanobssub,
                 method = "REML" ,
                 family =gaussian(link='log'))

#aic check
mod1$aic
mod2$aic
mod3$aic
mod4$aic

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)


#final model
modgam<- mgcv::gam(mean ~ s(yr)+region+codgr,
                   sp=c(10),
                   data= szmeanobssub,
                   method = "REML" ,
                   family =gaussian(link='log'))

#assign GAM model to the specie
modgam_wah=modgam

#diagnostic
summary(modgam)

par(mfrow=c(2,2))
gam.check(modgam)# model check

#residual 
#homocedasticity
lmtest::bptest(modgam)

#normality
shapiro.test(resid(modgam))

# Check overall concurvity
round(concurvity(modgam, full = TRUE),2)

plot.gam(modgam, # plot partial effects
         rug=TRUE,       #puts the x variable
         #select = c(1), #select some partial effects
         pages = 1,
         all.terms = TRUE,
         residuals = TRUE,
         #cex=2,
         #lwd=2,
         se=TRUE,
         #shade = TRUE,
         #shade.col = "lightblue", #Shading intervals
         shift = coef(modgam)[1]) #shift the scale by


#----------------------
#predictions "year"
newdataexpand= expand.grid(yr=seq(1950,2022,1),
                           region=szmeanobssub$region,
                           codgr=szmeanobssub$codgr)


pred<- data.frame(predict(modgam, newdata = newdataexpand,
                          type = "response", se.fit=TRUE))

pred$up<-pred$fit+qnorm(0.975)*pred$se.fit #confidence intervals
pred$lw<-pred$fit-qnorm(0.975)*pred$se.fit

gam=cbind(type='GAM',
          codsp=sp,
          family=modgam$family$family,
          link=modgam$family$link,
          formula=paste0(c(paste(modgam$formula)[2],
                           paste(modgam$formula)[1],
                           paste(modgam$formula)[3]),
                         collapse = ' '),
          newdataexpand,
          pred,
          aic=modgam$aic,
          r.sq=summary(modgam)$r.sq,
          dev.ex=summary(modgam)$dev.expl)


# Extract mean of all levels
szmeanpred <- gam %>%
  dplyr::group_by(type, codsp, family, link, formula, yr) %>%
  dplyr::summarize(across(c(fit, se.fit, up, lw, aic, r.sq, dev.ex), ~ round(mean(.), 2)))

#creating data frame
szmeanpred_wah=szmeanpred #size mean predicted for all factors combined
#----------------------------------------------

#partial effects 
newdata= szmeanobssub   #desired partial effect

pareffect<- data.frame(predict(modgam, newdata = newdata,
                               type = "terms", se.fit=TRUE)) #terms gets the effects

#intervals for the partial effect
pareffect$up.codgr<-pareffect$fit.codgr+qnorm(0.975)*pareffect$se.fit.codgr #confidence intervals
pareffect$lw.codgr<-pareffect$fit.codgr-qnorm(0.975)*pareffect$se.fit.codgr
pareffect$codgr= szmeanobssub$codgr

pareffect$up.region<-pareffect$fit.region+qnorm(0.975)*pareffect$se.fit.region #confidence intervals
pareffect$lw.region<-pareffect$fit.region-qnorm(0.975)*pareffect$se.fit.region
pareffect$region= szmeanobssub$region


#extract mean of all levels (gear)
pareffectgear= pareffect[grepl("codgr", names(pareffect))]
pareffectgear <- pareffectgear %>%
  dplyr::group_by(codgr) %>%
  dplyr::summarize(across(c(fit.codgr, se.fit.codgr, up.codgr, lw.codgr), ~ round(mean(.), 2))) %>%
  dplyr::rename(effect = codgr, 
    fit = fit.codgr, 
    se.fit = se.fit.codgr, 
    up = up.codgr, 
    lw = lw.codgr) %>%
  dplyr::mutate(codsp = sp, var = "Gear")


#extract mean of all levels (Region)
pareffectregion= pareffect[grepl("region", names(pareffect))]
pareffectregion <- pareffectregion %>%
  dplyr::group_by(region) %>%
  dplyr::summarize(across(c(fit.region, se.fit.region, up.region, lw.region), ~ round(mean(.), 2))) %>%
  dplyr::rename(effect = region, 
    fit = fit.region, 
    se.fit = se.fit.region, 
    up = up.region, 
    lw = lw.region) %>%
  dplyr::mutate(codsp = sp, var = "Region")

pareffect=rbind(pareffectregion,pareffectgear)

pareffect_wah=pareffect


#----------------------------------
#simulating frequency distributions
#----------------------------------

#estimating FL distributions from reconstructed mean length
flloop=data.frame(yr=NULL,codsp=NULL,fl=NULL,source=NULL) #empty data frame
sp=sp  #list of species for loop

#coefficient of variation for the corresponding distribution 
cv=mean(dist$cv[dist$codsp==sp & dist$dist==distribution]) #the shape is like mean for gamma, so the cv is corrupted

for (j in 1:length(sp)) {  #loop through all species
  #codes
  years=1950:2022
  n=100 #number of simulate individuals
  
  for (i in 1:length(years)) { #loop through each year of available information
    
    l= data.frame(yr=rep(years[i],n),
                  codsp=rep(sp[j],n),
                  fl=rnorm(n,  #distribution 
                           mean=szmeanpred$fit[szmeanpred$yr==years[i] & szmeanpred$codsp==sp], 
                           sd=szmeanpred$fit[szmeanpred$yr==years[i] & szmeanpred$codsp==sp]*cv),
                  source=rep(as.character('Pred'),n))
    
    flloop=rbind(flloop,l)
  }
}


# size frequency predicted data frame
szfreqpred= flloop
szfreqpred_wah=szfreqpred
#-------------------------------------


#-----------------------------------------
#Observed SPR vs Reconstructed SPR (LBSPR)
#-----------------------------------------
mypars <- new("LB_pars",verbose=F) #blank parameters object
mypars@Species <- sp
mypars@Linf <- pars$linf[pars$codsp==sp]
mypars@L50 <-  pars$l50[pars$codsp==sp] 
mypars@L95 <- 1.1*pars$l50[pars$codsp==sp] 
mypars@MK <- pars$m[pars$codsp==sp]/pars$k[pars$codsp==sp] 
mypars@L_units <- "cm"


#SPR from observed length frequencies 
#getting size classes and counts
szfreqobssub= szfreqobs[szfreqobs$codsp==sp,]

szloop=data.frame(codsp=NULL,yr=NULL,lmids=NULL,n=NULL)

years=unique(szfreqobssub$yr[szfreqobssub$codsp==sp])

for (i in years) {
  h=hist(szfreqobssub$fl[szfreqobssub$yr==i],
         breaks=seq(from=min(szfreqobssub$fl[szfreqobssub$codsp==sp])-5,
                    to=max(szfreqobssub$fl[szfreqobssub$codsp==sp])+60, by=5), #binwidth
         xlab = 'Fork length (cm)',
         main = paste(j))
  freq=data.frame(codsp=sp,yr=i,lmids=h$mids,n=h$counts)    
  szloop=rbind(szloop,freq)
}

#LBSPR lengths format
szfreqobssub =szloop %>% #removing first row 
  tidyr::spread(., yr, n) %>%  #years in columns
  `rownames<-`(.[,2]) %>%   #years in rownames
  dplyr::select(-codsp,-lmids) 

#observed length frequencies object
mylengthsobs <- new("LB_lengths",LB_pars=mypars,dataType="freq",verbose=F) #blank length object

mylengthsobs@LMids<- as.numeric(rownames(szfreqobssub)) #mid points
mylengthsobs@LData<- as.matrix(szfreqobssub)    #frequency table
mylengthsobs@L_units<- "cm"  #units
mylengthsobs@Years<- as.numeric(colnames(szfreqobssub)) #list of years
mylengthsobs@NYears<- length(as.numeric(colnames(szfreqobssub))) #number of years

#plotSize(mylengthsobs)
fitlbspr= LBSPRfit(mypars, mylengthsobs,verbose = FALSE)
#plotSize(fitlbspr)
#plotEsts(fitlbspr)

lbsprobs= data.frame(yr=fitlbspr@Years,  #years 
                     codsp=sp,
                     rawsl50=fitlbspr@SL50, #raw sl50
                     sl50sd=sqrt(fitlbspr@Vars[,1]),  #sl50 sd
                     sl50= fitlbspr@Ests[,1],   #smooth sl50
                     rawsl95=fitlbspr@SL95,
                     sl95sd=sqrt(fitlbspr@Vars[,2]),
                     sl95= fitlbspr@Ests[,2],
                     rawfm=fitlbspr@FM,
                     fmsd=sqrt(fitlbspr@Vars[,3]),
                     fm= fitlbspr@Ests[,3],
                     rawspr=fitlbspr@SPR,
                     sprsd=sqrt(fitlbspr@Vars[,4]),
                     spr= fitlbspr@Ests[,4],
                     source="Obs") #observed length frequency

#----------------------------------------------------------

#SPR from Reconstructed length frequencies 
#getting size classes and counts
szfreqpredsub=szfreqpred[szfreqpred$codsp==sp,]
#szfreqrecsub= szfreqrecsub[szfreqrecsub$yr!= 2017,] #removing year with no data

szloop=data.frame(codsp=NULL,yr=NULL,lmids=NULL,n=NULL)

years=unique(szfreqpredsub$yr[szfreqpredsub$codsp==sp])

for (i in years) {
  h=hist(szfreqpredsub$fl[szfreqpredsub$yr==i],
         breaks=seq(from=min(szfreqpredsub$fl[szfreqpredsub$codsp==sp])-5,
                    to=max(szfreqpredsub$fl[szfreqpredsub$codsp==sp])+60, by=5), #binwidth
         xlab = 'Fork length (cm)',
         main = paste(j))
  freq=data.frame(codsp=sp,yr=i,lmids=round(h$mids,1),n=h$counts)    
  szloop=rbind(szloop,freq)
}

#LBSPR lengths format
szfreqpredsub =szloop %>% #removing first row 
  tidyr::spread(., yr, n) %>%  #years in columns
  `rownames<-`(.[,2]) %>%   #years in rownames
  dplyr::select(-codsp,-lmids) 

#observed length frequencies object
mylengthspred <- new("LB_lengths",LB_pars=mypars,dataType="freq",verbose=F) #blank length object
mylengthspred@LMids<- as.numeric(rownames(szfreqpredsub)) #mid points
mylengthspred@LData<- as.matrix(szfreqpredsub)    #frequency table
mylengthspred@L_units<- "cm"  #units
mylengthspred@Years<- as.numeric(colnames(szfreqpredsub)) #list of years
mylengthspred@NYears<- length(as.numeric(colnames(szfreqpredsub))) #number of years

#plotSize(mylengthspred)

fitlbspr= LBSPRfit(mypars, mylengthspred,verbose = FALSE)
#plotSize(fitlbspr)
#plotEsts(fitlbspr)


lbsprpred= data.frame(yr=fitlbspr@Years,  #years 
                      codsp=sp,
                      rawsl50=fitlbspr@SL50, #raw sl50
                      sl50sd=sqrt(fitlbspr@Vars[,1]),  #sl50 sd
                      sl50= fitlbspr@Ests[,1],   #smooth sl50
                      rawsl95=fitlbspr@SL95,
                      sl95sd=sqrt(fitlbspr@Vars[,2]),
                      sl95= fitlbspr@Ests[,2],
                      rawfm=fitlbspr@FM,
                      fmsd=sqrt(fitlbspr@Vars[,3]),
                      fm= fitlbspr@Ests[,3],
                      rawspr=fitlbspr@SPR,
                      sprsd=sqrt(fitlbspr@Vars[,4]),
                      spr= fitlbspr@Ests[,4],
                      source="Pred") #Predicted length frequency


#binding spr results from observed and reconstructed lengths (data frame)
lbspr= rbind(lbsprobs,lbsprpred)
lbspr_wah= lbspr #general object for plotting


#------------------------------------------------------
#length frequency reconstructed vs LBSPR Pop Simulation 
#------------------------------------------------------
szfreqsim=data.frame(yr=NULL,codsp=NULL,lmids=NULL,nsim=NULL,npred=NULL,res=NULL,rmse=NULL,
                     smd=NULL,ks=NULL,t=NULL,x2=NULL) #empty frame for simulated frequencies

#getting subset for current specie of predicted sizes 
szfreqpredsub=szfreqpred[szfreqpred$codsp==sp,]


for (j in 1:length(sp)) {  #loop through all species codes
  
  years=1950:2022
  
  for (i in 1:length(years)) { #loop through each year of available information
    
    
    #------------------------------------------------------
    #biological parameters for LBSPR simulation
    mypars <- new("LB_pars",verbose=F) #blank parameters object for simulation
    mypars@Species <- sp
    mypars@Linf <- pars$linf[pars$codsp==sp]
    mypars@L50 <-  pars$l50[pars$codsp==sp] 
    mypars@L95 <- 1.1*pars$l50[pars$codsp==sp] 
    mypars@MK <- pars$m[pars$codsp==sp]/pars$k[pars$codsp==sp] #
    mypars@L_units <- "cm"
    
    #Get selectivity and mortality parameters (SL50,SL95, FM) from LBSPR estimation
    mypars@SL50 <- lbsprpred$rawsl50[lbsprpred$yr==years[i]] 
    mypars@SL95 <- lbsprpred$rawsl95[lbsprpred$yr==years[i]]
    mypars@FM <-   lbsprpred$rawfm[lbsprpred$yr==years[i]]
    mypars@BinMax<- 1.4*pars$linf[pars$codsp==sp]
    mypars@BinMin<-  0
    mypars@BinWidth <- 5
    
    #LBSPR simulation
    mysim <- LBSPRsim(mypars)
    
    #normalize (min-max) function to make frequencies comparable
    normalize <- function(x, na.rm = TRUE) {
      return((x- min(x)) /(max(x)-min(x)))
    }
    
    
    lsim= data.frame(yr=rep(years[i], length(mysim@pLPop[,1])),
                     codsp=rep(sp[j],length(mysim@pLPop[,1])),
                     lmids=mysim@pLPop[,1],  
                     nsim=normalize(mysim@pLPop[,5]), #surviving in terms of probability
                     npred=NA)
    
    #-----------------------------------------------------
    # transforming predicted frequencies in length classes
    #we need to cut intervals to match with the simulated sizes
    #creating intervals between min and max length bins 
    lpred= data.frame(fl=szfreqpredsub$fl[szfreqpredsub$codsp==sp  
                                          & szfreqpredsub$yr==years[i]])
    freqtable <- data.frame(table(cut(lpred$fl,
                                      seq(min(mysim@pLPop[,1]),max(mysim@pLPop[,1])+5,5))))
    
    intervals=as.character(freqtable$Var1) #character vector
    intervals=gsub("[()]", "", intervals)       #remove parenthesis
    intervals=sapply(strsplit(intervals, split=',', fixed=TRUE), `[`, 1) #get the first class
    freqtable$Var1=intervals
    
    #normalized predicted frequency
    lsim$npred= normalize(freqtable$Freq)
    #------------------------------------------------------
    
    #Residual length predicted and simulated frequencies
    lsim$res= lsim$nsim-lsim$npred
    lsim$res[lsim$npred==0]=NA #avoiding Residual length in classes with no data
    
    #Root Mean Squared Error- RMSE
    rmse= function(x, y, na.rm = TRUE) {
      return(sqrt(sum((x- y)^2) /length(x)))
    }
    
    lsim$rmse= rmse(lsim$npred,lsim$nsim)
    
    #standardized mean reserence (SMD)
    smd= function(x, y, na.rm = TRUE) {
      return( abs(mean(x)- mean(y)) / sqrt((sd(x)^2+sd(y)^2)/2) )
    }
    
    lsim$smd= smd(lsim$npred,lsim$nsim)
    
    
    #Kolmogorov-smirnov test (comparison of cumulative distributions)
    #we need to convert probabilities in sample distributions to perform KS test
    lpredfreq=lpred$fl #predicted length distribution
    
    #simulated frequency distribution 
    lsimfreq= data.frame(lmids=mysim@pLPop[,1],n=round(mysim@pLPop[,5]*100))
    lsimfreq= lsimfreq %>%
      tidyr::uncount(n)
    lsimfreq=lsimfreq$lmids
    
    # Kolmogorov smirnov test to compare if the samples come from the same probability distributions
    lsim$ks=suppressWarnings(ks.test(lpredfreq,lsimfreq)$p.value) #kolmogorov-smirnov test
    
    #t-student test to compare if the means of distributions are the same
    lsim$t=suppressWarnings(t.test(lpredfreq,lsimfreq)$p.value)
    
    # chi-square goodness-of-fit test
    lsim$x2=suppressWarnings(chisq.test(sample(lpredfreq, 94,replace = FALSE), sample(lsimfreq, 94,replace = FALSE))$p.value)
    
    #bind into size frequency simulate data frame
    szfreqsim=rbind(szfreqsim,lsim)
  }
}
#get comparison of frequencies data frame 
szfreqsim_wah= szfreqsim
#--------------------------------------


#-------------
#plot section
#-------------
#plot GAM mean size prediction
p50<-ggplot(data=szmeanpred,aes(x=yr,y=fit))+
  geom_ribbon(data=dplyr::filter(szmeanpred,codsp==sp),
              aes(ymin = lw, ymax = up,linetype='Pred'), 
              col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(szmeanpred,codsp==sp),
            aes(x=yr,y=fit,linetype='Pred'),col="gray60",size=2)+
  #geom_line(data=dplyr::filter(gam,codsp==sp),aes(x=yr,y=fit,col=codgr,linetype=region),size=1.2)+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=yr,y=mean),size=2,col="gray40")+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='')+
  scale_y_continuous() +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p50

ggsave(paste(c("Length_model_prediction_Small_Tunas_",sp,".png"),collapse=''),plot=p50, device = "png", units = "cm",
       width = 24, height = 14)

# plot GAM partial effects
p51<-ggplot(dplyr::filter(pareffect,codsp==sp),
            aes(x=effect,y=fit,col=var),size=2) +
  #geom_line(aes(x=effect,y=up, col="Gear",group=1)) +
  geom_point(aes(x=effect,y=fit,col=var),size=3.5)+
  geom_errorbar(aes(ymin=lw, ymax=up,col=var), width=0.5,
                position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Effect",x="Factor",col='',fill="")+
  scale_y_continuous() +
  scale_x_discrete() +
  scale_colour_manual(values=c("gray60","gray30")) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p51
 
ggsave(paste(c("Effects_model_Small_Tunas_",sp,".png"),collapse=''),plot= p51, device = "png", units = "cm",
       width = 24, height = 14)


#plot GAM boxplot size freq prediction
p52<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqpred,codsp==sp), 
               aes(x=factor(yr), y=fl),
               col="gray40",
               fill="gray60",
               alpha=0.6,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=factor(yr),y=mean),size=2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Fork length (cm)",col="Region",linetype='',fill='Gear')+
  scale_y_continuous() +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p52

ggsave(paste(c("Length_pred_boxplot_Small_Tunas_",sp,".png"),collapse = ''),plot=p52, device = "png", units = "cm",
       width = 24, height = 15)


#plot half violin for the predicted frequencies
p53<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(szfreqpred,codsp==sp), 
                             aes(x=factor(yr),y=fl),
                             fill="gray20",
                             col="gray20",
                             alpha=0.7,
                             side = "r")+ 
  #show.legend = TRUE,
  #position = "identity",
  #fill=,
  #nudge = 0,
  #scale = 'area',
  #na.rm = TRUE,
  #lwd=1,
  #trim = TRUE)+
  #facet_wrap(~codsp,scales= "free_y")+
  #draw_quantiles = c(0.5))+
  geom_point(data=dplyr::filter(szmeanobs,codsp==sp),
             aes(x=factor(yr),y=mean),size=2)+
  # #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p53

ggsave(paste(c("Length_pred_density_Small_Tunas_",sp,".png"),collapse = ''),plot=p53, device = "png", units = "cm",
       width = 35, height = 12)


#plot SPR trajectories comparison (observed vs Predicted)
p54<-ggplot(dplyr::filter(lbspr,codsp==sp), 
            aes(x=yr,y=rawspr,fill=source,col=source),size=2) + 
  #geom_line(size=1,aes(fill=source)) +
  geom_point(size=2)+
  geom_smooth(size=1.2,method = 'loess' ,formula = 'y ~ x',se =FALSE,span = 0.4)+
  geom_errorbar(aes(ymin=ifelse(rawspr-sprsd < 0, 0, rawspr-sprsd),
                    ymax=ifelse(rawspr+sprsd > 1, 1, rawspr+sprsd)), width=1.5,
                position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Spawning Potential Ratio (SPR)",x="Year",col='',fill="")+
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks = seq(1950, 2023, 30),limits = c(1950, 2023)) +
  scale_colour_manual(values = c("gray60","gray40")) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p54

ggsave(paste(c("SPR_Small_Tunas_",sp,".png"),collapse=''),plot=p54, device = "png", units = "cm",
       width = 24, height = 14)



#plot GAM boxplot size freq prediction
p55<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqsim,codsp==sp), 
               aes(group=as.character(yr),x=as.character(yr), y=res),
               col="gray40",
               fill="gray60",
               alpha=0.4,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_smooth(data=dplyr::filter(szfreqsim,codsp==sp),
              method = "loess",span=0.1,col="gray20", se = FALSE,
              aes(group=1,x=as.character(yr), y=res,linetype="LOESS"),size=1.5)+
  geom_text(data=dplyr::filter(data.frame(szfreqsim),codsp==sp), 
            aes(x="2014", y=0.9*max(res,na.rm = TRUE), 
                label=paste("RMSE=",round(mean(rmse,na.rm = TRUE),2))))+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Length residual",col="",linetype='',fill='')+
  scale_y_continuous() +
  scale_fill_manual(values = c("blue", "red"))+
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p55

ggsave(paste(c("Length_res_boxplot_Small_Tunas_",sp,".png"),collapse = ''),plot=p55, device = "png", units = "cm",
       width = 24, height = 15)
#x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x#


#Combining final data frames

#write the size means observed (programs and literature)
#Write Final data frames in csv files
outfile0  <- "szmeanobs.csv"
write.table(szmeanobs, file = outfile0, append =FALSE,dec=".",sep = ",",
            row.names = FALSE) 


#gam models + mean size predicted for all factors combined data frame
szmeanpred=rbind(szmeanpred_blf,szmeanpred_brs,szmeanpred_dol,
                 szmeanpred_fri,szmeanpred_kgm,szmeanpred_lta,
                 szmeanpred_wah)

#Write Final data frames in csv files
outfile1  <- "szmeanpred.csv"
write.table(szmeanpred, file = outfile1, append =FALSE,dec=".",sep = ",",
            row.names = FALSE) 

#partial effects of the models for all species
pareffect=rbind(pareffect_blf,pareffect_brs,pareffect_dol,
                pareffect_fri,pareffect_kgm,#pareffect_lta, #LTA does not have partial effects
                pareffect_wah)

#Write Final data frames in csv files
outfile2  <- "pareffect.csv"
write.table(pareffect, file = outfile2, append =FALSE,dec=".",sep = ",",
            row.names = FALSE) 


#size frequency for all species combined
szfreqpred=rbind(szfreqpred_blf,szfreqpred_brs,szfreqpred_dol,
                 szfreqpred_fri,szfreqpred_kgm,szfreqpred_lta,
                 szfreqpred_wah)

#Write Final data frames in csv files
outfile3  <- "szfreqpred.csv"
write.table(szfreqpred, file = outfile3, append =FALSE,dec=".",sep = ",",
            row.names = FALSE) 


#LBSPR metrics for all species data frame
lbspr=rbind(lbspr_blf,lbspr_brs,lbspr_dol,
            lbspr_fri,lbspr_kgm,lbspr_lta,
            lbspr_wah)

#Write Final data frames in csv files
outfile4  <- "lbspr.csv"
write.table(lbspr, file = outfile4, append =FALSE,dec=".",sep = ",",
            row.names = FALSE) 

#Predicted and simulated frequencies comparison
szfreqsim=rbind(szfreqsim_blf,szfreqsim_brs,szfreqsim_dol,
                szfreqsim_fri,szfreqsim_kgm,szfreqsim_lta,
                szfreqsim_wah)

#Write Final data frames in csv files
outfile5  <- "szfreqsim.csv"
write.table(szfreqsim, file = outfile5, append =FALSE,dec=".",sep = ",",
            row.names = FALSE) 


#combine predicted frequency size with observed size
szfreqpred$codgr= "UN"
szfreqpred$region= "ALL"
sizefreqobs_pred=rbind(szfreqobs,szfreqpred)

#Write Final data frames in csv files
outfile6  <- "sizefreqobs_pred.csv"
write.table(sizefreqobs_pred, file = outfile6, append =FALSE,dec=".",sep = ",",
            row.names = FALSE) 

#----------------------------------------
# Final Plots with all species 
#----------------------------------------

p56<-ggplot(data=dplyr::filter(szmeanpred,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                       "KGM","LTA","WAH")),aes(x=yr,y=fit))+
  geom_ribbon(data=dplyr::filter(szmeanpred,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                              "KGM","LTA","WAH")),
              aes(ymin = lw, ymax = up,linetype='Pred'), 
              col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(szmeanpred,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                            "KGM","LTA","WAH")),
            aes(x=yr,y=fit,linetype='Pred'),col="gray60",size=2)+
  geom_point(data=dplyr::filter(szmeanobs,codsp %in% c("BLF",'BRS', "DOL","FRI",
                                                             "KGM","LTA","WAH")),
             aes(x=yr,y=mean,shape="Obs"),size=2,col="gray40")+
  facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='',shape='')+
  scale_y_continuous() +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p56

ggsave("Length_model_prediction_Small_Tunas_ALLSP.png",plot=p56, device = "png", units = "cm",
       width = 30, height = 18)

# plot GAM partial effects
# Define the levels for "Gear" and "Region"
gear_levels <- c("GN", "HL", "LL", "RR", "TR", "UN", "BB", "PS")
region_levels <- c("N", "NE", "S", "SE", "ALL")

# Combine the levels for "Gear" and "Region"
combined_levels <- c(gear_levels, region_levels)

# Reorder the levels of the `effect` column based on the combined_levels vector
pareffect$effect <- factor(pareffect$effect, levels = combined_levels, ordered = TRUE)

# Ensure that the `var` variable is correct
pareffect$var <- factor(pareffect$var, levels = c("Gear", "Region"))


p57<-ggplot(dplyr::filter(pareffect,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                  "KGM","LTA","WAH")),
            aes(x=effect,y=fit,col=factor(var)),size=2) +
  #geom_line(aes(x=effect,y=up, col="Gear",group=1)) +
  geom_point(aes(x=effect,y=fit,col=var),size=2.5)+
  geom_errorbar(aes(ymin=lw, ymax=up,col=var), width=0.5,
                position="identity",size=1.2)+
  facet_wrap(~codsp,scales = "free_y")+
  labs(y="Effect",x="Factor level",col='',fill="")+
  scale_y_continuous() +
  scale_x_discrete() +
  scale_colour_manual(values=c("gray60","gray40")) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p57

ggsave("Effects_model_Small_Tunas_ALLSP.png",plot= p57, device = "png", units = "cm",
       width = 32, height = 10)


#plot GAM boxplot size freq prediction
p58<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqpred,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                         "KGM","LTA","WAH")), 
               aes(x=factor(yr), y=fl,fill="Sim"),
               #col="gray40",
               #fill="gray60",
               alpha=0.6,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_point(data=dplyr::filter(szmeanobs,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                  "KGM","LTA","WAH")),
             aes(x=factor(yr),y=mean,col="Obs"),size=2)+
  facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Fork length (cm)",col="",linetype='',fill='')+
  scale_fill_manual(values = c("gray70"))+
  scale_color_manual(values = c("gray10"))+
  scale_y_continuous() +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p58

ggsave("Length_pred_boxplot_Small_Tunas_ALLSP.png",plot=p58,device="png",units="cm",
       width = 30, height = 18)


#plot GAM boxplot size prediction + model prediction

#removing extreme values
szmeanobssub=szmeanobs
szmeanobssub$mean[szmeanobssub$codsp=='BLF' & szmeanobssub$mean<50] =NA
szmeanobssub$mean[szmeanobssub$codsp=='DOL' & szmeanobssub$mean<75] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean>100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='KGM' & szmeanobssub$mean<60] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean<100] =NA
szmeanobssub$mean[szmeanobssub$codsp=='WAH' & szmeanobssub$mean>180] =NA

p59<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqpred,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                          "KGM","LTA","WAH")), 
               aes(x=factor(yr), y=fl,fill="Sim"),
               #col="gray40",
               #fill="gray60",
               alpha=0.6,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_ribbon(data=dplyr::filter(szmeanpred,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                         "KGM","LTA","WAH")),
              aes(x=factor(yr),y=fit,ymin = lw, ymax = up,group=1,linetype="Pred"),
              color="gray80",fill="gray30",alpha=0.4)+
  geom_line(data=dplyr::filter(szmeanpred,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                       "KGM","LTA","WAH")),
            aes(x=factor(yr),y=fit,group=1,linetype="Pred"),size=2,col="gray20")+
  geom_point(data=dplyr::filter(szmeanobssub,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                       "KGM","LTA","WAH")),
             aes(x=factor(yr),y=mean,group=1,shape="Obs"),size=2)+
  facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Fork length (cm)",linetype="",shape="",col="",fill="")+
  scale_color_manual(values = c("gray40",col="gray10"))+
  scale_fill_manual(values = c("gray60"))+
  scale_y_continuous() +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p59

ggsave("Length_pred_boxplot_Model_prediction_Small_Tunas_ALLSP.png",plot=p59,device="png",units="cm",
       width = 30, height = 18)


#plot half violin for the predicted frequencies
p60<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(szfreqpred,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                                        "KGM","LTA","WAH")), 
                             aes(x=factor(yr),y=fl),
                             fill="gray20",
                             col="gray20",
                             alpha=0.7,
                             side = "r")+ 
  #show.legend = TRUE,
  #position = "identity",
  #fill=,
  #nudge = 0,
  #scale = 'area',
  #na.rm = TRUE,
  #lwd=1,
  #trim = TRUE)+
  facet_wrap(~codsp,scales= "free_y")+
  #draw_quantiles = c(0.5))+
  geom_point(data=dplyr::filter(szmeanobs,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                        "KGM","LTA","WAH")),
             aes(x=factor(yr),y=mean),size=2)+
  # #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p60

ggsave("Length_pred_density_Small_Tunas_ALLSP.png",plot=p60, device = "png", units = "cm",
       width = 35, height = 12)


#Plot residual vs fitted values of GAM models
resvsfit= data.frame(res= c(modgam_blf$residuals,modgam_brs$residuals,modgam_dol$residuals,
                     modgam_fri$residuals,modgam_kgm$residuals,modgam_lta$residuals,
                     modgam_wah$residuals),
                     fit= c(modgam_blf$fitted.values,modgam_brs$fitted.values,modgam_dol$fitted.values,
                            modgam_fri$fitted.values,modgam_kgm$fitted.values,modgam_lta$fitted.values,
                            modgam_wah$fitted.values),
                     codsp=c(rep('BLF',length(modgam_blf$residuals)),
                             rep('BRS',length(modgam_brs$residuals)),
                             rep('DOL',length(modgam_dol$residuals)),
                             rep('FRI',length(modgam_fri$residuals)),
                             rep('KGM',length(modgam_kgm$residuals)),
                             rep('LTA',length(modgam_lta$residuals)),
                             rep('WAH',length(modgam_wah$residuals))))

p61<-ggplot(dplyr::filter(resvsfit,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                             "KGM","LTA","WAH")), 
            aes(x=fit,y=res),size=1.2) + 
  #geom_line(size=1,aes(fill=source)) +
  geom_point(size=1.5)+
  geom_smooth(size=1.2,method = 'loess' ,formula = 'y ~ x',se =FALSE,span = 2,col="gray60")+
  facet_wrap(~codsp,scales = "free")+
  labs(y="Residuals",x="Fitted values",col='',fill="")+
  geom_hline(yintercept=0, linetype=2)+
  #scale_x_continuous(breaks = seq(1950, 2023, 30),limits = c(1950, 2023)) +
  #scale_colour_manual(values = c("gray60","gray40")) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.5,0.05,0.05),"mm"))
p61
  
ggsave("Residuals_vs_Fitted_Small_Tunas_ALLSP.png",plot=p61, device = "png", units = "cm",
       width = 30, height = 18)



#plot SPR trajectories comparison (observed vs Predicted)
p62<-ggplot(dplyr::filter(lbspr,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                 "KGM","LTA","WAH")), 
            aes(x=yr,y=rawspr,fill=source,col=source),size=1.2) + 
  #geom_line(size=1,aes(fill=source)) +
  geom_point(size=1.5)+
  geom_smooth(size=1.2,method = 'loess' ,formula = 'y ~ x',se =FALSE,span = 0.4)+
  geom_errorbar(aes(ymin=ifelse(rawspr-sprsd < 0, 0, rawspr-sprsd),
                    ymax=ifelse(rawspr+sprsd > 1, 1, rawspr+sprsd)), width=1.5,
                position="identity",size=1.2)+
  facet_wrap(~codsp,scales = "free_y")+
  labs(y="Spawning Potential Ratio (SPR)",x="Year",col='',fill="")+
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks = seq(1950, 2023, 30),limits = c(1950, 2023)) +
  scale_colour_manual(values = c("gray60","gray40")) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p62

ggsave("SPR_Small_Tunas_ALLSP.png",plot=p62, device = "png", units = "cm",
       width = 30, height = 18)



#plot GAM boxplot size freq prediction
rmsetext <- szfreqsim %>% 
  dplyr::group_by(codsp) %>%
  dplyr::summarize(rmse = round(mean(rmse), 2)) %>%
  dplyr::mutate(text = paste("RMSE=", rmse))

p63<- ggplot() +
  geom_boxplot(data=dplyr::filter(szfreqsim,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                     "KGM","LTA","WAH")), 
               aes(group=as.character(yr),x=as.character(yr), y=res),
               col="gray40",
               fill="gray60",
               alpha=0.4,
               #show.legend = TRUE)+
               position = "identity",
               lwd=0.75,
               outlier.shape=8,
               outlier.size = 0.8,
               outlier.alpha = 0.5)+
  geom_smooth(data=dplyr::filter(szfreqsim,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                   "KGM","LTA","WAH")),
              method = "loess",span=0.2,col="gray20", se = FALSE,
              aes(group=1,x=as.character(yr), y=res),size=1.5)+
  geom_text(data=rmsetext, 
            aes(x="2012", y=0.98*max(szfreqsim$res,na.rm = TRUE), 
                          label=text))+
  facet_wrap(~codsp,scales = "free_y")+
  labs(x="Year",y="Length residual",col="",linetype='',fill='')+
  scale_color_manual(values = c("gray20"))+
  scale_fill_manual(values = c("gray30"))+
  scale_y_continuous() +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2022,10))))+
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,1.3,0.05,0.05),"mm"))
p63

ggsave("Length_res_boxplot_Small_Tunas_ALLSP.png",plot=p63, device = "png", units = "cm",
       width = 30, height = 18)

#x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x#
# end of reconstructions
#x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x#


