#####################################
# Small Tunas Length Reconstruction
# by: Matheus Lourenco (Silva, M.L.S)
#####################################

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

#install.packages("TropFishR")
library(TropFishR) #catch-curve estimates

# the development version from github LBSPR
#install.packages("devtools")
#devtools::install_github("AdrianHordyk/LBSPR",force = TRUE)
library(LBSPR)

#-------------  Defining some important functions ----------------#



#------ function to create length classes efficiently-----#
lf<- function(x, by=5){ 
  
  h = hist(x, breaks = seq(from = min(x) - by, to = max(x) + by, by = by), plot = FALSE)
  lf_dat<- data.frame(mids= h$mids,counts=h$counts)
  
  return(lf_dat)
}

#-- (min-max) function to make frequencies comparable--#
norm <- function(x, na.rm = TRUE) {
  
  return((x- min(x)) /(max(x)-min(x)))
}

#---- Root Mean Squared Error- RMSE function -----#
rmse= function(x, y, na.rm = TRUE) {
  
  rmse=sqrt(sum((x- y)^2) /length(x))
  
  return(rmse)
}


#----- Relative fishing mortality function ------ #
fm <- function(n, x, by, linf, k, l50, l95, mk, year, lmean, lstart, method, sp, units) {
  
  if (method == "Beverton and Holt") {
    zk = (linf - lmean) / (lmean - lstart)
    fm = (zk / mk) - 1  # ZK = MK * (1 + FM)
    
    # verify if viable estimates are obtained
    if (zk < mk | zk < 0 | zk>10) {
      message("Impossible to estimate with Beverton and Holt, trying LBSPR")
        method <- "LBSPR"
    } else {
      return(fm)
    }
  }
  
  if (method == "LBSPR") {
    require(LBSPR)
    
    # Biological parameters
    mypars <- new("LB_pars", verbose = F)
    mypars@Species <- sp
    mypars@Linf <- linf
    mypars@L50 <- l50
    mypars@L95 <- l95
    mypars@MK <- mk
    mypars@L_units <- units
    
    # length data
    h <- lf(x, by = by)
    maxbin <- max(h$mids)
    
    # Avoid maxbin < linf
    if (maxbin < linf) {
      bin_width <- diff(h$mids)[1]  
      new_bins <- seq(maxbin + bin_width, linf + bin_width, by = bin_width)
      
      mids <- c(h$mids, new_bins)
      counts <- c(h$counts, rep(0, length(new_bins)))
      h <- data.frame(mids = mids, counts = counts)
      
      # Ordering
      ordering <- order(h$mids)
      h$mids <- h$mids[ordering]
      h$counts <- h$counts[ordering]
    }
    
    # length data
    mylengthsobs <- new("LB_lengths", LB_pars = mypars, dataType = "freq", verbose = F)
    mylengthsobs@LMids <- as.numeric(h$mids)
    mylengthsobs@LData <- as.matrix(h$counts)
    mylengthsobs@L_units <- units
    mylengthsobs@Years <- as.numeric(year)
    mylengthsobs@NYears <- length(as.numeric(year))
    
    # Fit data 
    fit <- LBSPRfit(mypars, mylengthsobs, verbose = FALSE)
    
    # Exploitation
    fm <- fit@FM
    return(fm)
  }
}
    
    
    
# Selectivity from the selection ogive - catch curve (TropfishR)
# Selectivity from the LBSPR fit function
sel <- function(x, linf, k, l50, l95, mk, t0, by_values, year, sp, units, method) {
  
  if (method == "Catch curve") {
    
    require(TropFishR)
    
    for (by in by_values) {
      
      h <- tryCatch({
        lf(x, by)
      }, error = function(e) {
        return(NULL)
      })
      
      if (is.null(h)) {
        next # if is not correct, test the next by
      }
      
      list <- list(
        midLengths = round(h$mids, 1),
        Linf = linf,
        K = k,
        t0 = t0,
        catch = h$counts
      )
      
      sl <- tryCatch({
        TropFishR::catchCurve(list, auto = TRUE, calc_ogive = TRUE)
      }, error = function(e) {
        return(NULL)
      })
      
      if (!is.null(sl)) {
        sl50 <- sl$L50
        
        # Test if selectivity is outside the limits
        if (sl50 < h$mids[1] || sl50 > h$mids[length(h$mids)]) {
          message("Impossible to estimate with catch-curve, trying LBSPR.")
          return(sel(x, linf, k, l50, l95, mk, t0, by_values, year, sp, units, method = "LBSPR"))
        }
        
        return(sl50) # Return the best sl50
      }
    }
    
    # If the model fails for all by_values, add extra lengths
    max_length <- max(x) 
    extra_lengths <- seq(max_length + 2, max_length + 10, by = 2) 
    x_extended <- c(x, extra_lengths)
    
    h <- lf(x_extended, by_values[1])
    
    list <- list(
      midLengths = round(h$mids, 1),
      Linf = linf,
      K = k,
      t0 = t0,
      catch = h$counts
    )
    
    sl <- tryCatch({
      TropFishR::catchCurve(list, auto = TRUE, calc_ogive = TRUE)
    }, error = function(e) {
      message("Unreasonable values of by, even after extending the length range. Estimating selectivity via LBSPR.")
      return(sel(x, linf, k, l50, l95, mk, t0, by_values, year, sp, units, method = "LBSPR"))
    })
    
    return(sl)
  }
  
  if (method == "LBSPR") {
    
    require(LBSPR)
    
    # Biological parameters
    mypars <- new("LB_pars", verbose = F)
    mypars@Species <- sp
    mypars@Linf <- linf
    mypars@L50 <- l50
    mypars@L95 <- l95
    mypars@MK <- mk
    mypars@L_units <- units
    
    # Length data
    h <- lf(x, by = by_values[1])
    maxbin <- max(h$mids)
    
    # Avoiding maxbin < linf
    if (maxbin < linf) {
      
      bin_width <- diff(h$mids)[1]  # Calculate the bin width
      new_bins <- seq(maxbin + bin_width, linf + bin_width, by = bin_width)
      
      mids <- c(h$mids, new_bins)
      counts <- c(h$counts, rep(0, length(new_bins)))
      h <- data.frame(mids = mids, counts = counts)
      
      # Ordering
      ordering <- order(h$mids)
      h$mids <- h$mids[ordering]
      h$counts <- h$counts[ordering]
    }
    
    # Length data for LBSPR
    mylengthsobs <- new("LB_lengths", LB_pars = mypars, dataType = "freq", verbose = F)
    mylengthsobs@LMids <- as.numeric(h$mids)
    mylengthsobs@LData <- as.matrix(h$counts) 
    mylengthsobs@L_units <- units  
    mylengthsobs@Years <- as.numeric(year)  
    mylengthsobs@NYears <- length(as.numeric(year)) 
    
    # Fit the model
    fit <- LBSPRfit(mypars, mylengthsobs, verbose = FALSE)
    
    # Selectivity parameters
    sl50 <- fit@SL50
    
    # verify if sl50 is above the linf
    if (sl50 < h$mids[1] || sl50 > h$mids[length(h$mids)] || sl50 > linf) {
      message("sl50 outside the mids range or greater than Linf. Adjusting sl50 to Linf - by.")
      sl50 <- linf - by_values[1]
    }
    
    return(sl50)
  }
}

#--------- Fit distributions function -----------------#
distfit<- function(len) { #len(yr,codsp,codgr,fl)
  require(fitdistrplus)
  
  dist=data.frame(yr=NULL,codsp=NULL,dist=NULL,codgr=NULL, 
    mean=NULL,shape=NULL,cv=NULL,aic=NULL) 
  
  sp= unique(len$codsp)
  years=unique(years)
  
  for (i in sp) {
    
    for (j in years) {
      
      #frequency for each year  
      sz= len$fl[len$codsp==i & len$yr==j]
      gear= unique(len$codgr[len$codsp==i & len$yr==j])[1]
      
      if (length(sz)>5 & (length(unique(sz))==1)==FALSE) {
        
        #fitting distributions
        fit_1 <- fitdistrplus::fitdist(sz, "norm")
        fit_2<- fitdistrplus::fitdist(sz, "invgauss",start=list(mean=mean(sz),shape= 1))
        fit_3 <- fitdistrplus::fitdist(sz, "gamma")
        
        dat=data.frame(
          yr=j,
          codsp=i,
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
  return(dist) 
}  

#----------- Distribution Simulation Function --------------#

distsim<- function(mean, cv, years, n, sp, distribution, source) {
  require(actuar)
  
  flloop=data.frame(yr=NULL,codsp=NULL,fl=NULL,source=NULL) 
  
  for (j in 1:length(sp)) { 
    
    years=years
    n=n 
    
    for (i in 1:length(years)) { 
      
      l= data.frame(
        yr=rep(years[i],n),
        codsp=rep(sp[j],n),
        
        fl=if(distribution=="norm"){rnorm(n=n, mean = mean[i],sd = mean[i]*cv) 
        } else if(distribution=="invgauss"){actuar::rinvgauss(n=n, mean = mean[i],shape= mean[i]*cv) 
        } else if(distribution=="gamma") {rgamma(n=n, shape=mean[i],scale =mean[i]*cv)},
        
        source=rep(as.character(paste(source)),n))
      
      flloop=rbind(flloop,l)
    }
  }
  return(flloop)
}


#-------- LBSPR simulation Function -----------#
lbsprsim<- function(linf,l50,l95,mk,units,sl50,sl95,fm,minbin,maxbin,binwidth,sp) {
  
    require(LBSPR)
    #biological parameters
    mypars <- new("LB_pars",verbose=F)
    mypars@Species <- sp
    mypars@Linf <- linf
    mypars@L50 <-  l50 
    mypars@L95 <-  l95 
    mypars@MK <-  mk 
    mypars@L_units <- units
  
    #Exploitation parameters
    mypars@SL50 <- sl50
    mypars@SL95 <- sl95
    mypars@FM <-   fm
    mypars@BinMin<- minbin
    mypars@BinMax<- maxbin
    mypars@BinWidth <-binwidth
  
    #LBSPR simulation
    mysim <-LBSPR:: LBSPRsim(mypars)
  
  return(mysim)
}




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
freq_obs=rbind(sz1sub,sz2sub)
freq_obs= freq_obs[order(freq_obs$yr),] #order by year

unique(freq_obs$codsp) #species list

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
mean_obs = freq_obs %>%
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
                     codsp= unique(mean_obs$codsp[mean_obs$codsp!=""]), 
                     codgr="UN",source='Literature',
                     mean=NA, 
                     sd=NA,cv=NA,min=NA,max=NA)

#Final observed mean size 
mean_obs= mean_obs %>%
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

p1<-ggplot2::ggplot(freq_obs, aes(x=fl, fill=codsp)) + 
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
sp<- unique(freq_obs$codsp)
#colours=heat.colors(length(sp),alpha = 0.5,rev = TRUE)
colours=rainbow(length(sp),alpha = 0.65,rev = FALSE)

#colours=c('tomato1','deepskyblue4','springgreen','firebrick',
#          'deepskyblue','black','dimgray', 'gold')

d=list(BLF= density(freq_obs$fl[freq_obs$codsp=='BLF']),
     FRI= density(freq_obs$fl[freq_obs$codsp=='FRI']),
     WAH= density(freq_obs$fl[freq_obs$codsp=='WAH']),
     BRS= density(freq_obs$fl[freq_obs$codsp=='BRS']),
     KGM= density(freq_obs$fl[freq_obs$codsp=='KGM']),
     CER= density(freq_obs$fl[freq_obs$codsp=='CER']),
     LTA= density(freq_obs$fl[freq_obs$codsp=='LTA']),
     DOL= density(freq_obs$fl[freq_obs$codsp=='DOL']))


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
  gghalves::geom_half_violin(data=freq_obs, 
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
  geom_point(data=dplyr::filter(mean_obs,codsp!=""),
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
  gghalves::geom_half_violin(data=freq_obs, 
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
  geom_point(data=dplyr::filter(mean_obs,codsp!=""),
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
  geom_point(data=dplyr::filter(mean_obs,codsp!="" & codsp!="CER"),
             aes(x=factor(yr),y=mean,col=source),size=1.5)+
  geom_boxplot(data=dplyr::filter(freq_obs,codsp!="" & codsp!="CER"), 
    aes(x=factor(yr), y=fl,fill=source),
    #show.legend = TRUE)+
    #position = "identity",
    lwd=0.55,
    outlier.shape=8,
    outlier.size = 0.5,
    outlier.alpha = 0.5)+
  geom_point(data=dplyr::filter(mean_obs,codsp!="" & codsp!="CER"),
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

p6<-ggplot(dplyr::filter(mean_obs,codsp!=""), 
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
p7<-ggplot(freq_obs, 
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
p8<-ggplot(freq_obs, 
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
p9<-ggplot(freq_obs, 
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
p10<-ggplot(freq_obs, 
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
p11<-ggplot(dplyr::filter(mean_obs,codsp!=""), 
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
p12<-ggplot(dplyr::filter(mean_obs,codsp!=""), 
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
p13<-ggplot(dplyr::filter(mean_obs,codsp!=""), 
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
p14<-ggplot(dplyr::filter(mean_obs,codsp!=""), 
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
years=1950:2025
#-------------------------------------
#find the best distribution for specie
#-------------------------------------
dist<-distfit(len = data.frame(
  yr=freq_obs$yr[freq_obs$codsp==sp],
  codsp=freq_obs$codsp[freq_obs$codsp==sp],
  codgr=freq_obs$codgr[freq_obs$codsp==sp],
  fl=freq_obs$fl[freq_obs$codsp==sp]))

dist=dist %>%
  dplyr::group_by(yr,codsp,codgr) %>% 
  slice_min(aic, n = 1)  %>%
  dplyr::arrange(dist) 

#assign dist for the current species
dist_blf=dist

#Table with the best distribution
table(dist$codsp,dist$dist)  
distribution='gamma'

#coefficient of variation for the corresponding distribution 
cv=mean(dist$cv[dist$codsp==sp & dist$dist==distribution])

#-------------------------------------------------------
#predict 1950 length when using the slope of the lengths
#-------------------------------------------------------

#assuming the 1950 mean length as the prediction of a linear model
lm_data<- data.frame(length=mean_obs$mean[mean_obs$codsp==sp],year=mean_obs$yr[mean_obs$codsp==sp])
lm_data<-lm_data[lm_data$length>50,]

#fiting the model
lm_model<- lm(length~year, data = lm_data[-1,])
#Parameters... On average the mean length is being reduced 0.10 cm year by year
a<- coef(lm_model)[1]
b<- coef(lm_model)[2]

#using the fitted parameters to estimate what would be the length in 1950
newdata<- data.frame(year=1950:2025)
lm_predict_10<- as.data.frame(predict.lm(lm_model, newdata = newdata,
                                        se.fit = TRUE,interval="prediction",level = c(0.1)))
lm_predict_20<- as.data.frame(predict.lm(lm_model, newdata = newdata,
                              se.fit = TRUE,interval="prediction",level = c(0.2)))
lm_predict_30<- as.data.frame(predict.lm(lm_model, newdata = newdata,
                      se.fit = TRUE,interval="prediction",level = c(0.3)))

#data frame with predictions for different percentile levels
lm_predict<- data.frame(yr=newdata$year,
                        fit=lm_predict_10$fit.fit, 
                        lw_10=lm_predict_10$fit.lwr,
                        up_10=lm_predict_10$fit.upr,
                        lw_20=lm_predict_20$fit.lwr,
                        up_20=lm_predict_20$fit.upr,
                        lw_30=lm_predict_30$fit.lwr,
                        up_30=lm_predict_30$fit.upr,
                        codsp=sp)

#lm model coefficients and predicitions for the current species
lm_model_blf<- data.frame(formula=paste(deparse(lm_model$call), collapse = ""),
                          intercept=a,
                          slope=b,
                          codsp=sp)
lm_predict_blf<- lm_predict

# controlling variables
lm_predict <- lm_predict %>%
  mutate(confidence_level = case_when(
    yr == 1950 & lw_10 == lw_10 ~ "10%",
    yr == 1950 & lw_20 == lw_20 ~ "20%",
    yr == 1950 & lw_30 == lw_30 ~ "30%"
  ))

# visual representation of the linear model
#-------------
#plot section
#-------------
p15 <- ggplot(lm_predict, aes(x = yr, y = fit)) +
  geom_line(size = 1, color = "grey20") +  # main line prediction
  # Lines for the percentiles
  geom_line(aes(y = lw_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = up_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = lw_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = up_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = lw_30, color = "30%", linetype = "30%"), size = 0.8) +
  geom_line(aes(y = up_30, color = "30%", linetype = "30%"), size = 0.8) +
    # Observed data
  geom_point(data = lm_data, aes(x = year, y = length), size = 2, color = "grey30") +
  # Predicted points for 1950 
  geom_point(aes(x = yr[1], y = lw_10[1], fill = "10%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_10[1], fill = "10%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = lw_20[1], fill = "20%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_20[1], fill = "20%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = lw_30[1], fill = "30%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_30[1], fill = "30%"), size = 3, shape = 21) +
  labs(y = 'Mean Length (cm)', x = "Year", color = "Percentile", fill = "Percentile", linetype = "Percentile") +
  # Type and line colors
  scale_color_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_fill_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_linetype_manual(values = c("10%" = "dashed", "20%" = "dotted", "30%" = "dotdash")) +
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p15
#saving graph...
ggsave(paste0("Linear_model_prediction_Small_Tunas_", sp, ".png"), plot = p15, device = "png", units = "cm", width = 24, height = 16)
#----------------------------------

#selecting variables .. subseting
mean_obs_sub= mean_obs[,c("yr",'region','codsp','codgr','mean')]

#removing extreme values
mean_obs_sub$mean[mean_obs_sub$codsp=='BLF' & mean_obs_sub$mean<50] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='DOL' & mean_obs_sub$mean<75] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='KGM' & mean_obs_sub$mean>100] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='KGM' & mean_obs_sub$mean<60] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean<100] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean>180] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='FRI' & mean_obs_sub$mean>45] =NA

#get only non-NA values and selecting the current species
mean_obs_sub=mean_obs_sub[is.na(mean_obs_sub$mean)==FALSE & mean_obs_sub$codsp==sp,]

#grouping factors most likely (Fredou et al., 2021)
mean_obs_sub$region[mean_obs_sub$region=='S']='SE'
mean_obs_sub$codgr[mean_obs_sub$codgr=='BB']='UN'

#transforming columns in to factors
mean_obs_sub$region= factor(mean_obs_sub$region)
mean_obs_sub$codgr= factor(mean_obs_sub$codgr)
mean_obs_sub$codsp= factor(mean_obs_sub$codsp)


# model GAM selection----------------
#list of values to be tested
tests<- names(lm_predict[, grepl("fit|lw|up", names(lm_predict))])
#filling an empty 1950 with NA 
mean_obs_sub<- rbind(data.frame(yr=1950, region="SE", codsp=sp,codgr="UN",mean=NA),
                  mean_obs_sub)

# Definindo os vetores de variáveis de resposta e explicativas
resp_vars <- c("mean")
exp_vars_1 <- c("s(yr)")
exp_vars_2 <- c("region")
exp_vars_3 <- c("codgr")

#criando a formula geral
formulas <-list(as.formula(paste(resp_vars, "~", paste(c(exp_vars_1), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_2), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_3), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_2,exp_vars_3), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_3,exp_vars_2), collapse = "+"))))

formulas<- as.character(formulas)

# Data frame com as distribuicoes de probabilidade e funcoes de ligacao a serem testadas
family=distribution #testing only for the best distribution
link = ("log") #log fixed

#empty frame to store results
gam=data.frame(
    test=NULL,
    type=NULL,
    codsp=NULL,
    family=NULL,
    link=NULL,
    formula=NULL,
    newdataexpand=NULL,
    pred=NULL,
    aic=NULL,
    r.sq=NULL,
    dev.ex=NULL)


     #-------------------- Fitting models  -------------------#
     for (i in tests) {
    
       #1- first sensibility for different 1950 starting values--
       test_1950<- lm_predict[1, grepl(i, names(lm_predict))]
       
       #replacing with the current 1950 mean length test
       mean_obs_sub$mean[mean_obs_sub$yr==1950]<- test_1950
    
       #2- Term selection within the sensibility----
       # Loop through the formulas 
       for (j in formulas) {
          
         family_name <- as.character(family) #these are fixed
         link_name <- as.character(link)
         
         #Creating switch object
         family_obj <- switch(family_name,
          norm = gaussian(link = link_name),  
          invgauss = inverse.gaussian(link = link_name),
          gamma= Gamma(link=link_name))
         
         #formula j to be tested
         formula <- as.formula(j)
          
         #fit the model
         gam_model <- mgcv::gam(as.formula(formula),
                                  data= mean_obs_sub,
                                  method = "REML" ,
                                  family =family_obj)
        #taking the outputs
       gam=rbind(gam, data.frame(
       test=i,
       type="GAM",
       codsp=sp,
       family=gam_model$family$family,
       link=gam_model$family$link,
       formula=paste0(c(paste(gam_model$formula)[2],
       paste(gam_model$formula)[1],
       paste(gam_model$formula)[3]),collapse = ' '),
       aic=gam_model$aic,
       r.sq=summary(gam_model)$r.sq,
      dev.ex=summary(gam_model)$dev.expl))
     }
    }#-------------------- end of loop -----------------------#

#assign the tests to the species
gam_blf<-gam

#replacing the 1950 mean length with the best tested value
best_model<- gam[which.min(gam$aic),]; dplyr::tibble(best_model)
mean_obs_sub$mean[mean_obs_sub$yr==1950]<- lm_predict[1, grepl(best_model$test, names(lm_predict))]

#Final model (lowest AIC) 
gam_model<- mgcv::gam(as.formula(best_model$formula), 
                   sp=0.01,
                   data= mean_obs_sub,
                   method = "REML" ,
                   family =family_obj)
#assign GAM model to the specie
gam_model_blf=gam_model

#diagnostic
summary(gam_model)

par(mfrow=c(2,2))
mgcv::gam.check(gam_model)# model check

#residual 
#homocedasticity
lmtest::bptest(gam_model)

#normality
shapiro.test(resid(gam_model))

# Check overall concurvity
round(mgcv::concurvity(gam_model, full = TRUE),2)
#----------------------
#predictions "year"
newdataexpand= expand.grid(yr=seq(years[1],years[length(years)],1),
                          region=mean_obs_sub$region,
                          codgr=mean_obs_sub$codgr)

pred<- data.frame(predict(gam_model, newdata = newdataexpand,
                          type = "response", se.fit=TRUE))

pred$up<-pred$fit+qnorm(0.975)*pred$se.fit #confidence intervals
pred$lw<-pred$fit-qnorm(0.975)*pred$se.fit

#take the predictions
mean_pred=cbind(type='GAM',
          codsp=sp,
          yr=newdataexpand$yr,
          mean=NA,
          newdata=newdataexpand,
          pred = round(pred$fit, 2),
          lw = round(pred$lw, 2),
          up = round(pred$up, 2))

#extract mean of all levels
mean_pred = mean_pred %>%
  group_by(type, yr, codsp) %>%
  summarize(across(c('pred', 'lw', 'up'), ~ round(mean(.), 2)), .groups = 'drop') %>%
  #binding the observed mean lenth data and the other factors
  dplyr::left_join(mean_obs_sub[, c("yr", "mean","region","codgr")], by = "yr")

#assign the mean predictions for the species
mean_pred_blf=mean_pred #size mean predicted for all factors combined

#plotting the GAM's prediction and the observed data
#-------------
#plot section
#-------------
#plot GAM mean size prediction
p16<-ggplot(data=mean_pred,aes(x=yr,y=pred))+
  geom_ribbon(data=dplyr::filter(mean_pred,codsp==sp),
    aes(ymin = lw, ymax = up,linetype='Pred'), 
    col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(mean_pred,codsp==sp),
    aes(x=yr,y=pred,linetype='Pred'),col="gray60",size=2)+
  #geom_line(data=dplyr::filter(gam,codsp==sp),aes(x=yr,y=fit,col=codgr,linetype=region),size=1.2)+
  geom_point(data=dplyr::filter(mean_pred,codsp==sp),
    aes(x=yr,y=mean),size=2,col="gray40")+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=yr[1],y=mean[1]),pch=8,size=2,col="gray40")+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='')+
  scale_y_continuous() +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p16

ggsave(paste(c("Length_model_prediction_Small_Tunas_",sp,".png"),collapse=''), plot = p16, device = "png", units = "cm",
  width = 24, height = 14)
#----------------------------------------------

#partial effects 
newdata= mean_obs_sub   #desired partial effect

pareffect<- data.frame(predict(gam_model, newdata = newdata,
                               type = "terms", se.fit=TRUE)) #terms gets the effects

#intervals for the partial effect
pareffect$up.codgr<-pareffect$fit.codgr+qnorm(0.975)*pareffect$se.fit.codgr #confidence intervals
pareffect$lw.codgr<-pareffect$fit.codgr-qnorm(0.975)*pareffect$se.fit.codgr
pareffect$codgr= mean_obs_sub$codgr

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

# plot GAM partial effects
p17<-ggplot(dplyr::filter(pareffect,codsp==sp), 
  aes(x=effect,y=fit,fill=var,col="Gear"),size=2) + 
  #geom_line(aes(x=effect,y=up, col="Gear",group=1)) +
  geom_point(aes(x=effect,y=fit,col="Gear"),size=3.5)+
  geom_errorbar(aes(ymin=lw, ymax=up,col="Gear"), width=0.5,
    position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Effect",x="Gear",col='',fill="")+
  scale_color_viridis_d()+
  scale_y_continuous() +
  scale_x_discrete() +
  #scale_colour_manual(values="gray60") +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p17

ggsave(paste(c("Effects_model_Small_Tunas_",sp,".png"),collapse=''), plot = p17, device = "png", units = "cm",
  width = 24, height = 14)


#----------------------------------
#simulating frequency distributions
#----------------------------------
dist_sim<-distsim(
  mean = mean_pred$pred[mean_pred$codsp==sp],
  cv= cv,
  years= unique(mean_pred$yr[mean_pred$codsp==sp]),
  n=100,
  sp=sp,
  distribution = distribution,
  source = "Dist_Sim") 

#assign the distribution simulated for the species
dist_sim_blf=dist_sim


#-------------
#plot section
#-------------
#Mean predicted (line) + Distribution simulation (boxplots)
p18 <- ggplot() +
  geom_boxplot(data = dplyr::filter(dist_sim, codsp == sp), 
    aes(x = factor(yr), y = fl, fill = codsp, col = codsp),
    alpha = 0.3,
    position = "identity",
    lwd = 0.75,
    outlier.shape = 8,
    outlier.size = 0.8,
    outlier.alpha = 0.5) +
  geom_ribbon(data = dplyr::filter(mean_pred, codsp == sp),
    aes(x = factor(yr), ymin = lw, ymax = up, group = 1, linetype = 'Pred'),
    color = "white", fill = "grey80", alpha = 0.4) +
  geom_line(data = dplyr::filter(mean_pred, codsp == sp),
    aes(x = factor(yr), y = pred, group = 1, linetype = 'Pred'), 
    col = "gray60", size = 2) +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr), y = mean, group = 1), size = 2) +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr[1]), y = mean[1]), pch = 8, size = 2, col = "deepskyblue") +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr[1]), y = mean[1]), size = 2, col = "deepskyblue") +
  labs(x = "Year", y = "Fork length (cm)", col = "", linetype = '', fill = '') +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  scale_y_continuous() +
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10)))) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p18

ggsave(paste(c("Distribution_Simulation_Boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p18, device = "png", units = "cm",
  width = 24, height = 15)
#-------------------------------------

#plot half violin for the predicted frequencies
p19<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(dist_sim,codsp==sp), 
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
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr),y=mean),size=2)+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr[1]),y=mean[1]),pch=8,size=2)+
  # #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2025,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p19

ggsave(paste(c("Distribution_Simulation_HalfViolin_Small_Tunas_",sp,".png"),collapse = ''), plot = p19, device = "png", units = "cm",
  width = 35, height = 12)
#--------------------------


#--------------------------------------------------------
#Comparing Distribution simulation with LBSPR simulation
#--------------------------------------------------------

#Life history data is necessary in this section
linf=pars$linf[pars$codsp==sp] #growth
k=pars$k[pars$codsp==sp]   #growth
t0=pars$t0[pars$codsp==sp]  #growth
l50=pars$l50[pars$codsp==sp] #maturity
l95= pars$l50[pars$codsp==sp]*1.1 #maturity
m=0.94 #mortality
mk=0.94/pars$k[pars$codsp==sp] #M=0.94 (Freire et al., 2005)  
by=5 #binwidth


                  #----- Comparing the Simulation Distributions -------#
                  #   Comparing with the LBSPR Simulation function     #
                  #----------------------------------------------------#              
#empty frame to store the simulated frequencies
comp_sim<-data.frame(yr=NULL,codsp=NULL,lmids=NULL,sim1=NULL,sim2=NULL,sel_method=NULL,
                                       fm_method=NULL,sl_fm_methods=NULL,sl50=NULL,f_m=NULL,res=NULL,rmse=NULL,x2=NULL)  

for (i in years) {
  
  # Sim1 is the distribution simulation from the present study
  sim1<-  lf(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], by=by)
  sim1$counts<- norm(sim1$counts)
  
  # Sim2 is the LBSPR simulation population 
  #Selectivity from Catch-curve  and Selectivity from the LBSPR fit 
  #F/M from the Beverton and Holt and F/M from the LBSPR fit 

    #Selectivity from Catch-curve and Selectivity from the Beverton and Holt-------------
    sim2_cc_bh<-lbsprsim(
      linf = linf,
      l50= l50,
      l95= l95,
      mk= mk,
      units = "cm",
      sl50= sel(x=c(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by-6,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
      sl95= 1.17*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
      fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"), method="Beverton and Holt",sp=sp, units="cm"), 
      minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
      maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
      binwidth = by,
      sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_cc_bh@SL50
      f_m= sim2_cc_bh@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_cc_bh@LMids,counts=sim2_cc_bh@pLCatch*100) %>% 
      tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_cc_bh<- data.frame(mids=sim2_cc_bh@LMids,counts=norm(sim2_cc_bh@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_cc_bh <- subset(sim2_cc_bh, mids >= min(sim1$mids) & mids <= max(sim1$mids))
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_cc_bh$counts,
      sel_method="CC", fm_method="BH",sl_fm_methods= paste("CC_BH"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_cc_bh$counts,
      rmse=rmse(sim1$counts,sim2_cc_bh$counts),
      x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
      sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
      comp_sim<-rbind(comp_sim,sim_dat)
  
      #Selectivity from Catch-curve and F/M from the LBSPR fit-------------
      sim2_cc_lb<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        sl95= 1.17*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"), method="LBSPR",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
      
      #store Selectivity and F/M estimated
      sl50=sim2_cc_lb@SL50
      f_m= sim2_cc_lb@FM
      
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_cc_lb@LMids,counts=sim2_cc_lb@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
      
      #storing LBSPR distribution simulation
      sim2_cc_lb<- data.frame(mids=sim2_cc_lb@LMids,counts=norm(sim2_cc_lb@pLCatch) )
      
      # Adjusting to keep the same mids interval
      sim2_cc_lb <- subset(sim2_cc_lb, mids >= min(sim1$mids) & mids <= max(sim1$mids))
      
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_cc_lb$counts,
        sel_method="CC", fm_method="LBSPR",sl_fm_methods= paste("CC_LB"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_cc_lb$counts,
        rmse=rmse(sim1$counts,sim2_cc_lb$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
          sample(freq_lbspr, 95,replace = FALSE))$p.value))
      
      comp_sim<-rbind(comp_sim,sim_dat)
  
  
      #Selectivity from LBSPR and F/M from Beverton and Holt------------
      sim2_lb_bh<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        sl95= 1.17*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),method="Beverton and Holt",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
      
      #store Selectivity and F/M estimated
      sl50=sim2_lb_bh@SL50
      f_m= sim2_lb_bh@FM
      
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_lb_bh@LMids,counts=sim2_lb_bh@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
      
      #storing LBSPR distribution simulation
      sim2_lb_bh<- data.frame(mids=sim2_lb_bh@LMids,counts=norm(sim2_lb_bh@pLCatch) )
      
      # Adjusting to keep the same mids interval
      sim2_lb_bh <- subset(sim2_lb_bh, mids >= min(sim1$mids) & mids <= max(sim1$mids))
      
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_lb_bh$counts,
        sel_method="LB", fm_method="BH",sl_fm_methods= paste("LB_BH"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_lb_bh$counts,
        rmse=rmse(sim1$counts,sim2_lb_bh$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
          sample(freq_lbspr, 95,replace = FALSE))$p.value))
      
      comp_sim<-rbind(comp_sim,sim_dat)
      
  
      #Selectivity from LBSPR and F/M from LBSPR -----------
      sim2_lb_lb<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        sl95= 1.17*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"), method="LBSPR",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
      
      #store Selectivity and F/M estimated
      sl50=sim2_lb_lb@SL50
      f_m= sim2_lb_lb@FM
      
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_lb_lb@LMids,counts=sim2_lb_lb@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
      
      #storing LBSPR distribution simulation
      sim2_lb_lb<- data.frame(mids=sim2_lb_lb@LMids,counts=norm(sim2_lb_lb@pLCatch) )
      
      # Adjusting to keep the same mids interval
      sim2_lb_lb <- subset(sim2_lb_lb, mids >= min(sim1$mids) & mids <= max(sim1$mids))
      
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_lb_lb$counts,
        sel_method="LB", fm_method="LB",sl_fm_methods= paste("LB_LB"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_lb_lb$counts,
        rmse=rmse(sim1$counts,sim2_lb_lb$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
          sample(freq_lbspr, 95,replace = FALSE))$p.value))
      
      comp_sim<-rbind(comp_sim,sim_dat)
}

#assing the comparisons for the current species
comp_sim_blf<-comp_sim

#-------------
#plot section
#-------------
#Length residual comparison 
p20 <- ggplot() +
  geom_boxplot(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(group = interaction(factor(yr), factor(sl_fm_methods)), 
      x = factor(yr), 
      y = res, 
      fill = factor(sl_fm_methods)),  # fill to methods
    alpha = 1,
    position = position_dodge(width = 0.7), 
    width = 0.9,
    lwd = 0.1,  
    outlier.shape = 8,
    outlier.size = 0.7,
    outlier.alpha = 0.5) +
  geom_hline(yintercept =0,linetype = "dashed" )+
  geom_smooth(data = dplyr::filter(comp_sim, codsp == sp),
    method = "loess", span = 0.2, se = FALSE,
    aes(group = 1, x = factor(yr), y = res,color="Loess"), size = 1.5) +
  geom_text(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = "2015", y = 0.9 * max(res, na.rm = TRUE), 
      label = paste("RMSE=", round(mean(rmse, na.rm = TRUE), 2))),
    inherit.aes = FALSE) +
  labs(x = "Year", y = "Length residual", col = "", linetype = '', fill = '') +
  scale_y_continuous() +
  scale_fill_viridis_d(direction = -1)+
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10)))) +
  theme_classic(base_size = 14) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p20

ggsave(paste(c("Length_Residual_Boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p20, device = "png", units = "cm",
  width = 25, height = 17)
#------------------------

#Scaterplot comparing the distributions (Dist sim vs LBSPR sim)
p21 <- ggplot() +
  geom_point(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = sim2, y = sim1, color = sl_fm_methods), 
    size = 2, alpha = 0.7) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", size = 1) +
  geom_smooth(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = sim2, y = sim1, color = sl_fm_methods, fill = sl_fm_methods), 
    method = "loess", se = FALSE, size = 1.2, alpha = 0.3) +  
  scale_color_viridis_d(direction = -1) +  
  labs(x = "LBSPR Distribution Simulation",y= "Distribution Simulation",col = "Method", shape = "Method", fill = "Method") +
  theme_classic(base_size = 14) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p21

ggsave(paste(c("Distributions_Comparison_Small_Tunas_",sp,".png"),collapse = ''), plot = p21, device = "png", units = "cm",
  width = 24, height = 15)
#x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x#





#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
# models by species
# BRS model - Length Reconstruction 
#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
sp='BRS'
years=1950:2025

#-------------------------------------
#find the best distribution for specie
#-------------------------------------
dist<-distfit(len = data.frame(
  yr=freq_obs$yr[freq_obs$codsp==sp],
  codsp=freq_obs$codsp[freq_obs$codsp==sp],
  codgr=freq_obs$codgr[freq_obs$codsp==sp],
  fl=freq_obs$fl[freq_obs$codsp==sp]))

dist=dist %>%
  dplyr::group_by(yr,codsp,codgr) %>% 
  slice_min(aic, n = 1)  %>%
  dplyr::arrange(dist) 

#assign dist for the current species
dist_brs=dist

#Table with the best distribution
table(dist$codsp,dist$dist)  
distribution='norm'

#coefficient of variation for the corresponding distribution 
cv=mean(dist$cv[dist$codsp==sp & dist$dist==distribution])

#-------------------------------------------------------
#predict 1950 length when using the slope of the lengths
#-------------------------------------------------------

#assuming the 1950 mean length as the prediction of a linear model
lm_data<- data.frame(length=mean_obs$mean[mean_obs$codsp==sp],year=mean_obs$yr[mean_obs$codsp==sp])


#fiting the model
lm_model<- lm(length~year, data = lm_data[-1,])
#Parameters... On average the mean length is being reduced 0.10 cm year by year
a<- coef(lm_model)[1]
b<- coef(lm_model)[2]

#using the fitted parameters to estimate what would be the length in 1950
newdata<- data.frame(year=1950:2025)
lm_predict_10<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.1)))
lm_predict_20<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.2)))
lm_predict_30<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.3)))

#data frame with predictions for different percentile levels
lm_predict<- data.frame(yr=newdata$year,
  fit=lm_predict_10$fit.fit, 
  lw_10=lm_predict_10$fit.lwr,
  up_10=lm_predict_10$fit.upr,
  lw_20=lm_predict_20$fit.lwr,
  up_20=lm_predict_20$fit.upr,
  lw_30=lm_predict_30$fit.lwr,
  up_30=lm_predict_30$fit.upr,
  codsp=sp)

#lm model coefficients and predicitions for the current species
lm_model_brs<- data.frame(formula=paste(deparse(lm_model$call), collapse = ""),
  intercept=a,
  slope=b,
  codsp=sp)
lm_predict_brs<- lm_predict

# controlling variables
lm_predict <- lm_predict %>%
  mutate(confidence_level = case_when(
    yr == 1950 & lw_10 == lw_10 ~ "10%",
    yr == 1950 & lw_20 == lw_20 ~ "20%",
    yr == 1950 & lw_30 == lw_30 ~ "30%"
  ))

# visual representation of the linear model
#-------------
#plot section
#-------------
p22 <- ggplot(lm_predict, aes(x = yr, y = fit)) +
  geom_line(size = 1, color = "grey20") +  # main line prediction
  # Lines for the percentiles
  geom_line(aes(y = lw_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = up_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = lw_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = up_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = lw_30, color = "30%", linetype = "30%"), size = 0.8) +
  geom_line(aes(y = up_30, color = "30%", linetype = "30%"), size = 0.8) +
  # Observed data
  geom_point(data = lm_data, aes(x = year, y = length), size = 2, color = "grey30") +
  # Predicted points for 1950 
  geom_point(aes(x = yr[1], y = lw_10[1], fill = "10%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_10[1], fill = "10%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = lw_20[1], fill = "20%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_20[1], fill = "20%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = lw_30[1], fill = "30%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_30[1], fill = "30%"), size = 3, shape = 21) +
  labs(y = 'Mean Length (cm)', x = "Year", color = "Percentile", fill = "Percentile", linetype = "Percentile") +
  # Type and line colors
  scale_color_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_fill_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_linetype_manual(values = c("10%" = "dashed", "20%" = "dotted", "30%" = "dotdash")) +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p22
#saving graph...
ggsave(paste0("Linear_model_prediction_Small_Tunas_", sp, ".png"), plot = p22, device = "png", units = "cm", width = 16, height = 10)
#----------------------------------

#selecting variables .. subseting
mean_obs_sub= mean_obs[,c("yr",'region','codsp','codgr','mean')]

#removing extreme values
mean_obs_sub$mean[mean_obs_sub$codsp=='BLF' & mean_obs_sub$mean<50] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='DOL' & mean_obs_sub$mean<75] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='KGM' & mean_obs_sub$mean>100] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='KGM' & mean_obs_sub$mean<60] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean<100] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean>180] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='FRI' & mean_obs_sub$mean>45] =NA

#get only non-NA values and selecting the current species
mean_obs_sub=mean_obs_sub[is.na(mean_obs_sub$mean)==FALSE & mean_obs_sub$codsp==sp,]

#transforming columns in to factors
mean_obs_sub$region= factor(mean_obs_sub$region)
mean_obs_sub$codgr= factor(mean_obs_sub$codgr)
mean_obs_sub$codsp= factor(mean_obs_sub$codsp)


# model GAM selection----------------
#list of values to be tested
tests<- names(lm_predict[, grepl("fit|lw|up", names(lm_predict))])
#filling an empty 1950 with NA 
mean_obs_sub<- rbind(data.frame(yr=1950, region="NE", codsp=sp,codgr="UN",mean=NA),
  mean_obs_sub)

# Definindo os vetores de variáveis de resposta e explicativas
resp_vars <- c("mean")
exp_vars_1 <- c("s(yr)")
exp_vars_2 <- c("region")
exp_vars_3 <- c("codgr")

#criando a formula geral
formulas <-list(as.formula(paste(resp_vars, "~", paste(c(exp_vars_1), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_2), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_3), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_2,exp_vars_3), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_3,exp_vars_2), collapse = "+"))))

formulas<- as.character(formulas)

# Data frame com as distribuicoes de probabilidade e funcoes de ligacao a serem testadas
family=distribution #testing only for the best distribution
link = ("log") #log fixed

#empty frame to store results
gam=data.frame(
  test=NULL,
  type=NULL,
  codsp=NULL,
  family=NULL,
  link=NULL,
  formula=NULL,
  newdataexpand=NULL,
  pred=NULL,
  aic=NULL,
  r.sq=NULL,
  dev.ex=NULL)


        #-------------------- Fitting models  -------------------#
        for (i in tests) {
  
          #1- first sensibility for different 1950 starting values--
          test_1950<- lm_predict[1, grepl(i, names(lm_predict))]
  
           #replacing with the current 1950 mean length test
           mean_obs_sub$mean[mean_obs_sub$yr==1950]<- test_1950
  
           #2- Term selection within the sensibility----
           # Loop through the formulas 
           for (j in formulas) {
    
            family_name <- as.character(family) #these are fixed
            link_name <- as.character(link)
    
            #Creating switch object
            family_obj <- switch(family_name,
            norm = gaussian(link = link_name),  
            invgauss = inverse.gaussian(link = link_name),
            gamma= Gamma(link=link_name))
    
            #formula j to be tested
            formula <- as.formula(j)
    
            #fit the model
            gam_model <- mgcv::gam(as.formula(formula),
            data= mean_obs_sub,
            method = "REML" ,
            family =family_obj)
    
            #taking the outputs
            gam=rbind(gam, data.frame(
            test=i,
            type="GAM",
            codsp=sp,
            family=gam_model$family$family,
            link=gam_model$family$link,
            formula=paste0(c(paste(gam_model$formula)[2],
            paste(gam_model$formula)[1],
          paste(gam_model$formula)[3]),collapse = ' '),
          aic=gam_model$aic,
         r.sq=summary(gam_model)$r.sq,
         dev.ex=summary(gam_model)$dev.expl))
         }
       }
#-------------------- end of loop -----------------------#

#assign the tests to the species
gam_brs<-gam

#replacing the 1950 mean length with the best tested value
best_model<- gam[which.min(gam$aic),]; dplyr::tibble(best_model)
mean_obs_sub$mean[mean_obs_sub$yr==1950]<- lm_predict[1, grepl(best_model$test, names(lm_predict))]

#Final model (lowest AIC) 
gam_model<- mgcv::gam(as.formula(best_model$formula), 
  #sp=0.01,
  data= mean_obs_sub,
  method = "REML" ,
  family =family_obj)

#assign GAM model to the specie
gam_model_brs=gam_model

#diagnostic
summary(gam_model)

par(mfrow=c(2,2))
mgcv::gam.check(gam_model)# model check

#residual 
#homocedasticity
lmtest::bptest(gam_model)

#normality
shapiro.test(resid(gam_model))

# Check overall concurvity
round(mgcv::concurvity(gam_model, full = TRUE),2)
#----------------------
#predictions "year"
newdataexpand= expand.grid(yr=seq(years[1],years[length(years)],1),
  region=mean_obs_sub$region,
  codgr=mean_obs_sub$codgr)

pred<- data.frame(predict(gam_model, newdata = newdataexpand,
  type = "response", se.fit=TRUE))

pred$up<-pred$fit+qnorm(0.975)*pred$se.fit #confidence intervals
pred$lw<-pred$fit-qnorm(0.975)*pred$se.fit

#take the predictions
mean_pred=cbind(type='GAM',
  codsp=sp,
  yr=newdataexpand$yr,
  mean=NA,
  newdata=newdataexpand,
  pred = round(pred$fit, 2),
  lw = round(pred$lw, 2),
  up = round(pred$up, 2))

#extract mean of all levels
mean_pred = mean_pred %>%
  group_by(type, yr, codsp) %>%
  summarize(across(c('pred', 'lw', 'up'), ~ round(mean(.), 2)), .groups = 'drop') %>%
  #binding the observed mean lenth data and the other factors
  dplyr::left_join(mean_obs_sub[, c("yr", "mean","region","codgr")], by = "yr")

#assign the mean predictions for the species
mean_pred_brs=mean_pred #size mean predicted for all factors combined

#plotting the GAM's prediction and the observed data
#-------------
#plot section
#-------------
#plot GAM mean size prediction
p23<-ggplot(data=mean_pred,aes(x=yr,y=pred))+
  geom_ribbon(data=dplyr::filter(mean_pred,codsp==sp),
    aes(ymin = lw, ymax = up,linetype='Pred'), 
    col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(mean_pred,codsp==sp),
    aes(x=yr,y=pred,linetype='Pred'),col="gray60",size=2)+
  #geom_line(data=dplyr::filter(gam,codsp==sp),aes(x=yr,y=fit,col=codgr,linetype=region),size=1.2)+
  geom_point(data=dplyr::filter(mean_pred,codsp==sp),
    aes(x=yr,y=mean),size=2,col="gray40")+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=yr[1],y=mean[1]),pch=8,size=2,col="gray40")+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='')+
  scale_y_continuous() +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p23

ggsave(paste(c("Length_model_prediction_Small_Tunas_",sp,".png"),collapse=''), plot = p23, device = "png", units = "cm",
  width = 24, height = 14)
#----------------------------------------------

# #partial effects 
# newdata= mean_obs_sub   #desired partial effect
# 
# pareffect<- data.frame(predict(gam_model, newdata = newdata,
#   type = "terms", se.fit=TRUE)) #terms gets the effects
# 
# #intervals for the partial effect
# pareffect$up.codgr<-pareffect$fit.codgr+qnorm(0.975)*pareffect$se.fit.codgr #confidence intervals
# pareffect$lw.codgr<-pareffect$fit.codgr-qnorm(0.975)*pareffect$se.fit.codgr
# pareffect$codgr= mean_obs_sub$codgr
# 
# # Extract mean of all levels
# pareffect = pareffect[grepl("codgr", names(pareffect))]
# pareffect = pareffect %>%  # select only the desired partial effect
#   group_by(codgr) %>%
#   summarize(across(c("fit.codgr", "se.fit.codgr", "up.codgr", "lw.codgr"),
#     ~ round(mean(.), 2))) %>%
#   rename(effect = codgr,
#     fit = fit.codgr, 
#     se.fit = se.fit.codgr,
#     up = up.codgr,
#     lw = lw.codgr) %>%
#   mutate(codsp = sp, var = "Gear")
# 
# pareffect_blf = pareffect
# 
# # plot GAM partial effects
# p17<-ggplot(dplyr::filter(pareffect,codsp==sp), 
#   aes(x=effect,y=fit,fill=var,col="Gear"),size=2) + 
#   #geom_line(aes(x=effect,y=up, col="Gear",group=1)) +
#   geom_point(aes(x=effect,y=fit,col="Gear"),size=3.5)+
#   geom_errorbar(aes(ymin=lw, ymax=up,col="Gear"), width=0.5,
#     position="identity",size=1.2)+
#   #facet_wrap(~codsp,scales = "free_y")+
#   labs(y="Effect",x="Gear",col='',fill="")+
#   scale_color_viridis_d()+
#   scale_y_continuous() +
#   scale_x_discrete() +
#   #scale_colour_manual(values="gray60") +
#   theme_classic(base_size = 14) %+replace%
#   theme(strip.background=element_blank(),
#     plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
# p17
# 
# ggsave(paste(c("Effects_model_Small_Tunas_",sp,".png"),collapse=''), plot = p17, device = "png", units = "cm",
#   width = 24, height = 14)


#----------------------------------
#simulating frequency distributions
#----------------------------------
dist_sim<-distsim(
  mean = mean_pred$pred[mean_pred$codsp==sp],
  cv= cv,
  years= unique(mean_pred$yr[mean_pred$codsp==sp]),
  n=100,
  sp=sp,
  distribution = distribution,
  source = "Dist_Sim") 

#assign the distribution simulated for the species
dist_sim_brs=dist_sim


#-------------
#plot section
#-------------
#Mean predicted (line) + Distribution simulation (boxplots)
p24 <- ggplot() +
  geom_boxplot(data = dplyr::filter(dist_sim, codsp == sp), 
    aes(x = factor(yr), y = fl, fill = codsp, col = codsp),
    alpha = 0.3,
    position = "identity",
    lwd = 0.75,
    outlier.shape = 8,
    outlier.size = 0.8,
    outlier.alpha = 0.5) +
  geom_ribbon(data = dplyr::filter(mean_pred, codsp == sp),
    aes(x = factor(yr), ymin = lw, ymax = up, group = 1, linetype = 'Pred'),
    color = "white", fill = "grey80", alpha = 0.4) +
  geom_line(data = dplyr::filter(mean_pred, codsp == sp),
    aes(x = factor(yr), y = pred, group = 1, linetype = 'Pred'), 
    col = "gray60", size = 2) +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr), y = mean, group = 1), size = 2) +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr[1]), y = mean[1]), pch = 8, size = 2, col = "deepskyblue") +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr[1]), y = mean[1]), size = 2, col = "deepskyblue") +
  labs(x = "Year", y = "Fork length (cm)", col = "", linetype = '', fill = '') +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  scale_y_continuous() +
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10)))) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p24

ggsave(paste(c("Distribution_Simulation_Boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p24, device = "png", units = "cm",
  width = 24, height = 15)
#-------------------------------------

#plot half violin for the predicted frequencies
p25<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(dist_sim,codsp==sp), 
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
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr),y=mean),size=2)+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr[1]),y=mean[1]),pch=8,size=2)+
  # #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2025,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p25

ggsave(paste(c("Distribution_Simulation_HalfViolin_Small_Tunas_",sp,".png"),collapse = ''), plot = p25, device = "png", units = "cm",
  width = 35, height = 12)
#--------------------------


#--------------------------------------------------------
#Comparing Distribution simulation with LBSPR simulation
#--------------------------------------------------------
#Life history data is necessary in this section
linf=pars$linf[pars$codsp==sp] #growth
k=0.16   #~k=0.16(Nobrega, 2002) #growth
t0=pars$t0[pars$codsp==sp]  #growth
l50=pars$l50[pars$codsp==sp] #maturity
l95= pars$l50[pars$codsp==sp]*1.1 #maturity
m=0.36 ##m=0.36, 
mk=0.36/0.16 #m=0.36, ~k=0.16(Nobrega, 2002)
by=5 #binwidth


#----- Comparing the Simulation Distributions -------#
#   Comparing with the LBSPR Simulation function     #
#----------------------------------------------------#              
#empty frame to store the simulated frequencies
comp_sim<-data.frame(yr=NULL,codsp=NULL,lmids=NULL,sim1=NULL,sim2=NULL,sel_method=NULL,
  fm_method=NULL,sl_fm_methods=NULL,sl50=NULL,f_m=NULL,res=NULL,rmse=NULL,x2=NULL)  

for (i in years) {
  
    # Sim1 is the distribution simulation from the present study
    sim1<-  lf(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], by=by)
    sim1$counts<- norm(sim1$counts)
  
    # Sim2 is the LBSPR simulation population 
    #Selectivity from Catch-curve  and Selectivity from the LBSPR fit 
    #F/M from the Beverton and Holt and F/M from the LBSPR fit 
  
    #Selectivity from Catch-curve and Selectivity from the Beverton and Holt-------------
    sim2_cc_bh<-lbsprsim(
      linf = linf,
      l50= l50,
      l95= l95,
      mk= mk,
      units = "cm",
      sl50= sel(x=c(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by-6,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
      sl95= 1.26*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
      fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"), method="Beverton and Holt",sp=sp, units="cm"), 
      minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
      maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
      binwidth = by,
      sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_cc_bh@SL50
      f_m= sim2_cc_bh@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_cc_bh@LMids,counts=sim2_cc_bh@pLCatch*100) %>% 
                  tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_cc_bh<- data.frame(mids=sim2_cc_bh@LMids,counts=norm(sim2_cc_bh@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_cc_bh <- subset(sim2_cc_bh, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_cc_bh$counts,
        sel_method="CC", fm_method="BH",sl_fm_methods= paste("CC_BH"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_cc_bh$counts,
        rmse=rmse(sim1$counts,sim2_cc_bh$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 90,replace = FALSE), 
        sample(freq_lbspr, 90,replace = FALSE))$p.value))
  
      comp_sim<-rbind(comp_sim,sim_dat)
  
      #Selectivity from Catch-curve and F/M from the LBSPR fit-------------
      sim2_cc_lb<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        sl95= 1.26*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"), method="LBSPR",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_cc_lb@SL50
      f_m= sim2_cc_lb@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_cc_lb@LMids,counts=sim2_cc_lb@pLCatch*100) %>% 
      tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_cc_lb<- data.frame(mids=sim2_cc_lb@LMids,counts=norm(sim2_cc_lb@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_cc_lb <- subset(sim2_cc_lb, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_cc_lb$counts,
        sel_method="CC", fm_method="LBSPR",sl_fm_methods= paste("CC_LB"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_cc_lb$counts,
        rmse=rmse(sim1$counts,sim2_cc_lb$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
        sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
      comp_sim<-rbind(comp_sim,sim_dat)
  
  
      #Selectivity from LBSPR and F/M from Beverton and Holt------------
      sim2_lb_bh<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        sl95= 1.26*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),method="Beverton and Holt",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_lb_bh@SL50
      f_m= sim2_lb_bh@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_lb_bh@LMids,counts=sim2_lb_bh@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_lb_bh<- data.frame(mids=sim2_lb_bh@LMids,counts=norm(sim2_lb_bh@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_lb_bh <- subset(sim2_lb_bh, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_lb_bh$counts,
        sel_method="LB", fm_method="BH",sl_fm_methods= paste("LB_BH"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_lb_bh$counts,
        rmse=rmse(sim1$counts,sim2_lb_bh$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
        sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
      comp_sim<-rbind(comp_sim,sim_dat)
  
  
      #Selectivity from LBSPR and F/M from LBSPR -----------
      sim2_lb_lb<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        sl95= 1.26*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"), method="LBSPR",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_lb_lb@SL50
      f_m= sim2_lb_lb@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_lb_lb@LMids,counts=sim2_lb_lb@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_lb_lb<- data.frame(mids=sim2_lb_lb@LMids,counts=norm(sim2_lb_lb@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_lb_lb <- subset(sim2_lb_lb, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_lb_lb$counts,
        sel_method="LB", fm_method="LB",sl_fm_methods= paste("LB_LB"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_lb_lb$counts,
        rmse=rmse(sim1$counts,sim2_lb_lb$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
        sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
  comp_sim<-rbind(comp_sim,sim_dat)
}
#assing the comparisons for the current species
comp_sim_brs<-comp_sim

#-------------
#plot section
#-------------
#Length residual comparison 
p26 <- ggplot() +
  geom_boxplot(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(group = interaction(factor(yr), factor(sl_fm_methods)), 
      x = factor(yr), 
      y = res, 
      fill = factor(sl_fm_methods)),  # fill to methods
    alpha = 1,
    position = position_dodge(width = 0.7), 
    width = 0.9,
    lwd = 0.1,  
    outlier.shape = 8,
    outlier.size = 0.7,
    outlier.alpha = 0.5) +
  geom_hline(yintercept =0,linetype = "dashed" )+
  geom_smooth(data = dplyr::filter(comp_sim, codsp == sp),
    method = "loess", span = 0.2, se = FALSE,
    aes(group = 1, x = factor(yr), y = res,color="Loess"), size = 1.5) +
  geom_text(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = "2015", y = 0.9 * max(res, na.rm = TRUE), 
      label = paste("RMSE=", round(mean(rmse, na.rm = TRUE), 2))),
    inherit.aes = FALSE) +
  labs(x = "Year", y = "Length residual", col = "", linetype = '', fill = '') +
  scale_y_continuous() +
  scale_fill_viridis_d(direction = -1)+
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10)))) +
  theme_classic(base_size = 14) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p26

ggsave(paste(c("Length_Residual_Boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p26, device = "png", units = "cm",
  width = 25, height = 17)
#------------------------

#Scaterplot comparing the distributions (Dist sim vs LBSPR sim)
p27 <- ggplot() +
  geom_point(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = sim2, y = sim1, color = sl_fm_methods), 
    size = 2, alpha = 0.7) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", size = 1) +
  geom_smooth(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = sim2, y = sim1, color = sl_fm_methods, fill = sl_fm_methods), 
    method = "loess", se = FALSE, size = 1.2, alpha = 0.3) +  
  scale_color_viridis_d(direction = -1) +  
  labs(x = "LBSPR Distribution Simulation",y= "Distribution Simulation",col = "Method", shape = "Method", fill = "Method") +
  theme_classic(base_size = 14) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p27

ggsave(paste(c("Distributions_Comparison_Small_Tunas_",sp,".png"),collapse = ''), plot = p27, device = "png", units = "cm",
  width = 24, height = 15)
#-----------------------------------------------------------------------------------------------------------


#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
# models by species
# DOL model - Length Reconstruction 
#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
sp='DOL'
years=1950:2025

#-------------------------------------
#find the best distribution for specie
#-------------------------------------
dist<-distfit(len = data.frame(
  yr=freq_obs$yr[freq_obs$codsp==sp],
  codsp=freq_obs$codsp[freq_obs$codsp==sp],
  codgr=freq_obs$codgr[freq_obs$codsp==sp],
  fl=freq_obs$fl[freq_obs$codsp==sp]))

dist=dist %>%
  dplyr::group_by(yr,codsp,codgr) %>% 
  slice_min(aic, n = 1)  %>%
  dplyr::arrange(dist) 

#assign dist for the current species
dist_dol=dist

#Table with the best distribution
table(dist$codsp,dist$dist)  
distribution='invgauss'

#coefficient of variation for the corresponding distribution 
cv=mean(dist$cv[dist$codsp==sp & dist$dist==distribution])

#-------------------------------------------------------
#predict 1950 length when using the slope of the lengths
#-------------------------------------------------------

#assuming the 1950 mean length as the prediction of a linear model
lm_data<- data.frame(length=mean_obs$mean[mean_obs$codsp==sp],year=mean_obs$yr[mean_obs$codsp==sp])
lm_data<- lm_data[lm_data$length>75,] #removing extreme values

#fiting the model
lm_model<- lm(length~year, data = lm_data[-1,])
#Parameters... On average the mean length is being reduced 0.10 cm year by year
a<- coef(lm_model)[1]
b<- coef(lm_model)[2]

#using the fitted parameters to estimate what would be the length in 1950
newdata<- data.frame(year=1950:2025)
lm_predict_10<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.1)))
lm_predict_20<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.2)))
lm_predict_30<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.3)))

#data frame with predictions for different percentile levels
lm_predict<- data.frame(yr=newdata$year,
  fit=lm_predict_10$fit.fit, 
  lw_10=lm_predict_10$fit.lwr,
  up_10=lm_predict_10$fit.upr,
  lw_20=lm_predict_20$fit.lwr,
  up_20=lm_predict_20$fit.upr,
  lw_30=lm_predict_30$fit.lwr,
  up_30=lm_predict_30$fit.upr,
  codsp=sp)

#lm model coefficients and predicitions for the current species
lm_model_dol<- data.frame(formula=paste(deparse(lm_model$call), collapse = ""),
  intercept=a,
  slope=b,
  codsp=sp)
lm_predict_dol<- lm_predict

# controlling variables
lm_predict <- lm_predict %>%
  mutate(confidence_level = case_when(
    yr == 1950 & lw_10 == lw_10 ~ "10%",
    yr == 1950 & lw_20 == lw_20 ~ "20%",
    yr == 1950 & lw_30 == lw_30 ~ "30%"
  ))

# visual representation of the linear model
#-------------
#plot section
#-------------
p28 <- ggplot(lm_predict, aes(x = yr, y = fit)) +
  geom_line(size = 1, color = "grey20") +  # main line prediction
  # Lines for the percentiles
  geom_line(aes(y = lw_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = up_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = lw_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = up_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = lw_30, color = "30%", linetype = "30%"), size = 0.8) +
  geom_line(aes(y = up_30, color = "30%", linetype = "30%"), size = 0.8) +
  # Observed data
  geom_point(data = lm_data, aes(x = year, y = length), size = 2, color = "grey30") +
  # Predicted points for 1950 
  geom_point(aes(x = yr[1], y = lw_10[1], fill = "10%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_10[1], fill = "10%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = lw_20[1], fill = "20%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_20[1], fill = "20%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = lw_30[1], fill = "30%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_30[1], fill = "30%"), size = 3, shape = 21) +
  labs(y = 'Mean Length (cm)', x = "Year", color = "Percentile", fill = "Percentile", linetype = "Percentile") +
  # Type and line colors
  scale_color_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_fill_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_linetype_manual(values = c("10%" = "dashed", "20%" = "dotted", "30%" = "dotdash")) +
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p28
#saving graph...
ggsave(paste0("Linear_model_prediction_Small_Tunas_", sp, ".png"), plot = p28, device = "png", units = "cm", width = 24, height = 16)
#----------------------------------

#selecting variables .. subseting
mean_obs_sub= mean_obs[,c("yr",'region','codsp','codgr','mean')]

#removing extreme values
mean_obs_sub$mean[mean_obs_sub$codsp=='BLF' & mean_obs_sub$mean<50] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='DOL' & mean_obs_sub$mean<75] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='KGM' & mean_obs_sub$mean>100] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='KGM' & mean_obs_sub$mean<60] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean<100] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean>180] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='FRI' & mean_obs_sub$mean>45] =NA

#get only non-NA values and selecting the current species
mean_obs_sub=mean_obs_sub[is.na(mean_obs_sub$mean)==FALSE & mean_obs_sub$codsp==sp,]

#grouping factors most likely (Fredou et al., 2021)
mean_obs_sub$region[mean_obs_sub$region=='ALL']='NE'
mean_obs_sub$codgr[mean_obs_sub$codgr=='UN']='TR'

#transforming columns in to factors
mean_obs_sub$region= factor(mean_obs_sub$region)
mean_obs_sub$codgr= factor(mean_obs_sub$codgr)
mean_obs_sub$codsp= factor(mean_obs_sub$codsp)


# model GAM selection----------------
#list of values to be tested
tests<- names(lm_predict[, grepl("fit|lw|up", names(lm_predict))])
#filling an empty 1950 with NA 
mean_obs_sub<- rbind(data.frame(yr=1950, region="NE", codsp=sp,codgr="UN",mean=NA),
  mean_obs_sub)

# Definindo os vetores de variáveis de resposta e explicativas
resp_vars <- c("mean")
exp_vars_1 <- c("s(yr)")
exp_vars_2 <- c("region")
exp_vars_3 <- c("codgr")

#criando a formula geral
formulas <-list(as.formula(paste(resp_vars, "~", paste(c(exp_vars_1), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_2), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_3), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_2,exp_vars_3), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_3,exp_vars_2), collapse = "+"))))

formulas<- as.character(formulas)

# Data frame com as distribuicoes de probabilidade e funcoes de ligacao a serem testadas
family=distribution #testing only for the best distribution
link = ("log") #log fixed

#empty frame to store results
gam=data.frame(
  test=NULL,
  type=NULL,
  codsp=NULL,
  family=NULL,
  link=NULL,
  formula=NULL,
  newdataexpand=NULL,
  pred=NULL,
  aic=NULL,
  r.sq=NULL,
  dev.ex=NULL)


#-------------------- Fitting models  -------------------#
for (i in tests) {
  
    #1- first sensibility for different 1950 starting values--
    test_1950<- lm_predict[1, grepl(i, names(lm_predict))]
  
      #replacing with the current 1950 mean length test
      mean_obs_sub$mean[mean_obs_sub$yr==1950]<- test_1950
  
      #2- Term selection within the sensibility----
      # Loop through the formulas 
      for (j in formulas) {
    
      family_name <- as.character(family) #these are fixed
      link_name <- as.character(link)
    
      #Creating switch object
      family_obj <- switch(family_name,
        norm = gaussian(link = link_name),  
        invgauss = inverse.gaussian(link = link_name),
        gamma= Gamma(link=link_name))
    
      #formula j to be tested
      formula <- as.formula(j)
    
      #fit the model
      gam_model <- mgcv::gam(as.formula(formula),
        data= mean_obs_sub,
        method = "REML" ,
        family =family_obj)
    
      #taking the outputs
      gam=rbind(gam, data.frame(
        test=i,
        type="GAM",
        codsp=sp,
        family=gam_model$family$family,
        link=gam_model$family$link,
        formula=paste0(c(paste(gam_model$formula)[2],
        paste(gam_model$formula)[1],
        paste(gam_model$formula)[3]),collapse = ' '),
      aic=gam_model$aic,
      r.sq=summary(gam_model)$r.sq,
 dev.ex=summary(gam_model)$dev.expl))
  }
}
#-------------------- end of loop -----------------------#

#assign the tests to the species
gam_dol<-gam

#replacing the 1950 mean length with the best tested value
best_model<- gam[which.min(gam$aic),]; dplyr::tibble(best_model)
mean_obs_sub$mean[mean_obs_sub$yr==1950]<- lm_predict[1, grepl(best_model$test, names(lm_predict))]

#Final model (lowest AIC) 
gam_model<- mgcv::gam(as.formula(best_model$formula), 
  #sp=0.0001,
  data= mean_obs_sub,
  method = "REML" ,
  family =family_obj)

#assign GAM model to the specie
gam_model_dol=gam_model

#diagnostic
summary(gam_model)

par(mfrow=c(2,2))
mgcv::gam.check(gam_model)# model check

#residual 
#homocedasticity
lmtest::bptest(gam_model)

#normality
shapiro.test(resid(gam_model))

# Check overall concurvity
round(mgcv::concurvity(gam_model, full = TRUE),2)
#----------------------
#predictions "year"
newdataexpand= expand.grid(yr=seq(years[1],years[length(years)],1),
  region=mean_obs_sub$region,
  codgr=mean_obs_sub$codgr)

pred<- data.frame(predict(gam_model, newdata = newdataexpand,
  type = "response", se.fit=TRUE))

pred$up<-pred$fit+qnorm(0.975)*pred$se.fit #confidence intervals
pred$lw<-pred$fit-qnorm(0.975)*pred$se.fit

#take the predictions
mean_pred=cbind(type='GAM',
  codsp=sp,
  yr=newdataexpand$yr,
  mean=NA,
  newdata=newdataexpand,
  pred = round(pred$fit, 2),
  lw = round(pred$lw, 2),
  up = round(pred$up, 2))

#extract mean of all levels
mean_pred = mean_pred %>%
  group_by(type, yr, codsp) %>%
  summarize(across(c('pred', 'lw', 'up'), ~ round(mean(.), 2)), .groups = 'drop') %>%
  #binding the observed mean lenth data and the other factors
  dplyr::left_join(mean_obs_sub[, c("yr", "mean","region","codgr")], by = "yr")

#assign the mean predictions for the species
mean_pred_dol=mean_pred #size mean predicted for all factors combined

#plotting the GAM's prediction and the observed data
#-------------
#plot section
#-------------
#plot GAM mean size prediction
p29<-ggplot(data=mean_pred,aes(x=yr,y=pred))+
  geom_ribbon(data=dplyr::filter(mean_pred,codsp==sp),
    aes(ymin = lw, ymax = up,linetype='Pred'), 
    col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(mean_pred,codsp==sp),
    aes(x=yr,y=pred,linetype='Pred'),col="gray60",size=2)+
  #geom_line(data=dplyr::filter(gam,codsp==sp),aes(x=yr,y=fit,col=codgr,linetype=region),size=1.2)+
  geom_point(data=dplyr::filter(mean_pred,codsp==sp),
    aes(x=yr,y=mean),size=2,col="gray40")+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=yr[1],y=mean[1]),pch=8,size=2,col="gray40")+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='')+
  scale_y_continuous() +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p29

ggsave(paste(c("Length_model_prediction_Small_Tunas_",sp,".png"),collapse=''), plot = p29, device = "png", units = "cm",
  width = 24, height = 14)
#----------------------------------------------

#partial effects
newdata= mean_obs_sub   #desired partial effect

pareffect<- data.frame(predict(gam_model, newdata = newdata,
  type = "terms", se.fit=TRUE)) #terms gets the effects

#intervals for the partial effect
pareffect$up.codgr<-pareffect$fit.codgr+qnorm(0.975)*pareffect$se.fit.codgr #confidence intervals
pareffect$lw.codgr<-pareffect$fit.codgr-qnorm(0.975)*pareffect$se.fit.codgr
pareffect$codgr= mean_obs_sub$codgr

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

pareffect_dol = pareffect

# plot GAM partial effects
p30<-ggplot(dplyr::filter(pareffect,codsp==sp),
  aes(x=effect,y=fit,fill=var,col="Gear"),size=2) +
  #geom_line(aes(x=effect,y=up, col="Gear",group=1)) +
  geom_point(aes(x=effect,y=fit,col="Gear"),size=3.5)+
  geom_errorbar(aes(ymin=lw, ymax=up,col="Gear"), width=0.5,
    position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Effect",x="Gear",col='',fill="")+
  scale_color_viridis_d()+
  scale_y_continuous() +
  scale_x_discrete() +
  #scale_colour_manual(values="gray60") +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p30

ggsave(paste(c("Effects_model_Small_Tunas_",sp,".png"),collapse=''), plot = p30, device = "png", units = "cm",
  width = 24, height = 14)


#----------------------------------
#simulating frequency distributions
#----------------------------------
dist_sim<-distsim(
  mean = mean_pred$pred[mean_pred$codsp==sp],
  cv= cv,
  years= unique(mean_pred$yr[mean_pred$codsp==sp]),
  n=100,
  sp=sp,
  distribution = distribution,
  source = "Dist_Sim") 

#assign the distribution simulated for the species
dist_sim_dol=dist_sim


#-------------
#plot section
#-------------
#Mean predicted (line) + Distribution simulation (boxplots)
p31 <- ggplot() +
  geom_boxplot(data = dplyr::filter(dist_sim, codsp == sp), 
    aes(x = factor(yr), y = fl, fill = codsp, col = codsp),
    alpha = 0.3,
    position = "identity",
    lwd = 0.75,
    outlier.shape = 8,
    outlier.size = 0.8,
    outlier.alpha = 0.5) +
  geom_ribbon(data = dplyr::filter(mean_pred, codsp == sp),
    aes(x = factor(yr), ymin = lw, ymax = up, group = 1, linetype = 'Pred'),
    color = "white", fill = "grey80", alpha = 0.4) +
  geom_line(data = dplyr::filter(mean_pred, codsp == sp),
    aes(x = factor(yr), y = pred, group = 1, linetype = 'Pred'), 
    col = "gray60", size = 2) +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr), y = mean, group = 1), size = 2) +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr[1]), y = mean[1]), pch = 8, size = 2, col = "deepskyblue") +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr[1]), y = mean[1]),  size = 2, col = "deepskyblue") +
  labs(x = "Year", y = "Fork length (cm)", col = "", linetype = '', fill = '') +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  scale_y_continuous() +
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10)))) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p31

ggsave(paste(c("Distribution_Simulation_Boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p31, device = "png", units = "cm",
  width = 24, height = 15)
#-------------------------------------

#plot half violin for the predicted frequencies
p32<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(dist_sim,codsp==sp), 
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
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr),y=mean),size=2)+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr[1]),y=mean[1]),pch=8,size=2)+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr[1]),y=mean[1]),size=2)+
  # #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2025,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p32

ggsave(paste(c("Distribution_Simulation_HalfViolin_Small_Tunas_",sp,".png"),collapse = ''), plot = p32, device = "png", units = "cm",
  width = 35, height = 12)
#--------------------------


#--------------------------------------------------------
#Comparing Distribution simulation with LBSPR simulation
#--------------------------------------------------------
#Life history data is necessary in this section
linf=pars$linf[pars$codsp==sp] #growth
k=pars$k[pars$codsp==sp]   #~k=0.16(Nobrega, 2002) #growth
t0=pars$t0[pars$codsp==sp]  #growth
l50=pars$l50[pars$codsp==sp] #maturity
l95= pars$l50[pars$codsp==sp]*1.1 #maturity
m=pars$m[pars$codsp==sp] ##m=0.36, 
mk=m/k #M/K
by=10 #binwidth


                      #----- Comparing the Simulation Distributions -------#
                      #   Comparing with the LBSPR Simulation function     #
                      #----------------------------------------------------#              
#empty frame to store the simulated frequencies
comp_sim<-data.frame(yr=NULL,codsp=NULL,lmids=NULL,sim1=NULL,sim2=NULL,sel_method=NULL,
  fm_method=NULL,sl_fm_methods=NULL,sl50=NULL,f_m=NULL,res=NULL,rmse=NULL,x2=NULL)  

for (i in years) {
  
      # Sim1 is the distribution simulation from the present study
      sim1<-  lf(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], by=by)
      sim1$counts<- norm(sim1$counts)
  
      # Sim2 is the LBSPR simulation population 
      #Selectivity from Catch-curve  and Selectivity from the LBSPR fit 
      #F/M from the Beverton and Holt and F/M from the LBSPR fit 
  
      #Selectivity from Catch-curve and Selectivity from the Beverton and Holt-------------
      sim2_cc_bh<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=c(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by-6,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        sl95= 1.2*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"), method="Beverton and Holt",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_cc_bh@SL50
      f_m= sim2_cc_bh@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_cc_bh@LMids,counts=sim2_cc_bh@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_cc_bh<- data.frame(mids=sim2_cc_bh@LMids,counts=norm(sim2_cc_bh@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_cc_bh <- subset(sim2_cc_bh, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_cc_bh$counts,
        sel_method="CC", fm_method="BH",sl_fm_methods= paste("CC_BH"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_cc_bh$counts,
        rmse=rmse(sim1$counts,sim2_cc_bh$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
            sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
      comp_sim<-rbind(comp_sim,sim_dat)
  
      #Selectivity from Catch-curve and F/M from the LBSPR fit-------------
      sim2_cc_lb<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        sl95= 1.2*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"), method="LBSPR",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_cc_lb@SL50
      f_m= sim2_cc_lb@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_cc_lb@LMids,counts=sim2_cc_lb@pLCatch*100) %>% 
      tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_cc_lb<- data.frame(mids=sim2_cc_lb@LMids,counts=norm(sim2_cc_lb@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_cc_lb <- subset(sim2_cc_lb, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_cc_lb$counts,
        sel_method="CC", fm_method="LBSPR",sl_fm_methods= paste("CC_LB"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_cc_lb$counts,
        rmse=rmse(sim1$counts,sim2_cc_lb$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
            sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
      comp_sim<-rbind(comp_sim,sim_dat)
  
  
      #Selectivity from LBSPR and F/M from Beverton and Holt------------
      sim2_lb_bh<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        sl95= 1.2*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),method="Beverton and Holt",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_lb_bh@SL50
      f_m= sim2_lb_bh@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_lb_bh@LMids,counts=sim2_lb_bh@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_lb_bh<- data.frame(mids=sim2_lb_bh@LMids,counts=norm(sim2_lb_bh@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_lb_bh <- subset(sim2_lb_bh, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_lb_bh$counts,
        sel_method="LB", fm_method="BH",sl_fm_methods= paste("LB_BH"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_lb_bh$counts,
        rmse=rmse(sim1$counts,sim2_lb_bh$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
            sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
      comp_sim<-rbind(comp_sim,sim_dat)
  
  
      #Selectivity from LBSPR and F/M from LBSPR -----------
      sim2_lb_lb<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        sl95= 1.2*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"), method="LBSPR",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_lb_lb@SL50
      f_m= sim2_lb_lb@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_lb_lb@LMids,counts=sim2_lb_lb@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_lb_lb<- data.frame(mids=sim2_lb_lb@LMids,counts=norm(sim2_lb_lb@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_lb_lb <- subset(sim2_lb_lb, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_lb_lb$counts,
        sel_method="LB", fm_method="LB",sl_fm_methods= paste("LB_LB"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_lb_lb$counts,
        rmse=rmse(sim1$counts,sim2_lb_lb$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
          sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
  comp_sim<-rbind(comp_sim,sim_dat)
}
#assing the comparisons for the current species
comp_sim_dol<-comp_sim

#-------------
#plot section
#-------------
#Length residual comparison 
p33 <- ggplot() +
  geom_boxplot(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(group = interaction(factor(yr), factor(sl_fm_methods)), 
      x = factor(yr), 
      y = res, 
      fill = factor(sl_fm_methods)),  # fill to methods
    alpha = 1,
    position = position_dodge(width = 0.7), 
    width = 0.9,
    lwd = 0.1,  
    outlier.shape = 8,
    outlier.size = 0.7,
    outlier.alpha = 0.5) +
  geom_hline(yintercept =0,linetype = "dashed" )+
  geom_smooth(data = dplyr::filter(comp_sim, codsp == sp),
    method = "loess", span = 0.2, se = FALSE,
    aes(group = 1, x = factor(yr), y = res,color="Loess"), size = 1.5) +
  geom_text(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = "2015", y = 0.9 * max(res, na.rm = TRUE), 
      label = paste("RMSE=", round(mean(rmse, na.rm = TRUE), 2))),
    inherit.aes = FALSE) +
  labs(x = "Year", y = "Length residual", col = "", linetype = '', fill = '') +
  scale_y_continuous() +
  scale_fill_viridis_d(direction = -1)+
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10)))) +
  theme_classic(base_size = 14) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p33

ggsave(paste(c("Length_Residual_Boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p33, device = "png", units = "cm",
  width = 25, height = 17)
#------------------------

#Scaterplot comparing the distributions (Dist sim vs LBSPR sim)
p34 <- ggplot() +
  geom_point(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = sim2, y = sim1, color = sl_fm_methods), 
    size = 2, alpha = 0.7) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", size = 1) +
  geom_smooth(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = sim2, y = sim1, color = sl_fm_methods, fill = sl_fm_methods), 
    method = "loess", se = FALSE, size = 1.2, alpha = 0.3) +  
  scale_color_viridis_d(direction = -1) +  
  labs(x = "LBSPR Distribution Simulation",y= "Distribution Simulation",col = "Method", shape = "Method", fill = "Method") +
  theme_classic(base_size = 14) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p34

ggsave(paste(c("Distributions_Comparison_Small_Tunas_",sp,".png"),collapse = ''), plot = p34, device = "png", units = "cm",
  width = 24, height = 15)
#---------------------------------------------------------------------------------------------------------------------------




#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
# models by species
# FRI model - Length Reconstruction 
#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
sp='FRI'
years=1950:2025

#-------------------------------------
#find the best distribution for specie
#-------------------------------------
dist<-distfit(len = data.frame(
  yr=freq_obs$yr[freq_obs$codsp==sp],
  codsp=freq_obs$codsp[freq_obs$codsp==sp],
  codgr=freq_obs$codgr[freq_obs$codsp==sp],
  fl=freq_obs$fl[freq_obs$codsp==sp]))

dist=dist %>%
  dplyr::group_by(yr,codsp,codgr) %>% 
  slice_min(aic, n = 1)  %>%
  dplyr::arrange(dist) 

#assign dist for the current species
dist_fri=dist

#Table with the best distribution
table(dist$codsp,dist$dist)  
distribution='norm'

#coefficient of variation for the corresponding distribution 
cv=mean(dist$cv[dist$codsp==sp & dist$dist==distribution])

#-------------------------------------------------------
#predict 1950 length when using the slope of the lengths
#-------------------------------------------------------

#assuming the 1950 mean length as the prediction of a linear model
lm_data<- data.frame(length=mean_obs$mean[mean_obs$codsp==sp],year=mean_obs$yr[mean_obs$codsp==sp])

#fiting the model
lm_model<- lm(length~year, data = lm_data[-1,])
#Parameters... On average the mean length is being reduced 0.10 cm year by year
a<- coef(lm_model)[1]
b<- coef(lm_model)[2]

#using the fitted parameters to estimate what would be the length in 1950
newdata<- data.frame(year=1950:2025)
lm_predict_10<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.1)))
lm_predict_20<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.2)))
lm_predict_30<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.3)))

#data frame with predictions for different percentile levels
lm_predict<- data.frame(yr=newdata$year,
  fit=lm_predict_10$fit.fit, 
  lw_10=lm_predict_10$fit.lwr,
  up_10=lm_predict_10$fit.upr,
  lw_20=lm_predict_20$fit.lwr,
  up_20=lm_predict_20$fit.upr,
  lw_30=lm_predict_30$fit.lwr,
  up_30=lm_predict_30$fit.upr,
  codsp=sp)

#lm model coefficients and predicitions for the current species
lm_model_fri<- data.frame(formula=paste(deparse(lm_model$call), collapse = ""),
  intercept=a,
  slope=b,
  codsp=sp)
lm_predict_fri<- lm_predict

# controlling variables
lm_predict <- lm_predict %>%
  mutate(confidence_level = case_when(
    yr == 1950 & lw_10 == lw_10 ~ "10%",
    yr == 1950 & lw_20 == lw_20 ~ "20%",
    yr == 1950 & lw_30 == lw_30 ~ "30%"
  ))

# visual representation of the linear model
#-------------
#plot section
#-------------
p35 <- ggplot(lm_predict, aes(x = yr, y = fit)) +
  geom_line(size = 1, color = "grey20") +  # main line prediction
  # Lines for the percentiles
  geom_line(aes(y = lw_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = up_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = lw_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = up_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = lw_30, color = "30%", linetype = "30%"), size = 0.8) +
  geom_line(aes(y = up_30, color = "30%", linetype = "30%"), size = 0.8) +
  # Observed data
  geom_point(data = lm_data, aes(x = year, y = length), size = 2, color = "grey30") +
  # Predicted points for 1950 
  geom_point(aes(x = yr[1], y = lw_10[1], fill = "10%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_10[1], fill = "10%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = lw_20[1], fill = "20%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_20[1], fill = "20%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = lw_30[1], fill = "30%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_30[1], fill = "30%"), size = 3, shape = 21) +
  labs(y = 'Mean Length (cm)', x = "Year", color = "Percentile", fill = "Percentile", linetype = "Percentile") +
  # Type and line colors
  scale_color_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_fill_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_linetype_manual(values = c("10%" = "dashed", "20%" = "dotted", "30%" = "dotdash")) +
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p35
#saving graph...
ggsave(paste0("Linear_model_prediction_Small_Tunas_", sp, ".png"), plot = p35, device = "png", units = "cm", width = 24, height = 16)
#----------------------------------

#selecting variables .. subseting
mean_obs_sub= mean_obs[,c("yr",'region','codsp','codgr','mean')]

#removing extreme values
mean_obs_sub$mean[mean_obs_sub$codsp=='BLF' & mean_obs_sub$mean<50] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='DOL' & mean_obs_sub$mean<75] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='KGM' & mean_obs_sub$mean>100] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='KGM' & mean_obs_sub$mean<60] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean<100] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean>180] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='FRI' & mean_obs_sub$mean>45] =NA

#get only non-NA values and selecting the current species
mean_obs_sub=mean_obs_sub[is.na(mean_obs_sub$mean)==FALSE & mean_obs_sub$codsp==sp,]

#grouping factors most likely (Fredou et al., 2021)
mean_obs_sub$region[mean_obs_sub$region=='S']='SE'

#transforming columns in to factors
mean_obs_sub$region= factor(mean_obs_sub$region)
mean_obs_sub$codgr= factor(mean_obs_sub$codgr)
mean_obs_sub$codsp= factor(mean_obs_sub$codsp)


# model GAM selection----------------
#list of values to be tested
tests<- names(lm_predict[, grepl("fit|lw|up", names(lm_predict))])
#filling an empty 1950 with NA 
mean_obs_sub<- rbind(data.frame(yr=1950, region="NE", codsp=sp,codgr="UN",mean=NA),
  mean_obs_sub)

# Definindo os vetores de variáveis de resposta e explicativas
resp_vars <- c("mean")
exp_vars_1 <- c("s(yr)")
exp_vars_2 <- c("region")
exp_vars_3 <- c("codgr")

#criando a formula geral
formulas <-list(as.formula(paste(resp_vars, "~", paste(c(exp_vars_1), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_2), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_3), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_2,exp_vars_3), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_3,exp_vars_2), collapse = "+"))))

formulas<- as.character(formulas)

# Data frame com as distribuicoes de probabilidade e funcoes de ligacao a serem testadas
family=distribution #testing only for the best distribution
link = ("log") #log fixed

#empty frame to store results
gam=data.frame(
  test=NULL,
  type=NULL,
  codsp=NULL,
  family=NULL,
  link=NULL,
  formula=NULL,
  newdataexpand=NULL,
  pred=NULL,
  aic=NULL,
  r.sq=NULL,
  dev.ex=NULL)


#-------------------- Fitting models  -------------------#
for (i in tests) {
  
  #1- first sensibility for different 1950 starting values--
  test_1950<- lm_predict[1, grepl(i, names(lm_predict))]
  
  #replacing with the current 1950 mean length test
  mean_obs_sub$mean[mean_obs_sub$yr==1950]<- test_1950
  
  #2- Term selection within the sensibility----
  # Loop through the formulas 
  for (j in formulas) {
    
    family_name <- as.character(family) #these are fixed
    link_name <- as.character(link)
    
    #Creating switch object
    family_obj <- switch(family_name,
      norm = gaussian(link = link_name),  
      invgauss = inverse.gaussian(link = link_name),
      gamma= Gamma(link=link_name))
    
    #formula j to be tested
    formula <- as.formula(j)
    
    #fit the model
    gam_model <- mgcv::gam(as.formula(formula),
      data= mean_obs_sub,
      method = "REML" ,
      family =family_obj)
    
    #taking the outputs
    gam=rbind(gam, data.frame(
      test=i,
      type="GAM",
      codsp=sp,
      family=gam_model$family$family,
      link=gam_model$family$link,
      formula=paste0(c(paste(gam_model$formula)[2],
        paste(gam_model$formula)[1],
        paste(gam_model$formula)[3]),collapse = ' '),
      aic=gam_model$aic,
      r.sq=summary(gam_model)$r.sq,
      dev.ex=summary(gam_model)$dev.expl))
  }
}
#-------------------- end of loop -----------------------#

#assign the tests to the species
gam_fri<-gam

#replacing the 1950 mean length with the best tested value
best_model<- gam[which.min(gam$aic),]; dplyr::tibble(best_model)
mean_obs_sub$mean[mean_obs_sub$yr==1950]<- lm_predict[1, grepl(best_model$test, names(lm_predict))]

#Final model (lowest AIC) 
gam_model<- mgcv::gam(as.formula(best_model$formula), 
  #sp=0.01,
  data= mean_obs_sub,
  method = "REML" ,
  family =family_obj)

#assign GAM model to the specie
gam_model_fri=gam_model

#diagnostic
summary(gam_model)

par(mfrow=c(2,2))
mgcv::gam.check(gam_model)# model check

#residual 
#homocedasticity
lmtest::bptest(gam_model)

#normality
shapiro.test(resid(gam_model))

# Check overall concurvity
round(mgcv::concurvity(gam_model, full = TRUE),2)
#----------------------
#predictions "year"
newdataexpand= expand.grid(yr=seq(years[1],years[length(years)],1),
  region=mean_obs_sub$region,
  codgr=mean_obs_sub$codgr)

pred<- data.frame(predict(gam_model, newdata = newdataexpand,
  type = "response", se.fit=TRUE))

pred$up<-pred$fit+qnorm(0.975)*pred$se.fit #confidence intervals
pred$lw<-pred$fit-qnorm(0.975)*pred$se.fit

#take the predictions
mean_pred=cbind(type='GAM',
  codsp=sp,
  yr=newdataexpand$yr,
  mean=NA,
  newdata=newdataexpand,
  pred = round(pred$fit, 2),
  lw = round(pred$lw, 2),
  up = round(pred$up, 2))

#extract mean of all levels
mean_pred = mean_pred %>%
  group_by(type, yr, codsp) %>%
  summarize(across(c('pred', 'lw', 'up'), ~ round(mean(.), 2)), .groups = 'drop') %>%
  #binding the observed mean lenth data and the other factors
  dplyr::left_join(mean_obs_sub[, c("yr", "mean","region","codgr")], by = "yr")

#assign the mean predictions for the species
mean_pred_fri=mean_pred #size mean predicted for all factors combined

#plotting the GAM's prediction and the observed data
#-------------
#plot section
#-------------
#plot GAM mean size prediction
p36<-ggplot(data=mean_pred,aes(x=yr,y=pred))+
  geom_ribbon(data=dplyr::filter(mean_pred,codsp==sp),
    aes(ymin = lw, ymax = up,linetype='Pred'), 
    col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(mean_pred,codsp==sp),
    aes(x=yr,y=pred,linetype='Pred'),col="gray60",size=2)+
  #geom_line(data=dplyr::filter(gam,codsp==sp),aes(x=yr,y=fit,col=codgr,linetype=region),size=1.2)+
  geom_point(data=dplyr::filter(mean_pred,codsp==sp),
    aes(x=yr,y=mean),size=2,col="gray40")+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=yr[1],y=mean[1]),pch=8,size=2,col="gray40")+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='')+
  scale_y_continuous() +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p36

ggsave(paste(c("Length_model_prediction_Small_Tunas_",sp,".png"),collapse=''), plot = p36, device = "png", units = "cm",
  width = 24, height = 14)
#----------------------------------------------

#partial effects
newdata= mean_obs_sub   #desired partial effect

pareffect<- data.frame(predict(gam_model, newdata = newdata,
  type = "terms", se.fit=TRUE)) #terms gets the effects

#intervals for the partial effect
pareffect$up.codgr<-pareffect$fit.codgr+qnorm(0.975)*pareffect$se.fit.codgr #confidence intervals
pareffect$lw.codgr<-pareffect$fit.codgr-qnorm(0.975)*pareffect$se.fit.codgr
pareffect$codgr= mean_obs_sub$codgr

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

pareffect_fri = pareffect

# plot GAM partial effects
p37<-ggplot(dplyr::filter(pareffect,codsp==sp),
  aes(x=effect,y=fit,fill=var,col="Gear"),size=2) +
  #geom_line(aes(x=effect,y=up, col="Gear",group=1)) +
  geom_point(aes(x=effect,y=fit,col="Gear"),size=3.5)+
  geom_errorbar(aes(ymin=lw, ymax=up,col="Gear"), width=0.5,
    position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Effect",x="Gear",col='',fill="")+
  scale_color_viridis_d()+
  scale_y_continuous() +
  scale_x_discrete() +
  #scale_colour_manual(values="gray60") +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p37

ggsave(paste(c("Effects_model_Small_Tunas_",sp,".png"),collapse=''), plot = p37, device = "png", units = "cm",
  width = 24, height = 14)


#----------------------------------
#simulating frequency distributions
#----------------------------------
dist_sim<-distsim(
  mean = mean_pred$pred[mean_pred$codsp==sp],
  cv= cv,
  years= unique(mean_pred$yr[mean_pred$codsp==sp]),
  n=100,
  sp=sp,
  distribution = distribution,
  source = "Dist_Sim") 

#assign the distribution simulated for the species
dist_sim_fri=dist_sim


#-------------
#plot section
#-------------
#Mean predicted (line) + Distribution simulation (boxplots)
p38 <- ggplot() +
  geom_boxplot(data = dplyr::filter(dist_sim, codsp == sp), 
    aes(x = factor(yr), y = fl, fill = codsp, col = codsp),
    alpha = 0.3,
    position = "identity",
    lwd = 0.75,
    outlier.shape = 8,
    outlier.size = 0.8,
    outlier.alpha = 0.5) +
  geom_ribbon(data = dplyr::filter(mean_pred, codsp == sp),
    aes(x = factor(yr), ymin = lw, ymax = up, group = 1, linetype = 'Pred'),
    color = "white", fill = "grey80", alpha = 0.4) +
  geom_line(data = dplyr::filter(mean_pred, codsp == sp),
    aes(x = factor(yr), y = pred, group = 1, linetype = 'Pred'), 
    col = "gray60", size = 2) +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr), y = mean, group = 1), size = 2) +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr[1]), y = mean[1]), pch = 8, size = 2, col = "deepskyblue") +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr[1]), y = mean[1]),  size = 2, col = "deepskyblue") +
  labs(x = "Year", y = "Fork length (cm)", col = "", linetype = '', fill = '') +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  scale_y_continuous() +
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10)))) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p38

ggsave(paste(c("Distribution_Simulation_Boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p38, device = "png", units = "cm",
  width = 24, height = 15)
#-------------------------------------

#plot half violin for the predicted frequencies
p39<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(dist_sim,codsp==sp), 
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
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr),y=mean),size=2)+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr[1]),y=mean[1]),pch=8,size=2)+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr[1]),y=mean[1]),size=2)+
  # #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2025,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p39

ggsave(paste(c("Distribution_Simulation_HalfViolin_Small_Tunas_",sp,".png"),collapse = ''), plot = p39, device = "png", units = "cm",
  width = 35, height = 12)
#--------------------------


#--------------------------------------------------------
#Comparing Distribution simulation with LBSPR simulation
#--------------------------------------------------------
#Life history data is necessary in this section
linf=pars$linf[pars$codsp==sp] #growth
k=pars$k[pars$codsp==sp]   #~k=0.16(Nobrega, 2002) #growth
t0=pars$t0[pars$codsp==sp]  #growth
l50=pars$l50[pars$codsp==sp] #maturity
l95= pars$l50[pars$codsp==sp]*1.1 #maturity
m=pars$m[pars$codsp==sp] # 
mk=m/k #M/K
by=5 #binwidth


                        #----- Comparing the Simulation Distributions -------#
                        #   Comparing with the LBSPR Simulation function     #
                        #----------------------------------------------------#              
#empty frame to store the simulated frequencies
comp_sim<-data.frame(yr=NULL,codsp=NULL,lmids=NULL,sim1=NULL,sim2=NULL,sel_method=NULL,
  fm_method=NULL,sl_fm_methods=NULL,sl50=NULL,f_m=NULL,res=NULL,rmse=NULL,x2=NULL)  

for (i in years) {
  
      # Sim1 is the distribution simulation from the present study
      sim1<-  lf(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], by=by)
      sim1$counts<- norm(sim1$counts)
  
      # Sim2 is the LBSPR simulation population 
      #Selectivity from Catch-curve  and Selectivity from the LBSPR fit 
      #F/M from the Beverton and Holt and F/M from the LBSPR fit 
  
      #Selectivity from Catch-curve and Selectivity from the Beverton and Holt-------------
      sim2_cc_bh<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=c(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by-6,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        sl95= 1.1*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"), method="Beverton and Holt",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_cc_bh@SL50
      f_m= sim2_cc_bh@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_cc_bh@LMids,counts=sim2_cc_bh@pLCatch*100) %>% 
       tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_cc_bh<- data.frame(mids=sim2_cc_bh@LMids,counts=norm(sim2_cc_bh@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_cc_bh <- subset(sim2_cc_bh, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_cc_bh$counts,
        sel_method="CC", fm_method="BH",sl_fm_methods= paste("CC_BH"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_cc_bh$counts,
        rmse=rmse(sim1$counts,sim2_cc_bh$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
           sample(freq_lbspr, 95,replace = FALSE))$p.value))
    
      comp_sim<-rbind(comp_sim,sim_dat)
  
      #Selectivity from Catch-curve and F/M from the LBSPR fit-------------
      sim2_cc_lb<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        sl95= 1.1*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"), method="LBSPR",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_cc_lb@SL50
      f_m= sim2_cc_lb@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_cc_lb@LMids,counts=sim2_cc_lb@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_cc_lb<- data.frame(mids=sim2_cc_lb@LMids,counts=norm(sim2_cc_lb@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_cc_lb <- subset(sim2_cc_lb, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_cc_lb$counts,
        sel_method="CC", fm_method="LBSPR",sl_fm_methods= paste("CC_LB"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_cc_lb$counts,
        rmse=rmse(sim1$counts,sim2_cc_lb$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
            sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
      comp_sim<-rbind(comp_sim,sim_dat)
  
  
      #Selectivity from LBSPR and F/M from Beverton and Holt------------
      sim2_lb_bh<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        sl95= 1.1*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),method="Beverton and Holt",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_lb_bh@SL50
      f_m= sim2_lb_bh@FM
      
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_lb_bh@LMids,counts=sim2_lb_bh@pLCatch*100) %>% 
       tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_lb_bh<- data.frame(mids=sim2_lb_bh@LMids,counts=norm(sim2_lb_bh@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_lb_bh <- subset(sim2_lb_bh, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_lb_bh$counts,
        sel_method="LB", fm_method="BH",sl_fm_methods= paste("LB_BH"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_lb_bh$counts,
        rmse=rmse(sim1$counts,sim2_lb_bh$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
            sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
      comp_sim<-rbind(comp_sim,sim_dat)
  
  
      #Selectivity from LBSPR and F/M from LBSPR -----------
      sim2_lb_lb<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        sl95= 1.1*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"), method="LBSPR",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_lb_lb@SL50
      f_m= sim2_lb_lb@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_lb_lb@LMids,counts=sim2_lb_lb@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_lb_lb<- data.frame(mids=sim2_lb_lb@LMids,counts=norm(sim2_lb_lb@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_lb_lb <- subset(sim2_lb_lb, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_lb_lb$counts,
        sel_method="LB", fm_method="LB",sl_fm_methods= paste("LB_LB"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_lb_lb$counts,
        rmse=rmse(sim1$counts,sim2_lb_lb$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
           sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
  comp_sim<-rbind(comp_sim,sim_dat)
}
#assing the comparisons for the current species
comp_sim_fri<-comp_sim

#-------------
#plot section
#-------------
#Length residual comparison 
p40 <- ggplot() +
  geom_boxplot(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(group = interaction(factor(yr), factor(sl_fm_methods)), 
      x = factor(yr), 
      y = res, 
      fill = factor(sl_fm_methods)),  # fill to methods
    alpha = 1,
    position = position_dodge(width = 0.7), 
    width = 0.9,
    lwd = 0.1,  
    outlier.shape = 8,
    outlier.size = 0.7,
    outlier.alpha = 0.5) +
  geom_hline(yintercept =0,linetype = "dashed" )+
  geom_smooth(data = dplyr::filter(comp_sim, codsp == sp),
    method = "loess", span = 0.2, se = FALSE,
    aes(group = 1, x = factor(yr), y = res,color="Loess"), size = 1.5) +
  geom_text(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = "2015", y = 0.9 * max(res, na.rm = TRUE), 
      label = paste("RMSE=", round(mean(rmse, na.rm = TRUE), 2))),
    inherit.aes = FALSE) +
  labs(x = "Year", y = "Length residual", col = "", linetype = '', fill = '') +
  scale_y_continuous() +
  scale_fill_viridis_d(direction = -1)+
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10)))) +
  theme_classic(base_size = 14) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p40

ggsave(paste(c("Length_Residual_Boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p40, device = "png", units = "cm",
  width = 25, height = 17)
#------------------------

#Scaterplot comparing the distributions (Dist sim vs LBSPR sim)
p41 <- ggplot() +
  geom_point(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = sim2, y = sim1, color = sl_fm_methods), 
    size = 2, alpha = 0.7) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", size = 1) +
  geom_smooth(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = sim2, y = sim1, color = sl_fm_methods, fill = sl_fm_methods), 
    method = "loess", se = FALSE, size = 1.2, alpha = 0.3) +  
  scale_color_viridis_d(direction = -1) +  
  labs(x = "LBSPR Distribution Simulation",y= "Distribution Simulation",col = "Method", shape = "Method", fill = "Method") +
  theme_classic(base_size = 14) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p41

ggsave(paste(c("Distributions_Comparison_Small_Tunas_",sp,".png"),collapse = ''), plot = p41, device = "png", units = "cm",
  width = 24, height = 15)
#---------------------------------------------------------------------------------------------------------------------------





#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
# models by species
# KGM model - Length Reconstruction 
#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
sp='KGM'
years=1950:2025

#-------------------------------------
#find the best distribution for specie
#-------------------------------------
dist<-distfit(len = data.frame(
  yr=freq_obs$yr[freq_obs$codsp==sp],
  codsp=freq_obs$codsp[freq_obs$codsp==sp],
  codgr=freq_obs$codgr[freq_obs$codsp==sp],
  fl=freq_obs$fl[freq_obs$codsp==sp]))

dist=dist %>%
  dplyr::group_by(yr,codsp,codgr) %>% 
  slice_min(aic, n = 1)  %>%
  dplyr::arrange(dist) 

#assign dist for the current species
dist_kgm=dist

#Table with the best distribution
table(dist$codsp,dist$dist)  
distribution='invgauss'

#coefficient of variation for the corresponding distribution 
cv=mean(dist$cv[dist$codsp==sp & dist$dist==distribution])

#-------------------------------------------------------
#predict 1950 length when using the slope of the lengths
#-------------------------------------------------------

#assuming the 1950 mean length as the prediction of a linear model
lm_data<- data.frame(length=mean_obs$mean[mean_obs$codsp==sp],year=mean_obs$yr[mean_obs$codsp==sp])
lm_data<- lm_data[lm_data$length<93 & lm_data$length>60 & lm_data$year!=1969,]

#fiting the model
lm_model<- lm(length~year, data = lm_data[-1,])
#Parameters... On average the mean length is being reduced 0.10 cm year by year
a<- coef(lm_model)[1]
b<- coef(lm_model)[2]

#using the fitted parameters to estimate what would be the length in 1950
newdata<- data.frame(year=1950:2025)
lm_predict_10<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.1)))
lm_predict_20<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.2)))
lm_predict_30<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.3)))

#data frame with predictions for different percentile levels
lm_predict<- data.frame(yr=newdata$year,
  fit=lm_predict_10$fit.fit, 
  lw_10=lm_predict_10$fit.lwr,
  up_10=lm_predict_10$fit.upr,
  lw_20=lm_predict_20$fit.lwr,
  up_20=lm_predict_20$fit.upr,
  lw_30=lm_predict_30$fit.lwr,
  up_30=lm_predict_30$fit.upr,
  codsp=sp)

#lm model coefficients and predicitions for the current species
lm_model_kgm<- data.frame(formula=paste(deparse(lm_model$call), collapse = ""),
  intercept=a,
  slope=b,
  codsp=sp)
lm_predict_kgm<- lm_predict

# controlling variables
lm_predict <- lm_predict %>%
  mutate(confidence_level = case_when(
    yr == 1950 & lw_10 == lw_10 ~ "10%",
    yr == 1950 & lw_20 == lw_20 ~ "20%",
    yr == 1950 & lw_30 == lw_30 ~ "30%"
  ))

# visual representation of the linear model
#-------------
#plot section
#-------------
p42 <- ggplot(lm_predict, aes(x = yr, y = fit)) +
  geom_line(size = 1, color = "grey20") +  # main line prediction
  # Lines for the percentiles
  geom_line(aes(y = lw_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = up_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = lw_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = up_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = lw_30, color = "30%", linetype = "30%"), size = 0.8) +
  geom_line(aes(y = up_30, color = "30%", linetype = "30%"), size = 0.8) +
  # Observed data
  geom_point(data = lm_data, aes(x = year, y = length), size = 2, color = "grey30") +
  # Predicted points for 1950 
  geom_point(aes(x = yr[1], y = lw_10[1], fill = "10%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_10[1], fill = "10%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = lw_20[1], fill = "20%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_20[1], fill = "20%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = lw_30[1], fill = "30%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_30[1], fill = "30%"), size = 3, shape = 21) +
  labs(y = 'Mean Length (cm)', x = "Year", color = "Percentile", fill = "Percentile", linetype = "Percentile") +
  # Type and line colors
  scale_color_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_fill_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_linetype_manual(values = c("10%" = "dashed", "20%" = "dotted", "30%" = "dotdash")) +
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p42
#saving graph...
ggsave(paste0("Linear_model_prediction_Small_Tunas_", sp, ".png"), plot = p42, device = "png", units = "cm", width = 24, height = 16)
#----------------------------------

#selecting variables .. subseting
mean_obs_sub= mean_obs[,c("yr",'region','codsp','codgr','mean')]

#removing extreme values
mean_obs_sub$mean[mean_obs_sub$codsp=='BLF' & mean_obs_sub$mean<50] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='DOL' & mean_obs_sub$mean<75] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='KGM' & mean_obs_sub$mean>100] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='KGM' & mean_obs_sub$mean<60] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean<100] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean>180] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='FRI' & mean_obs_sub$mean>45] =NA

#get only non-NA values and selecting the current species
mean_obs_sub=mean_obs_sub[is.na(mean_obs_sub$mean)==FALSE & mean_obs_sub$codsp==sp,]

#grouping factors most likely (Fredou et al., 2021)
#mean_obs_sub$region[mean_obs_sub$region=='S']='SE'

#transforming columns in to factors
mean_obs_sub$region= factor(mean_obs_sub$region)
mean_obs_sub$codgr= factor(mean_obs_sub$codgr)
mean_obs_sub$codsp= factor(mean_obs_sub$codsp)


# model GAM selection----------------
#list of values to be tested
tests<- names(lm_predict[, grepl("fit|lw|up", names(lm_predict))])
#filling an empty 1950 with NA 
mean_obs_sub<- rbind(data.frame(yr=1950, region="NE", codsp=sp,codgr="UN",mean=NA),
  mean_obs_sub)

# Definindo os vetores de variáveis de resposta e explicativas
resp_vars <- c("mean")
exp_vars_1 <- c("s(yr)")
exp_vars_2 <- c("region")
exp_vars_3 <- c("codgr")

#criando a formula geral
formulas <-list(as.formula(paste(resp_vars, "~", paste(c(exp_vars_1), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_2), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_3), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_2,exp_vars_3), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_3,exp_vars_2), collapse = "+"))))

formulas<- as.character(formulas)

# Data frame com as distribuicoes de probabilidade e funcoes de ligacao a serem testadas
family=distribution #testing only for the best distribution
link = ("log") #log fixed

#empty frame to store results
gam=data.frame(
  test=NULL,
  type=NULL,
  codsp=NULL,
  family=NULL,
  link=NULL,
  formula=NULL,
  newdataexpand=NULL,
  pred=NULL,
  aic=NULL,
  r.sq=NULL,
  dev.ex=NULL)


#-------------------- Fitting models  -------------------#
for (i in tests) {
  
  #1- first sensibility for different 1950 starting values--
  test_1950<- lm_predict[1, grepl(i, names(lm_predict))]
  
  #replacing with the current 1950 mean length test
  mean_obs_sub$mean[mean_obs_sub$yr==1950]<- test_1950
  
  #2- Term selection within the sensibility----
  # Loop through the formulas 
  for (j in formulas) {
    
    family_name <- as.character(family) #these are fixed
    link_name <- as.character(link)
    
    #Creating switch object
    family_obj <- switch(family_name,
      norm = gaussian(link = link_name),  
      invgauss = inverse.gaussian(link = link_name),
      gamma= Gamma(link=link_name))
    
    #formula j to be tested
    formula <- as.formula(j)
    
    #fit the model
    gam_model <- mgcv::gam(as.formula(formula),
      data= mean_obs_sub,
      method = "REML" ,
      family =family_obj)
    
    #taking the outputs
    gam=rbind(gam, data.frame(
      test=i,
      type="GAM",
      codsp=sp,
      family=gam_model$family$family,
      link=gam_model$family$link,
      formula=paste0(c(paste(gam_model$formula)[2],
        paste(gam_model$formula)[1],
        paste(gam_model$formula)[3]),collapse = ' '),
      aic=gam_model$aic,
      r.sq=summary(gam_model)$r.sq,
      dev.ex=summary(gam_model)$dev.expl))
  }
}
#-------------------- end of loop -----------------------#

#assign the tests to the species
gam_kgm<-gam

#replacing the 1950 mean length with the best tested value
best_model<- gam[which.min(gam$aic),]; dplyr::tibble(best_model)
mean_obs_sub$mean[mean_obs_sub$yr==1950]<-95 

#Final model (lowest AIC) 
gam_model<- mgcv::gam(as.formula(best_model$formula), 
  sp=0.001,
  data= mean_obs_sub,
  method = "REML" ,
  family =family_obj)

#assign GAM model to the specie
gam_model_kgm=gam_model

#diagnostic
summary(gam_model)

par(mfrow=c(2,2))
mgcv::gam.check(gam_model)# model check

#residual 
#homocedasticity
lmtest::bptest(gam_model)

#normality
shapiro.test(resid(gam_model))

# Check overall concurvity
round(mgcv::concurvity(gam_model, full = TRUE),2)
#----------------------
#predictions "year"
newdataexpand= expand.grid(yr=seq(years[1],years[length(years)],1),
  region=mean_obs_sub$region,
  codgr=mean_obs_sub$codgr)

pred<- data.frame(predict(gam_model, newdata = newdataexpand,
  type = "response", se.fit=TRUE))

pred$up<-pred$fit+qnorm(0.975)*pred$se.fit #confidence intervals
pred$lw<-pred$fit-qnorm(0.975)*pred$se.fit

#take the predictions
mean_pred=cbind(type='GAM',
  codsp=sp,
  yr=newdataexpand$yr,
  mean=NA,
  newdata=newdataexpand,
  pred = round(pred$fit, 2),
  lw = round(pred$lw, 2),
  up = round(pred$up, 2))

#extract mean of all levels
mean_pred = mean_pred %>%
  group_by(type, yr, codsp) %>%
  summarize(across(c('pred', 'lw', 'up'), ~ round(mean(.), 2)), .groups = 'drop') %>%
  #binding the observed mean lenth data and the other factors
  dplyr::left_join(mean_obs_sub[, c("yr", "mean","region","codgr")], by = "yr")

#assign the mean predictions for the species
mean_pred_kgm=mean_pred #size mean predicted for all factors combined

#plotting the GAM's prediction and the observed data
#-------------
#plot section
#-------------
#plot GAM mean size prediction
p43<-ggplot(data=mean_pred,aes(x=yr,y=pred))+
  geom_ribbon(data=dplyr::filter(mean_pred,codsp==sp),
    aes(ymin = lw, ymax = up,linetype='Pred'), 
    col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(mean_pred,codsp==sp),
    aes(x=yr,y=pred,linetype='Pred'),col="gray60",size=2)+
  #geom_line(data=dplyr::filter(gam,codsp==sp),aes(x=yr,y=fit,col=codgr,linetype=region),size=1.2)+
  geom_point(data=dplyr::filter(mean_pred,codsp==sp),
    aes(x=yr,y=mean),size=2,col="gray40")+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=yr[1],y=mean[1]),pch=8,size=2,col="gray40")+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='')+
  scale_y_continuous() +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p43

ggsave(paste(c("Length_model_prediction_Small_Tunas_",sp,".png"),collapse=''), plot = p43, device = "png", units = "cm",
  width = 24, height = 14)
#----------------------------------------------

#partial effects
newdata= mean_obs_sub   #desired partial effect

pareffect<- data.frame(predict(gam_model, newdata = newdata,
  type = "terms", se.fit=TRUE)) #terms gets the effects

#intervals for the partial effect
pareffect$up.codgr<-pareffect$fit.codgr+qnorm(0.975)*pareffect$se.fit.codgr #confidence intervals
pareffect$lw.codgr<-pareffect$fit.codgr-qnorm(0.975)*pareffect$se.fit.codgr
pareffect$codgr= mean_obs_sub$codgr


pareffect$up.region<-pareffect$fit.region+qnorm(0.975)*pareffect$se.fit.region #confidence intervals
pareffect$lw.region<-pareffect$fit.region-qnorm(0.975)*pareffect$se.fit.region
pareffect$region= mean_obs_sub$region

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


# plot GAM partial effects
p44<-ggplot(dplyr::filter(pareffect,codsp==sp),
  aes(x=effect,y=fit,fill=var,col=var),size=2) +
  #geom_line(aes(x=effect,y=up, col="Gear",group=1)) +
  geom_point(aes(x=effect,y=fit,col=var),size=3.5)+
  geom_errorbar(aes(ymin=lw, ymax=up,col=var), width=0.5,
    position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Effect",x="Gear",col='',fill="")+
  scale_color_viridis_d()+
  scale_y_continuous() +
  scale_x_discrete() +
  #scale_colour_manual(values="gray60") +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p44

ggsave(paste(c("Effects_model_Small_Tunas_",sp,".png"),collapse=''), plot = p44, device = "png", units = "cm",
  width = 24, height = 14)


#----------------------------------
#simulating frequency distributions
#----------------------------------
dist_sim<-distsim(
  mean = mean_pred$pred[mean_pred$codsp==sp],
  cv= cv,
  years= unique(mean_pred$yr[mean_pred$codsp==sp]),
  n=100,
  sp=sp,
  distribution = distribution,
  source = "Dist_Sim") 

#assign the distribution simulated for the species
dist_sim_kgm=dist_sim


#-------------
#plot section
#-------------
#Mean predicted (line) + Distribution simulation (boxplots)
p45 <- ggplot() +
  geom_boxplot(data = dplyr::filter(dist_sim, codsp == sp), 
    aes(x = factor(yr), y = fl, fill = codsp, col = codsp),
    alpha = 0.3,
    position = "identity",
    lwd = 0.75,
    outlier.shape = 8,
    outlier.size = 0.8,
    outlier.alpha = 0.5) +
  geom_ribbon(data = dplyr::filter(mean_pred, codsp == sp),
    aes(x = factor(yr), ymin = lw, ymax = up, group = 1, linetype = 'Pred'),
    color = "white", fill = "grey80", alpha = 0.4) +
  geom_line(data = dplyr::filter(mean_pred, codsp == sp),
    aes(x = factor(yr), y = pred, group = 1, linetype = 'Pred'), 
    col = "gray60", size = 2) +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr), y = mean, group = 1), size = 2) +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr[1]), y = mean[1]), pch = 8, size = 2, col = "deepskyblue") +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr[1]), y = mean[1]),  size = 2, col = "deepskyblue") +
  labs(x = "Year", y = "Fork length (cm)", col = "", linetype = '', fill = '') +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  scale_y_continuous() +
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10)))) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p45

ggsave(paste(c("Distribution_Simulation_Boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p45, device = "png", units = "cm",
  width = 24, height = 15)
#-------------------------------------

#plot half violin for the predicted frequencies
p46<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(dist_sim,codsp==sp), 
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
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr),y=mean),size=2)+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr[1]),y=mean[1]),pch=8,size=2)+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr[1]),y=mean[1]),size=2)+
  # #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2025,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p46

ggsave(paste(c("Distribution_Simulation_HalfViolin_Small_Tunas_",sp,".png"),collapse = ''), plot = p46, device = "png", units = "cm",
  width = 35, height = 12)
#--------------------------


#--------------------------------------------------------
#Comparing Distribution simulation with LBSPR simulation
#--------------------------------------------------------
#Life history data is necessary in this section
linf=pars$linf[pars$codsp==sp] #growth
k=pars$k[pars$codsp==sp]   #~k=0.16(Nobrega, 2002) #growth
t0=pars$t0[pars$codsp==sp]  #growth
l50=pars$l50[pars$codsp==sp] #maturity
l95= pars$l50[pars$codsp==sp]*1.1 #maturity
m=pars$m[pars$codsp==sp] # 
mk=m/k #M/K
by=5 #binwidth


#----- Comparing the Simulation Distributions -------#
#   Comparing with the LBSPR Simulation function     #
#----------------------------------------------------#              
#empty frame to store the simulated frequencies
comp_sim<-data.frame(yr=NULL,codsp=NULL,lmids=NULL,sim1=NULL,sim2=NULL,sel_method=NULL,
  fm_method=NULL,sl_fm_methods=NULL,sl50=NULL,f_m=NULL,res=NULL,rmse=NULL,x2=NULL)  

for (i in years) {
  
      # Sim1 is the distribution simulation from the present study
      sim1<-  lf(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], by=by)
      sim1$counts<- norm(sim1$counts)
  
      # Sim2 is the LBSPR simulation population 
      #Selectivity from Catch-curve  and Selectivity from the LBSPR fit 
      #F/M from the Beverton and Holt and F/M from the LBSPR fit 
  
      #Selectivity from Catch-curve and Selectivity from the Beverton and Holt-------------
      sim2_cc_bh<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=c(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by-6,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        sl95= 1.14*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"), method="Beverton and Holt",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_cc_bh@SL50
      f_m= sim2_cc_bh@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_cc_bh@LMids,counts=sim2_cc_bh@pLCatch*100) %>% 
          tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_cc_bh<- data.frame(mids=sim2_cc_bh@LMids,counts=norm(sim2_cc_bh@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_cc_bh <- subset(sim2_cc_bh, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_cc_bh$counts,
        sel_method="CC", fm_method="BH",sl_fm_methods= paste("CC_BH"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_cc_bh$counts,
        rmse=rmse(sim1$counts,sim2_cc_bh$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
            sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
      comp_sim<-rbind(comp_sim,sim_dat)
  
      #Selectivity from Catch-curve and F/M from the LBSPR fit-------------
      sim2_cc_lb<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        sl95= 1.14*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"), method="LBSPR",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_cc_lb@SL50
      f_m= sim2_cc_lb@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_cc_lb@LMids,counts=sim2_cc_lb@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_cc_lb<- data.frame(mids=sim2_cc_lb@LMids,counts=norm(sim2_cc_lb@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_cc_lb <- subset(sim2_cc_lb, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_cc_lb$counts,
        sel_method="CC", fm_method="LBSPR",sl_fm_methods= paste("CC_LB"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_cc_lb$counts,
        rmse=rmse(sim1$counts,sim2_cc_lb$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
            sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
      comp_sim<-rbind(comp_sim,sim_dat)
  
  
      #Selectivity from LBSPR and F/M from Beverton and Holt------------
      sim2_lb_bh<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        sl95= 1.14*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),method="Beverton and Holt",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_lb_bh@SL50
      f_m= sim2_lb_bh@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_lb_bh@LMids,counts=sim2_lb_bh@pLCatch*100) %>% 
          tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_lb_bh<- data.frame(mids=sim2_lb_bh@LMids,counts=norm(sim2_lb_bh@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_lb_bh <- subset(sim2_lb_bh, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_lb_bh$counts,
        sel_method="LB", fm_method="BH",sl_fm_methods= paste("LB_BH"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_lb_bh$counts,
        rmse=rmse(sim1$counts,sim2_lb_bh$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
            sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
      comp_sim<-rbind(comp_sim,sim_dat)
  
  
      #Selectivity from LBSPR and F/M from LBSPR -----------
      sim2_lb_lb<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        sl95= 1.14*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"), method="LBSPR",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_lb_lb@SL50
      f_m= sim2_lb_lb@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_lb_lb@LMids,counts=sim2_lb_lb@pLCatch*100) %>% 
          tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_lb_lb<- data.frame(mids=sim2_lb_lb@LMids,counts=norm(sim2_lb_lb@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_lb_lb <- subset(sim2_lb_lb, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_lb_lb$counts,
        sel_method="LB", fm_method="LB",sl_fm_methods= paste("LB_LB"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_lb_lb$counts,
        rmse=rmse(sim1$counts,sim2_lb_lb$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
            sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
 comp_sim<-rbind(comp_sim,sim_dat)
}
#assing the comparisons for the current species
comp_sim_kgm<-comp_sim

#-------------
#plot section
#-------------
#Length residual comparison 
p47 <- ggplot() +
  geom_boxplot(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(group = interaction(factor(yr), factor(sl_fm_methods)), 
      x = factor(yr), 
      y = res, 
      fill = factor(sl_fm_methods)),  # fill to methods
    alpha = 1,
    position = position_dodge(width = 0.7), 
    width = 0.9,
    lwd = 0.1,  
    outlier.shape = 8,
    outlier.size = 0.7,
    outlier.alpha = 0.5) +
  geom_hline(yintercept =0,linetype = "dashed" )+
  geom_smooth(data = dplyr::filter(comp_sim, codsp == sp),
    method = "loess", span = 0.2, se = FALSE,
    aes(group = 1, x = factor(yr), y = res,color="Loess"), size = 1.5) +
  geom_text(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = "2015", y = 0.9 * max(res, na.rm = TRUE), 
      label = paste("RMSE=", round(mean(rmse, na.rm = TRUE), 2))),
    inherit.aes = FALSE) +
  labs(x = "Year", y = "Length residual", col = "", linetype = '', fill = '') +
  scale_y_continuous() +
  scale_fill_viridis_d(direction = -1)+
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10)))) +
  theme_classic(base_size = 14) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p47

ggsave(paste(c("Length_Residual_Boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p47, device = "png", units = "cm",
  width = 25, height = 17)
#------------------------

#Scaterplot comparing the distributions (Dist sim vs LBSPR sim)
p48 <- ggplot() +
  geom_point(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = sim2, y = sim1, color = sl_fm_methods), 
    size = 2, alpha = 0.7) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", size = 1) +
  geom_smooth(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = sim2, y = sim1, color = sl_fm_methods, fill = sl_fm_methods), 
    method = "loess", se = FALSE, size = 1.2, alpha = 0.3) +  
  scale_color_viridis_d(direction = -1) +  
  labs(x = "LBSPR Distribution Simulation",y= "Distribution Simulation",col = "Method", shape = "Method", fill = "Method") +
  theme_classic(base_size = 14) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p48

ggsave(paste(c("Distributions_Comparison_Small_Tunas_",sp,".png"),collapse = ''), plot = p48, device = "png", units = "cm",
  width = 24, height = 15)
#---------------------------------------------------------------------------------------------------------------------------




#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
# models by species
# LTA model - Length Reconstruction 
#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
sp='LTA'
years=1950:2025

#-------------------------------------
#find the best distribution for specie
#-------------------------------------
dist<-distfit(len = data.frame(
  yr=freq_obs$yr[freq_obs$codsp==sp],
  codsp=freq_obs$codsp[freq_obs$codsp==sp],
  codgr=freq_obs$codgr[freq_obs$codsp==sp],
  fl=freq_obs$fl[freq_obs$codsp==sp]))

dist=dist %>%
  dplyr::group_by(yr,codsp,codgr) %>% 
  slice_min(aic, n = 1)  %>%
  dplyr::arrange(dist) 

#assign dist for the current species
dist_lta=dist

#Table with the best distribution
table(dist$codsp,dist$dist)  
distribution='invgauss'

#coefficient of variation for the corresponding distribution 
cv=mean(dist$cv[dist$codsp==sp & dist$dist==distribution])

#-------------------------------------------------------
#predict 1950 length when using the slope of the lengths
#-------------------------------------------------------

#assuming the 1950 mean length as the prediction of a linear model
lm_data<- data.frame(length=mean_obs$mean[mean_obs$codsp==sp],year=mean_obs$yr[mean_obs$codsp==sp])

#fiting the model
lm_model<- lm(length~year, data = lm_data[-1,])
#Parameters... On average the mean length is being reduced 0.10 cm year by year
a<- coef(lm_model)[1]
b<- coef(lm_model)[2]

#using the fitted parameters to estimate what would be the length in 1950
newdata<- data.frame(year=1950:2025)
lm_predict_10<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.1)))
lm_predict_20<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.2)))
lm_predict_30<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.3)))

#data frame with predictions for different percentile levels
lm_predict<- data.frame(yr=newdata$year,
  fit=lm_predict_10$fit.fit, 
  lw_10=lm_predict_10$fit.lwr,
  up_10=lm_predict_10$fit.upr,
  lw_20=lm_predict_20$fit.lwr,
  up_20=lm_predict_20$fit.upr,
  lw_30=lm_predict_30$fit.lwr,
  up_30=lm_predict_30$fit.upr,
  codsp=sp)

#lm model coefficients and predicitions for the current species
lm_model_lta<- data.frame(formula=paste(deparse(lm_model$call), collapse = ""),
  intercept=a,
  slope=b,
  codsp=sp)
lm_predict_lta<- lm_predict

# controlling variables
lm_predict <- lm_predict %>%
  mutate(confidence_level = case_when(
    yr == 1950 & lw_10 == lw_10 ~ "10%",
    yr == 1950 & lw_20 == lw_20 ~ "20%",
    yr == 1950 & lw_30 == lw_30 ~ "30%"
  ))

# visual representation of the linear model
#-------------
#plot section
#-------------
p49 <- ggplot(lm_predict, aes(x = yr, y = fit)) +
  geom_line(size = 1, color = "grey20") +  # main line prediction
  # Lines for the percentiles
  geom_line(aes(y = lw_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = up_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = lw_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = up_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = lw_30, color = "30%", linetype = "30%"), size = 0.8) +
  geom_line(aes(y = up_30, color = "30%", linetype = "30%"), size = 0.8) +
  # Observed data
  geom_point(data = lm_data, aes(x = year, y = length), size = 2, color = "grey30") +
  # Predicted points for 1950 
  geom_point(aes(x = yr[1], y = lw_10[1], fill = "10%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_10[1], fill = "10%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = lw_20[1], fill = "20%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_20[1], fill = "20%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = lw_30[1], fill = "30%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_30[1], fill = "30%"), size = 3, shape = 21) +
  labs(y = 'Mean Length (cm)', x = "Year", color = "Percentile", fill = "Percentile", linetype = "Percentile") +
  # Type and line colors
  scale_color_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_fill_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_linetype_manual(values = c("10%" = "dashed", "20%" = "dotted", "30%" = "dotdash")) +
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p49
#saving graph...
ggsave(paste0("Linear_model_prediction_Small_Tunas_", sp, ".png"), plot = p42, device = "png", units = "cm", width = 24, height = 16)
#----------------------------------

#selecting variables .. subseting
mean_obs_sub= mean_obs[,c("yr",'region','codsp','codgr','mean')]

#removing extreme values
mean_obs_sub$mean[mean_obs_sub$codsp=='BLF' & mean_obs_sub$mean<50] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='DOL' & mean_obs_sub$mean<75] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='KGM' & mean_obs_sub$mean>100] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='KGM' & mean_obs_sub$mean<60] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean<100] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean>180] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='FRI' & mean_obs_sub$mean>45] =NA

#get only non-NA values and selecting the current species
mean_obs_sub=mean_obs_sub[is.na(mean_obs_sub$mean)==FALSE & mean_obs_sub$codsp==sp,]

#transforming columns in to factors
mean_obs_sub$region= factor(mean_obs_sub$region)
mean_obs_sub$codgr= factor(mean_obs_sub$codgr)
mean_obs_sub$codsp= factor(mean_obs_sub$codsp)


# model GAM selection----------------
#list of values to be tested
tests<- names(lm_predict[, grepl("fit|lw|up", names(lm_predict))])
#filling an empty 1950 with NA 
mean_obs_sub<- rbind(data.frame(yr=1950, region="NE", codsp=sp,codgr="UN",mean=NA),
  mean_obs_sub)

# Definindo os vetores de variáveis de resposta e explicativas
resp_vars <- c("mean")
exp_vars_1 <- c("s(yr)")
exp_vars_2 <- c("region")
exp_vars_3 <- c("codgr")

#criando a formula geral
formulas <-list(as.formula(paste(resp_vars, "~", paste(c(exp_vars_1), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_2), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_3), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_2,exp_vars_3), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_3,exp_vars_2), collapse = "+"))))

formulas<- as.character(formulas)

# Data frame com as distribuicoes de probabilidade e funcoes de ligacao a serem testadas
family=distribution #testing only for the best distribution
link = ("log") #log fixed

#empty frame to store results
gam=data.frame(
  test=NULL,
  type=NULL,
  codsp=NULL,
  family=NULL,
  link=NULL,
  formula=NULL,
  newdataexpand=NULL,
  pred=NULL,
  aic=NULL,
  r.sq=NULL,
  dev.ex=NULL)


#-------------------- Fitting models  -------------------#
for (i in tests) {
  
  #1- first sensibility for different 1950 starting values--
  test_1950<- lm_predict[1, grepl(i, names(lm_predict))]
  
  #replacing with the current 1950 mean length test
  mean_obs_sub$mean[mean_obs_sub$yr==1950]<- test_1950
  
  #2- Term selection within the sensibility----
  # Loop through the formulas 
  for (j in formulas) {
    
      family_name <- as.character(family) #these are fixed
      link_name <- as.character(link)
    
      #Creating switch object
      family_obj <- switch(family_name,
        norm = gaussian(link = link_name),  
        invgauss = inverse.gaussian(link = link_name),
        gamma= Gamma(link=link_name))
    
      #formula j to be tested
      formula <- as.formula(j)
    
      #fit the model
      gam_model <- mgcv::gam(as.formula(formula),
        data= mean_obs_sub,
        method = "REML" ,
        family =family_obj)
    
      #taking the outputs
      gam=rbind(gam, data.frame(
        test=i,
        type="GAM",
        codsp=sp,
        family=gam_model$family$family,
        link=gam_model$family$link,
        formula=paste0(c(paste(gam_model$formula)[2],
          paste(gam_model$formula)[1],
          paste(gam_model$formula)[3]),collapse = ' '),
        aic=gam_model$aic,
        r.sq=summary(gam_model)$r.sq,
      dev.ex=summary(gam_model)$dev.expl))
  }
}
#-------------------- end of loop -----------------------#

#assign the tests to the species
gam_lta<-gam

#replacing the 1950 mean length with the best tested value
best_model<- gam[which.min(gam$aic),]; dplyr::tibble(best_model)
mean_obs_sub$mean[mean_obs_sub$yr==1950]<-lm_predict[1, grepl(best_model$test, names(lm_predict))]

#Final model (lowest AIC) 
gam_model<- mgcv::gam(as.formula(best_model$formula), 
  #sp=0.001,
  data= mean_obs_sub,
  method = "REML" ,
  family =family_obj)

#assign GAM model to the specie
gam_model_lta=gam_model

#diagnostic
summary(gam_model)

par(mfrow=c(2,2))
mgcv::gam.check(gam_model)# model check

#residual 
#homocedasticity
lmtest::bptest(gam_model)

#normality
shapiro.test(resid(gam_model))

# Check overall concurvity
round(mgcv::concurvity(gam_model, full = TRUE),2)
#----------------------
#predictions "year"
newdataexpand= expand.grid(yr=seq(years[1],years[length(years)],1),
  region=mean_obs_sub$region,
  codgr=mean_obs_sub$codgr)

pred<- data.frame(predict(gam_model, newdata = newdataexpand,
  type = "response", se.fit=TRUE))

pred$up<-pred$fit+qnorm(0.975)*pred$se.fit #confidence intervals
pred$lw<-pred$fit-qnorm(0.975)*pred$se.fit

#take the predictions
mean_pred=cbind(type='GAM',
  codsp=sp,
  yr=newdataexpand$yr,
  mean=NA,
  newdata=newdataexpand,
  pred = round(pred$fit, 2),
  lw = round(pred$lw, 2),
  up = round(pred$up, 2))

#extract mean of all levels
mean_pred = mean_pred %>%
  group_by(type, yr, codsp) %>%
  summarize(across(c('pred', 'lw', 'up'), ~ round(mean(.), 2)), .groups = 'drop') %>%
  #binding the observed mean lenth data and the other factors
  dplyr::left_join(mean_obs_sub[, c("yr", "mean","region","codgr")], by = "yr")

#assign the mean predictions for the species
mean_pred_lta=mean_pred #size mean predicted for all factors combined

#plotting the GAM's prediction and the observed data
#-------------
#plot section
#-------------
#plot GAM mean size prediction
p50<-ggplot(data=mean_pred,aes(x=yr,y=pred))+
  geom_ribbon(data=dplyr::filter(mean_pred,codsp==sp),
    aes(ymin = lw, ymax = up,linetype='Pred'), 
    col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(mean_pred,codsp==sp),
    aes(x=yr,y=pred,linetype='Pred'),col="gray60",size=2)+
  #geom_line(data=dplyr::filter(gam,codsp==sp),aes(x=yr,y=fit,col=codgr,linetype=region),size=1.2)+
  geom_point(data=dplyr::filter(mean_pred,codsp==sp),
    aes(x=yr,y=mean),size=2,col="gray40")+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=yr[1],y=mean[1]),pch=8,size=2,col="gray40")+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='')+
  scale_y_continuous() +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p50

ggsave(paste(c("Length_model_prediction_Small_Tunas_",sp,".png"),collapse=''), plot = p50, device = "png", units = "cm",
  width = 24, height = 14)
#----------------------------------------------

#partial effects
newdata= mean_obs_sub   #desired partial effect

pareffect<- data.frame(predict(gam_model, newdata = newdata,
  type = "terms", se.fit=TRUE)) #terms gets the effects

#intervals for the partial effect
pareffect$up.region<-pareffect$fit.region+qnorm(0.975)*pareffect$se.fit.region #confidence intervals
pareffect$lw.region<-pareffect$fit.region-qnorm(0.975)*pareffect$se.fit.region
pareffect$region= mean_obs_sub$region

#extract mean of all levels (Region)
pareffect= pareffect[grepl("region", names(pareffect))]
pareffect <- pareffect %>%
  dplyr::group_by(region) %>%
  dplyr::summarize(across(c(fit.region, se.fit.region, up.region, lw.region), ~ round(mean(.), 2))) %>%
  dplyr::rename(effect = region, 
    fit = fit.region, 
    se.fit = se.fit.region, 
    up = up.region, 
    lw = lw.region) %>%
  dplyr::mutate(codsp = sp, var = "Region")

pareffect_lta=pareffect


# plot GAM partial effects
p51<-ggplot(dplyr::filter(pareffect,codsp==sp),
  aes(x=effect,y=fit,fill=var,col=var),size=2) +
  #geom_line(aes(x=effect,y=up, col="Gear",group=1)) +
  geom_point(aes(x=effect,y=fit,col=var),size=3.5)+
  geom_errorbar(aes(ymin=lw, ymax=up,col=var), width=0.5,
    position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Effect",x="Gear",col='',fill="")+
  scale_color_viridis_d()+
  scale_y_continuous() +
  scale_x_discrete() +
  #scale_colour_manual(values="gray60") +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p51

ggsave(paste(c("Effects_model_Small_Tunas_",sp,".png"),collapse=''), plot = p51, device = "png", units = "cm",
  width = 24, height = 14)


#----------------------------------
#simulating frequency distributions
#----------------------------------
dist_sim<-distsim(
  mean = mean_pred$pred[mean_pred$codsp==sp],
  cv= cv,
  years= unique(mean_pred$yr[mean_pred$codsp==sp]),
  n=100,
  sp=sp,
  distribution = distribution,
  source = "Dist_Sim") 

#assign the distribution simulated for the species
dist_sim_lta=dist_sim


#-------------
#plot section
#-------------
#Mean predicted (line) + Distribution simulation (boxplots)
p52 <- ggplot() +
  geom_boxplot(data = dplyr::filter(dist_sim, codsp == sp), 
    aes(x = factor(yr), y = fl, fill = codsp, col = codsp),
    alpha = 0.3,
    position = "identity",
    lwd = 0.75,
    outlier.shape = 8,
    outlier.size = 0.8,
    outlier.alpha = 0.5) +
  geom_ribbon(data = dplyr::filter(mean_pred, codsp == sp),
    aes(x = factor(yr), ymin = lw, ymax = up, group = 1, linetype = 'Pred'),
    color = "white", fill = "grey80", alpha = 0.4) +
  geom_line(data = dplyr::filter(mean_pred, codsp == sp),
    aes(x = factor(yr), y = pred, group = 1, linetype = 'Pred'), 
    col = "gray60", size = 2) +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr), y = mean, group = 1), size = 2) +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr[1]), y = mean[1]), pch = 8, size = 2, col = "deepskyblue") +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr[1]), y = mean[1]),  size = 2, col = "deepskyblue") +
  labs(x = "Year", y = "Fork length (cm)", col = "", linetype = '', fill = '') +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  scale_y_continuous() +
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10)))) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p52

ggsave(paste(c("Distribution_Simulation_Boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p52, device = "png", units = "cm",
  width = 24, height = 15)
#-------------------------------------

#plot half violin for the predicted frequencies
p53<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(dist_sim,codsp==sp), 
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
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr),y=mean),size=2)+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr[1]),y=mean[1]),pch=8,size=2)+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr[1]),y=mean[1]),size=2)+
  # #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2025,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p53

ggsave(paste(c("Distribution_Simulation_HalfViolin_Small_Tunas_",sp,".png"),collapse = ''), plot = p53, device = "png", units = "cm",
  width = 35, height = 12)
#--------------------------


#--------------------------------------------------------
#Comparing Distribution simulation with LBSPR simulation
#--------------------------------------------------------
#Life history data is necessary in this section
linf=pars$linf[pars$codsp==sp] #growth
k=pars$k[pars$codsp==sp]   #~k=0.16(Nobrega, 2002) #growth
t0=pars$t0[pars$codsp==sp]  #growth
l50=pars$l50[pars$codsp==sp] #maturity
l95= pars$l50[pars$codsp==sp]*1.1 #maturity
m=pars$m[pars$codsp==sp] # 
mk=m/k #M/K
by=5 #binwidth


                        #----- Comparing the Simulation Distributions -------#
                        #   Comparing with the LBSPR Simulation function     #
                        #----------------------------------------------------#              
#empty frame to store the simulated frequencies
comp_sim<-data.frame(yr=NULL,codsp=NULL,lmids=NULL,sim1=NULL,sim2=NULL,sel_method=NULL,
  fm_method=NULL,sl_fm_methods=NULL,sl50=NULL,f_m=NULL,res=NULL,rmse=NULL,x2=NULL)  

for (i in years) {
  
      # Sim1 is the distribution simulation from the present study
      sim1<-  lf(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], by=by)
      sim1$counts<- norm(sim1$counts)
  
      # Sim2 is the LBSPR simulation population 
      #Selectivity from Catch-curve  and Selectivity from the LBSPR fit 
      #F/M from the Beverton and Holt and F/M from the LBSPR fit 
  
      #Selectivity from Catch-curve and Selectivity from the Beverton and Holt-------------
      sim2_cc_bh<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=c(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        sl95= 1.15*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"), method="Beverton and Holt",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_cc_bh@SL50
      f_m= sim2_cc_bh@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_cc_bh@LMids,counts=sim2_cc_bh@pLCatch*100) %>% 
          tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_cc_bh<- data.frame(mids=sim2_cc_bh@LMids,counts=norm(sim2_cc_bh@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_cc_bh <- subset(sim2_cc_bh, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_cc_bh$counts,
        sel_method="CC", fm_method="BH",sl_fm_methods= paste("CC_BH"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_cc_bh$counts,
        rmse=rmse(sim1$counts,sim2_cc_bh$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
            sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
      comp_sim<-rbind(comp_sim,sim_dat)
  
      #Selectivity from Catch-curve and F/M from the LBSPR fit-------------
      sim2_cc_lb<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        sl95= 1.15*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"), method="LBSPR",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_cc_lb@SL50
      f_m= sim2_cc_lb@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_cc_lb@LMids,counts=sim2_cc_lb@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_cc_lb<- data.frame(mids=sim2_cc_lb@LMids,counts=norm(sim2_cc_lb@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_cc_lb <- subset(sim2_cc_lb, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_cc_lb$counts,
        sel_method="CC", fm_method="LBSPR",sl_fm_methods= paste("CC_LB"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_cc_lb$counts,
        rmse=rmse(sim1$counts,sim2_cc_lb$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
        sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
      comp_sim<-rbind(comp_sim,sim_dat)
  
  
      #Selectivity from LBSPR and F/M from Beverton and Holt------------
      sim2_lb_bh<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        sl95= 1.15*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),method="Beverton and Holt",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_lb_bh@SL50
      f_m= sim2_lb_bh@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_lb_bh@LMids,counts=sim2_lb_bh@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_lb_bh<- data.frame(mids=sim2_lb_bh@LMids,counts=norm(sim2_lb_bh@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_lb_bh <- subset(sim2_lb_bh, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
     #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_lb_bh$counts,
        sel_method="LB", fm_method="BH",sl_fm_methods= paste("LB_BH"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_lb_bh$counts,
        rmse=rmse(sim1$counts,sim2_lb_bh$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
            sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
      comp_sim<-rbind(comp_sim,sim_dat)
  
  
      #Selectivity from LBSPR and F/M from LBSPR -----------
      sim2_lb_lb<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        sl95= 1.15*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"), method="LBSPR",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
  
      #store Selectivity and F/M estimated
      sl50=sim2_lb_lb@SL50
      f_m= sim2_lb_lb@FM
  
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_lb_lb@LMids,counts=sim2_lb_lb@pLCatch*100) %>% 
          tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
  
      #storing LBSPR distribution simulation
      sim2_lb_lb<- data.frame(mids=sim2_lb_lb@LMids,counts=norm(sim2_lb_lb@pLCatch) )
  
      # Adjusting to keep the same mids interval
      sim2_lb_lb <- subset(sim2_lb_lb, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
  
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_lb_lb$counts,
        sel_method="LB", fm_method="LB",sl_fm_methods= paste("LB_LB"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_lb_lb$counts,
        rmse=rmse(sim1$counts,sim2_lb_lb$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
            sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
  comp_sim<-rbind(comp_sim,sim_dat)
}
#assing the comparisons for the current species
comp_sim_lta<-comp_sim

#-------------
#plot section
#-------------
#Length residual comparison 
p54 <- ggplot() +
  geom_boxplot(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(group = interaction(factor(yr), factor(sl_fm_methods)), 
      x = factor(yr), 
      y = res, 
      fill = factor(sl_fm_methods)),  # fill to methods
    alpha = 1,
    position = position_dodge(width = 0.7), 
    width = 0.9,
    lwd = 0.1,  
    outlier.shape = 8,
    outlier.size = 0.7,
    outlier.alpha = 0.5) +
  geom_hline(yintercept =0,linetype = "dashed" )+
  geom_smooth(data = dplyr::filter(comp_sim, codsp == sp),
    method = "loess", span = 0.2, se = FALSE,
    aes(group = 1, x = factor(yr), y = res,color="Loess"), size = 1.5) +
  geom_text(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = "2015", y = 0.9 * max(res, na.rm = TRUE), 
      label = paste("RMSE=", round(mean(rmse, na.rm = TRUE), 2))),
    inherit.aes = FALSE) +
  labs(x = "Year", y = "Length residual", col = "", linetype = '', fill = '') +
  scale_y_continuous() +
  scale_fill_viridis_d(direction = -1)+
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10)))) +
  theme_classic(base_size = 14) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p54

ggsave(paste(c("Length_Residual_Boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p54, device = "png", units = "cm",
  width = 25, height = 17)
#------------------------

#Scaterplot comparing the distributions (Dist sim vs LBSPR sim)
p55 <- ggplot() +
  geom_point(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = sim2, y = sim1, color = sl_fm_methods), 
    size = 2, alpha = 0.7) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", size = 1) +
  geom_smooth(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = sim2, y = sim1, color = sl_fm_methods, fill = sl_fm_methods), 
    method = "loess", se = FALSE, size = 1.2, alpha = 0.3) +  
  scale_color_viridis_d(direction = -1) +  
  labs(x = "LBSPR Distribution Simulation",y= "Distribution Simulation",col = "Method", shape = "Method", fill = "Method") +
  theme_classic(base_size = 14) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p55

ggsave(paste(c("Distributions_Comparison_Small_Tunas_",sp,".png"),collapse = ''), plot = p55, device = "png", units = "cm",
  width = 24, height = 15)
#---------------------------------------------------------------------------------------------------------------------------






#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
# models by species
# WAH model - Length Reconstruction 
#X##X##X##X##X##X##X##X##X##X#X##X##X##X##X#
sp='WAH'
years=1950:2025

#-------------------------------------
#find the best distribution for specie
#-------------------------------------
dist<-distfit(len = data.frame(
  yr=freq_obs$yr[freq_obs$codsp==sp],
  codsp=freq_obs$codsp[freq_obs$codsp==sp],
  codgr=freq_obs$codgr[freq_obs$codsp==sp],
  fl=freq_obs$fl[freq_obs$codsp==sp]))

dist=dist %>%
  dplyr::group_by(yr,codsp,codgr) %>% 
  slice_min(aic, n = 1)  %>%
  dplyr::arrange(dist) 

#assign dist for the current species
dist_wah=dist

#Table with the best distribution
table(dist$codsp,dist$dist)  
distribution='norm'

#coefficient of variation for the corresponding distribution 
cv=mean(dist$cv[dist$codsp==sp & dist$dist==distribution])

#-------------------------------------------------------
#predict 1950 length when using the slope of the lengths
#-------------------------------------------------------

#assuming the 1950 mean length as the prediction of a linear model
lm_data<- data.frame(length=mean_obs$mean[mean_obs$codsp==sp],year=mean_obs$yr[mean_obs$codsp==sp])
lm_data<- lm_data[lm_data$length>100 & lm_data$length<180,]
#mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean<100] =NA
#mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean>180] =NA


#fiting the model
lm_model<- lm(length~year, data = lm_data[-1,])
#Parameters... On average the mean length is being reduced 0.10 cm year by year
a<- coef(lm_model)[1]
b<- coef(lm_model)[2]

#using the fitted parameters to estimate what would be the length in 1950
newdata<- data.frame(year=1950:2025)
lm_predict_10<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.1)))
lm_predict_20<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.2)))
lm_predict_30<- as.data.frame(predict.lm(lm_model, newdata = newdata,
  se.fit = TRUE,interval="prediction",level = c(0.3)))

#data frame with predictions for different percentile levels
lm_predict<- data.frame(yr=newdata$year,
  fit=lm_predict_10$fit.fit, 
  lw_10=lm_predict_10$fit.lwr,
  up_10=lm_predict_10$fit.upr,
  lw_20=lm_predict_20$fit.lwr,
  up_20=lm_predict_20$fit.upr,
  lw_30=lm_predict_30$fit.lwr,
  up_30=lm_predict_30$fit.upr,
  codsp=sp)

#lm model coefficients and predicitions for the current species
lm_model_wah<- data.frame(formula=paste(deparse(lm_model$call), collapse = ""),
  intercept=a,
  slope=b,
  codsp=sp)
lm_predict_wah<- lm_predict

# controlling variables
lm_predict <- lm_predict %>%
  mutate(confidence_level = case_when(
    yr == 1950 & lw_10 == lw_10 ~ "10%",
    yr == 1950 & lw_20 == lw_20 ~ "20%",
    yr == 1950 & lw_30 == lw_30 ~ "30%"
  ))

# visual representation of the linear model
#-------------
#plot section
#-------------
p56 <- ggplot(lm_predict, aes(x = yr, y = fit)) +
  geom_line(size = 1, color = "grey20") +  # main line prediction
  # Lines for the percentiles
  geom_line(aes(y = lw_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = up_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = lw_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = up_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = lw_30, color = "30%", linetype = "30%"), size = 0.8) +
  geom_line(aes(y = up_30, color = "30%", linetype = "30%"), size = 0.8) +
  # Observed data
  geom_point(data = lm_data, aes(x = year, y = length), size = 2, color = "grey30") +
  # Predicted points for 1950 
  geom_point(aes(x = yr[1], y = lw_10[1], fill = "10%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_10[1], fill = "10%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = lw_20[1], fill = "20%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_20[1], fill = "20%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = lw_30[1], fill = "30%"), size = 3, shape = 21) +
  geom_point(aes(x = yr[1], y = up_30[1], fill = "30%"), size = 3, shape = 21) +
  labs(y = 'Mean Length (cm)', x = "Year", color = "Percentile", fill = "Percentile", linetype = "Percentile") +
  # Type and line colors
  scale_color_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_fill_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_linetype_manual(values = c("10%" = "dashed", "20%" = "dotted", "30%" = "dotdash")) +
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p56
#saving graph...
ggsave(paste0("Linear_model_prediction_Small_Tunas_", sp, ".png"), plot = p56, device = "png", units = "cm", width = 24, height = 16)
#----------------------------------

#selecting variables .. subseting
mean_obs_sub= mean_obs[,c("yr",'region','codsp','codgr','mean')]

#removing extreme values
mean_obs_sub$mean[mean_obs_sub$codsp=='BLF' & mean_obs_sub$mean<50] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='DOL' & mean_obs_sub$mean<75] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='KGM' & mean_obs_sub$mean>100] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='KGM' & mean_obs_sub$mean<60] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean<100] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='WAH' & mean_obs_sub$mean>180] =NA
mean_obs_sub$mean[mean_obs_sub$codsp=='FRI' & mean_obs_sub$mean>45] =NA

#get only non-NA values and selecting the current species
mean_obs_sub=mean_obs_sub[is.na(mean_obs_sub$mean)==FALSE & mean_obs_sub$codsp==sp,]

#transforming columns in to factors
mean_obs_sub$region= factor(mean_obs_sub$region)
mean_obs_sub$codgr= factor(mean_obs_sub$codgr)
mean_obs_sub$codsp= factor(mean_obs_sub$codsp)


# model GAM selection----------------
#list of values to be tested
tests<- names(lm_predict[, grepl("fit|lw|up", names(lm_predict))])
#filling an empty 1950 with NA 
mean_obs_sub<- rbind(data.frame(yr=1950, region="NE", codsp=sp,codgr="UN",mean=NA),
  mean_obs_sub)

# Definindo os vetores de variáveis de resposta e explicativas
resp_vars <- c("mean")
exp_vars_1 <- c("s(yr)")
exp_vars_2 <- c("region")
exp_vars_3 <- c("codgr")

#criando a formula geral
formulas <-list(as.formula(paste(resp_vars, "~", paste(c(exp_vars_1), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_2), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_3), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_2,exp_vars_3), collapse = "+"))),
  as.formula(paste(resp_vars, "~", paste(c(exp_vars_1,exp_vars_3,exp_vars_2), collapse = "+"))))

formulas<- as.character(formulas)

# Data frame com as distribuicoes de probabilidade e funcoes de ligacao a serem testadas
family=distribution #testing only for the best distribution
link = ("log") #log fixed

#empty frame to store results
gam=data.frame(
  test=NULL,
  type=NULL,
  codsp=NULL,
  family=NULL,
  link=NULL,
  formula=NULL,
  newdataexpand=NULL,
  pred=NULL,
  aic=NULL,
  r.sq=NULL,
  dev.ex=NULL)


#-------------------- Fitting models  -------------------#
for (i in tests) {
  
      #1- first sensibility for different 1950 starting values--
      test_1950<- lm_predict[1, grepl(i, names(lm_predict))]
  
      #replacing with the current 1950 mean length test
      mean_obs_sub$mean[mean_obs_sub$yr==1950]<- test_1950
  
      #2- Term selection within the sensibility----
      # Loop through the formulas 
      for (j in formulas) {
    
      family_name <- as.character(family) #these are fixed
      link_name <- as.character(link)
    
      #Creating switch object
      family_obj <- switch(family_name,
        norm = gaussian(link = link_name),  
        invgauss = inverse.gaussian(link = link_name),
        gamma= Gamma(link=link_name))
    
      #formula j to be tested
      formula <- as.formula(j)
    
      #fit the model
      gam_model <- mgcv::gam(as.formula(formula),
        data= mean_obs_sub,
        method = "REML" ,
        family =family_obj)
    
      #taking the outputs
      gam=rbind(gam, data.frame(
        test=i,
        type="GAM",
      codsp=sp,
        family=gam_model$family$family,
        link=gam_model$family$link,
        formula=paste0(c(paste(gam_model$formula)[2],
        paste(gam_model$formula)[1],
        paste(gam_model$formula)[3]),collapse = ' '),
        aic=gam_model$aic,
        r.sq=summary(gam_model)$r.sq,
    dev.ex=summary(gam_model)$dev.expl))
  }
}
#-------------------- end of loop -----------------------#

#assign the tests to the species
gam_wah<-gam

#replacing the 1950 mean length with the best tested value
best_model<- gam[which.min(gam$aic),]; dplyr::tibble(best_model)
mean_obs_sub$mean[mean_obs_sub$yr==1950]<-lm_predict[1, grepl(best_model$test, names(lm_predict))]

#Final model (lowest AIC) 
gam_model<- mgcv::gam(as.formula(best_model$formula), 
  sp=10,
  data= mean_obs_sub,
  method = "REML" ,
  family =family_obj)

#assign GAM model to the specie
gam_model_wah=gam_model

#diagnostic
summary(gam_model)

par(mfrow=c(2,2))
mgcv::gam.check(gam_model)# model check

#residual 
#homocedasticity
lmtest::bptest(gam_model)

#normality
shapiro.test(resid(gam_model))

# Check overall concurvity
round(mgcv::concurvity(gam_model, full = TRUE),2)
#----------------------
#predictions "year"
newdataexpand= expand.grid(yr=seq(years[1],years[length(years)],1),
  region=mean_obs_sub$region,
  codgr=mean_obs_sub$codgr)

pred<- data.frame(predict(gam_model, newdata = newdataexpand,
  type = "response", se.fit=TRUE))

pred$up<-pred$fit+qnorm(0.975)*pred$se.fit #confidence intervals
pred$lw<-pred$fit-qnorm(0.975)*pred$se.fit

#take the predictions
mean_pred=cbind(type='GAM',
  codsp=sp,
  yr=newdataexpand$yr,
  mean=NA,
  newdata=newdataexpand,
  pred = round(pred$fit, 2),
  lw = round(pred$lw, 2),
  up = round(pred$up, 2))

#extract mean of all levels
mean_pred = mean_pred %>%
  group_by(type, yr, codsp) %>%
  summarize(across(c('pred', 'lw', 'up'), ~ round(mean(.), 2)), .groups = 'drop') %>%
  #binding the observed mean lenth data and the other factors
  dplyr::left_join(mean_obs_sub[, c("yr", "mean","region","codgr")], by = "yr")

#assign the mean predictions for the species
mean_pred_wah=mean_pred #size mean predicted for all factors combined

#plotting the GAM's prediction and the observed data
#-------------
#plot section
#-------------
#plot GAM mean size prediction
p57<-ggplot(data=mean_pred,aes(x=yr,y=pred))+
  geom_ribbon(data=dplyr::filter(mean_pred,codsp==sp),
    aes(ymin = lw, ymax = up,linetype='Pred'), 
    col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(mean_pred,codsp==sp),
    aes(x=yr,y=pred,linetype='Pred'),col="gray60",size=2)+
  #geom_line(data=dplyr::filter(gam,codsp==sp),aes(x=yr,y=fit,col=codgr,linetype=region),size=1.2)+
  geom_point(data=dplyr::filter(mean_pred,codsp==sp),
    aes(x=yr,y=mean),size=2,col="gray40")+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=yr[1],y=mean[1]),pch=8,size=2,col="gray40")+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='')+
  scale_y_continuous() +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p57

ggsave(paste(c("Length_model_prediction_Small_Tunas_",sp,".png"),collapse=''), plot = p57, device = "png", units = "cm",
  width = 24, height = 14)
#----------------------------------------------

#partial effects
newdata= mean_obs_sub   #desired partial effect

pareffect<- data.frame(predict(gam_model, newdata = newdata,
  type = "terms", se.fit=TRUE)) #terms gets the effects

#intervals for the partial effect
pareffect$up.codgr<-pareffect$fit.codgr+qnorm(0.975)*pareffect$se.fit.codgr #confidence intervals
pareffect$lw.codgr<-pareffect$fit.codgr-qnorm(0.975)*pareffect$se.fit.codgr
pareffect$codgr= mean_obs_sub$codgr

pareffect$up.region<-pareffect$fit.region+qnorm(0.975)*pareffect$se.fit.region #confidence intervals
pareffect$lw.region<-pareffect$fit.region-qnorm(0.975)*pareffect$se.fit.region
pareffect$region= mean_obs_sub$region


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


# plot GAM partial effects
p58<-ggplot(dplyr::filter(pareffect,codsp==sp),
  aes(x=effect,y=fit,fill=var,col=var),size=2) +
  #geom_line(aes(x=effect,y=up, col="Gear",group=1)) +
  geom_point(aes(x=effect,y=fit,col=var),size=3.5)+
  geom_errorbar(aes(ymin=lw, ymax=up,col=var), width=0.5,
    position="identity",size=1.2)+
  #facet_wrap(~codsp,scales = "free_y")+
  labs(y="Effect",x="Gear",col='',fill="")+
  scale_color_viridis_d()+
  scale_y_continuous() +
  scale_x_discrete() +
  #scale_colour_manual(values="gray60") +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p58

ggsave(paste(c("Effects_model_Small_Tunas_",sp,".png"),collapse=''), plot = p58, device = "png", units = "cm",
  width = 24, height = 14)


#----------------------------------
#simulating frequency distributions
#----------------------------------
dist_sim<-distsim(
  mean = mean_pred$pred[mean_pred$codsp==sp],
  cv= cv,
  years= unique(mean_pred$yr[mean_pred$codsp==sp]),
  n=100,
  sp=sp,
  distribution = distribution,
  source = "Dist_Sim") 

#assign the distribution simulated for the species
dist_sim_wah=dist_sim


#-------------
#plot section
#-------------
#Mean predicted (line) + Distribution simulation (boxplots)
p59 <- ggplot() +
  geom_boxplot(data = dplyr::filter(dist_sim, codsp == sp), 
    aes(x = factor(yr), y = fl, fill = codsp, col = codsp),
    alpha = 0.3,
    position = "identity",
    lwd = 0.75,
    outlier.shape = 8,
    outlier.size = 0.8,
    outlier.alpha = 0.5) +
  geom_ribbon(data = dplyr::filter(mean_pred, codsp == sp),
    aes(x = factor(yr), ymin = lw, ymax = up, group = 1, linetype = 'Pred'),
    color = "white", fill = "grey80", alpha = 0.4) +
  geom_line(data = dplyr::filter(mean_pred, codsp == sp),
    aes(x = factor(yr), y = pred, group = 1, linetype = 'Pred'), 
    col = "gray60", size = 2) +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr), y = mean, group = 1), size = 2) +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr[1]), y = mean[1]), pch = 8, size = 2, col = "deepskyblue") +
  geom_point(data = dplyr::filter(mean_obs_sub, codsp == sp),
    aes(x = factor(yr[1]), y = mean[1]),  size = 2, col = "deepskyblue") +
  labs(x = "Year", y = "Fork length (cm)", col = "", linetype = '', fill = '') +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  scale_y_continuous() +
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10)))) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p59

ggsave(paste(c("Distribution_Simulation_Boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p59, device = "png", units = "cm",
  width = 24, height = 15)
#-------------------------------------

#plot half violin for the predicted frequencies
p60<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(dist_sim,codsp==sp), 
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
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr),y=mean),size=2)+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr[1]),y=mean[1]),pch=8,size=2)+
  geom_point(data=dplyr::filter(mean_obs_sub,codsp==sp),
    aes(x=factor(yr[1]),y=mean[1]),size=2)+
  # #geom_half_boxplot()+
  #geom_half_point(side = "l", 
  #                 show.legend = F)+
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop=FALSE,breaks=c(as.character(seq(1950,2025,20))))+
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p60

ggsave(paste(c("Distribution_Simulation_HalfViolin_Small_Tunas_",sp,".png"),collapse = ''), plot = p60, device = "png", units = "cm",
  width = 35, height = 12)
#--------------------------


#--------------------------------------------------------
#Comparing Distribution simulation with LBSPR simulation
#--------------------------------------------------------
#Life history data is necessary in this section
linf=pars$linf[pars$codsp==sp] #growth
k=pars$k[pars$codsp==sp]   #~k=0.16(Nobrega, 2002) #growth
t0=pars$t0[pars$codsp==sp]  #growth
l50=pars$l50[pars$codsp==sp] #maturity
l95= pars$l50[pars$codsp==sp]*1.1 #maturity
m=pars$m[pars$codsp==sp] # 
mk=m/k #M/K
by=20 #binwidth


                          #----- Comparing the Simulation Distributions -------#
                          #   Comparing with the LBSPR Simulation function     #
                          #----------------------------------------------------#              
#empty frame to store the simulated frequencies
comp_sim<-data.frame(yr=NULL,codsp=NULL,lmids=NULL,sim1=NULL,sim2=NULL,sel_method=NULL,
  fm_method=NULL,sl_fm_methods=NULL,sl50=NULL,f_m=NULL,res=NULL,rmse=NULL,x2=NULL)  

for (i in years) {
  
      # Sim1 is the distribution simulation from the present study
      sim1<-  lf(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], by=by)
      sim1$counts<- norm(sim1$counts)
      
      # Sim2 is the LBSPR simulation population 
      #Selectivity from Catch-curve  and Selectivity from the LBSPR fit 
      #F/M from the Beverton and Holt and F/M from the LBSPR fit 
      
      #Selectivity from Catch-curve and Selectivity from the Beverton and Holt-------------
      sim2_cc_bh<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=c(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by-6,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        sl95= 1.3*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"), method="Beverton and Holt",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
      
      #store Selectivity and F/M estimated
      sl50=sim2_cc_bh@SL50
      f_m= sim2_cc_bh@FM
      
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_cc_bh@LMids,counts=sim2_cc_bh@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
      
      #storing LBSPR distribution simulation
      sim2_cc_bh<- data.frame(mids=sim2_cc_bh@LMids,counts=norm(sim2_cc_bh@pLCatch) )
      
      # Adjusting to keep the same mids interval
      sim2_cc_bh <- subset(sim2_cc_bh, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
      
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_cc_bh$counts,
        sel_method="CC", fm_method="BH",sl_fm_methods= paste("CC_BH"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_cc_bh$counts,
        rmse=rmse(sim1$counts,sim2_cc_bh$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
          sample(freq_lbspr, 95,replace = FALSE))$p.value))
      
      comp_sim<-rbind(comp_sim,sim_dat)
      
      #Selectivity from Catch-curve and F/M from the LBSPR fit-------------
      sim2_cc_lb<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        sl95= 1.3*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by-4,by-5,by+1,by+2,by+3,by+4,by+5),year=i,sp=sp,units="cm",method = "Catch curve"), method="LBSPR",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
      
      #store Selectivity and F/M estimated
      sl50=sim2_cc_lb@SL50
      f_m= sim2_cc_lb@FM
      
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_cc_lb@LMids,counts=sim2_cc_lb@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
      
      #storing LBSPR distribution simulation
      sim2_cc_lb<- data.frame(mids=sim2_cc_lb@LMids,counts=norm(sim2_cc_lb@pLCatch) )
      
      # Adjusting to keep the same mids interval
      sim2_cc_lb <- subset(sim2_cc_lb, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
      
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_cc_lb$counts,
        sel_method="CC", fm_method="LBSPR",sl_fm_methods= paste("CC_LB"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_cc_lb$counts,
        rmse=rmse(sim1$counts,sim2_cc_lb$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
          sample(freq_lbspr, 95,replace = FALSE))$p.value))
      
      comp_sim<-rbind(comp_sim,sim_dat)
      
      
      #Selectivity from LBSPR and F/M from Beverton and Holt------------
      sim2_lb_bh<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        sl95= 1.3*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),method="Beverton and Holt",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
      
      #store Selectivity and F/M estimated
      sl50=sim2_lb_bh@SL50
      f_m= sim2_lb_bh@FM
      
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_lb_bh@LMids,counts=sim2_lb_bh@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
      
      #storing LBSPR distribution simulation
      sim2_lb_bh<- data.frame(mids=sim2_lb_bh@LMids,counts=norm(sim2_lb_bh@pLCatch) )
      
      # Adjusting to keep the same mids interval
      sim2_lb_bh <- subset(sim2_lb_bh, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
      
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_lb_bh$counts,
        sel_method="LB", fm_method="BH",sl_fm_methods= paste("LB_BH"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_lb_bh$counts,
        rmse=rmse(sim1$counts,sim2_lb_bh$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
          sample(freq_lbspr, 95,replace = FALSE))$p.value))
      
      comp_sim<-rbind(comp_sim,sim_dat)
      
      
      #Selectivity from LBSPR and F/M from LBSPR -----------
      sim2_lb_lb<-lbsprsim(
        linf = linf,
        l50= l50,
        l95= l95,
        mk= mk,
        units = "cm",
        sl50= sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        sl95= 1.3*sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"),
        fm= fm(n=length(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],by=by,linf=linf,k=k,l50=l50,l95=l95,mk=mk,year=i,lmean=mean(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp]),lstart=sel(x=dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp],linf=linf,k=k,l50=l50,l95=l95,mk=mk,t0=t0,by_values=c(by,by-1,by-2,by-3,by+1,by+2,by+3),year=i,sp=sp,units="cm",method = "LBSPR"), method="LBSPR",sp=sp, units="cm"), 
        minbin = min(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])-by,
        maxbin = ifelse(max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])<linf,linf+by,max(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp])+by),
        binwidth = by,
        sp=sp)
      
      #store Selectivity and F/M estimated
      sl50=sim2_lb_lb@SL50
      f_m= sim2_lb_lb@FM
      
      #store lbspr frequency data for Chi-Squared test
      freq_lbspr= data.frame(mids=sim2_lb_lb@LMids,counts=sim2_lb_lb@pLCatch*100) %>% 
        tidyr::uncount(round(counts)) 
      freq_lbspr<- freq_lbspr$mids
      
      #storing LBSPR distribution simulation
      sim2_lb_lb<- data.frame(mids=sim2_lb_lb@LMids,counts=norm(sim2_lb_lb@pLCatch) )
      
      # Adjusting to keep the same mids interval
      sim2_lb_lb <- subset(sim2_lb_lb, mids >= min(sim1$mids) & mids <= max(sim1$mids)+1)
      
      #binding the combination
      sim_dat<-data.frame(yr=i,codsp=sp, lmids=sim1$mids, sim1=sim1$counts, sim2=sim2_lb_lb$counts,
        sel_method="LB", fm_method="LB",sl_fm_methods= paste("LB_LB"),sl50=sl50,f_m=f_m,res= sim1$counts-sim2_lb_lb$counts,
        rmse=rmse(sim1$counts,sim2_lb_lb$counts),
        x2= suppressWarnings(chisq.test(sample(dist_sim$fl[dist_sim$yr==i & dist_sim$codsp==sp], 95,replace = FALSE), 
            sample(freq_lbspr, 95,replace = FALSE))$p.value))
  
  comp_sim<-rbind(comp_sim,sim_dat)
}
#assing the comparisons for the current species
comp_sim_wah<-comp_sim

#-------------
#plot section
#-------------
#Length residual comparison 
p61 <- ggplot() +
  geom_boxplot(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(group = interaction(factor(yr), factor(sl_fm_methods)), 
      x = factor(yr), 
      y = res, 
      fill = factor(sl_fm_methods)),  # fill to methods
    alpha = 1,
    position = position_dodge(width = 0.7), 
    width = 0.9,
    lwd = 0.1,  
    outlier.shape = 8,
    outlier.size = 0.7,
    outlier.alpha = 0.5) +
  geom_hline(yintercept =0,linetype = "dashed" )+
  geom_smooth(data = dplyr::filter(comp_sim, codsp == sp),
    method = "loess", span = 0.2, se = FALSE,
    aes(group = 1, x = factor(yr), y = res,color="Loess"), size = 1.5) +
  geom_text(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = "2015", y = 0.9 * max(res, na.rm = TRUE), 
      label = paste("RMSE=", round(mean(rmse, na.rm = TRUE), 2))),
    inherit.aes = FALSE) +
  labs(x = "Year", y = "Length residual", col = "", linetype = '', fill = '') +
  scale_y_continuous() +
  scale_fill_viridis_d(direction = -1)+
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10)))) +
  theme_classic(base_size = 14) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p61

ggsave(paste(c("Length_Residual_Boxplot_Small_Tunas_",sp,".png"),collapse = ''), plot = p61, device = "png", units = "cm",
  width = 25, height = 17)
#------------------------

#Scaterplot comparing the distributions (Dist sim vs LBSPR sim)
p62 <- ggplot() +
  geom_point(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = sim2, y = sim1, color = sl_fm_methods), 
    size = 2, alpha = 0.7) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", size = 1) +
  geom_smooth(data = dplyr::filter(comp_sim, codsp == sp), 
    aes(x = sim2, y = sim1, color = sl_fm_methods, fill = sl_fm_methods), 
    method = "loess", se = FALSE, size = 1.2, alpha = 0.3) +  
  scale_color_viridis_d(direction = -1) +  
  labs(x = "LBSPR Distribution Simulation",y= "Distribution Simulation",col = "Method", shape = "Method", fill = "Method") +
  theme_classic(base_size = 14) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p62

ggsave(paste(c("Distributions_Comparison_Small_Tunas_",sp,".png"),collapse = ''), plot = p62, device = "png", units = "cm",
  width = 24, height = 15)
#---------------------------------------------------------------------------------------------------------------------------
#x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x#


#Combining final data frames

#write the size means observed (programs and literature)
#Write Final data frames in csv files
outfile0  <- "mean_obs.csv"
write.table(mean_obs, file = outfile0, append =FALSE,dec=".",sep = ",",
            row.names = FALSE) 

#write the observved size frequency data (programs)
#Write Final data frames in csv files
outfile1  <- "freq_obs.csv"
write.table(freq_obs, file = outfile1, append =FALSE,dec=".",sep = ",",
  row.names = FALSE) 

#write the lm model (formula+intercept+slope)
#Write Final data frames in csv files
outfile2  <- "lm_model.csv"
lm_model=rbind(lm_model_blf,lm_model_brs,lm_model_dol,
              lm_model_fri,lm_model_kgm,lm_model_lta,lm_model_wah)
write.table(x=lm_model, 
  file = outfile2, append =FALSE,dec=".",sep = ",",
  row.names = FALSE) 

#write the lm model predictions (fit+10%,20%,30%)
#Write Final data frames in csv files
outfile3  <- "lm_predict.csv"
lm_predict=rbind(lm_predict_blf,lm_predict_brs,lm_predict_dol,
                lm_predict_fri,lm_predict_kgm,lm_predict_lta,lm_predict_wah)

write.table(x=lm_predict, 
  file = outfile3, append =FALSE,dec=".",sep = ",",
  row.names = FALSE) 

#write the Gam models (family+link+formula+ AIC+R2)
#Write Final data frames in csv files
outfile4  <- "gam.csv"
gam=rbind(gam_blf,gam_brs,gam_dol,
             gam_fri,gam_kgm,gam_lta,gam_wah)

write.table(x=gam, 
  file = outfile4, append =FALSE,dec=".",sep = ",",
  row.names = FALSE) 

#write the mean length predicted+ used points
#Write Final data frames in csv files
outfile5  <- "mean_pred.csv"
mean_pred=rbind(mean_pred_blf,mean_pred_brs,mean_pred_dol,
  mean_pred_fri,mean_pred_kgm,mean_pred_lta,mean_pred_wah)

write.table(x=mean_pred, 
  file = outfile5, append =FALSE,dec=".",sep = ",",
  row.names = FALSE) 

#write the partial effects of the GAM models
#Write Final data frames in csv files
outfile6  <- "pareffect.csv"
pareffect=rbind(pareffect_blf,pareffect_dol,
  pareffect_fri,pareffect_kgm,pareffect_lta,pareffect_wah)

write.table(x=pareffect, 
  file = outfile6, append =FALSE,dec=".",sep = ",",
  row.names = FALSE) 

#write the length distribution simulation (sizes)
#Write Final data frames in csv files
outfile7  <- "dist_sim.csv"
dist_sim=rbind(dist_sim_blf,dist_sim_brs,dist_sim_dol,
  dist_sim_fri,dist_sim_kgm,dist_sim_lta,dist_sim_wah)

write.table(x=dist_sim, 
  file = outfile7, append =FALSE,dec=".",sep = ",",
  row.names = FALSE) 

#write the distribution simulation comparisons (Distribution Simulation vs LBSPR sim)
#Write Final data frames in csv files
outfile8  <- "comp_sim.csv"
comp_sim=rbind(comp_sim_blf,comp_sim_brs,comp_sim_dol,
  comp_sim_fri,comp_sim_kgm,comp_sim_lta,comp_sim_wah)

write.table(x=comp_sim, 
  file = outfile8, append =FALSE,dec=".",sep = ",",
  row.names = FALSE) 


#----------------------------------------
# Final Plots with all species 
#----------------------------------------

#--- Linear model predictions + percentiles ---# 
p63 <- ggplot(data=dplyr::filter(lm_predict, codsp %in% c("BLF", 'BRS', "DOL", "FRI", "KGM", "LTA", "WAH")), aes(x = yr, y = fit)) +
  geom_line(size = 1, color = "grey20") +  # main line prediction
  # Lines for the percentiles
  geom_line(aes(y = lw_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = up_10, color = "10%", linetype = "10%"), size = 0.8) +
  geom_line(aes(y = lw_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = up_20, color = "20%", linetype = "20%"), size = 0.8) +
  geom_line(aes(y = lw_30, color = "30%", linetype = "30%"), size = 0.8) +
  geom_line(aes(y = up_30, color = "30%", linetype = "30%"), size = 0.8) +
  
  # Observed data with custom filters for each species
  geom_point(data = mean_obs %>%
      dplyr::filter(codsp %in% c("BLF", 'BRS', "DOL", "FRI", "KGM", "LTA", "WAH")) %>%
      dplyr::filter(case_when(
        codsp == "BLF" ~ mean > 50,     # For BLF, mean > 50
        codsp == "DOL" ~ mean > 75,     # For DOL, mean > 75
        codsp == "KGM" ~ mean > 75 & mean < 100, # For KGM, mean between 75 and 100
        codsp == "WAH" ~ mean > 100 & mean < 180, # For WAH, mean between 100 and 180
        TRUE ~ TRUE  # Keep other species as they are 
      )), 
    aes(x = yr, y = mean), size = 2, color = "grey30") +
  
  # Predicted points for 1950, specific to each species
  geom_point(data = dplyr::group_by(lm_predict, codsp) %>% dplyr::slice(1), 
    aes(x = yr, y = lw_10,color = "10%",fill="10%"), size = 3, shape = 21) +
  geom_point(data = dplyr::group_by(lm_predict, codsp) %>% dplyr::slice(1), 
    aes(x = yr, y = up_10,color = "10%",fill="10%"), size = 3, shape = 21) +
  geom_point(data = dplyr::group_by(lm_predict, codsp) %>% dplyr::slice(1), 
    aes(x = yr, y = lw_20,color = "20%",fill="20%"), size = 3, shape = 21) +
  geom_point(data = dplyr::group_by(lm_predict, codsp) %>% dplyr::slice(1), 
    aes(x = yr, y = up_20,color = "20%",fill="20%"), size = 3, shape = 21) +
  geom_point(data = dplyr::group_by(lm_predict, codsp) %>% dplyr::slice(1), 
    aes(x = yr, y = lw_30,color = "30%",fill="30%"), size = 3, shape = 21) +
  geom_point(data = dplyr::group_by(lm_predict, codsp) %>% dplyr::slice(1), 
    aes(x = yr, y = up_30,color = "30%",fill="30%"), size = 3, shape = 21) +
  
  labs(y = 'Mean Length (cm)', x = "Year", color = "Percentile", fill = "Percentile", linetype = "Percentile") +
  facet_wrap(.~codsp, scales = "free") +
  
  # Type and line colors
  scale_color_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_fill_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "30%" = "#009E73")) +
  scale_linetype_manual(values = c("10%" = "dashed", "20%" = "dotted", "30%" = "dotdash")) +
  
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p63
#saving graph...
ggsave("Linear_model_prediction_Small_Tunas_ALLSP.png", plot = p63, device = "png", units = "cm", width = 24, height = 16)
#---------------------------------------------------------------


#Plot residual vs fitted values of GAM models
resvsfit= data.frame(res= c(gam_model_blf$residuals,gam_model_brs$residuals,gam_model_dol$residuals,
  gam_model_fri$residuals,gam_model_kgm$residuals,gam_model_lta$residuals,
  gam_model_wah$residuals),
  fit= c(gam_model_blf$fitted.values,gam_model_brs$fitted.values,gam_model_dol$fitted.values,
    gam_model_fri$fitted.values,gam_model_kgm$fitted.values,gam_model_lta$fitted.values,
    gam_model_wah$fitted.values),
  codsp=c(rep('BLF',length(gam_model_blf$residuals)),
    rep('BRS',length(gam_model_brs$residuals)),
    rep('DOL',length(gam_model_dol$residuals)),
    rep('FRI',length(gam_model_fri$residuals)),
    rep('KGM',length(gam_model_kgm$residuals)),
    rep('LTA',length(gam_model_lta$residuals)),
    rep('WAH',length(gam_model_wah$residuals))))

p64<-ggplot(dplyr::filter(resvsfit,codsp %in% c("BLF",'BRS',"DOL","FRI",
  "KGM","LTA","WAH")), 
  aes(x=fit,y=res),size=1.2) + 
  #geom_line(size=1,aes(fill=source)) +
  geom_smooth(aes(color = "Loess"),size=1.7,method = 'loess' ,
    formula = 'y ~ x',se =FALSE,span = 2,show.legend = FALSE)+
  geom_point(size=1.5)+
  facet_wrap(~codsp,scales = "free")+
  labs(y="Residuals",x="Fitted values",col='',fill="")+
  geom_hline(yintercept=0, linetype=2)+
  #scale_x_continuous(breaks = seq(1950, 2023, 30),limits = c(1950, 2023)) +
  #scale_colour_manual(values = c("gray60","gray40")) +
  scale_fill_viridis_d(direction = -1)+
  theme_classic(base_size = 12) %+replace%
  theme(strip.background=element_blank(),
    plot.margin =unit(c(0.05,0.5,0.05,0.05),"mm"))
p64

ggsave("Residuals_vs_Fitted_Small_Tunas_ALLSP.png",plot=p64, device = "png", units = "cm",
  width = 30, height = 19)
#-------------------------------------------------------------------------


#--- GAM mean length predictions + confidence intervals+ observed mean length points ---#
p65<-ggplot(data=dplyr::filter(mean_pred,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                       "KGM","LTA","WAH")),aes(x=yr,y=pred))+
  geom_ribbon(data=dplyr::filter(mean_pred,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                              "KGM","LTA","WAH")),
              aes(ymin = lw, ymax = up,linetype='Pred'), 
              col="white",fill="grey80",alpha=0.4)+
  geom_line(data=dplyr::filter(mean_pred,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                            "KGM","LTA","WAH")),
            aes(x=yr,y=pred,linetype='Pred'),col="gray60",size=2)+
  geom_point(data=dplyr::filter(mean_pred,codsp %in% c("BLF",'BRS', "DOL","FRI",
                                                             "KGM","LTA","WAH")),
             aes(x=yr,y=mean,shape="Obs"),size=2,col="gray40")+
  facet_wrap(~codsp,scales = "free_y")+
  labs(y='Mean Length (cm)',x="Year",col='Gear',fill='',linetype='',size='',shape='')+
  scale_y_continuous() +
  theme_classic(base_size = 12) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p65

ggsave("Length_model_prediction_Small_Tunas_ALLSP.png",plot=p65, device = "png", units = "cm",
       width = 30, height = 18)
#---------------------------------------------------------------

# ---- plot GAM partial effects ---- #
# Define the levels for "Gear" and "Region"
gear_levels <- c("GN", "HL", "LL", "RR", "TR", "UN", "BB", "PS")
region_levels <- c("N", "NE", "S", "SE", "ALL")

# Combine the levels for "Gear" and "Region"
combined_levels <- c(gear_levels, region_levels)

# Reorder the levels of the `effect` column based on the combined_levels vector
pareffect$effect <- factor(pareffect$effect, levels = combined_levels, ordered = TRUE)

# Ensure that the `var` variable is correct
pareffect$var <- factor(pareffect$var, levels = c("Gear", "Region"))


p66<-ggplot(dplyr::filter(pareffect,codsp %in% c("BLF",'BRS',"DOL","FRI",
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
  scale_colour_viridis_d(,alpha = 0.6)+
  #scale_colour_manual(values=c("gray60","gray40")) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p66

ggsave("Effects_model_Small_Tunas_ALLSP.png",plot= p66, device = "png", units = "cm",
       width = 32, height = 10)
#---------------------------------------------

# Distribution simulation (boxplots)
p67 <- ggplot() +
  
  # distribution simulation
  geom_boxplot(data=dplyr::filter(dist_sim, codsp %in% c("BLF", 'BRS', "DOL", "FRI", "KGM", "LTA", "WAH")), 
    aes(x = factor(yr), y = fl),  
    alpha = 0.5,
    position = "identity",
    lwd = 0.75,
    outlier.shape = 8,
    outlier.size = 0.8,
    outlier.alpha = 0.5,
    fill = "#6A0DAD", 
    color = "gray25",
    show.legend = FALSE) +  
  
  # Observed data with custom filters for each species
  geom_point(data = mean_obs %>%
      dplyr::filter(codsp %in% c("BLF", 'BRS', "DOL", "FRI", "KGM", "LTA", "WAH")) %>%
      dplyr::filter(case_when(
        codsp == "BLF" ~ mean > 50,     
        codsp == "DOL" ~ mean > 75,     
        codsp == "KGM" ~ mean > 60 & mean < 100, 
        codsp == "WAH" ~ mean > 100 & mean < 180, 
        codsp == "FRI" ~ mean < 45,
        TRUE ~ TRUE  
      )), 
    aes(x = factor(yr), y = mean, group = 1), size = 2) +
  
  # Assign another color for the 1950 prediction
  geom_point(data = mean_pred %>%
      dplyr::filter(codsp %in% c("BLF", 'BRS', "DOL", "FRI", "KGM", "LTA", "WAH") & yr == 1950),
    aes(x = factor(yr), y = mean), pch = 8, size = 2, col = "deepskyblue") +
  
  # some additional graph parametrization
  facet_wrap(.~codsp, scales = "free") +
  labs(x = "Year", y = "Fork length (cm)", col = "", linetype = '', fill = '') +
  scale_y_continuous() +
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10))), 
    expand = expansion(mult = c(0.015, 0.015))) +  # Add margin to x-axis
  theme_classic(base_size = 14) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p67

ggsave("Distribution_Simulation_Boxplot_Small_Tunas_ALLSP.png", plot = p67, device = "png", units = "cm",
  width = 30, height = 19)
#------------------------------------------------------


# Mean predicted (line) + Distribution simulation (boxplots)
p68 <- ggplot() +
  
  # distribution simulation
  geom_boxplot(data=dplyr::filter(dist_sim, codsp %in% c("BLF", 'BRS', "DOL", "FRI", "KGM", "LTA", "WAH")), 
    aes(x = factor(yr), y = fl),  
    alpha = 0.7,
    position = "identity",
    lwd = 0.7,
    outlier.shape = 8,
    outlier.size = 0.8,
    outlier.alpha = 0.5,
    fill = "#6A0DAD", 
    color = "gray10", 
    show.legend = FALSE) +  
  
  # Gam mean prediction
  geom_ribbon(data = dplyr::filter(mean_pred, codsp %in% c("BLF", 'BRS', "DOL", "FRI", "KGM", "LTA", "WAH")),
    aes(x = factor(yr), ymin = lw, ymax = up, group = 1, linetype = 'Pred'),
    color = "white", fill = "grey80", alpha = 0.4, show.legend = FALSE) +
  
  geom_line(data = dplyr::filter(mean_pred, codsp %in% c("BLF", 'BRS', "DOL", "FRI", "KGM", "LTA", "WAH")),
    aes(x = factor(yr), y = pred, group = 1, linetype = 'Pred'),
    col = "gray60", size = 1.5, show.legend = FALSE) +
  
  # Observed data with custom filters for each species
  geom_point(data = mean_obs %>%
      dplyr::filter(codsp %in% c("BLF", 'BRS', "DOL", "FRI", "KGM", "LTA", "WAH")) %>%
      dplyr::filter(case_when(
        codsp == "BLF" ~ mean > 50,     
        codsp == "DOL" ~ mean > 75,     
        codsp == "KGM" ~ mean > 60 & mean < 100, 
        codsp == "WAH" ~ mean > 100 & mean < 180, 
        codsp == "FRI" ~ mean < 45,
        TRUE ~ TRUE  
      )), 
    aes(x = factor(yr), y = mean, group = 1), size = 2) +
  
  # Assign another color for the 1950 prediction
  geom_point(data = mean_pred %>%
      dplyr::filter(codsp %in% c("BLF", 'BRS', "DOL", "FRI", "KGM", "LTA", "WAH") & yr == 1950),
    aes(x = factor(yr), y = mean), pch = 8, size = 2, col = "deepskyblue") +
  
  # some additional graph parametrization
  facet_wrap(.~codsp, scales = "free") +
  labs(x = "Year", y = "Fork length (cm)", col = "", linetype = '', fill = '') +
  scale_y_continuous() +
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10))), 
                    expand = expansion(mult = c(0.015, 0.015))) +  # Add margin to x-axis
  theme_classic(base_size = 12) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p68

ggsave("Mean_Length_Distribution_Simulation_boxplot_Small_Tunas_ALLSP.png", plot = p68, device = "png", units = "cm",
                        width = 30, height = 19)
#------------------------------------------------------


#plot half violin for the Distribution Simulation
p69<- ggplot() +
  gghalves::geom_half_violin(data=dplyr::filter(dist_sim,codsp %in% c("BLF",'BRS',"DOL","FRI",
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
  geom_point(data=dplyr::filter(mean_pred,codsp %in% c("BLF",'BRS',"DOL","FRI",
                                                        "KGM","LTA","WAH")),
             aes(x=factor(yr),y=mean),size=2)+
  # Assign another color for the 1950 prediction
  geom_point(data = mean_pred %>%
      dplyr::filter(codsp %in% c("BLF", 'BRS', "DOL", "FRI", "KGM", "LTA", "WAH") & yr == 1950),
    aes(x = factor(yr), y = mean), pch = 8, size = 2, col = "deepskyblue") +
  # 
  labs(x="Year",y="Fork length (cm)")+
  #scale_y_continuous(limits = c(0.9*min(szfreqpred$fl[szfreqpred$codsp==sp]),
  #                              1.1*max(szfreqpred$fl[szfreqpred$codsp==sp]))) +
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10))), 
    expand = expansion(mult = c(0.015, 0.015))) + 
  scale_fill_viridis_d(direction = 1) +
  theme_classic(base_size = 14) %+replace%
  theme(strip.background=element_blank(),
        plot.margin =unit(c(0.05,0.05,0.05,0.05),"mm"))
p69

ggsave("Distribution_Simulation_HalfViolin_Small_Tunas_ALLSP.png",plot=p69, device = "png", units = "cm",
        width = 35, height = 12)
#------------------------------------------------------


#Length residual comparison (Residual and Loess)- All methods 
p70 <- ggplot() +
  
  # Boxplot
  geom_boxplot(data = dplyr::filter(comp_sim, codsp %in% c("BLF", 'BRS', "DOL", "FRI", "KGM", "LTA", "WAH")), 
    aes(x = factor(yr), y = res),  
    alpha = 1,
    position = position_dodge(width = 0.7), 
    width = 0.9,
    lwd = 0.7,  
    outlier.shape = 8,
    outlier.size = 0.5,
    outlier.alpha = 0.3,
    fill = "#6A0DAD", 
    color = "gray10") +
  
  # Line for zero residual
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # Loess smoothing
  geom_smooth(data = dplyr::filter(comp_sim, codsp %in% c("BLF", 'BRS', "DOL", "FRI", "KGM", "LTA", "WAH")),
    method = "loess", span = 0.2, se = FALSE,
    aes(group = 1, x = factor(yr), y = res, color = "Loess"), size = 1.5,show.legend = FALSE) +
  
  # Facet by species
  facet_wrap(.~codsp, scales = "free") +
  
  # RMSE text by species
  geom_text(data = comp_sim %>%
      dplyr::filter(codsp %in% c("BLF", 'BRS', "DOL", "FRI", "KGM", "LTA", "WAH")) %>%
      dplyr::group_by(codsp) %>%
      dplyr::summarise(rmse_mean = round(mean(rmse, na.rm = TRUE), 2), 
        res_max = max(res, na.rm = TRUE)),  # Capture max residual to set y position
    aes(x = "2015", y = 0.95 * res_max,  # Use max residual for y position
      label = paste("RMSE=", rmse_mean)),
    inherit.aes = FALSE) +

  # Labels and scales
  labs(x = "Year", y = "Length residual", col = "", linetype = '', fill = '') +
  scale_y_continuous() +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(drop = FALSE, breaks = c(as.character(seq(1950, 2025, 10)))) +
  
  # Theme
  theme_classic(base_size = 12) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"))

p70

ggsave("Length_Residual_Boxplot_Small_Tunas_ALLSP.png", plot = p70, device = "png", units = "cm",
      width = 28, height = 19)
#----------------------------------------------------------------



#Scaterplot comparing the distributions (Dist sim vs LBSPR sim)
# Calculating the R² for each specie
r2 <- comp_sim %>%
  filter(codsp %in% c("BLF", "BRS", "DOL", "FRI", "KGM", "LTA", "WAH")) %>%
  group_by(codsp) %>%
  summarize(r_squared = summary(lm(sim1 ~ sim2))$r.squared)



#  #Scaterplot comparing the distributions (Dist sim vs LBSPR sim)
p71 <- ggplot() +
  geom_point(data = dplyr::filter(comp_sim, codsp %in% c("BLF", "BRS", "DOL", "FRI", "KGM", "LTA", "WAH")), 
    aes(x = sim2, y = sim1, color = sl_fm_methods), 
    size = 2, alpha = 0.6) +  
  geom_smooth(data = dplyr::filter(comp_sim, codsp %in% c("BLF", "BRS", "DOL", "FRI", "KGM", "LTA", "WAH")), 
    aes(x = sim2, y = sim1, color = sl_fm_methods, fill = sl_fm_methods), 
    method = "loess", se = FALSE, size = 1.2, alpha = 0.3) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "tomato1", size = 1) +
  facet_wrap(.~codsp, scales = "free") +
  scale_color_viridis_d(direction = -1) +  
  labs(x = "LBSPR Distribution Simulation", y = "Distribution Simulation", col = "Method", shape = "Method", fill = "Method") +
  theme_classic(base_size = 10) %+replace% 
  theme(strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm")) +
  # add R2 data
  geom_text(data = r2, 
    aes(x = Inf, y = -Inf, label = paste("R² = ", round(r_squared, 2))),
    hjust = 1.1, vjust = -0.5, size = 3, color = "black", inherit.aes = FALSE)

p71

ggsave("Distributions_Comparison_Scatterplot_Small_Tunas_ALLSP.png", plot = p71, device = "png", units = "cm",
                                                 width = 23, height = 15)
#x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x#
# end of reconstructions
#x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x##x#x#