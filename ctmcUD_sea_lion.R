### Steller Sea Lion UD

library(crawl) # devtools::install_github("NMML/crawl") for latest version
library(ctmcmove)
library(gdistance)
library(lubridate)
library(raster)
library(dplyr)
library(purrr)
library(readr)


# 1) Load the data
#------------------------------------------------
ssl_data <- readr::read_csv("sea_lion_telemetry.csv") %>% mutate(GMT=mdy_hm(GMT)) # telemetry data
aleut_hab <- brick("habitat/aleut_habitat.grd", values=TRUE) %>% stack()   # habitat covariates
land <- 1.0*(aleut_hab$bathy>=0)

# Project data telemetry data
dat_toproj <- filter(ssl_data, !is.na(latitude))
coordinates(dat_toproj) <- ~longitude+latitude
proj4string(dat_toproj) <- CRS("+proj=longlat")
dat_toproj <- spTransform(dat_toproj, CRS("+init=epsg:3338")) %>% 
  as.data.frame() %>% rename(x=longitude, y=latitude) %>% 
  dplyr::select(Deploy_ID, GMT, x, y)
ssl_data = full_join(ssl_data, dat_toproj, by=c("Deploy_ID", "GMT"))

###################################
# Fit CTMC UD to first animal    ##
###################################
temp_dat <- dplyr::filter(ssl_data, Deploy_ID==15137) %>% arrange(GMT) #pull out one sea lion

#set Argos error & add a constraint parameter for estimating Argos error
# in CTCRW imputation model
temp_dat$argos_class = factor(temp_dat$argos_class, levels=c("3","2","1","0","A","B"))

# Define 'haul' variable for determining when the animal is hauled out. These will
# be removed in the CTMC analysis
temp_dat %>%  mutate(
                  haul = (
                    DryTime==1 & c(0,diff(DryTime==1)==0) |
                      DryTime==0 & c(0,diff(DryTime==1)==-1)
                    )
                  ) -> temp_dat

## Initial state values from CTCRW model in crawl (for imputation of raw telemetry)
initial <- list(a=c(temp_dat$x[1],0, temp_dat$y[1],0),
                P=diag(c(10000^2,5400^2,10000^2,5400^2)))
# Fixed parameter values for CTCRW imputation model
fixPar = c(log(250), log(500), log(1500), rep(NA,5), 0) 
# Lower bounds for parameter estimates
constr = list(lower=c(rep(log(1500),3),rep(-Inf,2)), upper=rep(Inf,5))
# Check that parameters are as expected
crawl::displayPar( mov.model=~1, err.model=list(x=~argos_class-1),data=temp_dat,
            activity=~I(1-DryTime),fixPar=fixPar) #with dry time

# 3) Run crawl model with dry (haul-out/activity) data
# --------------------------------------------------------------------
set.seed(123)
mi_fit <- crawl::crwMLE(
  mov.model=~1, err.model=list(x=~argos_class-1), activity=~I(1-DryTime),
  data=temp_dat, coord=c("x","y"), Time.name="GMT",
  initial.state=initial, fixPar=fixPar,
  constr=constr, #prior=ln.prior,
  method="L-BFGS-B",
  control=list(maxit=2000, trace=1, REPORT=10),
  initialSANN=list(maxit=250, temp=10, trace=1, REPORT=10))


# Create simulation object for multiple imputation
set.seed(123)
temp_dat %>% mutate(
  idx = (!haul & minute(GMT)%in%c(0, 20, 40))
  ) -> temp_dat

simObj <- crawl::crwSimulator(mi_fit, parIS = 0)

# 4) Fit CTMC Model
# -----------------------------------------------------------------
P = 10
glm_data = NULL
for(i in 1:P){
  samp <- crawl::crwPostIS(simObj, fullPost = FALSE) %>% pluck(1) %>% .[temp_dat$idx,c(1,3)] %>%  # pull out xy coords and time
    as.data.frame() %>% mutate(t=1:n())
  
  # To substantially improve processing time we suggest cropping your rasters to include only the cells that surround each
  # imputed path
  crop_lim = extent(c(range(samp$mu.x)+c(-5000,5000), range(samp$mu.y)+c(-5000,5000)))
  grad.stack_crop = stack(crop(aleut_hab, crop_lim))                     # raster stack of gradient covariates
  names(grad.stack_crop) = paste0(names(grad.stack_crop), "_grad")
  loc.stack_crop = stack(crop(aleut_hab, crop_lim))                      # raster stack of motility covariates
  names(loc.stack_crop) = paste0(names(loc.stack_crop), "_loc")
  land_crop = crop(land,crop_lim) # grad.stack_crop[[1]]==0              # barrier cells
  water_crop = crop(1-land,crop_lim)
  holes_crop = which(land_crop@data@values==1)                           # index of barrier cells
  trans_crop = gdistance::transition(water_crop, prod,4)                            # transition matrix 
  newpath <- crawl::fix_path(as.matrix(samp[,1:2]), samp[,3], land_crop, trans_crop) %>% 
    data.frame(.)  # path that does not cross barriers
  
  path <- list(xy=as.matrix(newpath[,1:2]), t=as.vector(newpath$time))
  ctmc <- ctmcmove::path2ctmc(xy=path$xy,t=path$t, rast=grad.stack_crop, zero.idx = holes_crop) # extract discrete path & cell residence times
  glm_data <- rbind(glm_data, ctmcmove::ctmc2glm(ctmc, loc.stack_crop, grad.stack_crop, 
                                    zero.idx=holes_crop)) 
  message(paste('i =',i,Sys.time()))
  message(nrow(glm_data))
} 

 write.csv(glm_data,'glm.data_SSL-15137.csv',row.names=FALSE)

#fit GLM
fit <- glm(z~bathy_grad+slope_grad+d2site_grad+d2shelf_grad+
             bathy_loc+slope_loc+d2site_loc+d2shelf_loc,
           weights=rep(1/P,nrow(glm_data)),family="poisson",
           offset=log(tau),data=glm_data)
summary(fit)

# 5) Calculate a UD from the CTMC output
# --------------------------------------------------------------

# Set up Rate matrix
holes <- which(raster::getValues(land)==1)
grad.stack <- loc.stack <- aleut_hab
names(grad.stack) = paste0(names(grad.stack), "_grad")
names(loc.stack) = paste0(names(loc.stack), "_loc")

# Use rate matrix to calculate the UD
R = ctmcmove::get.rate.matrix(fit, loc.stack, grad.stack, zero.idx=holes)  #grad.stack0 & loc.stack0 do NOT have NA values for land
pi = ctmcmove::get.UD(R, method = "limit", maxiter = 10000)
UD.rast=land
values(UD.rast) <- pi
plot(log(UD.rast))

# to view UD in ggplot
ud1df <- as.data.frame(raster::rasterToPoints(UD.rast)); names(ud1df) <- c('x','y','Density')
ud1 <-  ggplot(ud1df,aes(x/1000,y/1000)) + geom_raster(aes(fill=log(Density))) + coord_equal() +
  scale_fill_gradientn(colours=viridis(100))+
  theme_light()+
  labs(x='Easting',y='Northing',title='Sea Lion: 14809',fill='Density')+
  theme(axis.ticks.y = element_blank(), 
        axis.ticks.x = element_blank(),
        #legend.position=c(0.90,0.73),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        panel.grid = element_blank())

# 5.1) Procedure for estimating the variance of the UD can increase reps for real analysis
# library(mvtnorm)

# tmp = fit
# V = vcov(fit)
# b = coef(fit)
# M2 <- m <-  0*pi
# reps = 10
# for(i in 1:10){
#   tmp$coefficients = rmvnorm(1, b, V)
#   R_tmp = get.rate.matrix(tmp, loc.stack, grad.stack, zero.idx=holes) 
#   pi_tmp = get.UD(R_tmp, method = "limit", maxiter = 10000)
#   # pi_tmp = get.UD(R_tmp)
#   delta = pi_tmp - m
#   m = m + delta/i
#   delta2 = pi_tmp - m
#   M2 = M2 + delta*delta2
# }
# UD.se.rast = land
# values(UD.se.rast) = sqrt(M2/(reps-1))



