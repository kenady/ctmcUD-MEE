### Steller Sea Lion UD

library(crawl) # devtools::install_github("NMML/crawl") for latest version
library(ctmcmove)
library(lubridate)
library(dplyr)
library(purrr)
library(readr)


#1) Load the data
#------------------------------------------------
ssl_data = read_csv("sea_lion_telemetry.csv") %>% mutate(GMT=mdy_hm(GMT))
aleut_hab = brick("habitat/aleut_habitat.grd", values=TRUE) %>% stack()
land = 1.0*(aleut_hab$bathy>=0)

#2) Project data
dat_toproj = filter(ssl_data, !is.na(latitude))
coordinates(dat_toproj) = ~longitude+latitude
proj4string(dat_toproj) <- CRS("+proj=longlat")
dat_toproj <- spTransform(dat_toproj, CRS("+init=epsg:3338")) %>% 
  as.data.frame() %>% rename(x=longitude, y=latitude) %>% 
  dplyr::select(Deploy_ID, GMT, x, y)
ssl_data = full_join(ssl_data, dat_toproj, by=c("Deploy_ID", "GMT"))

###################################
# Fit CTMC UD to first animal    ##
###################################
temp_dat <- filter(ssl_data, Deploy_ID==15136) %>% arrange(GMT) #pull out one sea lion

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
displayPar( mov.model=~1, err.model=list(x=~argos_class-1),data=temp_dat,
            activity=~I(1-DryTime),fixPar=fixPar) #with dry time

# Run crawl model for SDR tags w/stop model
set.seed(123)
mi_fit <- crwMLE(
  mov.model=~1, err.model=list(x=~argos_class-1), activity=~I(1-DryTime),
  data=temp_dat, coord=c("x","y"), Time.name="GMT",
  initial.state=initial, fixPar=fixPar,
  constr=constr, #prior=ln.prior,
  method="L-BFGS-B",
  control=list(maxit=2000, trace=1, REPORT=10),
  initialSANN=list(maxit=250, temp=10, trace=1, REPORT=10))

# # Create predicted path
# predTimes <- seq(
#   ceiling_date(min(temp_dat$datetime), "hour"), 
#   floor_date(max(temp_dat$datetime),"hour"), 
#   by = "1 hour")
# 
# pred <- crawl::crwPredict(temp.fit, predTime=predTimes) %>% 
#   mutate(GMT = as.POSIXct1970(TimeNum)) %>% 
#   mutate(DryTime = zoo::na.locf(DryTime)) %>% 
#   mutate(haul = (
#     DryTime==1 & c(0,diff(DryTime==1)==0) |
#       DryTime==0 & c(0,diff(DryTime==1)==-1)
#   )
#   ) %>% 
#   select(TimeNum, GMT, mu.x, mu.y, x.y, y.y, haul, locType)
# 
# #Simulate
# set.seed(123)
# idx = !pred$haul & pred$locType=="p"
# simObj <- crwSimulator(temp.fit, method="quadrature", parIS = 0, predTime = predTimes)

## This is were multiple imputation loops start
#------------------------------------------------------------------------------

# #get path
# set.seed(123)
# samp <- crwPostIS(simObj, fullPost = FALSE)
# samp <- cbind(samp[[1]][,c(1,3)],t=samp[[3]]) #pull out xy coordinates and time
# samp = samp[idx,]
# samp[,3] = seq(0,by=1,length=nrow(samp))
# 
# crop_lim = extent(c(range(samp[,"mu.x"])+c(-5000,5000), range(samp[,"mu.y"])+c(-5000,5000)))
# 
# grad.stack_crop = stack(crop(grad.stack, crop_lim))
# loc.stack_crop = stack(crop(loc.stack, crop_lim))
# land_crop = crop(land,crop_lim) #grad.stack_crop[[1]]==0
# water_crop = crop(water,crop_lim)
# # holes_crop = which(land_crop@data@values==1)
# holes_crop <- which(getValues(land_crop)==1)
# 
# trans_crop = transition(water_crop, prod,4)
# newpath <- fix_path(samp[,1:2], samp[,3], land_crop, trans_crop) %>% data.frame(.)
# 
# P=10
# path <- list(xy=as.matrix(newpath[,1:2]), t=as.vector(newpath$time))
# ctmc=path2ctmc(xy=path$xy,t=path$t, rast=grad.stack_crop, zero.idx = holes_crop)
# glm.data = ctmc2glm(ctmc, loc.stack_crop, grad.stack_crop, zero.idx=holes_crop) #use rasters with holes for land (grad.stack & loc.stack)

# Create simulation object for multiple imputation
set.seed(123)
temp_dat %>% mutate(
  idx = (!haul & minute(GMT)%in%c(0, 20, 40))
  ) -> temp_dat
# idx = !ssl_data$haul & pred$locType=="p"
simObj <- crwSimulator(mi_fit, parIS = 0)


P = 10
glm_data = NULL
for(i in 1:P){
  samp <- crwPostIS(simObj, fullPost = FALSE) %>% pluck(1) %>% .[temp_dat$idx,c(1,3)] %>%  # pull out xy coords and time
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
  trans_crop = transition(water_crop, prod,4)                            # transition matrix 
  newpath <- fix_path(as.matrix(samp[,1:2]), samp[,3], land_crop, trans_crop) %>% data.frame(.)  # path that does not cross barriers
  
  path <- list(xy=as.matrix(newpath[,1:2]), t=as.vector(newpath$time))
  ctmc=path2ctmc(xy=path$xy,t=path$t, rast=grad.stack_crop, zero.idx = holes_crop) # extract discrete path & cell residence times
  glm_data=rbind(glm_data, ctmc2glm(ctmc, loc.stack_crop, grad.stack_crop, zero.idx=holes_crop)) 
  message(paste('i =',i,Sys.time()))
  message(nrow(glm_data))
} 

#fit GLM
fit <- glm(z~bathy_grad+slope_grad+d2site_grad+d2shelf_grad+
             bathy_loc+slope_loc+d2site_loc+d2shelf_loc,
           weights=rep(1/P,nrow(glm_data)),family="poisson",
           offset=log(tau),data=glm_data)
summary(fit)

# Set up Rate matrix
holes=which(land@data@values==1)
grad.stack <- loc.stack <- aleut_hab
names(grad.stack) = paste0(names(grad.stack), "_grad")
names(loc.stack) = paste0(names(loc.stack), "_loc")

R = get.rate.matrix(fit, loc.stack, grad.stack, zero.idx=holes)  #grad.stack0 & loc.stack0 do NOT have NA values for land
pi = get.UD(R, method = "limit", maxiter = 10000)
UD.rast=land
values(UD.rast) <- pi
plot(log(UD.rast))

# Procedure for estimating variance of UD
# can increase reps for real analysis
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



