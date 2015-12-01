
##### INDIVIDUAL TRANSITION FUNCTIONS########################################

####### OVERALL GERMINATION FUNCTION #######################################
## NextGerm is a function to move the seed dormancy distribution (psi) and the hydro thermal unit (HTU) distribuion one step foward in time
#
# Parameters  	- HTU - matrix giving accumulated hydrothermal time for each seed
#           		   	- psi- matrix giving psi (seed dormancy) values
#            			- Tb.d -  base temperature for germination
#            			- To - optimal temperature for germination
#            			- kT - dormancy increase in psi for each degree C above To
#            			- theta.HT - threshold for germination
#            			- time.step - numeric describing timestep in hours (defaults to 1)
#            			- temp - vector of temperatures at scale of timestep
#           			- moisture - vector of moistures at scale of timestep
#            			- is.wet - threshold for after-ripening WC.sat germ progress
#            			- psi.mean -  mean dormancy level of cohort at dispersal
#           			- dorm.breadth- initial difference between lowest and highest dormancy classes at dispersal
#            			- psi.min  - minimum dormancy possible
#            			- dsat- Number of days it takes a seed to go from psi 0 to psi -1
#						- psi.scalar - scalar for dormancy loss   ADD TO FUNCTIONS!!!!
#						-Tsl -
#						-psi.ll - lower moistuer limit for after ripening
#						-psi.ul - upper moisture limit for after ripening
#
#
# Returns list including psi and HTU (matrix of psi + matrix of  HTU for every cohort in columns),
#                        
# all one time-step later


nextGerm <- function(HTU,
		     		psi,
                     Tb.d,
                     To,
                     kT,
                     theta.HT,
                     Tsl,
                     psi.ll,
                     psi.ul,
                     time.step,
                     temp,
                     moisture,
                     is.wet,
                     n.cohort,
                     psi.mean,
                     dorm.breadth,
                     psi.min,
                     dsat,
                     psi.scalar) {
    
    if (!is.wet){   #If it is drier than the set threshold for wetness, accumulate afteripening according to AR.loss
        psi<-AR.loss(psi=psi, dsat=dsat, psi.min=psi.min, temp=temp, moisture=moisture, Tsl=Tsl, psi.ll=psi.ll, psi.ul=psi.ul, psi.scalar=psi.scalar)
    }  else {    # Otherwise use the function accum.HTU to generate a new matrix of HTT's based on the psi matrix  
        add.HTT <- accum.HTU(psi=psi,temp=temp,moisture=moisture,To=To,Tb.d=Tb.d,kT=kT,time.step=time.step)
        HTU <- HTU+add.HTT   
    }
    return(list(HTU=HTU, psi=psi))   
}



####### GERMINATION SUB-FUNCTIONS ##### 

#########hydrothermal unit accumulation############
# Calculates hydrothermal units accumulated in a single time step based on seed dormancy levels and environmental factors in that hour
#
# Parameters 	- psi
#            			- temp
#            			- moisture 
#            			- To - 
#            			- Tb.d 
#            			- kT
#            			- time.step 
#
# Returns matrix of new HTUs that occured  in a single timestep to be added on to the current accumulated matrix 
#
accum.HTU <- function(psi,  temp, moisture, To, Tb.d,kT,  time.step) {

    added.HTU <- matrix(data=rep(NA,length(psi)),nrow=nrow(psi),ncol(psi)) 
    
    added.HTU[temp<=To & moisture<psi] <- 0
    added.HTU[temp<=Tb.d] <- 0
         # This adds 0 progress toward germination if the base water potential of the fraction is higher than the ambient base water potential or the temperature is lower than the base temperature
    added.HTU[temp<=To & moisture>=psi & temp>Tb.d] <- accum.sub (moisture=moisture,psi=psi[temp<=To & moisture>=psi],temp=temp,Tb.d=Tb.d,time.step=time.step)
         # This equation adds progress towards germ if the temperature is in the sub optimal range
    added.HTU[temp>To & 0<(moisture-psi-kT*(temp-To))] <- accum.supra(moisture=moisture,psi=psi[temp>To & 0<(moisture-psi-kT*(temp-To))],
                                                                      temp=temp,Tb.d=Tb.d,To=To,kT=kT,time.step=time.step)
         # This equation adds progress toward germination in the supra optimal range
    added.HTU[(moisture-psi-kT*(temp-To))<0] <- 0
         # at really high temperatures (at which the base water potential is above the ambient water potential) there is no movement toward germination
    return(added.HTU)
}


###### sub and supra optimal HTU calculations from Alvarado et al paper#####
# sub optimal temp accumulation function
accum.sub <-  function(moisture,psi,temp,Tb.d,time.step) (moisture-psi)*(temp-Tb.d)*time.step   
# supra optimal temp accumulation function
accum.supra <-  function (moisture,psi,temp,Tb.d,To,kT,time.step) (moisture-psi-kT*(temp-To))*(To-Tb.d)*time.step

############# After-ripening functions/ Dormancy loss functions ##################
# function to calculate the AR.HTU accumulated in an hour
# this function is used by the AR.loss function to determine how afteripening progress translates into dormancy loss 
#
# Parameters	 - temp 
#            			- moisture 
#            			- Tsl 
#            			- psi.ll 
#            			- psi.ul 
#
# Return - after-ripening time that hour
#  remember  after ripening only happens if the environment is drier than the is.wet threshold!!!!!!

AR.HTU <- function (temp,moisture,Tsl,psi.ll,psi.ul) {
    if (moisture<psi.ll & temp > Tsl) x <- 0  #too dry to after-ripen
    if (moisture>=psi.ll & moisture<=psi.ul & temp > Tsl) x <- ((psi.ll-moisture)/(psi.ll-psi.ul))*(temp-Tsl) #as env is wetter AR rate increases.....
    if (moisture>psi.ul & temp > Tsl) x <- temp-Tsl # up to a plateau where moisture stops mattering
    if (temp <= Tsl) x <- 0  # if Temp is below the low threshold, After-ripening doesn't happen
    return(as.numeric(x))
}

##################  Dormancy loss ######################
#function that translates an hourly after-ripening amount into dormancy loss
# Parameters 	-dsat 
#						-psi
#						-psi.min
#						-psi.scalar
#						-all the other parameters are passed to the AR.HTU function and are listed above in the main germination function
#
# Returns - a matrix of new psi for all currently existing classes of seeds

AR.loss<-function(psi=psi,dsat=dsat,psi.min=psi.min,temp=temp,moisture=moisture,Tsl=Tsl,psi.ll=psi.ll,psi.ul=psi.ul,psi.scalar=psi.scalar) {
    AR.Saturate<-AR.HTU(temp=20,moisture=-200,Tsl=3,psi.ll=-350,psi.ul=-50)*24*dsat
    psi.temp<-psi-(AR.HTU(temp=temp,moisture=moisture,Tsl=Tsl,psi.ll=psi.ll,psi.ul=psi.ul)*(psi.scalar/AR.Saturate))
    psi.temp[which(psi.temp<psi.min,arr.ind=TRUE)] <- psi.min
    return(psi.temp)
}


####################################################################################

############OVERALL FLOWERING FUNCTION####################################### 
# NextFlower is a function that advances the vector of winter chilling accumation (WC) and photothermal unit (PTU) accumulation by one time step.
#
# Parameters 	- PTU - current vector of accumulated photothermal units
#            			- WC - current vector of winter chilling/vernalization accumulation
#            			- temp - temperature
#            			- day_hour
#            			- day_length
#            			- Tb.f - base temperature for flowering
#			 			- ds - critical short photoperiod
#						- ps - rate of development at critical short photoperiod
#			 			- dl - critical long photoperiod
#						- pl - rate of development at critical long photoperiod
#						- Fi - Inital floral repression
#						- Fu - floral repression ast WC.sat
#						- WC.sat- winter chilling saturation point
#						- Tv.min - temperature minimum for winter chilling
#						- Tv.max - temperture maximum for winter chilling
#						- kappa - parameter for shaper of winter chilling effectiveness
#						- omega - parameter for shaper of winter chilling effectiveness
#						- xi - parameter for shaper of winter chilling effectiveness

# Returns - a list containing the vector of new PTU and WC values

nextFlower <- function( PTU,
                        WC,
                        PTU.growth,
                        temp=4,
                        day_length=7,
                        day_hours=1,
                        Tb.f=3,ds=10,ps=0.626,dl=14,pl=1,WC.sat=960,Fi=0.598,Fu=0,
                        Tv.min=-3.5,Tv.max=6,
                        kappa=-5.17485563597551, omega=2.22568888596679, xi=.995904471970916
                        ) {
    
    PTU <- PTU+hourly.progress(temp=temp,day_length=day_length,day_hours=day_hours,
                               Tb.f=Tb.f,d=d,ds=ds,ps=ps,dl=dl,pl=pl,WC.sat=WC.sat,Fi=Fi,Fu=Fu,wc=WC)
    WC  <- WC+as.numeric(acc.vern(temp=temp,Tv.min=Tv.min, Tv.max=Tv.max, kappa=kappa, omega=omega, xi=xi))
    
    PTU.growth<-hourly.progress.grow(temp=temp,day_length=day_length,day_hours=day_hours,
                               Tb.f=Tb.f,d=d,ds=ds,ps=ps,dl=dl,pl=pl,WC.sat=WC.sat)

    return(list(PTU=PTU,WC=WC,PTU.growth=PTU.growth))

}



######################SUB FUNCTIONS FOR FLOWERING ############################
# The accumulation of photo-thermal units (PTU's) is determined by the muliplication of of photo and thermal components.
# PTU's are multiplied by a winter chilling component which takes into account how much winter chilling (FLC repression) has happened so far.
# Therefore MPTU (for a given hour) = Temp*Photo*Chilling

######### The thermal component ##############
# Function that defines the amount of thermal progress a plant accumulates in an hour. 
# This function only is used during the day-time
#
# Parameters - temperature
#            - day_hours
#            - Tb.f - base temp
#
# Returns Temp factor for that hour

Temp <- function(temp,day_hours,Tb.f) {
    x <- 0
    if(temp >= Tb.f) x <-day_hours*(temp-Tb.f)
    if (temp <  Tb.f ) x <- 0
    return (as.numeric(x))
}

######## The daylength component function########
# Function that defines the amount of photoperiod progress a plant accumulates in an hour.
#
# Parameters	- day_length
#			 			- ds - critical short day length
#			 			- ps - progress in short days
# 			 			- dl - critical long day length
#					 	- pl - progress in long days
#
# Returns daylength factor for that hour


Photo <- function(day_length, ds, ps, dl, pl) {
    d <- day_length; 
    if (d<ds) x <- ps
    if (ds <= d & dl> d) x <- ps+(pl-ps)*(d-ds)/(dl-ds)
    if (d>=dl) x <- pl
    return (as.numeric(x))
}

######### Winter chilling modifier function #########
#Function that describes how the current floral repression state depends on the amount of accumulated winter chilling
#
# Parameters 	- WC.sat - winter chilling saturation level, 
#			 			- Fi - initial floral repression
#			 			- Fu - floral repression at saturated vernalization level
# 			 			- wc - a vector of current WC levels for each flowering cohort
#			 
# Returns the vernalization modifier 

Chilling <- function(wc,WC.sat, Fi, Fu) {
    ifelse (wc<WC.sat, 1-Fi+(Fi-Fu)*(wc/WC.sat),1-Fu)
}

######## Winter chilling accumulation function######
# This evaluates how much vernalization progress happens in a given hour.
#
# Parameters 	- Tv.min - the minimum temperature for vernalization
#			 			- Tv.max - the maximum temperature for vernalization
# 			 			- kappa, omega, xi - parameters for gamma distribution for vernalization temperature effectiveness 
#
# Returns- vector of values to add onto the wc vector which keeps track of total winter chlling accumulation for the plant			 

acc.vern <- function(temp, Tv.min, Tv.max, kappa, omega, xi) {
  ifelse (temp>=Tv.min & temp<=Tv.max, exp(kappa)*((temp-Tv.min)^omega)*((Tv.max-temp)^xi),0)
}  # This acumulates hourly. 40 days at 4C 24 hours a day is 960 which is saturating vernalization for Col...

################# Putting the flowering functions together#######################
# MPTU accumulation function: Putting the three together Temp*Photo*WC 
# Passes all the arguments into the above compnent functions.
# Returns a vector of PTU values to be added to the PTU accumulation vectors 
 
hourly.progress <- function(temp,day_length,day_hours,Tb.f,d,ds,ps,dl,pl,WC.sat,Fi,Fu,wc) {
    as.numeric(Temp(temp=temp,day_hours=day_hours,Tb.f=Tb.f)*
               Photo(day_length=day_length,ds=ds,ps=ps,dl=dl,pl=pl)*Chilling(wc=wc,WC.sat=WC.sat,Fi=Fi,Fu=Fu))             
}

hourly.progress.grow <- function(temp,day_length,day_hours,Tb.f,d,ds,ps,dl,pl,WC.sat,Fi,Fu,wc) {
    as.numeric(Temp(temp=temp,day_hours=day_hours,Tb.f=Tb.f))             
}


###########OVERALL SEED RELEASE FUNCTION#######################################
# Seed Release is a function of thermal time accumulation. Based on amount of thermal time it takes 10% of seeds to mature in chambers at 3 constant temperatures
# Parameters 	- SDU - vector of current Seed development units for current seed cohorts
#			 			- temp - temperature in a given hour
# 			 			- seed.release.params ["Tb.d"] - base temperature for fruit maturation

nextSeed <- function(simple.SDT,SDU,seed.release.params,t,temp,moisture) {
     x <- 0
    if(temp >= seed.release.params["Tb.d"]) x <-(temp-seed.release.params["Tb.d"])/24
    if (temp < seed.release.params["Tb.d"] ) x <- 0
    SDU<- SDU + x
    return (SDU)
}


###############Demographic functions (survival and fecundity)##############
# For now, just returns 1
seed.survival <- function(seed.surv.params,temp,moisture) {
    surv <- seed.surv.params[1]
    return(surv)
}

# For now, just returns 1
rosette.survival <- function(rosette.surv.params,temp,moisture, ptu,size) {
    surv <- rosette.surv.params[1]
    return(surv)
}

# constant seed production rate if you survive
n.seeds.simple <- function(rosette.seed.params, size,PTU.growth,age) {
    seeds <- rosette.seed.params[1]
    return(seeds)
}


###############################################################################
###############################################################################
###############MAIN WORKHORSE FUNCTION-ITERATEPOP#########################

## Function to link together all these functions and iterate the population using environmental inputs
#
#
# Parameters 	- Tmax - total number of time-steps
#           			- store.all- if FALSE do not include the huge matrices that keep track of hydro/photo/thermal development each day
#                      - start.PTU - set intial PTU level for first generation
#                      - start.WC- set intial PTU level for first generation
#                      - start.SDU - set intial SDU level for first generation
#                      - start.size - initiates seedling size for growth tracking (not used currently)
#                      - dispersal.day - sets the day of the year for seeds to be dispersed
#                      - max.cohort.store - set the maximum size of big matrices
#                      - do.flower.only- If TRUE restricts the loop to calculate the response of flowering to the environment indpendent of other life stages
#                       -do.germ.only - If TRUE restricts the loop to calculate the response of germination to the environment indpendent of other life stages
#
# All other parameters were described in the functions above and are set to the defaults defined in the paper

iteratePop <- function(Tmax,
                       time.step,
                       #environment
                       temp,
                       moisture,
                       day_hours,
                       day_length,
                       #germination parameters
                       Tb.d=3,
                       To=22,
                       kT=0.12,
                       theta.HT=1000,
                       #indexes and definitions for germination matrices
                       dorm.breadth=1,
                       n.seed.classes=10,
                       psi.mean = 0,  #defines starting range psi
                       psi.min=-1,
                       psi.scalar=1,
                       #after-ripening parameters
                       Tsl=3,
                       psi.ll=-350,
                       psi.ul=-50,
                       thresh.wet=-5,
                       dsat=40,
                       n.cohort.start=1000, #starting pop size
                       #flowering parameters
                       thresh.PTU=2604,
                       Tb.f=3,ds=10,ps=0.626,dl=14,pl=1,WC.sat=960,Fi=0.598,Fu=0,
                       Tv.min=-3.5,Tv.max=6,
                       kappa=-5.17485563597551, omega=2.22568888596679, xi=.995904471970916,
                       #seed development paramenters
                       #simple.SDT=1,
                       thresh.SDU=352, 
                       seed.release.params=c(Tb.d=3),
                       #demographic parameters
                       seed.surv.params=c(1),
                       rosette.surv.params=c(1),
                       rosette.seed.params=c(1),
                       fecundity="simple",
                       #starting conditions 
                       start.PTU=0,
                       start.WC=0, 
                       start.SDU=0,
                       start.size=1.001,
                       dispersal.day=0,
                       #maximum storage of big matrices
                       max.cohort.store=1000,
                       #use only one out of two loops?
                       do.flower.only=FALSE,
                       do.germ.only=FALSE,
                       store.all=FALSE) {
    
    #vectors of storage for daily records
    nSeeds <- nRosette <- nFlower <- flowering <- germination <- seeding <-  rep(NA,Tmax/24)
 
    if (store.all) {
   	 	#matrix of storage for daily values in all cohorts
   	 	storePsi <- storeHTU <- array(dim=c(Tmax/24,n.seed.classes,max.cohort.store)) 
    	storeSize <- storePTU<-storePTU.growth <- storeWC <- storeSDU <-  matrix(NA,Tmax/24,max.cohort.store)
    }
    
    #index cohort from which start storage (will change as you approach max.cohort.stre)
	cohort.start.store1 <- cohort.start.store2 <- cohort.start.store3 <- 1
    #set up a matrix to store which seed classes make up each germination cohort germinate 
    storeSeedcohort<-matrix(0,n.seed.classes,max.cohort.store)
    store.n.cohort<-matrix(NA,n.seed.classes,max.cohort.store)
    #set up matrices of cohort trackers (columns = cohorts, rows = days)    
    cohort.seeds <- cohort.flower <- cohort.germ <- matrix(0,Tmax/24,max.cohort.store)
  
    
    count <- 0    # starts counting days at 0
    #values to be updated including matrix of psi, HTU, PTU, etc
    size <- n.germinants <- PTU <- WC<-PTU.growth <-dispersalday<- germday <- flowerday <- n.flower <- SDU <-seedday<-cohort.dormancy<- c() 
    n.cohort <- psi <- HTU <- matrix(0,n.seed.classes,1);

	
    	#initiate cohort of seeds
    	n.cohort[,1] <-store.n.cohort[,1] <-n.cohort.start*dnorm(seq(-3,3,length=n.seed.classes))/sum(dnorm(seq(-3,3,length=n.seed.classes)))   #initiate cohort of seeds
   	 	psi[,1] <-seq(from=psi.mean-(dorm.breadth/2),to=psi.mean+(dorm.breadth/2),length.out=n.seed.classes)
    	dispersalday<-count  

    #for flowering only, set up as if germinated today
    if (do.flower.only) {
        PTU <- start.PTU; WC<-start.WC; SDU <- start.SDU;  size <- start.size
        flowerday <- c(flowerday,0)
        seedday <- c(seedday,0)
        n.germinants <- n.cohort.start    # begin all the seeds as germinants on the same day
    } 
      
    #loop over times
    for (t in 1:Tmax) {
        ##advance HTU and ART according to the current env
        if(count<dispersal.day) {
        	if (t%%24==0){
                count <- count+1
                #print(count)
            }
        }
        
        else{
        if (!do.flower.only) { 
            tmpG <- nextGerm(HTU=HTU,
                             psi=psi,
                             Tb.d=Tb.d,
                             To=To,
                             kT=kT,
                             theta.HT=theta.HT,
                             Tsl=Tsl,
                             psi.ll=psi.ll,
                             psi.ul=psi.ul,
                             time.step=time.step,
                             temp=temp[t],
                             moisture=moisture[t],
                             is.wet=moisture[t]>thresh.wet,
                             n.cohort=n.cohort,
                             psi.mean=psi.mean,
                             psi.min=psi.min,
                             dsat=dsat,
                             psi.scalar=psi.scalar)
        
            HTU <- tmpG$HTU     
            psi <- tmpG$psi     
     
            ##assess germination; and store 
            if (t%%24==0){
                count <- count+1
              # print(count)
                germination[count] <- sum(n.cohort[which(HTU>theta.HT,arr.ind=T)])   # count up how many seeds germinated in the past 24 hours
                #print(n.cohort[which(HTU>theta.HT,arr.ind=T)])

                if (ncol(n.cohort)>1) {
                    cohort.germ[count,1:ncol(n.cohort)] <-  colSums(n.cohort*(HTU>theta.HT))    #cohort tracker
                } else {
                    cohort.germ[count,1:ncol(n.cohort)] <-  sum(n.cohort[which(HTU>theta.HT,arr.ind=T)]) 
                }
              
      		 if (store.all) {
                storeHTU[count,,c(1:length(tmpG$HTU[1,]))] <- ((HTU*n.cohort)/rowSums(n.cohort))[,cohort.start.store1:ncol(HTU)]
                                        # Record the HTU accumulation of each psi class of seeds. 
                storePsi[count,,c(1:length(tmpG$psi[1,]))] <- ((psi*n.cohort)/rowSums(n.cohort))[,cohort.start.store1:ncol(HTU)]              
               }
               
                storeSeedcohort[which(HTU>theta.HT,arr.ind=T)]<-count
                #Add a storage component for n.cohort no each day                
                n.cohort[which(HTU>theta.HT,arr.ind=T)] <- 0      # remove the seeds that germinated from the Pop being tracked
                HTU[which(HTU>theta.HT,arr.ind=T)] <- NA       # set the HTU trackers of those that germinated to NA
                psi[which(HTU>theta.HT,arr.ind=T)] <- NA                         #do psi as well?                     

                
                if (count>0) {
                    if (germination[count]>0) {
                        germday <- c(germday, count)  #store day germinated on 
                        flowerday <- c(flowerday,0)   # plan ahead for flowering 
                        seedday <- c(seedday,0)       # plan ahead for the day it produces seeds
                        cohort.dormancy<- c(cohort.dormancy,0)      
                        PTU <- c(PTU,start.PTU)       # set flowering PTU counter
                        WC <- c(WC,start.WC)             # set vernalization counter
                        size <- c(size,start.size)      # set the size counter (not currently in use)
                       PTU.growth<-c(PTU.growth,0)
                        SDU <- c(SDU,NA)              # add an entry to the Seed development counter in prep
                        n.germinants <- c(n.germinants,germination[count])   # add the number of new germinants to cohort tracker 
                        n.flower <- c(n.flower,0)     # Add a space to the flowering cohort tracker in preperation
                        
                        
                    }
                }

               
                	if (ncol(HTU)>=max.cohort.store) { 
                   	 	if (store.all) {
                   	 		storeHTU[,1:(max.cohort.store-1)] <- storeHTU[,2:max.cohort.store]
                    		storePsi[,1:(max.cohort.store-1)] <- storePsi[,2:max.cohort.store]
                    		storeHTU[,max.cohort.store] <- storePsi[,max.cohort.store] <- NA
                    	}
                    	cohort.start.store1 <- cohort.start.store1+1
                	}            	  
            }
        }
        
        ##advance PTU and WC according to the current env
        if (!do.germ.only) { 
            if (length(PTU)>0) {
                hour.in.day <- t%%24
                tmpF <- nextFlower(PTU=PTU,
                                   WC=WC, PTU.growth=PTU.growth,
                                   Tb.f=Tb.f,ds=ds,ps=ps,dl=dl,pl=pl,WC.sat=WC.sat,Fi=Fi,Fu=Fu,
                                   Tv.min=Tv.min,Tv.max=Tv.max,
                                   kappa=kappa, omega=omega, xi=xi,
                                   temp=temp[t],
                                   day_length=day_length[t],
                                   day_hours=day_hours[t])         
                PTU <- tmpF$PTU
                WC <- tmpF$WC
                PTU.growth[which(PTU>0)]<-PTU.growth[which(PTU>0)]+tmpF$PTU.growth
            }
            
            ##assess  flowering
            ## at daily intervals
            if (t%%24==0 & length(PTU)>0){      
                if (do.flower.only) count <- count+1 
                flowering[count] <- sum(n.germinants[which(PTU>thresh.PTU,arr.ind=T)], na.rm=T)
                                        # record the # of individuals that flowered during the 24 hours
                #cohort.flower[count,which(PTU>thresh.PTU,arr.ind=T)] <- n.germinants[which(PTU>thresh.PTU,arr.ind=T)] # Does not function right across a broader range of contexts than those used in the paper. Not used in summaries so currently just hashed out.
               
               	if (store.all) {
                	storePTU[count,c(1:length(PTU))] <- PTU[cohort.start.store2:length(PTU)]  # store the PTUs for each cohort being tracked
                	storeWC[count,c(1:length(WC))] <- WC[cohort.start.store2:length(WC)]   # store the vern accumulation for each cohort being tracked
                	storeSize[count,c(1:length(size))] <- size[cohort.start.store2:length(size)] # store the sizes of each of the plants
                	storePTU.growth[count,c(1:length(PTU.growth))] <- PTU.growth[cohort.start.store2:length(PTU.growth)] # store the sizes of each of the plants
                }                                                
               
                if (flowering[count]>0) {
                    flowerday[which(PTU>thresh.PTU,arr.ind=T)] <- count  # record the day that the cohorts flowered
                    n.flower[which(PTU>thresh.PTU,arr.ind=T)] <- n.germinants[which(PTU>thresh.PTU,arr.ind=T)] # transfer the number of individuals in the n.germ cohorts that flowered into their respective n.flower cohort
                    n.germinants[which(PTU>thresh.PTU,arr.ind=T)] <- NA   #get rid of the ones that flowered
                    SDU[which(PTU>thresh.PTU,arr.ind=T)] <- start.SDU  # start the seed development accumulation                
                    WC[PTU>thresh.PTU] <- NA    # Set WC accumulation of those that flowered to NA
                    PTU[PTU>thresh.PTU] <- NA   #Set PTUs of those that flowered to NA
                }


                if (length(PTU)>=max.cohort.store) { 
                    if (store.all) {
                    	storePTU[1:(max.cohort.store-1)] <- storePTU[2:max.cohort.store]
                    	storeWC[1:(max.cohort.store-1)] <- storeWC[2:max.cohort.store]
                    	storeSize[1:(max.cohort.store-1)] <- storeSize[2:max.cohort.store]
                    	storePTU[max.cohort.store] <- storeWC[max.cohort.store] <- storeSize[max.cohort.store] <- NA
                    }
                    cohort.start.store2 <- cohort.start.store2+1
                }
              
            }
        
            ## advance seed development according to current env
            if (length(SDU)>0) {
                tmpSDU<- nextSeed(simple.SDT=simple.SDT,
                                  SDU=SDU,seed.release.params=seed.release.params,
                                  t=t,temp=temp[t],moisture=moisture[t])
                SDU<-tmpSDU
                
            }
            ## assess seed maturation at daily intervals
            if (t%%24==0 & length(SDU)>0){
                #print(SDU)
                seeding[count]<-sum(n.flower[which(SDU>thresh.SDU,arr.ind=T)],na.rm=T)  # store how many plants dispersed seeds that day
               if (store.all) {
                	storeSDU[count,c(1:length(SDU))] <- SDU[cohort.start.store3:length(SDU)] # store the Seed Development time accumulated in last 24 hrs
               } 
               
                #FECUNDITY--determine how many seeds each of those plants will disperse
                if (fecundity=="simple") { 
                		new.seeds<-sum(n.flower[which(SDU>thresh.SDU,arr.ind=T)]* n.seeds.simple(rosette.seed.params=rosette.seed.params,
                		size=size[which(SDU>thresh.SDU,arr.ind=T)],
						PTU.growth=PTU.growth[which(SDU>thresh.SDU, arr.ind=T)],
                        age=flowerday[which(SDU>thresh.SDU,arr.ind=T)]-germday[which(SDU>thresh.SDU,arr.ind=T)]))   
                }		                                            
                
                if (seeding[count] > 0){                                     
                    #initiation of new cohort of seeds 
                    seedday[which(SDU>thresh.SDU,arr.ind=T)] <- count   # store on which day each germ cohort dispersed seeds
                    dispersalday<-c(dispersalday,count) #store what day the new seed cohort (can come from multiple germination days if seeding becomes synchronized.)
                    cohort.seeds[count,which(SDU>thresh.SDU,arr.ind=T)] <- n.germinants[which(SDU>thresh.SDU,arr.ind=T)] 
                    n.cohort <- cbind(n.cohort,rep(0,length(n.seed.classes)))   # add a new column for this cohort to the n.cohort
                    n.cohort[,ncol(n.cohort)] <-store.n.cohort[,ncol(n.cohort)]<- new.seeds*dnorm(seq(-3,3,length=n.seed.classes))/sum(dnorm(seq(-3,3,length=n.seed.classes))) # will need to modify this if want store.n.cohort to change with a max # of cohorts. Not sure this is necessary
                   
                   #define the dormancy of the seed cohort ********This is not yet structured so different mothers that have different histories that disperse seeds on the same day can have different dormancies... this will take a restructuring of the code. Right now combines all the plants that flower on the same day into a dispersal cohort.

                    # distribute new seeds in a normal distribution across the new cohort column                                                      
                    HTU <- cbind(HTU,rep(0,n.seed.classes))  # same as above but for adding a preperatory column to the HTU matrix
                    psi <- cbind(psi,rep(0,n.seed.classes)) # same as above but for adding a preperatory column to the psi matrix
                    psi[,ncol(n.cohort)] <-seq(from=psi.mean-(dorm.breadth/2),to=psi.mean+(dorm.breadth/2),length.out=n.seed.classes)  
                 	# determine the psi's that go along with each seed fraction with a given mean starting psi and step between psis
                    n.flower[SDU>thresh.SDU] <- NA  # get rid of the individuals that produced seeds
                    SDU[SDU>thresh.SDU] <- NA   # clear the Seed maturation tracking cell for those individuals
                           
                }
                
                if (length(SDU)>=max.cohort.store) { 
                    if (store.all) {
                    	storeSDU[1:(max.cohort.store-1)] <- storeSDU[2:max.cohort.store]
                    }
                    cohort.start.store3 <- cohort.start.store3+1
                }

            }  
          }      
        if (t%%24==0){
        #Assess survival of all life stages 
        # seeds
        surv <- seed.survival(seed.surv.params=seed.surv.params,
                                          temp=temp[t],moisture=moisture[t])
        #rosettes and flowering plants       
        surv.rosette <- rosette.survival(rosette.surv.params=rosette.surv.params,
                                                 temp=temp[t],moisture=moisture[t], ptu=PTU, size=size)                
        
        #Reduce each of the cohort tracking devices accordingly
        n.cohort <- n.cohort*surv       # death!! of seeds
        n.germinants <- n.germinants*surv.rosette # reduce rosettes according to rossette survival for the day                
        n.flower <- n.flower*surv.rosette # reduce flowering plants (right now is the same as rosette survival)
          
        nSeeds[count] <- sum(n.cohort)          # record number of seeds that day
        nRosette[count] <- sum(n.germinants, na.rm=T)  # record number of rosettes that day
        nFlower[count] <- sum(n.flower,na.rm=T)        # record the no of individuals currently flowering
      } 
     }
     print(count)                        
    }
    
    if(store.all){ 
      	return(list(germination=germination, flowering=flowering, seeding=seeding,
                cohort.flower=cohort.flower,cohort.seeds=cohort.seeds,cohort.germ=cohort.germ,
                dispersalday=dispersalday,germday=germday,flowerday=flowerday, seedday=seedday, cohort.dormancy=cohort.dormancy,             
                
                storePsi=storePsi, storeHTU=storeHTU,      
                storePTU=storePTU, storeWC=storeWC, storeSize=storeSize,
                storeSDU=storeSDU,storePTU.growth=storePTU.growth,
                
                 cohort.start.store1= cohort.start.store1, cohort.start.store2= cohort.start.store2, cohort.start.store3= cohort.start.store3,
                
                storeSeedcohort=storeSeedcohort,store.n.cohort=store.n.cohort,                                                                  
                nSeeds=nSeeds, nRosette=nRosette, nFlower=nFlower,
                n.germinants=n.germinants, n.flower=n.flower,n.cohort=n.cohort,
                psi=psi,HTU=HTU,PTU=PTU,WC=WC,SDU=SDU,size=size,PTU.growth=PTU.growth))
	}  
    
    else{
    return(list(germination=germination, flowering=flowering, seeding=seeding,
                cohort.flower=cohort.flower,cohort.seeds=cohort.seeds,cohort.germ=cohort.germ,
                dispersalday=dispersalday,germday=germday,flowerday=flowerday, seedday=seedday, cohort.dormancy=cohort.dormancy,                             
                
                 cohort.start.store1= cohort.start.store1, cohort.start.store2= cohort.start.store2, cohort.start.store3= cohort.start.store3,
                
                storeSeedcohort=storeSeedcohort,store.n.cohort=store.n.cohort,                                                                  
                nSeeds=nSeeds, nRosette=nRosette, nFlower=nFlower,
                n.germinants=n.germinants, n.flower=n.flower,n.cohort=n.cohort))
	}  
}

#########Explanation of the output from iteratePop#######################
#germination, flowering, and seeding each are a vector on the day scale with as many 
#spaces as days in the model run (or Tmax/24) length. This vector contains the 
# number of individuals which underwent that stage transition on a given day

# germday, flowerday, and seedday are vectors of equal length that store the day 
#that life history transitions of a given cohort of germinants occured. The vectors 
#correspond to one another so the first space in each is the first cohort. 
#A seed cohort is all those seeds that were produced on the same day. This is broken up because 
# seeds that germinate on the same day are pooled together into flowering cohorts. 

#n.germinants, n.flower, n.cohort all simply hold the number of individuals in a 
#life stage at one time. Once the individuals have transfered out is set to NA. 
#These are of variable length and do not match up cohort to cohort.

#storeSeedcohort is a matrix of size the number of seed classes by the potential number of cohorts that stores the day that seeds from each class and cohort germinate
#store.n.cohort is a matrix the same dimensions as storeSeedcohort which stores the number of seeds in each seed class and cohort...

#nSeeds,nRosettes,nFlower records the number of individuals that transitioned on a given day over the simulation
#n.germinants, n.flower, n.cohort keeps track of the number of individuals in each catogory at a particular time

# storage objects only included if store.all =TRUE
#storePTU, storeWC, storeSize, store SDU are all matrices 365 rows and 100 colums 
#which acts to record daily progress of each cohort toward each transition (so by default we
#store up to 1000 cohorts)

#storeHTU similar to storeARU etc... but with a rows that correspond to cohort fraction
#storePsi identical to storeHTU

#psi, HTU, PTU etc... give the last image of each of these tracking structures.
# size is the size at flowering and PTU growth is based on PTU accumulated 
#####################################################################

######## FUNCTION THAT SUMMARIZES THE OUTPUT FROM ITERATEPOP#######

whatistheLifecycle.output<-function(rc, first.year=15,last.year=40) {

first.day<-366+365*(first.year-1)
last.day<-last.year*365

# Create a matrix of the dispersal day for each seed cohort
seed.cohort.classes<-dim(rc$storeSeedcohort)[1]
disp.mat<-matrix(rep(rc$dispersalday,seed.cohort.classes),seed.cohort.classes,length(rc$dispersalday),byrow=TRUE)
Seedcohort<-rc$storeSeedcohort[,1:dim(disp.mat)[2]] ### Cut off Seedcohort structure based on the dispersal matrix
n.seedcohort<-rc$store.n.cohort[,1:dim(disp.mat)[2]] ### Cut off Cohort number storage structure based on the dispersal matrix

#For seed period you can directly calculate the length of time from dispersal to germination
seedperiod<-Seedcohort-disp.mat # This is a matrix of dimensions the # of seed classes and the # of cohorts with the time spent as a seed for each class of seed... 

# Because germinants are subsequently grouped by germination cohort. We need to backcalculate the vegetative and flowering lengths
abovegroundperiod<-rc$seedday[rc$seedday>0]-rc$germday[rc$seedday>0] #gives you the aboveground length for each GERMINATION cohort
vegperiod<-rc$flowerday[rc$seedday>0]-rc$germday[rc$seedday>0]
flowperiod<-rc$seedday[rc$seedday>0]-rc$flowerday[rc$seedday>0]
seeding<-rc$seedday[rc$seedday>0]

# Determine the germination day for each germinaton
germ.temp <-rc$germday[rc$seedday>0]#take the subset of germdays for which a vegetative period could be calculated

# We match each germ cohort length to seed classes via germination day to determine the generation times within in each cohort of seeds
nonseedperiod<-Seedcohort #create a structure to hold the nonseed, veggie, and flowering periods for each germ cohort 
veggieperiod<-Seedcohort
floweringperiod<-Seedcohort
seedingtime<-Seedcohort  # this allows you to cut off cohorts based on whether they completed their lifecycle with-in your timeframe

#replace the germination day that each seed group germinated with the vegetative length associated with that germination cohort
	for(s in 1:length(germ.temp)){
		nonseedperiod[which(nonseedperiod==germ.temp[s])]<-abovegroundperiod[s]
	}
nonseedperiod[which(nonseedperiod>max(germ.temp))]<-NA  # This is the matrix for the vegitative length
nonseedperiod[which(nonseedperiod==0)]<-NA
	
	for(s in 1:length(germ.temp)){
		veggieperiod[which(veggieperiod==germ.temp[s])]<-vegperiod[s]
	}
veggieperiod[which(veggieperiod>max(germ.temp))]<-NA  # This is the matrix for the vegitative length
veggieperiod[which(veggieperiod==0)]<-NA 

	for(s in 1:length(germ.temp)){
		floweringperiod[which(floweringperiod==germ.temp[s])]<-flowperiod[s]
	}
floweringperiod[which(floweringperiod>max(germ.temp))]<-NA  # This is the matrix for the vegitative length
floweringperiod[which(floweringperiod==0)]<-NA 

	for(s in 1:length(germ.temp)){
		seedingtime[which(seedingtime==germ.temp[s])]<-seeding[s]
	}
seedingtime[which(seedingtime==0)]<-NA 

totperiod<-nonseedperiod+seedperiod

#Cut off the end part for cohorts that did not complete their life cycle
totperiod<-totperiod[,which(!is.na(apply(floweringperiod,2,FUN=mean)))]
nonseedperiod<-nonseedperiod[,which(!is.na(apply(floweringperiod,2,FUN=mean)))]
seedperiod<-seedperiod[,which(!is.na(apply(floweringperiod,2,FUN=mean)))]
veggieperiod<-veggieperiod[,which(!is.na(apply(floweringperiod,2,FUN=mean)))]
floweringperiod<-floweringperiod[,which(!is.na(apply(floweringperiod,2,FUN=mean)))]
seedingtime<-seedingtime[,which(!is.na(apply(floweringperiod,2,FUN=mean)))]
n.seedcohort<-n.seedcohort[,which(!is.na(apply(floweringperiod,2,FUN=mean)))]
disp.mat<-disp.mat[,which(!is.na(apply(floweringperiod,2,FUN=mean)))]
dispersalday<-rc$dispersalday[which(!is.na(apply(floweringperiod,2,FUN=mean)))]%%365


#Cut off these results based on the time period you are summarizing over based on disp.mat and dispersaltime
cohort.seedingtimes<-apply(seedingtime,2,FUN=max)
totperiod<-totperiod[,intersect(which(rc$dispersalday>first.day),which(cohort.seedingtimes<last.day))]
nonseedperiod<-nonseedperiod[,intersect(which(rc$dispersalday>first.day),which(cohort.seedingtimes<last.day))]
seedperiod<-seedperiod[,intersect(which(rc$dispersalday>first.day),which(cohort.seedingtimes<last.day))]
veggieperiod<-veggieperiod[,intersect(which(rc$dispersalday>first.day),which(cohort.seedingtimes<last.day))]
floweringperiod<-floweringperiod[,intersect(which(rc$dispersalday>first.day),which(cohort.seedingtimes<last.day))]
seedingtime<-seedingtime[,intersect(which(rc$dispersalday>first.day),which(cohort.seedingtimes<last.day))]
n.seedcohort<-n.seedcohort[,intersect(which(rc$dispersalday>first.day),which(cohort.seedingtimes<last.day))]
disp.mat<-disp.mat[,intersect(which(rc$dispersalday>first.day),which(cohort.seedingtimes<last.day))]
dispersalday<-dispersalday[intersect(which(rc$dispersalday>first.day),which(cohort.seedingtimes<last.day))]

#store which month each dispersal event happened in....
months<- findInterval(dispersalday,c(seq(0,365,30.42)))

#Calculate the length averages for each cohort for each stage ***regardless of the # of individuals in each cohort ****
if(length(totperiod)>10){
totperiod.cohort<-colSums(n.seedcohort*totperiod)/colSums(n.seedcohort)
seedperiod.cohort<-colSums(n.seedcohort*seedperiod)/colSums(n.seedcohort)
veggieperiod.cohort<-colSums(n.seedcohort*veggieperiod)/colSums(n.seedcohort)
floweringperiod.cohort<-colSums(n.seedcohort*floweringperiod)/colSums(n.seedcohort)
nonseedperiod.cohort<-colSums(n.seedcohort*nonseedperiod)/colSums(n.seedcohort)
cohort.size<-colSums(n.seedcohort)
}

if(length(totperiod)<=10){
totperiod.cohort<-sum(n.seedcohort*totperiod)/sum(n.seedcohort)
seedperiod.cohort<-sum(n.seedcohort*seedperiod)/sum(n.seedcohort)
veggieperiod.cohort<-sum(n.seedcohort*veggieperiod)/sum(n.seedcohort)
floweringperiod.cohort<-sum(n.seedcohort*floweringperiod)/sum(n.seedcohort)
nonseedperiod.cohort<-sum(n.seedcohort*nonseedperiod)/sum(n.seedcohort)
cohort.size<-sum(n.seedcohort)	
}

avg.gen<-sum(n.seedcohort*totperiod)/sum(n.seedcohort)
avg.seed<-sum(n.seedcohort*seedperiod)/sum(n.seedcohort)
avg.veggie<-sum(n.seedcohort*veggieperiod)/sum(n.seedcohort)
avg.flowering<-sum(n.seedcohort*floweringperiod)/sum(n.seedcohort)
avg.nonseed<-sum(n.seedcohort*nonseedperiod)/sum(n.seedcohort)

# Reformat period data so there is a data point for every individual indicating its stage length
cohorts<-as.vector(round(n.seedcohort)) #vectorize and round off the number of seeds in each cohort 
floweringperiod.all<-rep(as.vector(floweringperiod),cohorts) #repeat each life length by the # of individuals in that cohort
veggieperiod.all<-rep(as.vector(veggieperiod),cohorts)
seedperiod.all<-rep(as.vector(seedperiod),cohorts)
nonseedperiod.all<-rep(as.vector(nonseedperiod),cohorts)
totperiod.all<-rep(as.vector(totperiod),cohorts)
dispersal.all<-rep(as.vector(disp.mat),cohorts)

 # A data frame holding all of the total life stage lengths  
  	period.all<-data.frame(rbind(cbind(seedperiod.all,rep(1,length(seedperiod.all)),rep("Seed",length(seedperiod.all))),
  		cbind(veggieperiod.all,rep(2,length(veggieperiod.all)),rep("Vegetative",length(seedperiod.all))),
  		cbind(floweringperiod.all,rep(3,length(floweringperiod.all)),rep("Reproductive",length(seedperiod.all))),
  		cbind(totperiod.all,rep(4,length(totperiod.all)),rep("Total",length(seedperiod.all)))))
	colnames(period.all)<-c("Length","Stage.numeric","Stage.name")
	period.all$Length<-as.numeric(as.character(period.all$Length))
	period.all$Stage.numeric<-as.numeric(as.character(period.all$Stage.numeric))
	period.all$Stage.name<-(factor(period.all$Stage.name,levels=c("Seed","Vegetative","Reproductive","Total")))

###### Summarize the cohort results by month for spaghetti plots #########
if(length(totperiod)<=10){
	months<-c(0,months)
	totperiod<-matrix(c(rep(0,20-length(totperiod)),totperiod),nrow=10,ncol=2)
	seedperiod<-matrix(c(rep(0,20-length(seedperiod)),totperiod),nrow=10,ncol=2)
	veggieperiod<-matrix(c(rep(0,20-length(veggieperiod)),totperiod),nrow=10,ncol=2)
	floweringperiod<-matrix(c(rep(0,20-length(floweringperiod)),totperiod),nrow=10,ncol=2)
	nonseedperiod<-matrix(c(rep(0,20-length(nonseedperiod)),totperiod),nrow=10,ncol=2)
	cohort.size<-c(0,cohort.size)
}

sumbymonth.tot<-numeric(0)
for(z in 1:seed.cohort.classes){
c.z<-tapply(totperiod[z,-1],months[-1],mean.NA)
c.z[which(c.z=="NaN")]<-NA
sumbymonth.tot<-rbind(sumbymonth.tot,c.z)
}
sumbymonth.seed<-numeric(0)
for(z in 1:seed.cohort.classes){
c.z<-tapply(seedperiod[z,-1],months[-1],mean.NA)
c.z[which(c.z=="NaN")]<-NA
sumbymonth.seed<-rbind(sumbymonth.seed,c.z)
}
sumbymonth.veggie<-numeric(0)
for(z in 1:seed.cohort.classes){
c.z<-tapply(veggieperiod[z,-1],months[-1],mean.NA)
c.z[which(c.z=="NaN")]<-NA
sumbymonth.veggie<-rbind(sumbymonth.veggie,c.z)
}
sumbymonth.flowering<-numeric(0)
for(z in 1:seed.cohort.classes){
c.z<-tapply(floweringperiod[z,-1],months[-1],mean.NA)
c.z[which(c.z=="NaN")]<-NA
sumbymonth.flowering<-rbind(sumbymonth.flowering,c.z)
}
sumbymonth.nonseed<-numeric(0)
for(z in 1:seed.cohort.classes){
c.z<-tapply(nonseedperiod[z,-1],months[-1],mean.NA)
c.z[which(c.z=="NaN")]<-NA
sumbymonth.nonseed<-rbind(sumbymonth.nonseed,c.z)
}
rownames(sumbymonth.tot)<-rownames(sumbymonth.veggie)<-rownames(sumbymonth.seed)<-rownames(sumbymonth.flowering)<-rownames(sumbymonth.nonseed)<-c(1:seed.cohort.classes)

# Create the summary of how many seed there are in each month
month.cohortsize.temp<-tapply(cohort.size[-1],months[-1],sum)
month.cohortsize<-rep(0,12)
names(month.cohortsize)<-1:12
month.cohortsize[as.numeric(names(month.cohortsize.temp))]<-month.cohortsize.temp
month.cohortsize.log<-findInterval(month.cohortsize,c(10,100,1000,10000,100000))

######## Summarize the timing of life stage transitions by day,week, and month   ######
   week<-c(1:52)
   weeks<- findInterval(c(1:365),c(seq(1,365-7,7)))
   month<- findInterval(c(1:365),c(seq(0,365,30.42)))
   
   germ.temp<-data.frame(cbind(rc$germination[first.day:last.day],rep(1:365,length(first.day:last.day)/365)))
   germ.days<-tapply(germ.temp[,1],germ.temp[,2],mean)        
   germ.weeks<-sapply(split(germ.days,weeks),sum,na.rm=T)
   germ.months<-sapply(split(germ.days,month),sum,na.rm=T)
   
   flow.temp<-data.frame(cbind(rc$flowering[first.day:last.day],rep(1:365,length(first.day:last.day)/365)))
   flow.days<-tapply(flow.temp[,1],flow.temp[,2],mean)        
	flow.weeks<-sapply(split(flow.days,weeks),sum,na.rm=T)
   	flow.months<-sapply(split(flow.days,month),sum,na.rm=T) 
   	
   	 seed.temp<-data.frame(cbind(rc$seeding[first.day:last.day],rep(1:365,length(first.day:last.day)/365)))
     seed.days<-tapply( seed.temp[,1], seed.temp[,2],mean)  
	seed.weeks<-sapply(split(seed.days,weeks),sum,na.rm=T)
   	seed.months<-sapply(split(seed.days,month),sum,na.rm=T) 
  
  germ.all=rc$germination[first.day:last.day]
  flow.all=rc$germination[first.day:last.day]
  seed.all=rc$germination[first.day:last.day]
  
  
  #### How many of each type of plant exists and any given time#######
 
   avg.year<-(last.day-first.day)/365
    times <- 1+c(first.day:last.day)%%365
    
    nSeeds.all <- rc$nSeeds[first.day:last.day]
    nRosette.all <- rc$nRosette[first.day:last.day]
    nFlower.all <- rc$nFlower[first.day:last.day]
    
    nSeeds.daily <- as.vector(sapply(split(nSeeds.all[1:(365* avg.year)],times[1:(365* avg.year)]),mean)) #get total numbers
    nSeeds.weekly<-sapply(split(nSeeds.daily,weeks),sum,na.rm=T)
    nSeeds.monthly<-sapply(split(nSeeds.daily,month),sum,na.rm=T)
    
    nRosette.daily <- as.vector(sapply(split(nRosette.all[1:(365* avg.year)],times[1:(365* avg.year)]),mean))
    nRosette.weekly<-sapply(split(nRosette.daily,weeks),sum,na.rm=T)
    nRosette.monthly<-sapply(split(nRosette.daily,month),sum,na.rm=T)
    
    nFlower.daily<- as.vector(sapply(split(nFlower.all[1:(365* avg.year)],times[1:(365* avg.year)]),mean))
    nFlower.weekly<-sapply(split(nFlower.daily,weeks),sum,na.rm=T)
    nFlower.monthly<-sapply(split(nFlower.daily,month),sum,na.rm=T)
    
   states.all <- rbind(nSeeds.daily,nRosette.daily,nFlower.daily) # Why do I have this?
   
   	seed.hist<-hist.line(seedperiod.all)
	veggie.hist<-hist.line(veggieperiod.all)
	flowering.hist<-hist.line(floweringperiod.all)
	tot.hist<-hist.line(totperiod.all)
	nonseed.hist<-hist.line(nonseedperiod.all)
	
	tot.hist.mid<-tot.hist$mid
	seed.hist.mid<-seed.hist$mid
	 veggie.hist.mid<-veggie.hist$mid
	 flowering.hist.mid<-flowering.hist$mid
	 nonseed.hist.mid<-nonseed.hist$mid	
	  
	  tot.hist.density<-tot.hist$density
	  seed.hist.density<-seed.hist$density
	   veggie.hist.density<-veggie.hist$density
	   flowering.hist.density<-flowering.hist$density
	   nonseed.hist.density<-nonseed.hist$density
                
return(list(

disp.mat=disp.mat,dispersalday=dispersalday, n.seedcohort=n.seedcohort,cohort.size=cohort.size, totperiod=totperiod, seedperiod=seedperiod, veggieperiod=veggieperiod, floweringperiod=floweringperiod, nonseedperiod=nonseedperiod, seed.cohort.classes=seed.cohort.classes,	months=months,
		
avg.gen=avg.gen,avg.seed=avg.seed,avg.flowering=avg.flowering,avg.veggie=avg.veggie,avg.nonseed=avg.nonseed,

totperiod.cohort=totperiod.cohort, seedperiod.cohort=seedperiod.cohort, nonseedperiod.cohort=nonseedperiod.cohort, veggieperiod.cohort=veggieperiod.cohort, floweringperiod.cohort=floweringperiod.cohort,
		
seedperiod.all=seedperiod.all, veggieperiod.all=veggieperiod.all, floweringperiod.all=floweringperiod.all, nonseedperiod.all=nonseedperiod.all, totperiod.all=totperiod.all, dispersal.all=dispersal.all, period.all=period.all,

month.cohortsize=month.cohortsize, month.cohortsize.log=month.cohortsize.log, 

sumbymonth.tot=sumbymonth.tot,sumbymonth.seed=sumbymonth.seed, sumbymonth.veggie=sumbymonth.veggie, sumbymonth.flowering=sumbymonth.flowering, sumbymonth.nonseed=sumbymonth.nonseed,
		
germ.days=germ.days ,flow.days=flow.days, seed.days=seed.days, germ.weeks=germ.weeks, flow.weeks=flow.weeks, seed.weeks=seed.weeks, germ.months=germ.months, flow.months=flow.months, seed.months=seed.months, germ.all=germ.all, flow.all=flow.all, flow.all=flow.all,
		
nSeeds.all=nSeeds.all, nRosette.all=nRosette.all, nFlower.all=nFlower.all,
nSeeds.daily=nSeeds.daily, nRosette.daily=nRosette.daily, nFlower.daily=nFlower.daily, nSeeds.weekly=nSeeds.weekly, nRosette.weekly=nRosette.weekly, nFlower.weekly=nFlower.weekly,
nSeeds.monthly=nSeeds.monthly, nRosette.monthly=nRosette.monthly, nFlower.monthly=nFlower.monthly,		
seed.hist=seed.hist,veggie.hist=veggie.hist,flowering.hist=flowering.hist,nonseed.hist=nonseed.hist,tot.hist=tot.hist, tot.hist.mid=tot.hist.mid,seed.hist.mid=seed.hist.mid, veggie.hist.mid=veggie.hist.mid, flowering.hist.mid=flowering.hist.mid,nonsed.hist.mid=nonseed.hist.mid, tot.hist.density=tot.hist.density,seed.hist.density=seed.hist.density, veggie.hist.density=veggie.hist.density, flowering.hist.density=flowering.hist.density,nonsed.hist.density=nonseed.hist.density
		))
}

################### Explanation of the output of whatistheLifecycle.output ###############
### all of these data frames are restricted to show only data from year 15 of the simulation to year 60 of the simulation.
#disp.mat - matrix of the dispersal day for each seed class (rows) x seed cohort (columns)
#dispersalday -vector the length of the number of cohorts giving the dispersal day for each cohort
#totperiod - matrix of the generation length of for each seed class (rows) x seed cohort (columns)
#seedperiod - matrix of the length of time spent as a seed for each seed class (rows) x seed cohort (columns) 
#veggieperiod - matrix of the length of time spent as a rosette for each seed class (rows) x seed cohort (columns) 
#floweringperiod - matrix of the length of time spent flowering for each seed class (rows) x seed cohort (columns) 
#nonseedperiod -matrix of the length of time spent as a rosette or flowering for each seed class (rows) x seed cohort (columns) 
# cohort.size - vector of length the # of seed cohorts indicating the # of individuals in each cohort.
#seed.cohort.classes - # of seed classes used in the run	
#n.seedcohort - matrix of # of individuals in each seed class (rows) x seed cohort (columns)	
#months - indication of which month each dispersal event happened in 

####### Specialized data frames ###############
##### X stands in for various life stages or transitions. tot is short for total and refers to a value measure across the entire life cycle of a plant such as total life-cycle length

#avg.X - average (generation, seed stage, vegetative stage, or flowering stage) length for all individuals across the whole simulation
# period.all - data frame holding all the life-stage lengths for each individual
#Xperiod.all- vector holding the life-stage lengths for each individual in the simulation
#Xperiod.cohort - same as above but averaged across cohorts so there is an average life stage length for each cohort not for each individual.
#month.cohortsize.X - vector of cohort sizes averaged per month
#sumbymonth.X - vector of average stage length (or generation length) for each seed class summarized for each month of the year

#Output that has to do with when particular life stage transitions occur
#X.all - number of germination events, flowering events, or dispersal events on every day of the simulation (vector the length of the # of days of your simulation)
#X.days - average # of germination, flowering, or dispersal events on each day of the year across the entire simulation
#X.weeks - average # of germination, flowering, or dispersal events on each week of the year across the entire simulation
#X.months - average # of germination, flowering, or dispersal events on each month of the year across the entire simulation

# Output that has to do with the # of plants in a given life stage a given time
#nX.all - number of seeds, rosettes, or flowering plants on every day of the simulation (vector the length of the # of days of your simulation) 
#nX.daily- number of seeds, rosettes, or flowering plants one each day of the year (vector of 365 days of the year) 
#.nX.weekly - number of seeds, rosettes, or flowering plants one each day of the year (vector of 52 weeks of the year) 
#.nX.monthly - number of seeds, rosettes, or flowering plants one each day of the year (vector of 12 months of the year) 

#Outputs that are designed for making histograms of stage lengths
#X.hist - all of the summarized output from the hist function in R for the given life stage length
#.hist.density - density output from the hist function in R for the given life stage length (paired with the mid values below to create a histogram)
#X.hist.mid - mid point values to pair with the density values above from the hist function in R for the given life stage length 