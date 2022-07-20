#title: "Quantitative genetic model to examine the effects of assortative mating on population recruitment and viability"
#authors: "Samuel May, Eric Ward, Jeffery Hard, and Kerry Naish"
#date: "07/05/2022"


##### Introduction #####  
# Here we provide code for an individual based model to simulate wild population 
# dynamics within a quantitative genetic framework for inheritance of 
# phenotypically correlated traits. This model can be somewhat computationally 
# intensive if population sizes get very large. We recommend running the 
# model using parallel computing - to assist with this we have coded the model 
# as a function Assortative_Mating_IBM() compatible with the foreach and 
# doParallel packages for parallel computing in R. We implemented this function
# in a High Performance Computing (HPC) environment to permit many iterations 
# and combinations of input parameters. The function generates one .csv file 
# per iteration from the internal object "gen_data", which is a pedigree (animal, dam, sire)
# with additional phenotypic data for each individual.

# Beneath the function is the code used to run the simulations detailed 
# in the text.

#### Load Data and Packages ####

library(dplyr)
library(TruncatedNormal)
library(rockchalk)
library(foreach)
library(doParallel)

#### Description of Function Arguments: ####

#iteration = number of iterations to run
#rho = phenotypic correlation between entry day and RLS
#variance_entry_day = phenotypic variance in entry day
#variance_RLS = phenotypic variance in RLS
#assortative = should mating occur randomly or assortatively? (T or F) 
#PopSize = allow the population size to vary or remain constant at Nc_initial? (variable or constant)... 
#G = number of generations to run for
#Nc_initial = number of individuals in the starting population
#OMEGA = shape of the fitness surface - approx 1-4 phenotypic standard deviations
#var_env_entry = Interannual variability in optimal entry day
#var_env_RLS = Interannual variability in optimal lifespan
#PATH = the directory to which output files should be written

#### Model Initialization ####

Assortative_Mating_IBM<-function(iteration = 1,
                                 rho = -0.3,
                                 MU_entry_day = 10,
                                 MU_RLS = 5,
                                 variance_entry_day = 20,
                                 variance_RLS = 20,
                                 assortative = T,
                                 PopSize = "variable",
                                 G = 1,
                                 OMEGA = 2,
                                 var_env_entry = 20,
                                 var_env_RLS = 1,
                                 Nc_initial = 500,
                                 PATH = "output/gen_data_"){

  
#Make Sigma (The phenotpyic VCV matrix, or "P matrix"):
#Sigma = Phenotypic Matrix = Additive Genetic VCV matrix + Residual (environmental) VCV matrix
  
Sigma<-rockchalk::lazyCov(Rho = rho, Sd = c(sqrt(variance_RLS),sqrt(variance_entry_day)))
rownames(Sigma)<-c("RLS","entry_day")
colnames(Sigma)<-c("RLS","entry_day")

Corr<-cov2cor(Sigma)

#Initiate the starting population (F0) from Sigma and population trait means. 
F0<-data.frame(id=1:Nc_initial,entry_day=NA,RLS=NA,sex=rep(c("m","f")),year=2022) #year is irrelevant and redundant with Generation... but could be used if desired.
u = c(MU_RLS,MU_entry_day) 

#Draw entry days for Nc_initial individuals
F0$entry_day<-TruncatedNormal::rtnorm(n=Nc_initial,mu=u[2],sd=sqrt(Sigma[2,2]),lb=0,ub=30) #30 day spawning season

#Draw RLS values for those individuals, conditional upon entry day such that the traits will be correlated
for(i in 1:nrow(F0)){
  cond_mean = u[1] + Corr[1,2]*sqrt(Sigma[1,1])/sqrt(Sigma[2,2]) * (F0[i,"entry_day"] - u[2])
  cond_sd = sqrt(Sigma[1,1]) * sqrt(1 - Corr[1,2])
  
  # RLS cannot be greater than the length of the spawning season (30) - entry day
  F0$RLS[i]<-TruncatedNormal::rtnorm(1,mu=cond_mean,sd=cond_sd,lb=0,ub=30-F0[i,"entry_day"]) 
}

#Round these trait values to the nearest day
F0$entry_day<-round(F0$entry_day)
F0$RLS<-round(F0$RLS)


#Initialize the data to keep track of for each individual within each generation (gen_data)
#The gen_data object will be be the eventual output of the model
#Add F0 information to this data frame
gen_data <- matrix(nrow=Nc_initial,ncol=9)
colnames(gen_data) <- c("ID", "G", "sire", "dam", "sex", "entry_day", "RLS","RS_exp","RS_obs") 
gen_data <- as_tibble(gen_data)
gen_data$ID<-seq(1:nrow(F0))
gen_data$G<-0
gen_data$sex<-F0$sex
gen_data$entry_day<-F0$entry_day
gen_data$RLS<-F0$RLS
gen_data$year<-F0$year


for (g in 0:G){ #Begin looping through generations
  
  #Set parents of current generation
  parents<-gen_data%>%filter(G==g)
   
  #### Fitness Estimation Module ####
  
  #Fitness surface from Lande and Arnold, determined by trait optima (theta values) and omega
  #Here, OMEGA is the number of standard deviations for the multivariate normal distribution
  #Which determines the strength of selection. We use mean theta values equal to the initial
  #trait mean in the population (MU_entry_day, MU_RLS) 
  
  #Draw a new theta value for each trait, for each generation. 
  #environmental variance (var_env) determines how much theta can vary from year to year
  THETA_entry_day <- rtnorm(n=1, mu = MU_entry_day, lb=1,ub=30, sd = sqrt(var_env_entry))
  THETA_RLS <- rtnorm(n=1, mu = MU_RLS, lb=1, ub=30-THETA_entry_day, sd = sqrt(var_env_RLS))
  
  #Fitness estimation of a two-trait surface, from Lande and Arnold
  parents$RS_exp<-exp((-1/(2*(1-(rho^2))*OMEGA))*(((parents$entry_day - (THETA_entry_day+0))/sqrt(variance_entry_day))^2 + 
                                     ((parents$RLS - (THETA_RLS+0))/sqrt(variance_RLS))^2 - 
                                     2*rho*((parents$entry_day - (THETA_entry_day+0))/sqrt(variance_entry_day))*((parents$RLS - 
                                     (THETA_RLS+0))/sqrt(variance_RLS))))

  #Draw Expected RS values from a Poisson distribution, using the 0-1 scale generated above as probabilities
  parents$RS_exp<-qpois(parents$RS_exp,lambda=2)
  
  
  #Size of the next generation is the sum of RS_exp values of females in the populations,
  #or we can hold the population size constant (Nc_initial is the carrying capacity)
  if(PopSize=="variable"){
    Nc<-round(parents%>%filter(sex=='f')%>%
                #group_by(spawning_location)%>%
                dplyr::summarize(Nc=sum(RS_exp))%>%
                dplyr::select(Nc)%>%as.matrix%>%c())
  }
  
  if(PopSize=="constant"){
    Nc<-Nc_initial
  }
  
  if(Nc<50){break} #If the population crashes below 50, assume extinction.
  
#### Reproduction Module ####
# Here we draw parents for the number of offspring in the next generation
  
#days_df provides one row for each day that each individual was available to mate  
  days_df<-as.data.frame(matrix(nrow=0,ncol=4)) 
  colnames(days_df)<-c('day','ID','sex')
  
  for (i in 1:nrow(parents)){
    ind_i<-parents$ID[i]	
    days<-seq(parents$entry_day[i],parents$entry_day[i]+parents$RLS[i])
    sex<-parents$sex[i]
    temp_df<-data.frame(day=days,ID=rep(ind_i,length(days)), sex=rep(sex,length(days)))
    days_df<-rbind(days_df,temp_df)
  }
  
#days_matrix is a large Nc x Nc matrix where data represent the number of days
#each pair of individuals overlapped in the stream - which are then turned into
#overlap weights
  
  days_matrix<-as.data.frame(matrix(nrow=nrow(parents),ncol=nrow(parents),data = 0))
  colnames(days_matrix)<-parents$ID
  rownames(days_matrix)<-parents$ID
  
  #If modeling an assortative mating system (recommended), calculate overlap weight (W[o i,j])
  #See Equation 1 in the text. 
  #W[o i,j] = 0 for same-sex pairs, self-self (diagonal), and individuals who do not overlap
  #If a possible mating could have occurred, W[o i,j] = 1. 
  #This makes the assumption that all pairs who could have mated, have an equal probability of mating
  #In reality, could weight this by the relative number of days individuals i and j overlapped,
  #and/or individuals who returned to the stream within x days of each other. |entry_day[i]-entry_day[j]| > x
  if(assortative==T){ 
    for (i in parents$ID){
      days<-days_df%>%filter(ID==i)%>%dplyr::select(day)%>%as.matrix()%>%c()
      overlaps<-days_df%>%filter(sex!=parents$sex[which(parents$ID==i)], #Opposite sex - same sex gets 0 probability of mating
                                 day%in%days)%>% #overlaps days - probability depends on number of days
        filter(!ID==i)%>% #remove ID i from list - can't mate with self
        group_by(ID)%>% 
        dplyr::summarize(n=n(),.groups="drop")
      days_matrix[as.character(overlaps$ID),as.character(i)]<- 1 #To weight by relative number days overlapped: <- overlaps$n/sum(overlaps$n)
      
    }
  }

  #Assign Parent Weights (W[i,j]) with overlap weight * fitness weight of dad * fitness weight of mum:
  
  #relative reproductive success values are not sex-specific.
  parents<-parents%>%mutate(Nc=Nc)
  parents<-parents%>%mutate(meanKexp=mean(RS_exp),fitness_weight=RS_exp/meanKexp)
  
  #Extract weights from days_matrix format to data frame format
  X<- which(upper.tri(days_matrix, diag = F), arr.ind = T)
  weights_df<-as.data.frame(cbind(X, days_matrix[X]))
  weights_df$row<-rownames(days_matrix)[weights_df$row]
  weights_df$col<-colnames(days_matrix)[weights_df$col]
  colnames(weights_df)<-c("parent1","parent2","overlap")
  
  
  #Get parent sexes
  weights_df$sex1<-parents$sex[match(weights_df$parent1,parents$ID)]
  weights_df$sex2<-parents$sex[match(weights_df$parent2,parents$ID)]
  weights_df<-weights_df%>%filter(sex1!=sex2) #get only opposite sexes
  
  #fix order of sexes so females in dam column, males in sire column
  weights_df[weights_df$sex1 == "f", c("parent1", "parent2","sex1","sex2")] <- weights_df[weights_df$sex1 == "f", c("parent2", "parent1","sex2","sex1")]
  colnames(weights_df)[1:2]<-c("sire","dam")
  
  #Add fitness weights (Relative reproductive success values)
  weights_df$dam_fitness<-parents$fitness_weight[match(weights_df$dam,parents$ID)]
  weights_df$sire_fitness<-parents$fitness_weight[match(weights_df$sire,parents$ID)]
  
  if(assortative==T){
    weights_df<-weights_df %>% mutate(weight=dam_fitness*sire_fitness*overlap)} #Overlap weights when assortative is true
  if(assortative==F){
    weights_df<-weights_df %>% mutate(weight=1) #No weighting for random mating  
  }
  


#Draw parents using weights
  
  #Initialize new data frame for offspring
  offspring<-data.frame(ID=seq(max(parents$ID)+1,max(parents$ID)+sum(Nc)),
                        dam=NA,
                        sire=NA,
                        gen=g+1)
  
  #Parents:
  offspring[,c(3,2)]<-weights_df[sample(x=1:nrow(weights_df),
                                 size=Nc,replace=T,
                                 prob=weights_df$weight),1:2]
  
  
  #### Inheritance Module ####
  #Here we assign traits to offspring
  
  #Get the trait values of their parents
  offspring$dam_return<-gen_data$entry_day[match(offspring$dam,gen_data$ID)]
  offspring$sire_return<-gen_data$entry_day[match(offspring$sire,gen_data$ID)]
  offspring$midparent_return<-rowMeans(offspring%>%dplyr::select(dam_return,sire_return))
  
  offspring$dam_RLS<-gen_data$RLS[match(offspring$dam,gen_data$ID)]
  offspring$sire_RLS<-gen_data$RLS[match(offspring$sire,gen_data$ID)]
  offspring$midparent_RLS<-rowMeans(offspring%>%dplyr::select(dam_RLS,sire_RLS))
  
  #Inheritance of correlated traits
  for(i in 1:nrow(offspring)){
    
      #Means equal midparent values
      u = c(offspring$midparent_RLS[i],offspring$midparent_return[i]) #midparent RLS , midparents entry_day
      
      #Draw entry day from truncated normal from 0 to 30
      offspring[i,"entry_day"]<-TruncatedNormal::rtnorm(n=1,mu=u[2],sd=sqrt(Sigma[2,2]),lb=0,ub=30) #30 day spawning season
      
      #RLS is correlated with entry day, so drawn from conditional mean and sd.
      cond_mean = u[1] + (Corr[1,2] * sqrt(Sigma[1,1])/sqrt(Sigma[2,2]) * (offspring[i,"entry_day"] - u[2]))
      cond_sd = sqrt(Sigma[1,1]) * sqrt(1 - Corr[1,2])
      
      #Individuals can't live past the last day of the season.
      offspring[i,"RLS"]<-TruncatedNormal::rtnorm(n=1,mu=cond_mean,sd=cond_sd,lb=0,ub=30-offspring[i,"entry_day"]) 
  }
  
  offspring$RLS<-round(offspring$RLS)
  offspring$entry_day<-round(offspring$entry_day)
  
  #Assuming a sex ratio of 1:1 and randomly draw sexes. 
  #Thus traits are completely independent of sex
  offspring$sex<-sample(c('m','f'),size = nrow(offspring),replace = T) 
  
  #Assign parents RS_obs
  parents$RS_obs<-0
  sires_RS_obs<-offspring%>%group_by(sire)%>%dplyr::summarize(RS_obs=n(),.groups="drop")%>%dplyr::rename(parent=sire) %>% as.data.frame()
  dams_RS_obs<-offspring%>%group_by(dam)%>%dplyr::summarize(RS_obs=n(),.groups="drop")%>%dplyr::rename(parent=dam)%>%as.data.frame()
  parents_RS_obs<-rbind(sires_RS_obs,dams_RS_obs)
  
  parents$RS_obs[match(parents_RS_obs$parent,parents$ID)]<-parents_RS_obs$RS_obs
  
  
 
  ############......Finish parents in gen_data################
  
  gen_data$RS_exp[which(gen_data$G==g)]<-parents$RS_exp
  gen_data$RS_obs[which(gen_data$G==g)]<-parents$RS_obs
  
  ############......Add offspring as next generation in gen_data################
  
  
  offspring$RS_obs<-NA
  offspring$RS_exp<-NA
  offspring$year<-max(gen_data$year)+4 #arbitrary generation length of 4 years. 
  
  gen_data_g<-offspring%>%dplyr::select(ID,gen,sire,dam,sex,entry_day,RLS,RS_exp,RS_obs,year)
  colnames(gen_data_g)<-colnames(gen_data)
  
  gen_data<-rbind(gen_data,gen_data_g)
  
  
} #end G Loop

#Add input parameters to gen data
gen_data$rho<-rho
gen_data$var_entry_day<-variance_entry_day
gen_data$var_RLS<-variance_RLS
gen_data$PopSize<-PopSize
gen_data$assortative<-assortative
gen_data$iter<-iteration
gen_data$var_env_entry = var_env_entry
gen_data$var_env_RLS = var_env_RLS
gen_data$OMEGA<- OMEGA


#At the end of each iteration, export gen_data to a .csv
#The title of this .csv contains all unique input parameters,
#so that files can easily read in and concatenated with rbind()

write.csv(gen_data,file=paste(PATH,OMEGA,rho,variance_entry_day,variance_RLS,
                              PopSize,assortative,
                              var_env_entry,var_env_RLS,
                              iteration,".csv",sep=""))

  
} #end function




#Register a cluster of cores to run the function on in parallel:
clust<-makeCluster(8)
registerDoParallel(cl = clust, cores = 8) #Run on 8 cores


#for each input parameter to vary, define a new parameter with desired values
#Here the new parameter is just the input parameter with a "z" on the end.

#For all combinations of parameters used in text:
foreach(iterz = 1:10, .packages = c("tidyverse","rockchalk","TruncatedNormal")) %:%
  foreach(rhoz = c(0,-0.3,-0.6)) %:%
  foreach(variance_entry_dayz = c(10,20,30)) %:%
  foreach(var_env_entryz = c(10,20,30)) %:%
  foreach(assortativez = c(TRUE,FALSE)) %dopar% {
  Assortative_Mating_IBM(G = 10,
                         iteration = iterz,
                         variance_entry_day = variance_entry_dayz,
                         variance_RLS = 20,
                         assortative = assortativez,
                         rho = rhoz,
                         var_env_entry = var_env_entryz)}

###Question 1: Assortative Mating
foreach(iterz = 1:100, .packages = c("dplyr","rockchalk","TruncatedNormal")) %:%
  foreach(variance_entry_dayz = c(10,20,30)) %:%
  foreach(assortativez = c(TRUE,FALSE)) %dopar% {
    Assortative_Mating_IBM(G = 10,
                           OMEGA = 2,
                           iteration = iterz,
                           variance_entry_day = variance_entry_dayz,
                           variance_RLS = 20,
                           assortative = assortativez,
                           rho = -0.3,
                           var_env_entry = 20,
                           var_env_RLS = 1,
                           PATH = "output/Assortative_Mating/gen_data_")}

###Question 2: Trait Correlations
foreach(iterz = 1:100, .packages = c("dplyr","rockchalk","TruncatedNormal")) %:%
  foreach(rhoz = c(0,-0.3,-0.6)) %dopar% {
    Assortative_Mating_IBM(G = 10,
                           OMEGA = 2,
                           iteration = iterz,
                           variance_entry_day = 20,
                           variance_RLS = 20,
                           rho = rhoz,
                           assortative = T,
                           var_env_entry = 20,
                           var_env_RLS = 1,
                           PATH = "output/Trait_Correlations/gen_data_")}

###Question 3: Environmental Variation
foreach(iterz = 1:100, .packages = c("dplyr","rockchalk","TruncatedNormal")) %:%
  foreach(var_env_entryz = c(10,20,30))  %dopar% {
    Assortative_Mating_IBM(G = 10,
                           OMEGA = 2,
                           iteration = iterz,
                           variance_entry_day = 20,
                           variance_RLS = 20,
                           rho = -0.3,
                           assortative = T,
                           var_env_entry = var_env_entryz,
                           var_env_RLS = 1,
                           PATH = "output/Environmental_Variation/gen_data_")}

###Question 4: Strength of Selection
foreach(iterz = 1:100, .packages = c("dplyr","rockchalk","TruncatedNormal")) %:%
  foreach(OMEGAz = c(1,2,3)) %dopar% {
    Assortative_Mating_IBM(G = 10,
                           OMEGA = OMEGAz,
                           iteration = iterz,
                           variance_entry_day = 20,
                           variance_RLS = 20,
                           rho = -0.3,
                           assortative = T,
                           var_env_entry = 20,
                           var_env_RLS = 1,
                           PATH = "output/Selection_Strength/gen_data_")}

stopImplicitCluster() #closes unused processor connections after running foreach
