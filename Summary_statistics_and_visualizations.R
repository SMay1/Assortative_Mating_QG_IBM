#title: "Sam_PopDy_Model"
#author: "Sam"
#date: "2/18/2020"

#Originally an Rmd file - transferred here to enable looping over multiple code chunks.


#### 1: Load Data and Packages ####

#Set Working Directory to Source File Location

library(tidyverse)
library(sjPlot)
library(AICcmodavg)
library(mgcv)
library(MASS)
library(ggpubr)
library(TruncatedNormal)
library(rockchalk)
#devtools::install_github("zeehio/facetscales")
library(facetscales)
library(RColorBrewer)


#Script to load simulation iterations, calculate summary statistics, and plot data.
#Summary statistics data provided. Individual model outputs are not provided, as they are very large.

###Because the output is very large, and because we are only using the summary statistics for each iteration
###Easier to calculate summary stats for each iter, and then combine after.
###If, in the future, one wants to load in all together - here are two different methods to do so:

#gen_data_all <- list.files(path = "../data/iters/",pattern = "*.csv") %>% 
#  map_df(~read_csv(paste("../data/iters/",.,sep="")))
 
#file_names = list.files(path = "../data/iters/",pattern="*.csv",full.names = T)
#gen_data_all<-do.call(rbind, lapply(file_names, read.csv, header = T))


#summary stats loop:
#Reads each output file, calculates summary stats, adds them to "summary_stats_all"

file_names = list.files(path = "../data/iters/",pattern="*.csv",full.names = T,include.dirs = F)
counter<-1
for(i in file_names){
  
  file_name<-substr(i,start = 28,nchar(i)-4)
  gen_data_i<-read.csv(i)

summary_stats_i<-gen_data_i %>% group_by(G,rho,var_entry_day,var_RLS,PopSize,assortative,iter,var_env_entry,var_env_RLS,OMEGA) %>% 
  dplyr::summarise(Nc = n(), #number of individuals per generation
                   meanK = mean(RS_obs), #mean RS
                   varK = var(RS_obs), #variance in RS
                   S = sum(RS_obs)/2, #sum of RS per gen/2 (size of next gen)
                   Ne = ((meanK*Nc)-2)/(meanK-1+(varK/meanK)), #Ne from Eq 5 of Waples et al. 2011 Inbreeding effective population size and parentage analysis without parents
                   Ne_N = Ne/Nc, #Ne/Nc ratio 
                   Nm = sum(sex=="m"),
                   Nf = sum(sex=="f"),
                   sex_ratio = Nm/Nf,
                   mean_entry_day = mean(entry_day),
                   sd_entry_day = sd(entry_day),
                   mean_RLS = mean(RLS),
                   sd_RLS = sd(RLS),
                   Nb = ((meanK*Nc)-2)/(meanK-1+(varK/meanK)), #Effective Number of Breeders
                   breeders = sum(RS_obs>0)) #Actual Number of Breeders

if(counter==1){summary_stats_all<-summary_stats_i}
if(counter>1){summary_stats_all<-rbind(summary_stats_all,summary_stats_i)}

counter<-counter+1
}

#save(list = "summary_stats_all",file = "../data/iters/summary_stats_all.Rdata")

#load("../data/iters/summary_stats_all.Rdata")
summary_stats_all<-read.csv("summary_stats_all.csv")


##Assortative Mating Plots
  Nc_plot_AM<-summary_stats_all %>% filter(sim=="AM") %>% 
  filter(rho==-0.3) %>% mutate(assortative=ifelse(assortative==T,"Assortative","Random")) %>% 
  group_by(var_entry_day,G,assortative) %>% dplyr::summarize(sd=sd(Nc),Nc=mean(Nc),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=Nc,color=as.factor(var_entry_day)))+
  geom_ribbon(aes(fill=as.factor(var_entry_day),x=G,ymin=Nc-error,ymax=Nc+error,group=var_entry_day),alpha=0.4,show.legend = F)+
  geom_point()+
  facet_grid(.~assortative)+
  theme_classic()+
  labs(x="Generation",color= "Variance in\nEntry Day")+
  scale_color_brewer(type="seq",palette = "Blues")+
  scale_fill_brewer(type="seq",palette="Blues")


Nc_plot_AM<-summary_stats_all %>% filter(sim=="AM") %>% 
  filter(rho==-0.3) %>% mutate(assortative=ifelse(assortative==T,"Assortative","Random")) %>% 
  group_by(var_entry_day,G,assortative) %>% dplyr::summarize(sd=sd(Nc),Nc=mean(Nc),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=Nc,color=as.factor(var_entry_day)))+
  geom_ribbon(aes(fill=as.factor(var_entry_day),x=G,ymin=Nc-error,ymax=Nc+error,group=var_entry_day),alpha=0.4,show.legend = F)+
  geom_point()+
  facet_grid(.~assortative)+
  theme_classic()+
  labs(x="Generation",color= "Variance in\nEntry Day")+
  scale_color_brewer(type="seq",palette = "Blues")+
  scale_fill_brewer(type="seq",palette="Blues")


Entry_day_plot_AM<-summary_stats_all %>% filter(sim=="AM") %>% 
  filter(rho==-0.3) %>% mutate(assortative=ifelse(assortative==T,"Assortative","Random")) %>% 
  group_by(var_entry_day,G,assortative) %>% dplyr::summarize(sd=sd(mean_entry_day),mean=mean(mean_entry_day),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=mean,color=as.factor(var_entry_day)))+
  geom_ribbon(aes(fill=var_entry_day,x=G,ymin=mean-error,ymax=mean+error,group=var_entry_day),alpha=0.3,show.legend = F)+
  geom_point()+
  facet_grid(.~assortative)+
  theme_classic()+
  labs(y = "Mean Entry Day",x="Generation",color= "Variance in\nEntry Day")

RLS_plot_AM<-summary_stats_all %>% filter(sim=="AM") %>% 
  filter(rho==-0.3) %>% mutate(assortative=ifelse(assortative==T,"Assortative","Random")) %>% 
  group_by(var_entry_day,G,assortative) %>% dplyr::summarize(sd=sd(mean_RLS),mean=mean(mean_RLS),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=mean,color=as.factor(var_entry_day)))+
  geom_ribbon(aes(fill=var_entry_day,x=G,ymin=mean-error,ymax=mean+error,group=var_entry_day),alpha=0.3,show.legend = F)+
  geom_point()+
  facet_grid(.~assortative)+
  theme_classic()+
  labs(y = "Mean RLS",x="Generation",color= "Variance in\nEntry Day")


###### Trait Correlation Plots (rho):

Nc_plot_TC<-summary_stats_all %>% filter(sim=="TC") %>% 
  mutate(assortative=ifelse(assortative==T,"Assortative","Random")) %>% 
  group_by(rho,G,assortative) %>% dplyr::summarize(sd=sd(Nc),Nc=mean(Nc),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=Nc,color=as.factor(rho)))+
  geom_ribbon(aes(fill=rho,x=G,ymin=Nc-error,ymax=Nc+error,group=rho),alpha=0.3,show.legend = F)+
  geom_point()+
  facet_grid(.~assortative)+
  theme_classic()+
  labs(x="Generation",color= "Rho")

Entry_day_plot_TC<-summary_stats_all %>% filter(sim=="TC") %>% 
  mutate(assortative=ifelse(assortative==T,"Assortative","Random")) %>% 
  group_by(rho,G,assortative) %>% dplyr::summarize(sd=sd(mean_entry_day),mean=mean(mean_entry_day),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=mean,color=as.factor(rho)))+
  geom_ribbon(aes(fill=rho,x=G,ymin=mean-error,ymax=mean+error,group=rho),alpha=0.3,show.legend = F)+
  geom_point()+
  facet_grid(.~assortative)+
  theme_classic()+
  labs(y = "Mean Entry Day",x="Generation",color= "Rho")

RLS_plot_TC<-summary_stats_all %>% filter(sim=="TC") %>% 
  mutate(assortative=ifelse(assortative==T,"Assortative","Random")) %>% 
  group_by(rho,G,assortative) %>% dplyr::summarize(sd=sd(mean_RLS),mean=mean(mean_RLS),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=mean,color=as.factor(rho)))+
  geom_ribbon(aes(fill=rho,x=G,ymin=mean-error,ymax=mean+error,group=rho),alpha=0.3,show.legend = F)+
  geom_point()+
  facet_grid(.~assortative)+
  theme_classic()+
  labs(y = "Mean RLS",x="Generation",color= "Rho")

###### Environmental Variability:


Nc_plot_EV<-summary_stats_all %>% filter(sim=="EV") %>% 
  mutate(assortative=ifelse(assortative==T,"Assortative","Random")) %>% 
  group_by(var_env_entry,G,assortative) %>% dplyr::summarize(sd=sd(Nc),Nc=mean(Nc),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=Nc,color=as.factor(var_env_entry)))+
  geom_ribbon(aes(fill=var_env_entry,x=G,ymin=Nc-error,ymax=Nc+error,group=var_env_entry),alpha=0.3,show.legend = F)+
  geom_point()+
  facet_grid(.~assortative)+
  theme_classic()+
  labs(x="Generation",color= "Variance in\nTrait Optima")

Entry_day_plot_EV<-summary_stats_all %>% filter(sim=="EV") %>% 
  mutate(assortative=ifelse(assortative==T,"Assortative","Random")) %>% 
  group_by(var_env_entry,G,assortative) %>% dplyr::summarize(sd=sd(mean_entry_day),mean=mean(mean_entry_day),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=mean,color=as.factor(var_env_entry)))+
  geom_ribbon(aes(fill=var_env_entry,x=G,ymin=mean-error,ymax=mean+error,group=var_env_entry),alpha=0.3,show.legend = F)+
  geom_point()+
  facet_grid(.~assortative)+
  theme_classic()+
  labs(y = "Mean Entry Day",x="Generation",color= "Variance in\nTrait Optima")

RLS_plot_EV<-summary_stats_all %>% filter(sim=="EV") %>% 
  mutate(assortative=ifelse(assortative==T,"Assortative","Random")) %>% 
  group_by(var_env_entry,G,assortative) %>% dplyr::summarize(sd=sd(mean_RLS),mean=mean(mean_RLS),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=mean,color=as.factor(var_env_entry)))+
  geom_ribbon(aes(fill=var_env_entry,x=G,ymin=mean-error,ymax=mean+error,group=var_env_entry),alpha=0.3,show.legend = F)+
  geom_point()+
  facet_grid(.~assortative)+
  theme_classic()+
  labs(y = "Mean RLS",x="Generation",color= "Variance in\nTrait Optima")

###### Selection Plots (omega):



Nc_plot_SS<-summary_stats_all %>% filter(sim=="SS") %>% 
  mutate(assortative=ifelse(assortative==T,"Assortative","Random")) %>% 
  group_by(OMEGA,G) %>% dplyr::summarize(sd=sd(Nc),Nc=mean(Nc),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=Nc,color=as.factor(OMEGA)))+
  geom_ribbon(aes(fill=OMEGA,x=G,ymin=Nc-error,ymax=Nc+error,group=OMEGA),alpha=0.3,show.legend = F)+
  geom_point()+
  theme_classic()+
  labs(x="Generation",color= "Strength of Selection")

Entry_day_plot_SS<-summary_stats_all %>% filter(sim=="SS") %>% 
  mutate(assortative=ifelse(assortative==T,"Assortative","Random")) %>% 
  group_by(OMEGA,G,assortative) %>% dplyr::summarize(sd=sd(mean_entry_day),mean=mean(mean_entry_day),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=mean,color=as.factor(OMEGA)))+
  geom_ribbon(aes(fill=OMEGA,x=G,ymin=mean-error,ymax=mean+error,group=OMEGA),alpha=0.3,show.legend = F)+
  geom_point()+
  theme_classic()+
  labs(y = "Mean Entry Day",x="Generation",color= "Strength of Selection")

RLS_plot_SS<-summary_stats_all %>% filter(sim=="SS") %>% 
  mutate(assortative=ifelse(assortative==T,"Assortative","Random")) %>% 
  group_by(OMEGA,G,assortative) %>% dplyr::summarize(sd=sd(mean_RLS),mean=mean(mean_RLS),error=qnorm(0.975)*sd/sqrt(n())) %>% 
  ggplot(aes(x=G,y=mean,color=as.factor(OMEGA)))+
  geom_ribbon(aes(fill=OMEGA,x=G,ymin=mean-error,ymax=mean+error,group=OMEGA),alpha=0.3,show.legend = F)+
  geom_point()+
  theme_classic()+
  labs(y = "Mean RLS",x="Generation",color= "Strength of Selection")



####Combination plots:

summary_stats_all<- summary_stats_all %>% mutate(breeders = breeders/Nc) %>% gather("parameter","value",Nc, Ne, Ne_N, breeders, mean_entry_day, mean_RLS) %>% 
  group_by(G,rho,var_entry_day,var_RLS,PopSize,assortative,var_env_entry,var_env_RLS,OMEGA,parameter,sim) %>% 
  dplyr::summarise(mean=mean(value,na.rm=T),sd=sd(value,na.rm=T),error = qnorm(0.975)*sd/sqrt(n())) %>% filter(G<11) %>% 
  mutate(parameter=ifelse(parameter=="Ne_N","Ne / Nc",
                          ifelse(parameter=="breeders","Breeding Proportion",
                                 ifelse(parameter=="mean_entry_day","Return Day",
                                        ifelse(parameter=="mean_RLS","RLS",parameter)))))


scales_y <- list(
  'Nc' = scale_y_continuous(limits = c(0, 5000)),
  'Ne / Nc' = scale_y_continuous(limits = c(0.4, 0.8)),
  'Return Day' = scale_y_continuous(limits = c(9, 11.5)),
  'RLS' = scale_y_continuous(limits = c(6,8)))

AM_plot<-summary_stats_all %>% mutate(assortative=ifelse(assortative==T,"Assortative","Random")) %>% 
  filter(sim=="AM") %>% filter(parameter%in%c("Nc","Ne / Nc","Return Day","RLS")) %>% #"mean_entry_day","mean_RLS"
  mutate(parameter=factor(parameter,levels = c("Nc","Ne / Nc","Return Day","RLS"))) %>% 
  ggplot(aes(x=G,y=mean,color=as.factor(var_entry_day)))+
  geom_ribbon(aes(fill=as.factor(var_entry_day),x=G,ymin=mean-error,ymax=mean+error,group=var_entry_day),alpha=0.4,show.legend=F)+
  geom_point()+
  facet_grid(parameter~assortative,scales="free_y")+
  theme_classic()+
  labs(x="Generation",y="Mean of 100 Model Iterations",color= "Variance in\nReturn Day")+
  scale_color_manual(values=brewer.pal(n = 6, "Blues")[3:6])+
  scale_fill_manual(values=brewer.pal(n = 6, "Blues")[3:6])+
  scale_y_continuous()+
  scale_x_continuous(breaks = c(0,5,10))+
  theme(legend.position = "top")

VED_plot<-summary_stats_all %>% #Variance in Entry day for assortative = T only
  filter(sim=="AM",assortative==T) %>% filter(parameter%in%c("Nc","Ne / Nc","Return Day","RLS")) %>% #"mean_entry_day","mean_RLS"
  mutate(parameter=factor(parameter,levels = c("Nc","Ne / Nc","Return Day","RLS"))) %>% 
  ggplot(aes(x=G,y=mean,color=as.factor(var_entry_day)))+
  geom_ribbon(aes(fill=as.factor(var_entry_day),x=G,ymin=mean-error,ymax=mean+error,group=var_entry_day),alpha=0.4,show.legend=F)+
  geom_point()+
  facet_grid_sc(rows=vars(parameter),scales=list(y=scales_y))+
  theme_classic()+
  labs(x="Generation",y="Mean of 100 Model Iterations",color= "Variance in\nReturn Day")+
  scale_color_manual(values=brewer.pal(n = 6, "Blues")[3:6])+
  scale_fill_manual(values=brewer.pal(n = 6, "Blues")[3:6])+
  scale_y_continuous()+
  scale_x_continuous(breaks = c(0,5,10))+
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8))


TC_plot<-summary_stats_all %>% 
  filter(sim=="TC") %>% filter(parameter%in%c("Nc","Ne / Nc","Return Day","RLS")) %>% #"mean_entry_day","mean_RLS"
  mutate(parameter=factor(parameter,levels = c("Nc","Ne / Nc","Return Day","RLS"))) %>% 
  ggplot(aes(x=G,y=mean,color=as.factor(rho)))+
  geom_ribbon(aes(fill=as.factor(rho),x=G,ymin=mean-error,ymax=mean+error,group=rho),alpha=0.4,show.legend=F)+
  geom_point()+
  #facet_grid(parameter~.,scales="free_y")+
  theme_classic()+
  labs(x="Generation",y="",color= "Correlation\nBetween Traits")+
  scale_color_manual(values=brewer.pal(n = 6, "Oranges")[3:6])+
  scale_fill_manual(values=brewer.pal(n = 6, "Oranges")[3:6])+
  scale_x_continuous(breaks = c(0,5,10))+
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8))+
  facet_grid_sc(rows=vars(parameter),scales=list(y=scales_y))

EV_plot<-summary_stats_all %>% 
  filter(sim=="EV") %>% filter(parameter%in%c("Nc","Ne / Nc","Return Day","RLS")) %>% #"mean_entry_day","mean_RLS"
  mutate(parameter=factor(parameter,levels = c("Nc","Ne / Nc","Return Day","RLS"))) %>% 
  ggplot(aes(x=G,y=mean,color=as.factor(var_env_entry)))+
  geom_ribbon(aes(fill=as.factor(var_env_entry),x=G,ymin=mean-error,ymax=mean+error,group=var_env_entry),alpha=0.4,show.legend=F)+
  geom_point()+
  #facet_grid(parameter~.,scales="free_y")+
  theme_classic()+
  labs(x="Generation",y="",color= "Variance in\nTrait Optimum")+
  scale_color_manual(values=brewer.pal(n = 6, "Greens")[3:6])+
  scale_fill_manual(values=brewer.pal(n = 6, "Greens")[3:6])+
  scale_x_continuous(breaks = c(0,5,10))+
  theme(legend.position = "top",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y=element_blank(),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8))+
  facet_grid_sc(rows=vars(parameter),scales=list(y=scales_y))

SS_plot<-summary_stats_all %>% 
  filter(sim=="SS") %>% filter(parameter%in%c("Nc","Ne / Nc","Return Day","RLS")) %>% #"mean_entry_day","mean_RLS"
  mutate(parameter=factor(parameter,levels = c("Nc","Ne / Nc","Return Day","RLS"))) %>% 
  ggplot(aes(x=G,y=mean,color=as.factor(OMEGA)))+
  geom_ribbon(aes(fill=as.factor(OMEGA),x=G,ymin=mean-error,ymax=mean+error,group=OMEGA),alpha=0.4,show.legend=F)+
  geom_point()+
  facet_grid(parameter~.,scales="free_y")+
  theme_classic()+
  labs(x="Generation",y="",color= "Selection\nIntensity")+
  scale_color_manual(values=brewer.pal(n = 6, "Purples")[3:6])+
  scale_fill_manual(values=brewer.pal(n = 6, "Purples")[3:6])+
  scale_y_continuous()+
  scale_x_continuous(breaks = c(0,5,10))+
  theme(legend.position = "top",
        axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8))+
  facet_grid_sc(rows=vars(parameter),scales=list(y=scales_y))

ggarrange(VED_plot,TC_plot,EV_plot,SS_plot,ncol=4,nrow=1,align = "h")




## Get stats from G10 for results section:
summary_stats_all %>% filter(G==10,
                             sim=="SS",
                             !parameter%in%c("Breeding Proportion","Ne")) %>% View()

summary_stats_all %>% filter(G==10,
                             sim=="AM",
                             !parameter%in%c("Breeding Proportion","Ne")) %>% 
  mutate(var_entry_day=factor(var_entry_day,levels=c("10","20","30"))) %>% 
  ggplot(aes(x=var_entry_day,y=mean,fill=assortative))+
  geom_col(position = position_dodge())+
  facet_grid(parameter~.,scales="free")+
  geom_errorbar(aes(ymax=mean+error,ymin=mean-error),width=0.2,position=position_dodge(width=0.9))+
  theme_classic()

scales_y <- list(
  'Nc' = scale_y_continuous(limits = c(0, 5000)),
  'Ne / Nc' = scale_y_continuous(limits = c(0.4,1)),
  'Return Day' = scale_y_continuous(limits = c(9, 14)),
  'RLS' = scale_y_continuous(limits = c(6.5,10)))


G10_AM<- summary_stats_all %>% filter(G==10,
                             sim=="AM",
                             !parameter%in%c("Breeding Proportion","Ne")) %>% 
  mutate(var_entry_day=factor(var_entry_day,levels=c("10","20","30")),
         assortative=ifelse(assortative==TRUE,"Assortative","Random")) %>% 
  ggplot(aes(x=var_entry_day,y=mean,pch=assortative,col=var_entry_day))+
  geom_point()+
  facet_grid(parameter~.,scales="free")+
  geom_errorbar(aes(ymax=mean+error,ymin=mean-error),width=0.2)+
  theme_classic()+
  theme(legend.margin=margin(),legend.box="vertical",
        legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.border = element_blank())+
  labs(y="Mean of 100 Model Iterations",x = expression(sigma["Return Day"]^2),
       pch="Mating System", col= expression(sigma["Return Day"]^2))+
  scale_color_manual(values=brewer.pal(n = 6, "Blues")[3:6])+
  facet_grid_sc(rows=vars(parameter),scales=list(y=scales_y))

G10_TC<- summary_stats_all %>% filter(G==10,
                             sim=="TC",
                             !parameter%in%c("Breeding Proportion","Ne")) %>% 
  mutate(rho=factor(rho,levels=c("-0.6","-0.3","0"))) %>%  
  ggplot(aes(x=rho,y=mean,col=rho))+
  geom_point()+
  facet_grid(parameter~.,scales="free")+
  geom_errorbar(aes(ymax=mean+error,ymin=mean-error),width=0.2)+
  theme_classic()+
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.border = element_blank())+
  labs(y="Mean of 100 Model Iterations",x = expression(rho),
       col=expression(rho))+
  scale_color_manual(values=brewer.pal(n = 6, "Oranges")[3:6])+
  facet_grid_sc(rows=vars(parameter),scales=list(y=scales_y))


G10_EV<-summary_stats_all %>% filter(G==10,
                             sim=="EV",
                             !parameter%in%c("Breeding Proportion","Ne")) %>% 
  mutate(var_env_entry =factor(var_env_entry ,levels=c("10","20","30"))) %>%  
  ggplot(aes(x=var_env_entry ,y=mean,color=var_env_entry))+
  geom_point()+
  facet_grid(parameter~.,scales="free")+
  geom_errorbar(aes(ymax=mean+error,ymin=mean-error),width=0.2)+
  theme_classic()+
  theme(legend.position = "top",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y=element_blank(),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.border = element_blank())+
  labs(y="Mean of 100 Model Iterations",x = expression(sigma[theta]^2),
       color=expression(sigma[theta]^2))+
  scale_color_manual(values=brewer.pal(n = 6, "Greens")[3:6])+
  facet_grid_sc(rows=vars(parameter),scales=list(y=scales_y))


G10_SS <- summary_stats_all %>% filter(G==10,
                             sim=="SS",
                             !parameter%in%c("Breeding Proportion","Ne")) %>% 
  mutate(OMEGA =factor(OMEGA ,levels=c("1","2","3"))) %>%  
  ggplot(aes(x=OMEGA ,y=mean,color=OMEGA))+
  geom_point()+
  facet_grid(parameter~.,scales="free")+
  geom_errorbar(aes(ymax=mean+error,ymin=mean-error),width=0.2)+
  theme_classic()+
  theme(legend.position = "top",
        axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.border = element_blank())+
  labs(y="Mean of 100 Model Iterations",x = expression(omega),
       col=expression(omega))+
  scale_color_manual(values=brewer.pal(n = 6, "Purples")[3:6])+
  facet_grid_sc(rows=vars(parameter),scales=list(y=scales_y))

ggarrange(G10_AM,G10_TC,G10_EV,G10_SS,nrow=1,ncol=4,align = "h")

