# initialization ----------------------------------------------------------
rm(list=ls()) #clear all variables

library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggforce)
library(dplyr)
library(survival)
library(survminer)
library(stringr)
library(janitor)
library(cowplot)
library(ggrepel)
library(ggstance)
library(grDevices)
library(gridExtra)
library(igraph)


#### Figure 1A - Plot temperature in each treatment ####

# Import temperature treatment data
temps <- read.csv("Temperature profiles.csv")

F1A <- ggplot(temps, aes(x=Days, y=Temperature, color=Treatment, fill=Treatment))+
  geom_line(aes(color=Treatment), linewidth=1)+
  labs(x="Days", y="Temperature (°C)")+
  scale_color_manual(values=c("#03cafc", "#f24455")) +
  theme_classic(base_size = 8)+
  theme(legend.position = c(0.8,0.3))+
  scale_x_continuous(breaks=seq(1,9,1)) +
  scale_y_continuous(breaks=seq(27,36,2)); F1A


#### Import survival data #####
survival <- read.csv("survival.csv")

##### Figure 1B - Plot by treatment for WT larvae only ####

### Keep only WT
WT_survival <- survival[survival$Parental_cross=="WT",]

### Model WT by treatment
WT_treatment_fit <- surv_fit(coxph(Surv(timepoint, mortality) ~ treatment, data=WT_survival), data=WT_survival)

F1B <- ggsurvplot(WT_treatment_fit, WT_survival, 
           facet.by=c("Parental_cross"),
           palette = c("#03cafc", "#f24455"),
           pval = FALSE,
           legend.title = "",
           conf.int=TRUE)+
  theme_classic(base_size=8)+
  theme(legend.position="none",
        strip.background=element_blank(),
        strip.text=element_blank())+
  labs(y="Survival Probability",x="Days")+
  annotate("text",x=1.5,y=0.05,label="p<0.0001",size=3)+
  xlab("Days")$plot; F1B

#### Make Figure 1 
quartz(w=4.3,h=2.1)
plot_grid(F1A,F1B,align="h",axis="tb",rel_widths = c(1,1),labels="AUTO",label_size=10)

## Test whether treatments differ for WT larvae (Model 1)
survdiff(Surv(timepoint, mortality) ~ treatment, data = WT_survival) #p = 2e-16 

##### Figure 2A - Plot survival by parental crosses and treatment, removing WT larvae ####

### Remove WT larvae
survival_no_WT <- survival%>%filter(Parental_cross!="WT")

### With WT larvae removed, plot treatment by parental cross (Model 2)
Par_treatment_fit_2 <- surv_fit(coxph(Surv(timepoint, mortality) ~ treatment, data=survival_no_WT), data=survival_no_WT)

update_geom_defaults("text",list(size=2))

F2A <- ggsurvplot(Par_treatment_fit_2, data=survival_no_WT, 
                 facet.by=c("Father","Mother"),
                 palette = c("#03cafc", "#f24455"),
                 pval = TRUE,
                 ggtheme=theme_classic(base_size=8))+
                 xlab("Days")+
                 ylab("Survival Probability")+
                 theme(legend.position="none"); F2A


## Overall, do treatments differ?
survdiff(Surv(timepoint, mortality) ~ treatment, data = survival_no_WT) # Yes. p = 2e-13 


##### Figure 2B - Plot number of larvae per larval family ####
# Import larval count data
larval_count <- read.csv("larval_count.csv")    
  
mean(larval_count$Total)
sd(larval_count$Total)
  
# Order crosses
larval_count$Cross = factor(larval_count$Cross,levels=c("WT","A1","B1","C1","D1","E1","F1",
                                                          "A2","C2","D2","E2","F2",
                                                          "A3","C3","D3",
                                                          "A4","B4","C4","D4","E4","F4",
                                                          "A5","D5","F5",
                                                          "A6","C6","D6","E6"))
  
#### Plot number of larvae produced for each parental cross
F2B <-  ggplot(larval_count, aes(x = Total, y = reorder(Cross, desc(Cross)))) +
  geom_point(color="black", size=0.5) +
  labs(x="N Larvae\nProduced", y="Parental Cross")+
  scale_x_continuous(breaks=c(0,500,1000))+
  theme_classic(base_size = 8)+
  annotate("text",y=27,x=700,label="self-cross",size=2,fontface="italic")+
  annotate("text",y=15,x=700,label="self-cross",size=2,fontface="italic")+
  annotate("text",y=10,x=700,label="self-cross",size=2,fontface="italic"); F2B

### Make Figure 2
quartz(w=5.152941, h=3.796078)
plot_grid(F2A,F2B, align="h",axis="tb",rel_widths = c(4,1),labels="AUTO",label_size=10)


##### Get main effects of Mother and Father and their interaction for heated treatment only, with failed cross and WT larvae removed ####

### Keep only heated treatment, Remove WT, and cross D1 
Heat_survival_no_WT <- survival%>%filter(Parental_cross!="WT", treatment=="Heated", Parental_cross!="D1")

### (Model 3)
Heat_survival_no_WT_fit <- coxph(Surv(timepoint, mortality) ~ Mother * Father + frailty(larval_ID), iter.max=2000, 
                                 data=Heat_survival_no_WT)  

cox.zph(Heat_survival_no_WT_fit) # Mother and Mother*Father significant, but not Father alone
ggforest(Heat_survival_no_WT_fit, Heat_survival_no_WT)


####### Figure 3A - Plot hazard ratio in heated treatment for larval families as a heatmap ####
### Keep WT larvae, remove ambient treatment and failed cross
Heated_survival <- survival %>% filter(treatment=="Heated", Parental_cross!="D1")

## Set WT as reference
Heated_survival$Mother <- factor(Heated_survival$Mother, levels = c("WT","1","2","3","4","5","6"))
Heated_survival$Father <- factor(Heated_survival$Father, levels = c("WT","A","B","C","D","E","F"))
Heated_survival$Parental_cross <- factor(Heated_survival$Parental_cross, levels = c("WT","A1","A2","A3","A4","A5","A6",
                                                                                    "B1","B4",
                                                                                    "C1","C2","C3","C4","C6",
                                                                                    "D2","D3","D4","D5","D6",
                                                                                    "E1","E2","E4","E6",
                                                                                    "F1","F2","F4","F5"))

## Model survival for larval family in heated treatment to extract hazard ratios (Model 4)
par_cross_heated <- coxph(Surv(timepoint, mortality) ~ Parental_cross + frailty(larval_ID), data=Heated_survival)
cox.zph(par_cross_heated)
ggforest(par_cross_heated, data=Heated_survival)

# Extract model stats from model
mod_par_cross_heated <- summary(par_cross_heated)$conf.int
mod_par_cross_heated <- cbind(Parent = rownames(mod_par_cross_heated), mod_par_cross_heated)
rownames(mod_par_cross_heated) <- 1:nrow(mod_par_cross_heated)
mod_par_cross_heated <- as.data.frame(mod_par_cross_heated)
colnames(mod_par_cross_heated)[colnames(mod_par_cross_heated) == 'exp(coef)'] <- 'exp_coef'
mod_par_cross_heated$exp_coef <- as.numeric(mod_par_cross_heated$exp_coef)
colnames(mod_par_cross_heated)[colnames(mod_par_cross_heated) == 'upper .95'] <- 'upper_95'
colnames(mod_par_cross_heated)[colnames(mod_par_cross_heated) == 'lower .95'] <- 'lower_95'
mod_par_cross_heated$upper_95 <- as.numeric(mod_par_cross_heated$upper_95)
mod_par_cross_heated$lower_95 <- as.numeric(mod_par_cross_heated$lower_95)

# Remove extra columns
mod_par_cross_heated <- mod_par_cross_heated[-c(3)]

# Make new column with mother and father identification
mod_par_cross_heated <- mod_par_cross_heated %>%mutate(across('Parent', str_replace, 'Parental_cross', ''))
mod_par_cross_heated <- mod_par_cross_heated %>% mutate(Mother = substr(mod_par_cross_heated$Parent,2,2))
mod_par_cross_heated <- mod_par_cross_heated %>% mutate(Father = substr(mod_par_cross_heated$Parent,1,1))

# Remove WT from data frame and add failed crosses in order to plot
mod_par_cross_heated_noWT <- mod_par_cross_heated %>% filter(Parent!="WT")%>%
  add_row(Parent=c("B2","B3","B5","B6","C5","D1","E3","E5","F3","F6"))%>%
  select(-Mother,-Father)%>%separate(Parent,into=c("Father","Mother"),sep=1,remove=FALSE)

F3A <-  ggplot(mod_par_cross_heated_noWT, aes(x=Mother, y=reorder(Father, desc(Father)), fill=exp_coef)) +
    geom_tile() +
    scale_fill_gradient2(low="#0b8c19",mid="#fdfffc",high="#d1261d", midpoint=1,na.value="lightgray") +
    labs(x="Mother", y="Father", fill="Hazard \nRatio") +
    theme_classic(base_size=8)+
    theme(legend.position="left"); F3A
  
##### Figure 3B - Plot hazard ratio of larval families for heated treatment ####

# Add WT back to data set in order to plot
WT_heat <- data.frame(Parent=c("WT"),
                        exp_coef=c(1),
                        lower_95=c(0),
                        upper_95=c(0),
                        Mother=("WT"),
                        Father=("WT"))
  
# Combine datasets
mod_par_cross_heated <- rbind(mod_par_cross_heated, WT_heat)
  
# Order crosses
mod_par_cross_heated$Parent = factor(mod_par_cross_heated$Parent,levels=c("WT","A1","B1","C1","D1","E1","F1",
                                                                            "A2","C2","D2","E2","F2",
                                                                            "A3","C3","D3",
                                                                            "A4","B4","C4","D4","E4","F4",
                                                                            "A5","D5","F5",
                                                                            "A6","C6","D6","E6"))
  
#### Plot hazard ratio in heated treatment
F3B <-  ggplot(mod_par_cross_heated, aes(x = exp_coef, y = reorder(Parent, desc(Parent)))) +
    geom_point(color="#f24455",size=0.7) +
    geom_errorbar(aes(xmin=lower_95, xmax=upper_95), width=0, color="#f24455") +
    geom_vline(xintercept = 1, lwd=0.5) +
    scale_x_continuous(trans="log10", breaks = c(0.01,0.05,0.2,1,5)) +
    labs(x="Hazard Ratio" ,y="Parental Cross")+
    theme_classic(base_size = 8)+
    annotate("text",x=0.01, y=26.5, label="More \ntolerant", size=2,hjust=0)+
    annotate("text",x=7, y=26.5, label="Less \ntolerant", size=2,hjust=1); F3B

### Make Figure 3
quartz(w=5.231372, h=3.043137)
plot_grid(F3A,F3B, align="h",axis="tb",rel_widths = c(2,1),labels="AUTO",label_size=10)
  
  

####### Figure 4A - Extract hazard ratios in heated treatment for mothers and fathers separately to plot hazard ratios by parent ####

## Model mothers in comparison to WT (Model 5)
  Mot_heated_fit <- coxph(Surv(timepoint, mortality) ~ Mother + frailty(larval_ID), data=Heated_survival)
  cox.zph(Mot_heated_fit)
  summary(Mot_heated_fit)
  ggforest(Mot_heated_fit, data=Heated_survival)
  
# Extract model stats from mother heated treatment model
  mod_Mot_heated <- summary(Mot_heated_fit)$conf.int
  mod_Mot_heated <- as.data.frame(mod_Mot_heated)
  mod_Mot_heated <- tibble::rownames_to_column(mod_Mot_heated, "Parent")
  mod_Mot_heated <- mod_Mot_heated[-c(3)]
  
  mod_Mot_heated <- mod_Mot_heated %>%mutate(across('Parent', str_replace, 'Mother', ''))
  
## Model fathers in comparison to WT (Model 6)
  Fat_heated_fit <- coxph(Surv(timepoint, mortality) ~ Father + frailty(larval_ID), data=Heated_survival)
  cox.zph(Fat_heated_fit)
  summary(Fat_heated_fit)
  ggforest(Fat_heated_fit, data=Heated_survival)
  
# Extract model stats from father heated treatment model
  mod_Fat_heated <- summary(Fat_heated_fit)$conf.int
  mod_Fat_heated <- as.data.frame(mod_Fat_heated)
  mod_Fat_heated <- tibble::rownames_to_column(mod_Fat_heated, "Parent")
  mod_Fat_heated <- mod_Fat_heated[-c(3)]
  
  mod_Fat_heated <- mod_Fat_heated %>%mutate(across('Parent', str_replace, 'Father', ''))
  
# Combine dataframes
  mod_Mot_Fat_heated <- rbind(mod_Mot_heated,mod_Fat_heated)

# Reorder dataframes
  mod_Mot_Fat_heated$Parent = factor(mod_Mot_Fat_heated$Parent,levels=c("1","2","3","4","5","6","A","B","C","D","E","F"))

# Add Mother/Father identifier column   
  mod_Mot_Fat_heated<-mod_Mot_Fat_heated%>%rownames_to_column(var="r")%>%mutate(r=as.numeric(r))%>%
    mutate(group=case_when(r<=6~"Mother",TRUE~"Father"))

# Order group column in order to plot    
  mod_Mot_Fat_heated$group = factor(mod_Mot_Fat_heated$group,levels=c("Mother","Father"))
  
## Make functions for below plot
  tag_facet1 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                         hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
    
    gb <- ggplot_build(p)
    lay <- gb$layout$layout
    tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
    p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                  vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
  }
  
  my_tag1 <- c("","Mother p<0.001 \nFather p=0.622 \nInteraction p=0.001")
  
  
  tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                         hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
    
    gb <- ggplot_build(p)
    lay <- gb$layout$layout
    tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
    p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                  vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
  }
  
  my_tag2 <- c("Less tolerant \nMore tolerant","")
  

  
#### Plot hazard ratio by parent
F4A <- ggplot(mod_Mot_Fat_heated, aes(x=Parent, y=`exp(coef)`))+
    geom_abline(intercept = 1, slope = 0) +
    geom_errorbar(aes(ymin=`lower .95`, ymax=`upper .95`), width=0, color="gray") +
    geom_point(size=2, color="navy") +
    labs(x="Parent", y="Hazard Ratio")+
    scale_y_continuous(breaks=seq(-2,4,1)) +
    theme_classic(base_size = 12) +
    facet_wrap(~group,scale="free_x"); F4A

F4A <-tag_facet1(F4A, 
              x = 0, y = 3.2, 
              vjust = 0, hjust = 0,
              open = "", close = "",
              fontface = 'italic',
              size = 2.5,
              tag_pool = my_tag1); F4A

F4A <-tag_facet1(F4A, 
                 x = 1.4, y = 0.82, 
                 vjust = 0, hjust = 0,
                 open = "", close = "",
                 fontface = 'italic',
                 size = 2.5,
                 tag_pool = my_tag2); F4A
  
  
############ Figure 4B -  Using hazard ratios, plot mothers and fathers against each other #######
  
colnames(mod_Mot_heated)[colnames(mod_Mot_heated) == 'exp(coef)'] <- 'coef_Mot'
colnames(mod_Fat_heated)[colnames(mod_Fat_heated) == 'exp(coef)'] <- 'coef_Fat'
colnames(mod_Mot_heated)[colnames(mod_Mot_heated) == 'Parent'] <- 'Parent_M'
colnames(mod_Fat_heated)[colnames(mod_Fat_heated) == 'Parent'] <- 'Parent_F'
  
mod_Mot_heated <- mod_Mot_heated[-c(3,4)]
mod_Fat_heated <- mod_Fat_heated[-c(3,4)]

mod_Mot_Fat_heated_wide <- cbind(mod_Mot_heated,mod_Fat_heated)
mod_Mot_Fat_heated_wide$Colony <- paste(mod_Mot_Fat_heated_wide$Parent_M,mod_Mot_Fat_heated_wide$Parent_F)


## When plotting mother against father, is there a relationship between hazard ratios? (Model 7)
mod <- lm(coef_Mot ~coef_Fat, data=mod_Mot_Fat_heated_wide)
summary(mod) ## No, Multiple R-squared=0.01203; p=0.8362

# Plot
F4B <- ggplot(mod_Mot_Fat_heated_wide, aes(x=coef_Mot, y=coef_Fat))+
  geom_smooth(method = "lm", se = FALSE, color="lightgray") +
  geom_point(size=2, color="navy") +
  geom_text(aes(label=Colony, x=coef_Mot+0, y=coef_Fat+.15), size=3) +
  labs(x="Hazard Ratio Mother", y="Hazard Ratio Father")+
  scale_x_continuous(limits = c(0,3)) +
  scale_y_continuous(limits = c(0,3)) +
  geom_abline(intercept = 0, slope = 1) +
  annotate("text", x=2.3, y=0.1, label="p = 0.836 \nr-squared = 0.012",size=3) +
  theme_classic()+
  theme(legend.position = "none"); F4B


### Make Figure 4
quartz(w=5.529412, h=3.184314)
plot_grid(F4A,F4B, align="h",axis="tb",rel_widths = c(1,1),labels="AUTO",label_size=10)


##### Figure 5 - Plot distribution of hazard ratios (Model 8) ####
haz_dist <-mod_par_cross_heated %>% filter(Parent!="WT")%>%
  mutate(log=log(exp_coef))
mean(haz_dist$log);sd(haz_dist$log)

set.seed(3839)
out<-as.data.frame(rnorm(5000,mean(haz_dist$log),sd(haz_dist$log)))%>%
  rename(dist=1)%>%mutate(exp_coef=exp(dist))

cdf<-ecdf(out$exp_coef)
intervals<-c(0.01,0.05,0.1,0.2,0.25,0.33,0.5)
cdf_out<-data.frame(interval=numeric(),cut_percent=numeric())
for (i in intervals){
  cut_percent<-cdf(i)
  cdf_out<-cdf_out%>%add_row(interval=i,cut_percent=cut_percent)%>%mutate(times=1/interval)
}

# Plot
quartz(w=2.7,h=2.7)
a<-ggplot(out,aes(exp_coef))+
  geom_vline(xintercept=0.25,linetype="dotted",color="lightgray")+
  geom_vline(xintercept=0.5,linetype="dotted",color="lightgray")+
  geom_hline(yintercept=0.0914,linetype="dotted",color="lightgray")+
  geom_hline(yintercept=0.247,linetype="dotted",color="lightgray")+
  stat_ecdf(geom = "step")+
  theme_classic(base_size=8)+
  xlab("Hazard Ratio")+
  ylab("Cumulative Proportion")+
  coord_cartesian(xlim=c(0,1),ylim=c(0,0.4))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  annotate("text", x=0.55,y=0.2,label="24.7% at 2x\nWT tolerance",angle=0,size=2.5,hjust=0)+
  annotate("text", x=0.3,y=0.05,label="9.1% at 4x\nWT tolerance",angle=0,size=2.5,hjust=0)+
  annotate("segment",x=0.25,xend=0.3,y=0.09,yend=0.05,color="red")+
  annotate("segment",x=0.5,xend=0.55,y=0.24,yend=0.2,color="red")+
  theme(legend.position="none");a

v2<-ggplot(out,aes(exp_coef))+
  geom_vline(xintercept=1,color="darkgray")+
  annotate("segment",x=0.25,xend=1.15,y=0.09,yend=0.19,color="red")+
  annotate("segment",x=0.5,xend=1.15,y=0.247,yend=0.347,color="red")+
  annotate("segment",x=-1,xend=0.25,y=0.09,yend=0.09,linetype="dotted",color="red")+
  annotate("segment",x=0.25,xend=0.25,y=0.09,yend=-1,linetype="dotted",color="red")+
  annotate("segment",x=-1,xend=0.5,y=0.247,yend=0.247,linetype="dotted",color="red")+
  annotate("segment",x=0.5,xend=0.5,y=0.247,yend=-1,linetype="dotted",color="red")+
  stat_ecdf(geom = "step",size=1)+
  theme_classic(base_size=8)+
  xlab("Hazard Ratio")+
  ylab("Cumulative Proportion")+
  coord_cartesian(xlim=c(0,4),ylim=c(0,0.9))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  scale_x_continuous(breaks=c(0,0.25,0.5,1,2,3,4))+
  annotate("text", x=1.2,y=0.35,label="24.7% at 2x WT tolerance",angle=0,size=2.5,hjust=0)+
  annotate("text", x=1.2,y=0.19,label="9.1% at 4x WT tolerance",angle=0,size=2.5,hjust=0)+
  geom_rug(data=haz_dist,aes(y=0,x=exp_coef,alpha=0.5),length=unit(0.1,'cm'),sides="b")+
  annotate("text",x=0.95,y=0.75,angle=90,label="WT Hazard Ratio",size=2.5,color="darkgray",vjust=0)+
  theme(legend.position="none",
        axis.text.x=element_text(angle=90,vjust=0.5));v2



####### ITS2 Analyses ###########
### Import data
ITS2_profiles_per_colony <- read.csv("ITS2_profiles_per_colony.csv")

##### Figure 6A - ITS2 type profile by colony ####

#Color type profiles
matrix_profiles <- c("C31.C17d.C21.C31.9.C31.5.C17i.C31.10.C21ac.C17"="#8ef6fa",
                       "C31.C17d.C21.C31a.C21ac.C17i.C17"="#27a4a8",
                       "D1.D4.D6.D17d.D1r.D17e.D17c"="#e9b8fc",
                       "D1.D4.D6.D1ab.D17d.D17j"="#de98fa",
                       "D1.D4.D1ab.D6.D4d"="#ae41d9",
                       "D1.D4.D6.D1ab.D3h"="#793594",
                       "D1.D6.D4.D1r"="#4f1d63")

# Plot
F6A <- ggplot(ITS2_profiles_per_colony, aes(x=colony, y=colony_profile_prop))+
  geom_col(aes(fill = Type_profile), colour="grey", size=0.005)+
  labs(x="Colony", y="ITS2 Type Profile Proportion")+
  scale_y_continuous(labels=function(Prop)Prop)+
  guides(color = guide_legend(override.aes = list(size=1.5)), fill=guide_legend(ncol=1, title.theme = element_text(angle = 90)))+
  theme_classic(base_size = 8)+
  scale_fill_manual(values=matrix_profiles, name="")+
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1),axis.ticks.x=element_blank(),legend.text = element_text(size=6),legend.position = "bottom",legend.justification = "left",legend.key.size = unit(0.25, "cm")); F6A


##### Figure 6B - Plot proportion of larvae alive in ambient and heated treatments on day 8 by direction of symbiosis ####
## Subset survival data to remove heated treatment, WT larvae, while only keeping the final time point
survival_day8 <- survival %>% filter(Parental_cross!="WT", timepoint=="8")

## Add new columns to survival_day8 dataframe to include symbiont data by parent
survival_day8$Mot_symbio <- survival_day8$Mother
survival_day8$Fat_symbio <- survival_day8$Father

survival_day8 <- survival_day8 %>%      # Rename colonies with symbiont hosted
  mutate(Mot_symbio=case_when(Mother=="1"~"Mother(Clad)",                 
                              Mother=="2"~"Mother(Duru)",
                              Mother=="3"~"Mother(Duru)",
                              Mother=="4"~"Mother(Clad)",
                              Mother=="5"~"Mother(Duru)",
                              Mother=="6"~"Mother(Duru)"))

survival_day8 <- survival_day8 %>%      # Rename colonies with symbiont hosted
  mutate(Fat_symbio=case_when(Father=="A"~"Father(Clad)",                 
                              Father=="B"~"Father(Duru)",
                              Father=="C"~"Father(Duru)",
                              Father=="D"~"Father(Clad)",
                              Father=="E"~"Father(Duru)",
                              Father=="F"~"Father(Duru)"))

# Concatenate Mot_symbio and Fat_symbio into one column to indicate direction of symbiosis
survival_day8$Cross_symbios <- str_c(survival_day8$Mot_symbio," - ",survival_day8$Fat_symbio)

# Determine proportion of larvae still alive in each larval family on day 8
survival_day8_prop <- survival_day8 %>% group_by(treatment, Parental_cross) %>%
                             summarise(dead_larvae=sum(mortality),
                                       n=n,
                                       Cross_symbios=Cross_symbios) %>%
                             distinct() %>%
                             mutate(prop=(n-dead_larvae)/n)

## Plot jitter plot
F6B <- ggplot(survival_day8_prop, aes(x=Cross_symbios, y=prop, fill=treatment))+
  geom_boxplot()+
  scale_fill_manual(values=c("#03cafc", "#f24455")) +
  labs(x="Direction of Symbiosis", y="Proportion of Larval Survival",fill="Treatment")+
  theme_classic(base_size = 8) +
  theme(axis.title.x = element_text(margin = unit(c(1, 0, 0, -4), "mm")), 
        legend.position = c(0.9,0.9),
        legend.title = element_text(size=5), 
        legend.text = element_text(size=4),
        legend.key.size = unit(0.3, 'cm'),
        axis.text.x=element_text(angle=45,hjust=1)); F6B

### Make Figure 6
quartz(w=5.701961, h=2.933333)
plot_grid(F6A,F6B, rel_widths = c(1,1),labels="AUTO",label_size=10)

# Test whether there are differences in the proportion of larvae alive on day 8 depending on direction of symbiosis and treatment (Model 9)
symbios_mod <- aov(prop~Cross_symbios*treatment, data=survival_day8_prop)
summary(symbios_mod)  
TukeyHSD(symbios_mod)

