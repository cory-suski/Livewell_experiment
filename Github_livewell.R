# Load library ------------------------------------------------------------

library(drc)
library(dplyr)
library(tidyverse)
library(lme4)
library(emmeans)
library(car)
library(MuMIn)
library(cowplot)
library(ggborderline)
library(ggprism)
library(openxlsx)


# load in data ------------------------------------------------------------

cortisol<- read_csv("cortisol_new.csv") %>%  mutate(treatment= case_when(treatment=="control1" ~"resting",
                                                                         treatment=="control2" ~"ambient",
                                                                         treatment=="control3" ~"exercise",
                                                                         treatment=="treatment1" ~"increasing",
                                                                         treatment=="treatment2" ~"oscillating",
                                                                         treatment=="treatment3" ~"ice"))

glucose<- read_csv("Glucose_all.csv", 
                   col_types = cols(`spec_ 1` = col_number(), 
                                    spec_2 = col_number(), mean.abs = col_number(), 
                                    cv = col_number(), `glucose.conc (mmol/L)` = col_number())) %>% 
  mutate(treatment= case_when(treatment=="control1" ~"resting",
                              treatment=="control2" ~"ambient",
                              treatment=="treatment3" ~"exercise",
                              treatment=="treatment1" ~"increasing",
                              treatment=="treatment2" ~"oscillating",
                              treatment=="treatment4" ~"ice"))
hemoglobin<- read_csv("Hemoglobin_all.csv", 
                      col_types = cols(od1 = col_number(), 
                                       od2 = col_number())) %>% mutate(treatment= case_when(treatment=="control1" ~"resting",
                                                                                            treatment=="control2" ~"ambient",
                                                                                            treatment=="treatment3" ~"exercise",
                                                                                            treatment=="treatment1" ~"increasing",
                                                                                            treatment=="treatment2" ~"oscillating",
                                                                                            treatment=="treatment4" ~"ice"))
lactate<- read_csv("Lactate_all.csv", 
                   col_types = cols(spec_1 = col_number(), 
                                    spec_2 = col_number(), mean.abs = col_number(), 
                                    cv = col_number(), `lactate (mmol/L)` = col_number()))%>% 
  rename(sample = "sample_id") %>% mutate(treatment= case_when(treatment=="control1" ~"resting",
                                                               treatment=="control2" ~"ambient",
                                                               treatment=="treatment3" ~"exercise",
                                                               treatment=="treatment1" ~"increasing",
                                                               treatment=="treatment2" ~"oscillating",
                                                               treatment=="treatment4" ~"ice"))
livewellexp <- read_csv("livewell_exp_data.csv", 
                        col_types = cols(start.t = col_time(format = "%H:%M:%S"), 
                                         end.t = col_time(format = "%H:%M:%S"))) %>% 
  rename(sample = "sample_id") 
glucose<-glucose %>% rename(sample = "sample_id")
vetmed<- read_csv("VetmedAll.csv")%>% 
  rename(sample = "sample_id") %>% mutate(treatment= case_when(treatment=="control1" ~"resting",
                                                               treatment=="control2" ~"ambient",
                                                               treatment=="treatment3" ~"exercise",
                                                               treatment=="treatment1" ~"increasing",
                                                               treatment=="treatment2" ~"oscillating",
                                                               treatment=="treatment4" ~"ice"))

# Fish size ---------------------------------------------------------------
lengthaov<-aov(length ~ treatment, data= livewellexp)
anova(lengthaov)
Anova(lengthaov)

weightaov<-aov(weight ~ treatment, data= livewellexp)
anova(weightaov)
Anova(weightaov)


##Cortisol ----------------------------------------------------------------
#cortisol is in pg/L
cortisol_exp <-merge(cortisol,livewellexp, by=intersect(x="sample", y="sample"))
cortisolfish<-cortisol_exp %>% select(sample,trueconcen,chase.t,livewell,start.t,end.t,length,weight,hemat.prop1,hemat.prop2,ramp1,ramp2,treatment)
cortisolfish<-cortisolfish[!(cortisolfish$sample %in% c('3H4','5H4')),]
cortisolfish<-cortisolfish%>% mutate(nanocort=(trueconcen*.001))
###Plot---------------------------------------------------

plotrename<- as_labeller(c("resting" ="Resting","increasing" ="Increasing","exercise" ="Exercise",
                           "ambient" ="Ambient","oscillating" ="Oscillating","ice" ="Ice"))
cortisolfish$treatment <- factor(cortisolfish$treatment, levels=c("resting","exercise" ,"ambient","increasing" ,
                                                                  "oscillating" ,"ice" ))
cortbar<-ggplot(data=cortisolfish, aes(x=treatment, y=nanocort, fill=treatment)) + geom_boxplot() +
  scale_fill_manual(name= "Zone", values = c("#3B3B3B", "#3B3B3B","#3B3B3B","#3B3B3B","#3B3B3B","#3B3B3B"))+
  labs(x=NULL, y = "Cortisol (pg/L)")+
  scale_x_discrete(labels=plotrename)+
  theme(
    legend.position = 'none',
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 12,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(angle = 25,vjust = 1, hjust=1,size = 12, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())
cortbar
#subset by treatments to get avg concentration
controlcort<-subset(cortisolfish, treatment == 'control1')
mean(controlcort$trueconcen)
control2cort<-subset(cortisolfish, treatment == 'control2')
mean(control2cort$trueconcen)
trt1cort<-subset(cortisolfish, treatment == 'treatment1')
mean(trt1cort$trueconcen)
trt2cort<-subset(cortisolfish, treatment == 'treatment2')
mean(trt2cort$trueconcen)
trt3cort<-subset(cortisolfish, treatment == 'control3')
mean(trt3cort$trueconcen)
trt4cort<-subset(cortisolfish, treatment == 'treatment3')
mean(trt4cort$trueconcen)
###analyses ----------------------------------------------------------------
#oneway anova
cortaov<-aov(nanocort~treatment,data=cortisolfish)
anova(cortaov)
Anova(cortaov)
par(mfrow=c(2,2))
plot(cortaov)

## Section 2.3.6.1 Zuur et al Mixed Effects Model Book - Model Validation ##

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(cortaov, add.smooth = FALSE, which = 1)
E <- resid(cortaov)
hist(E, xlab = "Residuals", main = "")
plot(cortaov$treatment, E, xlab = "treatment",
     ylab = "Residuals")
par(op)

summary(cortaov)
TukeyHSD(cortaov)
write.csv(Anova(cortaov), file = 'Cortisol.Anova.csv')
write.csv(anova(cortaov), file = 'Cortisol.table.csv')
#twoway anova

##Glucose -----------------------------------------------------------------
#glucose is in mmol/L
#remove individuals with high cv
glucose<-glucose %>% filter(cv < 12) 
glucosefish<-merge(livewellexp,glucose, by=intersect(x="sample", y="sample"))
glucosefish<-glucosefish %>% rename(glucose.conc="glucose.conc (mmol/L)")
glucosefish<-glucosefish%>% select(sample,chase.t,livewell,start.t,end.t,length,weight,hemat.prop1,hemat.prop2,ramp1,ramp2,glucose.conc,treatment)
glucosefish<-glucosefish[!(glucosefish$sample %in% c('3H4','5H4')),]
#determine if minutes of holding time influences treatments
glucosefish.t<-glucosefish%>% mutate(hold.t=(end.t-start.t))
glucosefish.t<-glucosefish.t[!(glucosefish.t$sample %in% c('1H1','1H2','1H3','1H4','1H5','1H6','1H7','1H8','1H9','1H10')),]

###Plot to investigate ---------------------------------------------------
plotrename<- as_labeller(c("resting" ="Resting","increasing" ="Increasing","exercise" ="Exercise",
                           "ambient" ="Ambient","oscillating" ="Oscillating","ice" ="Ice"))
glucosefish$treatment <- factor(glucosefish$treatment, levels=c("resting","exercise" ,"ambient","increasing" ,
                                                                "oscillating" ,"ice" ))
#asterisk different from treatment, + different than resting control
glucbar<-ggplot(data=glucosefish, aes(x=treatment, y=glucose.conc, fill=treatment)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name= "treatment", values = c("lightgray", "lightgray","lightgray","lightgray","lightgray","lightgray"))+
  labs(x=NULL, y = "Glucose (mmol/L)")+
  scale_x_discrete(labels=plotrename)+
  annotate("text", x =4.81, y = 21.5, label = "+", fontface = "bold", colour = "black", size = 12) +
  annotate("text", x =5.13, y = 21, label = "*", fontface = "bold", colour = "black", size = 15) +
  annotate("text", x =4.00, y = 17, label = "*", fontface = "bold", colour = "black", size = 15) +
  theme(
    legend.position='none',
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 18,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length.y =  unit(.50,"cm"),
    axis.line.y.left = element_line(size = 1.00, colour = "black"),
    axis.ticks.y =  element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    axis.line.x =  element_line(size = 1.00, colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x =  element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())+
  geom_point(size=2, position = position_jitter(width = .1))

glucbar

#subset by treatments to get avg concentration
controlgluc<-subset(glucosefish, treatment == 'control1')
mean(controlgluc$glucose.conc)
control2gluc<-subset(glucosefish, treatment == 'control2')
mean(control2gluc$glucose.conc)
trt1gluc<-subset(glucosefish, treatment == 'treatment1')
mean(trt1gluc$glucose.conc)
trt2gluc<-subset(glucosefish, treatment == 'treatment2')
mean(trt2gluc$glucose.conc)
trt3gluc<-subset(glucosefish, treatment == 'control3')
mean(trt3gluc$glucose.conc)
trt4gluc<-subset(glucosefish, treatment == 'treatment3')
mean(trt4gluc$glucose.conc)
###analyses ----------------------------------------------------------------

#check for holding time influence

gluc.t<-aov(glucose.conc~treatment+hold.t, data=glucosefish.t)
anova(gluc.t)
Anova(gluc.t)


glucaov<-aov(glucose.conc~treatment, data=glucosefish)

anova(glucaov)
Anova(glucaov)
par(mfrow=c(2,2))
plot(glucaov)

## Section 2.3.6.1 Zuur et al Mixed Effects Model Book - Model Validation ##

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(glucaov, add.smooth = FALSE, which = 1)
E <- resid(glucaov)
hist(E, xlab = "Residuals", main = "")
plot(glucaov$treatment, E, xlab = "treatment",
     ylab = "Residuals")
par(op)

TukeyHSD(glucaov)
write.csv(Anova(glucaov), file = 'Glucose.Anova.csv')
write.csv(anova(glucaov), file = 'Glucose.table.csv')

###analyses ----------------------------------------------------------------
hemogaov<-aov(hemog.conc~treatment, data=hemoglobinfish)
anova(hemogaov)
Anova(hemogaov)
par(mfrow = c(2, 2))
plot(lm(hemogaov))

TukeyHSD(hemogaov)
#two way anova

## lactate -----------------------------------------------------------------
#lactate is in mmol/l
#remove individuals with high cv
lactate<-lactate %>% filter(cv < 12) 
lactatefish<-merge(livewellexp,lactate, by=intersect(x="sample", y="sample"))
lactatefish<-lactatefish %>% rename(lactate.conc="lactate (mmol/L)")
lactatefish<-lactatefish%>% select(sample,chase.t,livewell,start.t,end.t,length,weight,hemat.prop1,hemat.prop2,ramp1,ramp2,lactate.conc,treatment)
lactatefish<-lactatefish[!(lactatefish$sample %in% c('3H4','5H4')),]
lactatefish.t<-lactatefish%>% mutate(hold.t=(end.t-start.t))
lactatefish.t<-lactatefish.t[!(lactatefish.t$sample %in% c('1H1','1H2','1H3','1H4','1H5','1H6','1H7','1H8','1H9','1H10')),]
###Plot to investigate ---------------------------------------------------

plotrename<- as_labeller(c("resting" ="Resting","increasing" ="Increasing","exercise" ="Exercise",
                           "ambient" ="Ambient","oscillating" ="Oscillating","ice" ="Ice"))
lactatefish$treatment <- factor(lactatefish$treatment, levels=c("resting","exercise" ,"ambient","increasing" ,
                                                                "oscillating" ,"ice" ))
lactbar<-ggplot(data=lactatefish, aes(x=treatment, y=lactate.conc, fill=treatment)) + geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name= "treatment", values = c("lightgray", "lightgray","lightgray","lightgray","lightgray","lightgray"))+
  labs(x=NULL, y = "Lactate (mmol/L)")+
  scale_x_discrete(labels=plotrename, expand = c(0, .5,0,0))+
  annotate("text", x =2.00, y =10, label = "+", fontface = "bold", colour = "black", size = 12) +
  annotate("text", x =5.0, y = 12, label = "+", fontface = "bold", colour = "black", size = 12) +
  theme(
    legend.position = 'none',
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 18,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(angle = 25,vjust = 1, hjust=1,size = 12, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())+
  geom_point(size=2, position = position_jitter(width = .1))
lactbar


#subset by treatments to get avg concentration
controllact<-subset(lactatefish, treatment == 'control1')
mean(controllact$lactate.conc)
control2lact<-subset(lactatefish, treatment == 'control2')
mean(control2lact$lactate.conc)
trt1lact<-subset(lactatefish, treatment == 'treatment1')
mean(trt1lact$lactate.conc)
trt2lact<-subset(lactatefish, treatment == 'treatment2')
mean(trt2lact$lactate.conc)
trt3lact<-subset(lactatefish, treatment == 'control3')
mean(trt3lact$lactate.conc)
trt4lact<-subset(lactatefish, treatment == 'treatment3')
mean(trt4lact$lactate.conc)
###analyses ----------------------------------------------------------------
#check for holding influence
lact.t<-aov(lactate.conc~treatment + hold.t, data=lactatefish.t)
anova(lact.t)
Anova(lact.t)


lactaov<-aov(lactate.conc~treatment, data=lactatefish)
anova(lactaov)
Anova(lactaov)
par(mfrow=c(2,2))
plot(lactaov)

## Section 2.3.6.1 Zuur et al Mixed Effects Model Book - Model Validation ##

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(lactaov, add.smooth = FALSE, which = 1)
E <- resid(lactaov)
hist(E, xlab = "Residuals", main = "")
plot(lactaov$treatment, E, xlab = "treatment",
     ylab = "Residuals")
par(op)

TukeyHSD(lactaov)
write.csv(Anova(lactaov), file = 'Lactate.Anova.csv')
write.csv(anova(lactaov), file = 'Lactate.table.csv')
##vetmed ------------------------------------------------------------------
#this df holds all analytes measured by the diagnostic laboratory at the veternary medicine school
livewellexp<-livewellexp%>% select(sample,chase.t,livewell,start.t,end.t,length,weight,hemat.prop1,hemat.prop2,ramp1,ramp2)
vetmedfish<-merge(livewellexp,vetmed, by=intersect(x="sample", y="sample"))
vetmedfish<-vetmedfish %>% rename(K.mmol.L="K.mmol/L")
vetmedfish<-vetmedfish[!(vetmedfish$sample %in% c('3H4','5H4')),]
vetmedfish.t<-vetmedfish%>% mutate(hold.t=(end.t-start.t))
vetmedfish.t<-vetmedfish.t[!(vetmedfish.t$sample %in% c('1H1','1H2','1H3','1H4','1H5','1H6','1H7','1H8','1H9','1H10')),]
###Plot to investigate ---------------------------------------------------

####bicarb ------------------------------------------------------------------

plotrename<- as_labeller(c("resting" ="Resting","increasing" ="Increasing","exercise" ="Exercise",
                           "ambient" ="Ambient","oscillating" ="Oscillating","ice" ="Ice"))
vetmedfish$treatment <- factor(vetmedfish$treatment, levels=c("resting","exercise" ,"ambient","increasing" ,
                                                              "oscillating" ,"ice" ))
bicarbbar<-ggplot(data=vetmedfish, aes(x=treatment, y=bicarb.mmol.L)) + geom_bar(stat='summary') +
  geom_errorbar(stat='summary', width=.2)+
  labs(x=NULL, y = "Bicarbonate (mmol/L)")+
  scale_x_discrete(labels=plotrename)+
  theme(
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(colour = "black", size = 15),
    legend.title = element_text(face = "bold", colour = "black", size = 18),
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 12,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(angle = 25,vjust = 1, hjust=1,size = 12, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())
bicarbbar


#subset by treatments to get avg concentration
controlbicarb<-subset(vetmedfish, treatment == 'control1')
mean(controlbicarb$bicarb.mmol.L)
control2bicarb<-subset(vetmedfish, treatment == 'control2')
mean(control2bicarb$bicarb.mmol.L)
trt1bicarb<-subset(vetmedfish, treatment == 'treatment1')
mean(trt1bicarb$bicarb.mmol.L)
trt2bicarb<-subset(vetmedfish, treatment == 'treatment2')
mean(trt2bicarb$bicarb.mmol.L)
trt3bicarb<-subset(vetmedfish, treatment == 'control3')
mean(trt3bicarb$bicarb.mmol.L)
trt4bicarb<-subset(vetmedfish, treatment == 'treatment3')
mean(trt4bicarb$bicarb.mmol.L)
###analyses ----------------------------------------------------------------

#holding check
bicarb.t<-aov(bicarb.mmol.L~treatment+hold.t, data=vetmedfish.t)
anova(bicarb.t)
Anova(bicarb.t)

bicarbaov<-aov(bicarb.mmol.L~treatment, data=vetmedfish)
anova(bicarbaov)
Anova(bicarbaov)
par(mfrow = c(2, 2))
plot(lm(bicarbaov))
TukeyHSD(bicarbaov)

write.csv(Anova(bicarbaov), file = 'bicarb.Anova.csv')
write.csv(anova(bicarbaov), file = 'bicarb.table.csv')


####Ca ------------------------------------------------------------------

plotrename<- as_labeller(c("resting" ="Resting","increasing" ="Increasing","exercise" ="Exercise",
                           "ambient" ="Ambient","oscillating" ="Oscillating","ice" ="Ice"))
vetmedfish$treatment <- factor(vetmedfish$treatment, levels=c("resting","exercise" ,"ambient","increasing" ,
                                                              "oscillating" ,"ice" ))
cabar<-ggplot(data=vetmedfish, aes(x=treatment, y=Ca.mg.dL)) + geom_bar(stat='summary') +
  annotate("text", x =3.00, y = 19, label = "*", fontface = "bold", colour = "black", size = 8) +
  
  geom_errorbar(stat='summary', width=.2)+
  labs(x=NULL, y = "Calcium (mg/dL)")+
  scale_x_discrete(labels=plotrename)+
  theme(
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(colour = "black", size = 15),
    legend.title = element_text(face = "bold", colour = "black", size = 18),
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 12,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(angle = 25,vjust = 1, hjust=1,size = 12, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())
cabar
###analyses ----------------------------------------------------------------

CA.t<-aov(Ca.mg.dL~treatment+hold.t, data=vetmedfish.t)
anova(CA.t)
Anova(CA.t)

Caaov<-aov(Ca.mg.dL~treatment, data=vetmedfish)
anova(Caaov)
Anova(Caaov)
par(mfrow = c(2, 2))
plot(lm(Caaov))
TukeyHSD(Caaov)

write.csv(Anova(Caaov), file = 'calcium.Anova.csv')
write.csv(anova(Caaov), file = 'calcium.table.csv')

#subset by treatments to get avg concentration
controlCa<-subset(vetmedfish, treatment == 'control1')
mean(controlCa$Ca.mg.dL)
control2Ca<-subset(vetmedfish, treatment == 'control2')
mean(control2Ca$Ca.mg.dL)
trt1Ca<-subset(vetmedfish, treatment == 'treatment1')
mean(trt1Ca$Ca.mg.dL)
trt2Ca<-subset(vetmedfish, treatment == 'treatment2')
mean(trt2Ca$Ca.mg.dL,na.rm = TRUE)
trt3Ca<-subset(vetmedfish, treatment == 'control3')
mean(trt3Ca$Ca.mg.dL,na.rm = TRUE)
trt4Ca<-subset(vetmedfish, treatment == 'treatment3')
mean(trt4Ca$Ca.mg.dL)
####Cl ------------------------------------------------------------------

plotrename<- as_labeller(c("resting" ="Resting","increasing" ="Increasing","exercise" ="Exercise",
                           "ambient" ="Ambient","oscillating" ="Oscillating","ice" ="Ice"))
vetmedfish$treatment <- factor(vetmedfish$treatment, levels=c("resting","exercise" ,"ambient","increasing" ,
                                                              "oscillating" ,"ice" ))
clbar<-ggplot(data=vetmedfish, aes(x=treatment, y=Cl.mmol.L)) + geom_bar(stat='summary') +
  geom_errorbar(stat='summary', width=.2)+
  labs(x=NULL, y = "Chloride (mmol/L)")+
  scale_x_discrete(labels=plotrename)+
  theme(
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(colour = "black", size = 15),
    legend.title = element_text(face = "bold", colour = "black", size = 18),
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 12,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(angle = 25,vjust = 1, hjust=1,size = 12, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())
clbar

#subset by treatments to get avg concentration
controlCl<-subset(vetmedfish, treatment == 'control1')
mean(controlCl$Cl.mmol.L)
control2Cl<-subset(vetmedfish, treatment == 'control2')
mean(control2Cl$Cl.mmol.L)
trt1Cl<-subset(vetmedfish, treatment == 'treatment1')
mean(trt1Cl$Cl.mmol.L)
trt2Cl<-subset(vetmedfish, treatment == 'treatment2')
mean(trt2Cl$Cl.mmol.L)
trt3Cl<-subset(vetmedfish, treatment == 'control3')
mean(trt3Cl$Cl.mmol.L)
trt4Cl<-subset(vetmedfish, treatment == 'treatment3')
mean(trt4Cl$Cl.mmol.L)
###analyses ----------------------------------------------------------------

Cl.t<-aov(Cl.mmol.L~treatment+hold.t, data=vetmedfish.t)
anova(Cl.t)
Anova(Cl.t)

Claov<-aov(Cl.mmol.L~treatment, data=vetmedfish)
anova(Claov)
Anova(Claov)
par(mfrow = c(2, 2))
plot(lm(Claov))

write.csv(Anova(Claov), file = 'chloride.Anova.csv')
write.csv(anova(Claov), file = 'chloride.table.csv')
#### Cholesterol------------------------------------------------------------------
plotrename<- as_labeller(c("resting" ="Resting","increasing" ="Increasing","exercise" ="Exercise",
                           "ambient" ="Ambient","oscillating" ="Oscillating","ice" ="Ice"))
vetmedfish$treatment <- factor(vetmedfish$treatment, levels=c("resting","exercise" ,"ambient","increasing" ,
                                                              "oscillating" ,"ice" ))
cholbar<-ggplot(data=vetmedfish, aes(x=treatment, y=t.chol.mg.dL, na.rm=TRUE)) + geom_bar(stat='summary') +
  geom_errorbar(stat='summary', width=.2)+
  labs(x=NULL, y = "Cholesterol (mg/dL)")+
  scale_x_discrete(labels=plotrename)+
  theme(
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(colour = "black", size = 15),
    legend.title = element_text(face = "bold", colour = "black", size = 18),
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 12,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(angle = 25,vjust = 1, hjust=1,size = 12, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())
cholbar


#subset by treatments to get avg concentration
controlChol<-subset(vetmedfish, treatment == 'control1')
mean(controlChol$t.chol.mg.dL,na.rm = TRUE)
control2Chol<-subset(vetmedfish, treatment == 'control2')
mean(control2Chol$t.chol.mg.dL)
trt1Chol<-subset(vetmedfish, treatment == 'treatment1')
mean(trt1Chol$t.chol.mg.dL)
trt2Chol<-subset(vetmedfish, treatment == 'treatment2')
mean(trt2Chol$t.chol.mg.dL,na.rm = TRUE)
trt3Chol<-subset(vetmedfish, treatment == 'control3')
mean(trt3Chol$t.chol.mg.dL,na.rm = TRUE)
trt4Chol<-subset(vetmedfish, treatment == 'treatment3')
mean(trt4Chol$t.chol.mg.dL,na.rm = TRUE)
###analyses ----------------------------------------------------------------
Chol.t<-aov(t.chol.mg.dL~treatment+hold.t, data=vetmedfish.t)
anova(Chol.t)
Anova(Chol.t)


cholaov<-aov(t.chol.mg.dL~treatment, data=vetmedfish)
anova(cholaov)
Anova(cholaov)
par(mfrow = c(2, 2))
plot(lm(cholaov))
TukeyHSD(cholaov)
write.csv(Anova(cholaov), file = 'cholesterol.Anova.csv')
write.csv(anova(cholaov), file = 'cholesterol.table.csv')

####K ------------------------------------------------------------------
plotrename<- as_labeller(c("resting" ="Resting","increasing" ="Increasing","exercise" ="Exercise",
                           "ambient" ="Ambient","oscillating" ="Oscillating","ice" ="Ice"))
vetmedfish$treatment <- factor(vetmedfish$treatment, levels=c("resting","exercise" ,"ambient","increasing" ,
                                                              "oscillating" ,"ice" ))
kbar<-ggplot(data=vetmedfish, aes(x=treatment, y=K.mmol.L,fill=treatment)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name= "treatment", values = c("lightgray", "lightgray","lightgray","lightgray","lightgray","lightgray"))+
  labs(x=NULL, y = "Potassium (mmol/L)")+
  scale_x_discrete(labels=plotrename)+
  annotate("text", x =2.00, y = 5, label = "+", fontface = "bold", colour = "black", size = 12) +
  annotate("text", x =3, y = 4.5, label = "*", fontface = "bold", colour = "black", size = 15) +
  annotate("text", x =4, y = 5.5, label = "+", fontface = "bold", colour = "black", size = 12) +
  annotate("text", x =5, y = 4.5, label = "+", fontface = "bold", colour = "black", size = 12) +
  theme(
    legend.position='none',
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 18,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length.y =  unit(.50,"cm"),
    axis.line.y.left = element_line(size = 1.00, colour = "black"),
    axis.ticks.y =  element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    axis.line.x =  element_line(size = 1.00, colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x =  element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())+
  geom_point(size=2, position = position_jitter(width = .1))

kbar

#subset by treatments to get avg concentration
controlK<-subset(vetmedfish, treatment == 'control1')
mean(controlK$K.mmol.L,na.rm = TRUE)
control2K<-subset(vetmedfish, treatment == 'control2')
mean(control2K$K.mmol.L)
trt1K<-subset(vetmedfish, treatment == 'treatment1')
mean(trt1K$K.mmol.L)
trt2K<-subset(vetmedfish, treatment == 'treatment2')
mean(trt2K$K.mmol.L,na.rm = TRUE)
trt3K<-subset(vetmedfish, treatment == 'control3')
mean(trt3K$K.mmol.L,na.rm = TRUE)
trt4K<-subset(vetmedfish, treatment == 'treatment3')
mean(trt4K$K.mmol.L,na.rm = TRUE)
###analyses ----------------------------------------------------------------
Ka.t<-aov(K.mmol.L~treatment+hold.t, data=vetmedfish.t)
anova(Ka.t)
Anova(Ka.t)

Kaov<-aov(K.mmol.L~treatment, data=vetmedfish)
anova(Kaov)
Anova(Kaov)
par(mfrow=c(2,2))
plot(Kaov)

## Section 2.3.6.1 Zuur et al Mixed Effects Model Book - Model Validation ##

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(Kaov, add.smooth = FALSE, which = 1)
E <- resid(Kaov)
hist(E, xlab = "Residuals", main = "")
plot(Kaov$treatment, E, xlab = "treatment",
     ylab = "Residuals")
par(op)

TukeyHSD(Kaov)

write.csv(Anova(Kaov), file = 'potassium.Anova.csv')
write.csv(anova(Kaov), file = 'potassium.table.csv')
####Na ------------------------------------------------------------------
plotrename<- as_labeller(c("resting" ="Resting","increasing" ="Increasing","exercise" ="Exercise",
                           "ambient" ="Ambient","oscillating" ="Oscillating","ice" ="Ice"))
vetmedfish$treatment <- factor(vetmedfish$treatment, levels=c("resting","exercise" ,"ambient","increasing" ,
                                                              "oscillating" ,"ice" ))
Nabar<-ggplot(data=vetmedfish, aes(x=treatment, y=Na.mmol.L, fill=treatment)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name= "treatment", values = c("lightgray", "lightgray","lightgray","lightgray","lightgray","lightgray"))+
  labs(x=NULL, y = "Sodium (mmol/L)")+
  scale_x_discrete(labels=plotrename)+
  annotate("text", x =2.00, y = 180, label = "+", fontface = "bold", colour = "black", size = 12) +
  annotate("text", x =2.81, y = 160, label = "+", fontface = "bold", colour = "black", size = 12) +
  annotate("text", x =3.13, y = 159, label = "*", fontface = "bold", colour = "black", size = 15) +
  annotate("text", x =4, y = 159, label = "*", fontface = "bold", colour = "black", size = 15) +
  annotate("text", x =4.81, y = 160, label = "+", fontface = "bold", colour = "black", size = 12) +
  annotate("text", x =5.13, y = 159, label = "*", fontface = "bold", colour = "black", size = 15) +
  annotate("text", x =5.81, y = 160, label = "+", fontface = "bold", colour = "black", size = 12) +
  annotate("text", x =6.13, y = 159, label = "*", fontface = "bold", colour = "black", size = 15) +
  theme(
    legend.position='none',
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 18,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length.y =  unit(.50,"cm"),
    axis.line.y.left = element_line(size = 1.00, colour = "black"),
    axis.ticks.y =  element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    axis.line.x =  element_line(size = 1.00, colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x =  element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())+
  geom_point(size=2, position = position_jitter(width = .1))


Nabar

#subset by treatments to get avg concentration
controlNa<-subset(vetmedfish, treatment == 'control1')
mean(controlNa$Na.mmol.L,na.rm = TRUE)
control2Na<-subset(vetmedfish, treatment == 'control2')
mean(control2Na$Na.mmol.L)
trt1Na<-subset(vetmedfish, treatment == 'treatment1')
mean(trt1Na$Na.mmol.L)
trt2Na<-subset(vetmedfish, treatment == 'treatment2')
mean(trt2Na$Na.mmol.L,na.rm = TRUE)
trt3Na<-subset(vetmedfish, treatment == 'control3')
mean(trt3Na$Na.mmol.L,na.rm = TRUE)
trt4Na<-subset(vetmedfish, treatment == 'treatment3')
mean(trt4Na$Na.mmol.L,na.rm = TRUE)

###analyses ----------------------------------------------------------------
Na.t<-aov(Na.mmol.L~treatment+hold.t, data=vetmedfish.t)
anova(Na.t)
Anova(Na.t)

Naaov<-aov(Na.mmol.L~treatment, data=vetmedfish)
anova(Naaov)
Anova(Naaov)
par(mfrow=c(2,2))
plot(Naaov)

## Section 2.3.6.1 Zuur et al Mixed Effects Model Book - Model Validation ##

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(Naaov, add.smooth = FALSE, which = 1)
E <- resid(Naaov)
hist(E, xlab = "Residuals", main = "")
plot(Kaov$treatment, E, xlab = "treatment",
     ylab = "Residuals")
par(op)
TukeyHSD(Naaov)

write.csv(Anova(Naaov), file = 'sodium.Anova.csv')
write.csv(anova(Naaov), file = 'sodium.table.csv')
####protein ------------------------------------------------------------------
plotrename<- as_labeller(c("resting" ="Resting","increasing" ="Increasing","exercise" ="Exercise",
                           "ambient" ="Ambient","oscillating" ="Oscillating","ice" ="Ice"))
vetmedfish$treatment <- factor(vetmedfish$treatment, levels=c("resting","exercise" ,"ambient","increasing" ,
                                                              "oscillating" ,"ice" ))
protbar<-ggplot(data=vetmedfish, aes(x=treatment, y=t.protien.g.dL, fill=treatment))+ 
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(name= "treatment", values = c("lightgray", "lightgray","lightgray","lightgray","lightgray","lightgray"))+
  labs(x=NULL, y = "Protein (g/dL)")+
  scale_x_discrete(labels=plotrename)+
  annotate("text", x =3, y = 5.6, label = "*", fontface = "bold", colour = "black", size = 15) +
  annotate("text", x =4, y = 5.6, label = "*", fontface = "bold", colour = "black", size = 15) +
  annotate("text", x =5, y = 5.6, label = "*", fontface = "bold", colour = "black", size = 15) +
  annotate("text", x =6, y = 5.6, label = "*", fontface = "bold", colour = "black", size = 15) +
  theme(
    legend.position = 'none',
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 18,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(angle = 25,vjust = 1, hjust=1,size = 12, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())+
  geom_point(size=2, position = position_jitter(width = .1))

protbar

###analyses ----------------------------------------------------------------
prot.t<-aov(t.protien.g.dL~treatment+hold.t, data=vetmedfish.t)
anova(prot.t)
Anova(prot.t)

proteinaov<-aov(t.protien.g.dL~treatment, data=vetmedfish)
anova(proteinaov)
Anova(proteinaov)
par(mfrow = c(2, 2))
plot(lm(proteinaov))
## Section 2.3.6.1 Zuur et al Mixed Effects Model Book - Model Validation ##

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(proteinaov, add.smooth = FALSE, which = 1)
E <- resid(proteinaov)
hist(E, xlab = "Residuals", main = "")
plot(proteinaov$treatment, E, xlab = "treatment",
     ylab = "Residuals")
par(op)

TukeyHSD(proteinaov)

write.csv(Anova(proteinaov), file = 'protien.Anova.csv')
write.csv(anova(proteinaov), file = 'protien.table.csv')


#Data from experiment day -----------------------------------------------------
livewellexp<-vetmedfish%>% mutate(avg.hemat=(hemat.prop1+hemat.prop2)/2)
livewellexp<-livewellexp%>% mutate(rampscore1=(ramp1/5))
livewellexp<-livewellexp%>% mutate(rampscore2=(ramp2/5))
livewellexp<-livewellexp%>% mutate(rampdiff=(rampscore2-rampscore1))
livewellexp.t<-livewellexp%>% mutate(hold.t=(end.t-start.t))
livewellexp.t<-livewellexp.t[!(livewellexp.t$sample %in% c('1H1','1H2','1H3','1H4','1H5','1H6','1H7','1H8','1H9','1H10')),]
livewellexp <-livewellexp%>% select(sample,chase.t,livewell,start.t,end.t,length,weight,avg.hemat,rampscore1, rampscore2,rampdiff,treatment)


##plot to investigate -----------------------------------------------------

####hemat ------------------------------------------------------------------
plotrename<- as_labeller(c("resting" ="Resting","increasing" ="Increasing","exercise" ="Exercise",
                           "ambient" ="Ambient","oscillating" ="Oscillating","ice" ="Ice"))
livewellexp$treatment <- factor(livewellexp$treatment, levels=c("resting","exercise" ,"ambient","increasing" ,
                                                                "oscillating" ,"ice" ))
hematbar<-ggplot(data=livewellexp, aes(x=treatment, y=avg.hemat, fill=treatment)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name= "treatment", values = c("lightgray", "lightgray","lightgray","lightgray","lightgray","lightgray"))+
  labs(x=NULL, y = "Hematocrit (%)")+
  scale_x_discrete(labels=plotrename)+
  annotate("text", x =6.00, y = .5, label = "*", fontface = "bold", colour = "black", size = 15) +
  
  theme(
    legend.position = 'none',
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 18,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(angle = 25,vjust = 1, hjust=1,size = 12, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())+
  geom_point(size=2, position = position_jitter(width = .1))

hematbar

#subset by treatments to get avg concentration
controlhemat<-subset(livewellexp, treatment == 'control1')
mean(controlhemat$avg.hemat,na.rm = TRUE)
control2hemat<-subset(livewellexp, treatment == 'control2')
mean(control2hemat$avg.hemat)
trt1hemat<-subset(livewellexp, treatment == 'treatment1')
mean(trt1hemat$avg.hemat)
trt2hemat<-subset(livewellexp, treatment == 'treatment2')
mean(trt2hemat$avg.hemat,na.rm = TRUE)
trt3hemat<-subset(livewellexp, treatment == 'control3')
mean(trt3hemat$avg.hemat,na.rm = TRUE)
trt4hemat<-subset(livewellexp, treatment == 'treatment3')
mean(trt4hemat$avg.hemat,na.rm = TRUE)
###analyses ----------------------------------------------------------------
#check for holding time influence
hemataov<-aov(avg.hemat~treatment+hold.t, data=livewellexp.t)

hemataov<-aov(avg.hemat~treatment, data=livewellexp)
anova(hemataov)
Anova(hemataov)
par(mfrow = c(2, 2))
plot(lm(hemataov))

## Section 2.3.6.1 Zuur et al Mixed Effects Model Book - Model Validation ##

op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(hemataov, add.smooth = FALSE, which = 1)
E <- resid(hemataov)
hist(E, xlab = "Residuals", main = "")
plot(hemataov$treatment, E, xlab = "treatment",
     ylab = "Residuals")
par(op)

TukeyHSD(hemataov)
write.csv(Anova(hemataov), file = 'hemat.Anova.csv')
write.csv(anova(hemataov), file = 'hemat.table.csv')
####ramp ------------------------------------------------------------------
plotrename<- as_labeller(c("resting" ="Resting","increasing" ="Increasing","exercise" ="Exercise",
                           "ambient" ="Ambient","oscillating" ="Oscillating","ice" ="Ice"))
livewellexp$treatment <- factor(livewellexp$treatment, levels=c("resting","exercise" ,"ambient","increasing" ,
                                                                "oscillating" ,"ice" ))
rampbar<-ggplot(data=rampdf, aes(x=treatment, y=rampdiff, fill=treatment)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name= "treatment", values = c("lightgray", "lightgray","lightgray","lightgray","lightgray","lightgray"))+
  labs(x=NULL, y = "Change in Reflex Impairment (%)")+
  scale_x_discrete(labels=plotrename)+
  theme(
    legend.position = 'none',
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 25,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 20),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(size = 20, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())+
  geom_point(size=3.5, position = position_jitter(width = .1), alpha=.5)

rampbar

#final and initial ramp
rampdfintitial<-livewellexp %>% select(rampscore1, sample, treatment, rampdiff) %>% mutate(initial_final="initial") %>% rename(rampscore = "rampscore1")
rampdffinal<-livewellexp %>% select(rampscore2, sample, treatment, rampdiff) %>% mutate(initial_final="final") %>% rename(rampscore = "rampscore2")
rampall<-rbind(rampdffinal,rampdfintitial)
plotrename<- as_labeller(c("resting" ="Resting","increasing" ="Increasing","exercise" ="Exercise",
                           "ambient" ="Ambient","oscillating" ="Oscillating","ice" ="Ice"))
rampall$initial_final <- factor(rampall$initial_final, levels=c("initial","final"))
rampfinalinitialplot<-ggplot(data=rampall, aes(x=treatment, y=rampscore, fill=initial_final, group_by = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name= "RAMP Score", values = c("black", "lightgray","black","lightgray","black","lightgray","black", "lightgray","black","lightgray","black","lightgray"))+
  labs(x=NULL, y = "Reflex Impairment Score (%)")+
  scale_x_discrete(labels=plotrename)+
  scale_y_continuous(limits = NULL, breaks = c(0.0,0.2,0.4))+
  theme(
    legend.position = 'right',
    legend.text = element_text(size=15),
    legend.title = element_text(size = 20),
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 25,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 20),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(size = 20, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())
rampfinalinitialplot

###analyses ----------------------------------------------------------------
#Add ramp fial
rampdfintitial<-livewellexp %>% select(rampscore1, sample, treatment, rampdiff) %>% mutate(initial_final="initial") %>% rename(rampscore = "rampscore1")

rampdf<-livewellexp %>% select(rampscore1,rampscore2, sample, treatment, rampdiff)

#check for influenc of hold time
rampdf.t<-livewellexp.t %>% select(rampscore1,rampscore2, sample, treatment, rampdiff, hold.t)

rampdf.t$transformed_response <- (rampdf.t$rampdiff + 0.4) / 0.8
ramplinear_mode.t <- glm(transformed_response ~ treatment * hold.t, family= binomial, data = rampdf.t)
anova(ramplinear_mode.t)
Anova(ramplinear_mode.t)

# Apply linear transformation to the response variable
rampdf$transformed_response <- (rampdf$rampdiff + 0.4) / 0.8
ramplinear_model <- glm(transformed_response ~ treatment, family= binomial, data = rampdf)

summary(ramplinear_model)
anova(ramplinear_model)
Anova(ramplinear_model)
ramplinear_modelcontrast <-emmeans(ramplinear_model, ~ treatment)
pwpp(ramplinear_modelcontrast)
ramplinear_modelcontrast<-pairs(ramplinear_modelcontrast)
ramplinear_modelcontrast
r_squared <- r.squaredGLMM(ramplinear_model)
r_squared

par(mfrow = c(2, 2))
plot(lm(ramplinear_model))

write.csv(Anova(ramplinear_model), file = 'ramplinear_model.Anova.csv')
write.csv(anova(ramplinear_model), file = 'ramplinear_model.csv')
write.xlsx(ramplinear_modelcontrast, 'C:/Users/Allison/Box/Hay_Allison/Data/Livewell Exp/R Directory') 

#ramp model for  final
rampfinalmodel <- glm(rampscore2 ~ treatment, family= binomial, data = rampdf)

summary(rampfinalmodel)
anova(rampfinalmodel)
Anova(rampfinalmodel)
rampfinalmodelcontrast <-emmeans(rampfinalmodel, ~ treatment)
pwpp(rampfinalmodelcontrast)
rampfinalmodelcontrast<-pairs(rampfinalmodelcontrast)
rampfinalmodelcontrast
r_squared <- r.squaredGLMM(rampfinalmodel)
r_squared

par(mfrow = c(2, 2))
plot(lm(rampfinalmodel))

write.csv(Anova(rampfinalmodel), file = 'rampfinalmodel.Anova.csv')
write.csv(anova(rampfinalmodel), file = 'rampfinalmodel.output.csv')
write.xlsx(ramplinear_modelcontrast, 'C:/Users/Allison/Box/Hay_Allison/Data/Livewell Exp/R Directory') 
# Combining plots ---------------------------------------------------------

#secondary response plot
plot_grid(kbar, Nabar,glucbar,lactbar,protbar,hematbar, labels=c('A','B','C','D','E','F'), label_size = 18,ncol = 3, align = 'v')


ggsave("bloodgrid.png",
       plot = last_plot(),
       device = NULL,
       path = NULL,
       scale = 1,
       width = 17,
       height = 15,
       units = c("in"),
       dpi = 300,
       limitsize = FALSE,
       bg = NULL)
#tertiary response plot

rampbar
ggsave("ramp.png",
       plot = last_plot(),
       device = NULL,
       path = NULL,
       scale = 1,
       width = 15,
       height = 13,
       units = c("in"),
       dpi = 300,
       limitsize = FALSE,
       bg = NULL)

rampboth
ggsave("rampboth.png",
       plot = last_plot(),
       device = NULL,
       path = NULL,
       scale = 1,
       width = 15,
       height = 13,
       units = c("in"),
       dpi = 300,
       limitsize = FALSE,
       bg = NULL)



# In-Text metrics ---------------------------------------------------------

view(livewellexp)
class(livewellexp$chase.t)
sd(livewellexp$length,na.rm=TRUE)
summary(livewellexp)


# Load in Temperature  ----------------------------------------------------
####control 2
cntrl2fish1<-read_csv("Control 2 livewell 1 2023-10-31 21_27_53 CDT (Data CDT).csv", 
                      col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
cntrl2fish2<-read_csv("Control 2 livewell 7 2023-10-31 21_30_07 CDT (Data CDT).csv", 
                      col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
cntrl2fish3<-read_csv("2H8 2023-11-01 16_20_14 CDT (Data CDT).csv", 
                      col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
control2temps<-rbind(cntrl2fish1,cntrl2fish2,cntrl2fish3) %>% rename(date= 'Date-Time (CDT)',temperature= 'Ch:1 - Temperature   (°C)')

####treatment 1
trt1fish1<-read_csv("Livewell 1 trtmnt 3 2023-10-31 21_32_14 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
trt1fish2<-read_csv("Livewell 2 trtmt 3 2023-10-31 21_28_52 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
trt1fish3<-read_csv("Livewell 5 trtmt 3 2023-10-31 21_30_57 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
trt1fish4<-read_csv("Livewell 7 trtmt 3 2023-10-31 21_33_07 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
treatment1temps<-rbind(trt1fish1,trt1fish2,trt1fish3,trt1fish4) %>% rename(date= 'Date-Time (CDT)',temperature= 'Ch:1 - Temperature   (°C)')
####treatment 2
trt2fish1<- read_csv("5H1 2023-11-01 22_50_21 CDT (Data CDT).csv", 
                     col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S"))) 
trt2fish2<-read_csv("5H2 2023-11-01 22_55_32 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
trt2fish3<-read_csv("5H3 2023-11-01 22_51_24 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
trt2fish4<-read_csv("5H5 2023-11-01 22_52_12 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
trt2fish5<-read_csv("5H6 2023-11-01 22_53_47 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
trt2fish6<-read_csv("5H7 2023-11-01 22_53_02 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
trt2fish7<-read_csv("5H8 2023-11-01 22_49_39 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
treatment2temps<-rbind(trt2fish1,trt2fish2,trt2fish3,trt2fish4,trt2fish5,trt2fish6,trt2fish7) %>% rename(date= 'Date-Time (CDT)',temperature= 'Ch:1 - Temperature   (°C)')

###treatment 3
trt3fish1<- read_csv("7H1 2023-11-01 14_49_41 CDT (Data CDT).csv", 
                     col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S"))) 
trt3fish2<-read_csv("7H2 2023-11-01 14_46_09 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
trt3fish3<-read_csv("7H3 2023-11-01 14_43_09 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
trt3fish4<-read_csv("7H5 2023-11-01 14_48_27 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
trt3fish5<-read_csv("7H6 2023-11-01 14_54_10 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
trt3fish6<-read_csv("7H7 2023-11-01 14_44_15 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
trt3fish7<-read_csv("7H8 2023-11-01 14_47_19 CDT (Data CDT).csv", 
                    col_types = cols(`Date-Time (CDT)` = col_datetime(format = "%m/%d/%Y %H:%M:%S")))
treatment3temps<-rbind(trt3fish1,trt3fish2,trt3fish3,trt3fish4,trt3fish5,trt3fish6,trt3fish7) %>% rename(date= 'Date-Time (CDT)',temperature= 'Ch:1 - Temperature   (°C)')

###EPilimnion
marina_temp_epi<- read_csv("Marina 2023-11-22 15_52_59 CST (Data CST).csv", 
                           col_types = cols(`Date-Time (CST/CDT)` = col_datetime(format = "%m/%d/%Y %H:%M"))) %>% 
  rename(Date='Date-Time (CST/CDT)',temperature='Ch:1 - Temperature   (°C)')
marina_temp_epi['Time'] <-format(as.POSIXct(marina_temp_epi$Date,format="%Y-%m-%d %H:%M"),"%H:%M:%S")
marina_temp_epi$Date <- format(as.POSIXct(marina_temp_epi$Date,format="%Y-%m-%d %H:%M"),"%Y-%m-%d")

marina_temp_epi <-marina_temp_epi %>%  mutate(Type = case_when(locations == 'Lake'~'Epi'))
#marina_temp_epi$time <- as.POSIXct(marina_temp_epi$Time, format = '%H:%M', tz = 'UTC')
marina_temp_epi$Date<-as.Date(marina_temp_epi$Date)
marina_temp_epi$Time<-as.POSIXct(marina_temp_epi$Time)

#marina_temp_epi$Timestamp<-as.times(marina_temp_epi$Timestamp)
class(marina_temp_epi$Time)

Epilimnion<-marina_temp_epi %>% filter(Date == "2023-05-06"| Date == "2023-06-10" | 
                                         Date == "2023-06-25" | Date == "2023-07-08")


Epilimnion <-Epilimnion %>%  mutate(Tournament = case_when(Date == "2023-05-06" ~ "tournament1", Date== "2023-06-10" ~ "tournament2",
                                                           Date== "2023-06-25" ~ "tournament3",Date== "2023-07-08" ~ "tournament4"))

Epilimnion<-Epilimnion %>% select(Type,temperature,Date, Time,Tournament) %>% rename(Temperature= 'temperature')
angler_all_df <- read_csv("angler_all_df.csv")%>% select(Type,Temperature,Date, Time,Tournament)
field_data_plot<- rbind(angler_all_df,Epilimnion)
#Field data plot

#########################Plot for paper
tournament1sub<-subset(field_data_plot, Tournament=='tournament1')
tournament1sub$Time<-as.POSIXct(tournament1sub$Time)
write.xlsx(tournament1sub,"C:/Users/Allison/Box/Hay_Allison/Data/Livewell Exp/R Directory")
tournament1sub <- read_csv("tournament1sub.csv", 
                           col_types = cols(Time = col_time(format = "%H:%M:%S")))
tournament1sub$Time<-as.POSIXct(tournament1sub$Time, tz = "UTC")
class(tournament1sub$Time)


p1<-tournament1sub %>%
  ggplot(aes(x=Time,
             y=Temperature, 
             group=Type,
             colour=Type))+
  geom_line(aes(linewidth=3,alpha=.5))+
  geom_borderline(linewidth=2, bordercolour = "black")+
  scale_color_manual(name= "Type", values = c("darkgray","darkgray", "darkgray","darkgray","darkgray","yellow", "darkgray","darkgray","darkgray","darkgray", "red"))+
  coord_cartesian(ylim = c(14,25))+
  scale_y_continuous(breaks = c(15,20,25),guide = "prism_minor")+
  labs(x=NULL, y="Temperature (°C)", title="Tournament 1") + 
  scale_x_datetime(date_labels = '%H:%M',breaks = '1 hour',limits =  as.POSIXct(c(ymd_hms("1970-01-01 07:00:00"),ymd_hms("1970-01-01 13:00:00")),
                                                                                format = "%Y-%M-%D %H:%M:%S"))+
  theme(
    legend.position = 'none',
    plot.background = element_blank(),
    plot.title = element_text(size = 25),
    axis.title.y =  element_text(size = 20,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 18),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank()) 
p1
###tournament 2
tournament2sub<-subset(field_data_plot, Tournament=='tournament2')
tournament2sub$Time<-as.POSIXct(tournament2sub$Time, tz = "UTC")
class(tournament2sub$Time)

p2<-tournament2sub %>%
  ggplot(aes(x=Time,
             y=Temperature, 
             group=Type,
             colour=Type))+
  geom_line(aes(linewidth=3,alpha=.5))+
  geom_borderline(linewidth=2, bordercolour = "black")+
  scale_color_manual(name= "Type", values = c("darkgray","darkgray", "darkgray","yellow","blue","darkgray", "darkgray","darkgray","red"))+
  coord_cartesian(ylim = c(10,30))+
  scale_y_continuous(breaks = c(10,15,20,25,30),guide = "prism_minor")+
  labs(x=NULL, y="Temperature (°C)", title="Tournament 2") + 
  scale_x_datetime(date_labels = '%H:%M',breaks = '1 hour',limits =  as.POSIXct(c(ymd_hms("1970-01-01 07:00:00"),ymd_hms("1970-01-01 13:00:00")),
                                                                                format = "%Y-%M-%D %H:%M:%S"))+
  theme(
    legend.position = 'none',
    plot.background = element_blank(),
    plot.title = element_text(size = 25),
    axis.title.y =  element_text(size = 20,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 18),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank()) 
p2
###tournament 3
tournament3sub<-subset(field_data_plot, Tournament=='tournament3')
tournament3sub$Time<-as.POSIXct(tournament3sub$Time, tz = "UTC")
class(tournament3sub$Time)

p3<-tournament3sub %>%
  ggplot(aes(x=Time,
             y=Temperature, 
             group=Type,
             colour=Type))+
  geom_line(aes(linewidth=3,alpha=.5))+
  geom_borderline(linewidth=2, bordercolour = "black")+
  scale_color_manual(name= "Type", values = c("yellow","yellow", "darkgray","yellow","darkgreen","yellow", "darkgray","red"))+
  coord_cartesian(ylim = c(25,32))+
  scale_y_continuous(breaks = c(25,30),guide = "prism_minor")+
  labs(x=NULL, y="Temperature (°C)", title="Tournament 3") + 
  scale_x_datetime(date_labels = '%H:%M',breaks = '1 hour',limits =  as.POSIXct(c(ymd_hms("1970-01-01 07:00:00"),ymd_hms("1970-01-01 13:00:00")),
                                                                                format = "%Y-%M-%D %H:%M:%S"))+
  theme(
    legend.position = 'none',
    plot.background = element_blank(),
    plot.title = element_text(size = 25),
    axis.title.y =  element_text(size = 20,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 18),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank()) 
p3
###tournament 4
tournament4sub<-subset(field_data_plot, Tournament=='tournament4')
tournament4sub$Time<-as.POSIXct(tournament4sub$Time, tz = "UTC")
class(tournament4sub$Time)

p4<-tournament4sub %>%
  ggplot(aes(x=Time,
             y=Temperature, 
             group=Type,
             colour=Type))+
  geom_line(aes(linewidth=3,alpha=.5))+
  geom_borderline(linewidth=2, bordercolour = "black")+
  scale_color_manual(name= "Type", values = c("darkgreen","darkgray", "blue","darkgray","darkgray","darkgreen", "darkgray","darkgray","darkgray","darkgray","red"))+
  coord_cartesian(ylim = c(18,30))+
  scale_y_continuous(breaks = c(20,25,30),guide = "prism_minor")+
  labs(x=NULL, y="Temperature (°C)", title="Tournament 4") + 
  scale_x_datetime(date_labels = '%H:%M',breaks = '1 hour',limits =  as.POSIXct(c(ymd_hms("1970-01-01 07:00:00"),ymd_hms("1970-01-01 13:00:00")),
                                                                                format = "%Y-%M-%D %H:%M:%S"))+
  theme(
    legend.position = 'none',
    plot.background = element_blank(),
    plot.title = element_text(size = 25),
    axis.title.y =  element_text(size = 20,face = "bold"),
    axis.title.x =  element_text(size = 15, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 18),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank()) 
p4
top<-plot_grid(p1,p2, ncol = 1)
bottom<-plot_grid(p3,p4,ncol = 1)
plot_grid(top,bottom, ncol = 1)
ggsave("fieldlivewell.png",
       plot = last_plot(),
       device = NULL,
       path = NULL,
       scale = 1,
       width = 27,
       height = 25,
       units = c("in"),
       dpi = 300,
       limitsize = FALSE,
       bg = NULL)
# plot treatment temps ----------------------------------------------------
lims <- as.POSIXct(strptime(c("2023-10-31 11:00","2023-10-31 21:00"), format = "%Y-%m-%d %H:%M"))    



p1<-treatment1temps%>% ggplot(aes(x=date,
                                  y=temperature, 
                                  group=id,
                                  colour=id)) +
  
  geom_line(aes(linewidth=2,alpha=.5))+
  geom_borderline(linewidth=2, bordercolour = "black")+
  scale_color_manual(name= "id", values = c("#3B3B3B", "#3B3B3B","#3B3B3B","#3B3B3B"))+
  scale_x_datetime(date_minor_breaks = "1 hour", date_breaks = "2 hours",
                   date_labels = "%H:%M",limits =  as.POSIXct(c("2023-10-31 06:00:00","2023-10-31 14:00:00"),
                                                              format = "%Y-%m-%d %H:%M:%S")  )+
  labs(y= "Temperature °C", x=NULL)+
  coord_cartesian(ylim = c(20,32))+
  scale_y_continuous(limits = NULL, breaks = c(20,24,28,32))+
  theme(
    legend.position = 'none',
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 25,face = "bold"),
    axis.title.x =  element_text(size = 10, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(size = 20, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())
p1


vtreatment1temps%>% ggplot+geom_point(aes(x=date,
                                          y=temperature, 
                                          group=id,
                                          colour=id),alpha=.9,size=1) + 
  scale_x_datetime(
    date_minor_breaks = "1 hour", date_breaks = "2 hours",
    date_labels = "%H:%M")+
  coord_cartesian(ylim = c(20,32))+
  scale_y_continuous(limits = NULL, breaks = c(20,24,28,32))



p2<-treatment2temps%>% ggplot(aes(x=date,
                                  y=temperature, 
                                  group=id,
                                  colour=id)) +
  
  geom_line(aes(linewidth=2,alpha=.5))+
  geom_borderline(linewidth=2, bordercolour = "black")+
  scale_color_manual(name= "id", values = c("#3B3B3B", "#3B3B3B","#3B3B3B","#3B3B3B","#3B3B3B","#3B3B3B","#3B3B3B"))+
  scale_x_datetime(date_minor_breaks = "1 hour", date_breaks = "2 hours",
                   date_labels = "%H:%M",limits =  as.POSIXct(c("2023-11-1 10:00:00","2023-11-1 17:00:00"),
                                                              format = "%Y-%m-%d %H:%M:%S")  )+
  labs(y= "Temperature °C", x=NULL)+
  coord_cartesian(ylim = c(18,27))+
  scale_y_continuous(limits = NULL, breaks = c(20,24,27))+
  theme(
    legend.position = 'none',
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 25,face = "bold"),
    axis.title.x =  element_text(size = 10, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(size = 20, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())

p3<-treatment3temps%>% ggplot(aes(x=date,
                                  y=temperature, 
                                  group=id,
                                  colour=id)) +
  
  geom_line(aes(linewidth=2,alpha=.5))+
  geom_borderline(linewidth=2, bordercolour = "black")+
  scale_color_manual(name= "id", values = c("#3B3B3B", "#3B3B3B","#3B3B3B","#3B3B3B","#3B3B3B","#3B3B3B","#3B3B3B"))+
  scale_x_datetime(date_minor_breaks = "1 hour", date_breaks = "2 hours",
                   date_labels = "%H:%M", limits =  as.POSIXct(c("2023-11-1 03:00:00","2023-11-1 09:00:00"),
                                                               format = "%Y-%m-%d %H:%M:%S")  )+
  labs(y= "Temperature °C", x=NULL)+
  coord_cartesian(ylim = c(15,21))+
  scale_y_continuous(limits = NULL, breaks = c(15,17,19,21))+
  theme(
    legend.position = 'none',
    plot.background = element_blank(),
    axis.title.y =  element_text(size = 25,face = "bold"),
    axis.title.x =  element_text(size = 10, face = "bold"),
    axis.ticks.length = unit(.50,"cm"),
    axis.ticks = element_line(size = 1.00, colour = "black"),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    axis.line = element_line(size = 1.00, color = "black"),
    axis.text.x = element_text(size = 20, face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor =  element_blank(),
    panel.grid.major =  element_blank(),
    panel.background =  element_blank())


plot_grid(p1,p2,p3, ncol = 1)
ggsave("treatmenttemp.png",
       plot = last_plot(),
       device = NULL,
       path = NULL,
       scale = 1,
       width = 22,
       height = 15,
       units = c("in"),
       dpi = 500,
       limitsize = FALSE,
       bg = NULL) 
plot_grid(top2,p3, labels = '','C', label_size = 20, ncol = 1, rel_heights = 2,1)