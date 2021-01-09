####----INTRODUCTION----#####
# Manuscript title: "Proprioceptive afferents differentially contribute to effortful perception of object heaviness and length"
# Manuscript authors: Madhur Mangalam, Nisarg Desai, and Damian G. Kelty-Stephen
# Code author: Nisarg Desai, desai054[at]umn[dot]edu
# Last update: March 21, 2020
# Status: Manuscript submitted for peer-review

rm(list=ls())

library(nlme)
library(lme4)
library(dplyr)      #%>%
library(emmeans)    #emmeans
library(DescTools)  #EtaSq
setwd("...to the directory containing Data.csv...")
DATA <- read.csv("Data.csv", TRUE, ",",na.strings=" ")
DATA$Wrist_Angle <- as.factor(DATA$Wrist_Angle)
DATA$Wrist_Angle <- relevel(DATA$Wrist_Angle, ref = "Neutral")
DATA$Exploratory_Kinematics <- as.factor(DATA$Exploratory_Kinematics)
DATA$Trial_Order <-as.factor((DATA$Trial_Order))
DATA$Torque <-ordered(DATA$Torque)
DATA$Object <- as.factor(DATA$Object)
DATA$Participant <- as.factor(DATA$Participant)

# Let's first check the distributions of the dependent variables

hist(DATA$Hperceived) # Looks somewhat Poisson
hist(DATA$Lperceived) # Looks somewhat Normal


shapiro.test(DATA$Hperceived) # Not normal

library(fitdistrplus)
fit = fitdist(DATA$Hperceived, "pois", discrete=TRUE)
gofstat(fit) # goodness-of-fit test for Poisson--Not Poisson

library(car)
library(MASS)

qqp(DATA$Hperceived, "norm") # Not normal
qqp(DATA$Hperceived, "lnorm") # Not lognormal

poisson <- fitdistr(DATA$Hperceived, "Poisson")
qqp(DATA$Hperceived, "pois", lambda = poisson$estimate) # Not Poisson

nbinom <- fitdistr(DATA$Hperceived, "Negative Binomial")
qqp(DATA$Hperceived, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]]) # Not Negative binomial

gamma <- fitdistr(DATA$Hperceived, "gamma")
qqp(DATA$Hperceived, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]]) # Not Gamma

# So Hperceived does not fit any of the standard distributions

# Now lets check Lperceived

qqp(DATA$Lperceived, "norm") # Some deviations from normality
qqp(DATA$Lperceived, "lnorm") # Not lognormal

poisson <- fitdistr(round(DATA$Lperceived), "Poisson")
qqp(DATA$Lperceived, "pois", lambda = poisson$estimate) # some deviations from Poisson

nbinom <- fitdistr(DATA$Lperceived, "Negative Binomial")
qqp(DATA$Lperceived, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]]) # Not negative binomial

gamma <- fitdistr(DATA$Lperceived, "gamma")
qqp(DATA$Lperceived, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]]) # Some deviations from gamma

shapiro.test(DATA$Lperceived) # Not normal

# Scaling the continuous variavles to mean 0 and var 1
DATA$sLogI1 <- scale(DATA$LogI1)
DATA$sLogI3 <- scale(DATA$LogI3)

library(plyr)
# Descriptive summaries of Hperceived
ddply(DATA, ~ Wrist_Angle * Exploratory_Kinematics, function(DATA) summary(DATA$Hperceived))
ddply(DATA, ~ Wrist_Angle * Exploratory_Kinematics, summarise, Hperceived.mean=mean(Hperceived), Hperceived.sd=sd(Hperceived))
ddply(DATA, ~ Object, summarise, Hperceived.mean=mean(Hperceived), Hperceived.sd=sd(Hperceived))


# Checking if Hperceived looks poisson for combinations of factors
hist(DATA[DATA$Wrist_Angle == "Neutral" & DATA$Exploratory_Kinematics == "0",]$Hperceived)
hist(DATA[DATA$Wrist_Angle == "Neutral" & DATA$Exploratory_Kinematics == "1",]$Hperceived)
hist(DATA[DATA$Wrist_Angle == "Neutral" & DATA$Exploratory_Kinematics == "2",]$Hperceived)
hist(DATA[DATA$Wrist_Angle == "Radial" & DATA$Exploratory_Kinematics == "0",]$Hperceived)
hist(DATA[DATA$Wrist_Angle == "Radial" & DATA$Exploratory_Kinematics == "1",]$Hperceived)
hist(DATA[DATA$Wrist_Angle == "Radial" & DATA$Exploratory_Kinematics == "2",]$Hperceived)
hist(DATA[DATA$Wrist_Angle == "Ulnar" & DATA$Exploratory_Kinematics == "0",]$Hperceived)
hist(DATA[DATA$Wrist_Angle == "Ulnar" & DATA$Exploratory_Kinematics == "1",]$Hperceived)
hist(DATA[DATA$Wrist_Angle == "Ulnar" & DATA$Exploratory_Kinematics == "2",]$Hperceived)

# Check if Hperceived differs by combinations of factors
boxplot(Hperceived ~ Wrist_Angle * Exploratory_Kinematics, data=DATA, xlab="WA EK", ylab="Hperceived") # boxplots

# Interaction plots
with(DATA, interaction.plot(Exploratory_Kinematics, Wrist_Angle, Hperceived, ylim=c(150, 200))) # interaction plot
with(DATA, interaction.plot(Torque, Wrist_Angle, Hperceived, ylim=c(100, 300))) # interaction plot
with(DATA, interaction.plot(Torque, Exploratory_Kinematics, Hperceived, ylim=c(100, 300))) # interaction plot

with(DATA, interaction.plot(Object, Torque, Hperceived, ylim=c(100, 300))) # interaction plot


library(ggplot2)

(prelim_plot <- ggplot(DATA, aes(x = LogI1, y = Hperceived)) +
    geom_point() +
    geom_smooth(method = "lm"))

(prelim_plot <- ggplot(DATA, aes(x = LogI3, y = Hperceived)) +
    geom_point() +
    geom_smooth(method = "lm"))



# So since we don't have any standard distributions fitting the Hperceived or Lperceived, 
# we can use Aligned Rank Transform, a non parametric technique that will incorporate 
# pseudoreplication due to many participants. However, ART only takes in categorical predictors. 

library(ARTool)

m1 <- art(Hperceived ~ Wrist_Angle*Exploratory_Kinematics*Object + (1|Participant), data = DATA)
m1.anova <- anova(m1)
shapiro.test(residuals(m1)) # normality? Even though it's non-parametric, it's still an ANOVA and residual normality should be met
qqnorm(residuals(m1)); qqline(residuals(m1)) # seems to not conform, reduces confidence in the results

print(m1.anova, verbose = T)
m1.anova$eta.sq.part = with(m1.anova, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))
m1.anova

# For Cohen's d
mlme1 <- lmer(Hperceived ~ Wrist_Angle*Exploratory_Kinematics*Object + (1|Participant), data = DATA)
lsmeans(mlme1, pairwise ~ Exploratory_Kinematics, adjust = "tukey")
lsmeans(mlme1, pairwise ~ Wrist_Angle, adjust = "tukey")
lsmeans(mlme1, pairwise ~ Object, adjust = "tukey")


# Pairwise contrasts
library(lsmeans) # for lsmeans
lsmeans(artlm(m1, "Exploratory_Kinematics"), pairwise ~ Exploratory_Kinematics)
lsmeans(artlm(m1, "Wrist_Angle"), pairwise ~ Wrist_Angle)
lsmeans(artlm(m1, "Object"), pairwise ~ Object)

#"interaction contrasts." see vignette("art-contrasts")
library(phia)
testInteractions(artlm(m1, "Wrist_Angle:Exploratory_Kinematics"), pairwise=c("Wrist_Angle", "Exploratory_Kinematics"), adjustment="holm")
testInteractions(artlm(m1, "Exploratory_Kinematics:Object"), pairwise=c("Exploratory_Kinematics", "Object"), adjustment="holm")
testInteractions(artlm(m1, "Wrist_Angle:Torque"), pairwise=c("Wrist_Angle", "Torque"), adjustment="holm")
# in the output, A-B : C-D is interpreted as a difference-of-differences, i.e., the difference 
# between (A-B | C) and (A-B | D). in words, is the difference between A and B significantly 
# different in condition C from condition D?

# With Object
m3 <- art(Hperceived ~ Object + (1|Participant), data = DATA)
anova(m3)

lsmeans(artlm(m3, "Object"), pairwise ~ Object)


                      
## Lperceived

# histograms for two factors
hist(DATA[DATA$Wrist_Angle == "Neutral" & DATA$Exploratory_Kinematics == "0",]$Lperceived)
hist(DATA[DATA$Wrist_Angle == "Neutral" & DATA$Exploratory_Kinematics == "1",]$Lperceived)
hist(DATA[DATA$Wrist_Angle == "Neutral" & DATA$Exploratory_Kinematics == "2",]$Lperceived)
hist(DATA[DATA$Wrist_Angle == "Radial" & DATA$Exploratory_Kinematics == "0",]$Lperceived)
hist(DATA[DATA$Wrist_Angle == "Radial" & DATA$Exploratory_Kinematics == "1",]$Lperceived)
hist(DATA[DATA$Wrist_Angle == "Radial" & DATA$Exploratory_Kinematics == "2",]$Lperceived)
hist(DATA[DATA$Wrist_Angle == "Ulnar" & DATA$Exploratory_Kinematics == "0",]$Lperceived)
hist(DATA[DATA$Wrist_Angle == "Ulnar" & DATA$Exploratory_Kinematics == "1",]$Lperceived)
hist(DATA[DATA$Wrist_Angle == "Ulnar" & DATA$Exploratory_Kinematics == "2",]$Lperceived)
                    
boxplot(Lperceived ~ Wrist_Angle * Exploratory_Kinematics, data=DATA, xlab="WA EK", ylab="Lperceived") # boxplots
with(DATA, interaction.plot(Exploratory_Kinematics, Wrist_Angle, Lperceived, ylim=c(45, 60))) # interaction plot
with(DATA, interaction.plot(Torque, Wrist_Angle, Lperceived, ylim=c(45, 60))) # interaction plot
with(DATA, interaction.plot(Torque, Exploratory_Kinematics, Lperceived, ylim=c(45, 60))) # interaction plot
with(DATA, interaction.plot(Object, Torque, Lperceived, ylim=c(40, 60)))

                    
library(ggplot2)
                    
(prelim_plot <- ggplot(DATA, aes(x = LogI1, y = Lperceived)) +
                        geom_point() +
                        geom_smooth(method = "lm"))
                    
(prelim_plot <- ggplot(DATA, aes(x = LogI3, y = Lperceived)) +
                        geom_point() +
                        geom_smooth(method = "lm"))
                    
m2 <- art(Lperceived ~ Exploratory_Kinematics*Wrist_Angle*Object + (1|Participant), data = DATA)
Resultm2 = anova(m2)
Result$part.eta.sq = with(Result, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))

# For Cohen's d
mlme2 <- lmer(Lperceived ~ Exploratory_Kinematics*Wrist_Angle*Object + (1|Participant), data = DATA)
lsmeans(mlme2, pairwise ~ Exploratory_Kinematics, adjust = "tukey")
lsmeans(mlme2, pairwise ~ Wrist_Angle, adjust = "tukey")
lsmeans(mlme2, pairwise ~ Object, adjust = "tukey")
                    
shapiro.test(residuals(m2)) # normality? Even though it's non-parametric, it's still an ANOVA and residual normality should be met
qqnorm(residuals(m2)); qqline(residuals(m2)) # seems to mostly conform

# Pairwise contrasts
library(lsmeans) # for lsmeans
lsmeans(artlm(m2, "Exploratory_Kinematics"), pairwise ~ Exploratory_Kinematics)
lsmeans(artlm(m2, "Wrist_Angle"), pairwise ~ Wrist_Angle)
lsmeans(artlm(m2, "Object"), pairwise ~ Object)

#"interaction contrasts." see vignette("art-contrasts")
library(phia)
testInteractions(artlm(m2, "Exploratory_Kinematics:Wrist_Angle"), pairwise=c("Exploratory_Kinematics", "Wrist_Angle"), adjustment="holm")
testInteractions(artlm(m2, "Exploratory_Kinematics:Torque"), pairwise=c("Exploratory_Kinematics", "Torque"), adjustment="holm")
testInteractions(artlm(m2, "Wrist_Angle:Torque"), pairwise=c("Wrist_Angle", "Torque"), adjustment="holm")


# Further analyses for LogI and torque interactions. I will do Poisson regression for Heaviness
# and Linear regression for Length because the "appear" to be close to those distributions.
# However, the results have to be taken with a pinch of salt as the assumptions are violated in all cases

with(DATA, interaction.plot(Object, Torque, Hperceived, ylim=c(0, 400))) # interaction plot
with(DATA, interaction.plot(Object, Torque, Lperceived, ylim=c(0, 100))) # interaction plot
with(DATA, interaction.plot(Object, Torque, LogI1, ylim=c(4, 6))) # interaction plot
with(DATA, interaction.plot(Object, Torque, LogI3, ylim=c(2, 5))) # interaction plot

m4 <- art(Lperceived ~ Object + (1|Participant), data = DATA)
anova(m4)

lsmeans(artlm(m4, "Object"), pairwise ~ Object)


############## MODEL SELECTIONS

modsel <- read.csv("Data_2.csv")

# Scale the variables
modsel$sMass <- scale(modsel$Mass)
modsel$sTorque <- scale(modsel$Torque)
modsel$sI <- scale(modsel$LogI1)

# Models for Lperceived
LP001<-lmer(Lperceived ~ (1|Participant), data = modsel, REML = F)
LP002<-lmer(Lperceived ~ sMass + (1|Participant), data = modsel, REML = F)
LP003<-lmer(Lperceived ~ sTorque + (1|Participant), data = modsel, REML = F)
LP004<-lmer(Lperceived ~ sI + (1|Participant), data = modsel, REML = F)
LP005<-lmer(Lperceived ~ sMass + sTorque + (1|Participant), data = modsel, REML = F)
LP006<-lmer(Lperceived ~ sMass + sI + (1|Participant), data = modsel, REML = F)
LP007<-lmer(Lperceived ~ sTorque + sI + (1|Participant), data = modsel, REML = F)
LP008<-lmer(Lperceived ~ sMass + sTorque + sI + (1|Participant), data = modsel, REML = F, na.action = "na.fail")
summary(LP008)

# Unscaled models # Skip next 9 lines if you want to use scaled
LP001<-lmer(Lperceived ~ (1|Participant), data = modsel, REML = F)
LP002<-lmer(Lperceived ~ Mass + (1|Participant), data = modsel, REML = F)
LP003<-lmer(Lperceived ~ Torque + (1|Participant), data = modsel, REML = F)
LP004<-lmer(Lperceived ~ LogI1 + (1|Participant), data = modsel, REML = F)
LP005<-lmer(Lperceived ~ Mass + Torque + (1|Participant), data = modsel, REML = F)
LP006<-lmer(Lperceived ~ Mass + LogI1 + (1|Participant), data = modsel, REML = F)
LP007<-lmer(Lperceived ~ Torque + LogI1 + (1|Participant), data = modsel, REML = F)
LP008<-lmer(Lperceived ~ Mass + Torque + LogI1 + (1|Participant), data = modsel, REML = F, na.action = "na.fail")
summary(LP008)



# Calculate the dispersion parameter
chat_lp <- deviance(LP008) / df.residual(LP008)

library(MuMIn)
# Create a model selection table by ranking models based on QAICc
ms_tab_lp <- dredge(LP008, rank = "QAICc", chat = chat_lp) 
ms_tab_lp
#write.csv(ms_tab, "Length models.csv")

# Perform model-averaging by ranking models based on QAICc
p <- model.avg(LP008, LP007, LP006, LP005, LP004, LP003, LP002, LP001, rank = QAICc, rank.args = list(chat = chat_lp)) # perform model averaging
summary(p)
confint(p)

# Models for Hoerceived
HP001<-lmer(Hperceived ~ (1|Participant), data = modsel, REML = F)
HP002<-lmer(Hperceived ~ sMass + (1|Participant), data = modsel, REML = F)
HP003<-lmer(Hperceived ~ sTorque + (1|Participant), data = modsel, REML = F)
HP004<-lmer(Hperceived ~ sI + (1|Participant), data = modsel, REML = F)
HP005<-lmer(Hperceived ~ sMass + sTorque + (1|Participant), data = modsel, REML = F)
HP006<-lmer(Hperceived ~ sMass + sI + (1|Participant), data = modsel, REML = F)
HP007<-lmer(Hperceived ~ sTorque + sI + (1|Participant), data = modsel, REML = F)
HP008<-lmer(Hperceived ~ sMass + sTorque + sI + (1|Participant), data = modsel, REML = F, na.action = "na.fail")

# Unscaled models # Skip next 9 lines if you want to use scaled
HP001<-lmer(Hperceived ~ (1|Participant), data = modsel, REML = F)
HP002<-lmer(Hperceived ~ Mass + (1|Participant), data = modsel, REML = F)
HP003<-lmer(Hperceived ~ Torque + (1|Participant), data = modsel, REML = F)
HP004<-lmer(Hperceived ~ LogI1 + (1|Participant), data = modsel, REML = F)
HP005<-lmer(Hperceived ~ Mass + Torque + (1|Participant), data = modsel, REML = F)
HP006<-lmer(Hperceived ~ Mass + LogI1 + (1|Participant), data = modsel, REML = F)
HP007<-lmer(Hperceived ~ Torque + LogI1 + (1|Participant), data = modsel, REML = F)
HP008<-lmer(Hperceived ~ Mass + Torque + LogI1 + (1|Participant), data = modsel, REML = F, na.action = "na.fail")



# Calculate the dispersion parameter
chat_hp <- deviance(HP008) / df.residual(HP008)

library(MuMIn)
# Create a model selection table by ranking models based on QAICc
ms_tab_hp <- dredge(HP008, rank = "QAICc", chat = chat_hp) 
ms_tab_hp
#write.csv(ms_tab, "Heaviness models.csv")

# Perform model-averaging by ranking models based on QAICc
p <- model.avg(HP008, HP007, HP006, HP005, HP004, HP003, HP002, HP001, rank = QAICc, rank.args = list(chat = chat_hp)) # perform model averaging
summary(p)
confint(p)
