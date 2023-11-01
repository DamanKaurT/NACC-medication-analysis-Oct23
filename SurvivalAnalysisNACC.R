################## All Meds #####################

rm(list = ls())
library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)
library(car)
library(coxphw)
library(survival)
library(ggfortify)
library(ggplot2)
library(survminer)
library(sjPlot)
library(broom)

File <- file.path("C:/Users/Daman Kaur/Desktop/Ulster/Research Projects/NACC drug classes analysis paper/Cox regression csv")
setwd(File)
nacc <- fread("MDallmedscdrcox.csv")
df <- nacc

##### only for anxiolytics ####
df <- df %>%
  group_by(NACCID) %>%
  mutate(check1=is.na(ANX1) & is.na(ANX2) & is.na(ANX3))
df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))
df <- df[,-(148)] 
df <- as.data.frame(df)
df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(ANX1, na.rm=TRUE))
df$ANX1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$ANX1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(ANX2, na.rm=TRUE))
df$ANX2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$ANX2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(ANX3, na.rm=TRUE))
df$ANX3[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$ANX3[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(148)] 
df <- as.data.frame(df)
df$ANX[(df$ANX1==1 | df$ANX2==1 | df$ANX3==1)] <- 1
df$ANX[(df$ANX1==0 & df$ANX2==0 & df$ANX3==0)] <- 0

##### only for antipsychotics######
df <- df %>%
  group_by(NACCID) %>%
  mutate(check1=is.na(APS1) & is.na(APS2) & is.na(APS3) & is.na(APS4) & is.na(APS5))
df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))
df <- df[,-(148)] 
df <- as.data.frame(df)

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(APS1, na.rm=TRUE))
df$APS1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$APS1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(APS2, na.rm=TRUE))
df$APS2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$APS2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(APS3, na.rm=TRUE))
df$APS3[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$APS3[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(APS4, na.rm=TRUE))
df$APS4[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$APS4[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(APS5, na.rm=TRUE))
df$APS5[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$APS5[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df[,-(148)] 
df <- as.data.frame(df)

df$PSYC[(df$APS1==1 | df$APS2==1 | df$APS3==1 | df$APS4==1 | df$APS5==1)] <- 1
df$PSYC[(df$APS1==0 & df$APS2==0 & df$APS3==0 & df$APS4==0 & df$APS5==0)] <- 0



#### only for lipid lowering drugs #####
df <- df %>%
  group_by(NACCID) %>%
  mutate(check1=is.na(HYPERC1) & is.na(HYPERC2))
df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))
df <- df[,-(148)] 
df <- as.data.frame(df)

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(HYPERC1, na.rm=TRUE))
df$HYPERC1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$HYPERC1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(HYPERC2, na.rm=TRUE))
df$HYPERC2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$HYPERC2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df[,-(148)] 
df <- as.data.frame(df)

df$HYPERC[(df$HYPERC1==1 | df$HYPERC2==1)] <- 1
df$HYPERC[(df$HYPERC1==0 & df$HYPERC2==0)] <- 0

##### only for Parkinson's drugs ####
df <- df %>%
  group_by(NACCID) %>%
  mutate(check1=is.na(PARK1) & is.na(PARK2))
df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))
df <- df[,-(148)] 
df <- as.data.frame(df)

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(PARK1, na.rm=TRUE))
df$PARK1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$PARK1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(PARK2, na.rm=TRUE))
df$PARK2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$PARK2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df[,-(148)] 
df <- as.data.frame(df)

df$PARK[(df$PARK2==1 | df$PARK2==1)] <- 1
df$PARK[(df$PARK2==0 & df$PARK2==0)] <- 0

##### continue here ####
df <- as.data.frame(df)

df$SEX[df$SEX==2] <- 0  # female 0, male 1
df$RACE[df$RACE==50] <- NA

xy <- df %>%
  mutate(date_visit = as.Date(Date, "%Y/%m/%d"))
xy <- xy %>% 
  group_by(NACCID) %>%
  summarise(first_visit = min(date_visit), 
            last_visit = max(date_visit), 
            exposure = last_visit - first_visit)
xy$exposure <- as.numeric(xy$exposure)

df <- df[df$NEW_VISITNUM == 1, ] 

new <- left_join(df,xy)
df <- as.data.frame(new)

df <- subset(df, NACCAGE>=40,
             select=c(NACCID:exposure))
df$years <- df$exposure / 365.25
df <- as.data.frame(df)
#df <- df %>% filter(SEX == 1)

# for analyzing only AD cases
#df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
#df <- df[df$Diagnosis != 1, ]
#df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% 
 replace_with_na(replace = list(NACCAANX = 9, NACCAC = 9, NACCADEP = 9, NACCADMD = 9,
                                 NACCAHTN = 9, NACCAPSY = 9, NACCDBMD = 9, NACCLIPL = 9, 
                                 NACCNSD = 9, NACCPDMD = 9))


#df.baseDHH <- df[,c(1,8:10,14,26,108,109,111:113,115,118,121:123,127,130:137,148,152)]
#df.baseDHH <- df.baseDHH[,c(1:6,16,17:27)] # select relevant drug class

df.baseDHH <- df[,c(1,8:10,14,26,108,109,111:113,115,118,121:123,127,130:137,151)]
df.baseDHH <- df.baseDHH[,c(1:6,10,17:26)] # select relevant drug class

df.baseDHH.compl <- df.baseDHH[complete.cases(df.baseDHH), ]
lg.baseDHH <- df.baseDHH.compl


# calculate ipw (case weights) using logistic regression
#1 AANX
logmod <- glm(NACCAANX ~ ANX +
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod, exponentiate = TRUE)
drug_ipw <- augment_columns(logmod, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCAANX/.fitted) + ((1-NACCAANX)/(1-.fitted)))
summary(drug_ipw$ipw) ### TRUNCATE ONLY IF NECESSARY

anx <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ANX,
  lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ANX, Diagnosis))))
write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
    write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}
write.excel(anx)

drugSMw <- coxphw(Surv(years, Diagnosis) ~ NACCAANX + ANX +
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, #trunc.weights = 0.99, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw$ipw) 
summary(drugSMw)

#2 AC
logmod2 <- glm(NACCAC ~ 
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod2, exponentiate = TRUE)
drug_ipw2 <- augment_columns(logmod2, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCAC/.fitted) + ((1-NACCAC)/(1-.fitted)))
summary(drug_ipw2$ipw) ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~ NACCAC + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, #trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw2$ipw) 
summary(drugSMw)


#3 ADEP
logmod3 <- glm(NACCADEP ~ 
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod3, exponentiate = TRUE)
drug_ipw3 <- augment_columns(logmod3, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCADEP/.fitted) + ((1-NACCADEP)/(1-.fitted)))
summary(drug_ipw3$ipw) ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~ NACCADEP + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, #trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw3$ipw) 
summary(drugSMw)

#4 AHTN
logmod4 <- glm(NACCAHTN ~ 
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod4, exponentiate = TRUE)
drug_ipw4 <- augment_columns(logmod4, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCAHTN/.fitted) + ((1-NACCAHTN)/(1-.fitted)))
summary(drug_ipw4$ipw) ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~ NACCAHTN + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, #trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw4$ipw) 
summary(drugSMw)


#5 APSY
logmod5 <- glm(NACCAPSY ~ PSYC +
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod5, exponentiate = TRUE)
drug_ipw5 <- augment_columns(logmod5, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCAPSY/.fitted) + ((1-NACCAPSY)/(1-.fitted)))
summary(drug_ipw5$ipw) ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~ NACCAPSY + PSYC +
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, #trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw5$ipw) 
summary(drugSMw)

psyc <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$PSYC,
  lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(PSYC, Diagnosis))))
write.excel(psyc)

#6 DBMD
logmod6 <- glm(NACCDBMD ~ 
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod6, exponentiate = TRUE)
drug_ipw6 <- augment_columns(logmod6, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCDBMD/.fitted) + ((1-NACCDBMD)/(1-.fitted)))
summary(drug_ipw6$ipw) ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~ NACCDBMD + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw6$ipw) 
summary(drugSMw)



#7 LIPL
logmod7 <- glm(NACCLIPL ~ HYPERC +
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT +  DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod7, exponentiate = TRUE)
drug_ipw7 <- augment_columns(logmod7, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCLIPL/.fitted) + ((1-NACCLIPL)/(1-.fitted)))
summary(drug_ipw7$ipw) ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~ NACCLIPL + HYPERC +
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD,# trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw7$ipw) 
summary(drugSMw)

hyc <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HYPERC,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HYPERC, Diagnosis))))
write.excel(hyc)


#8 NSD
logmod8 <- glm(NACCNSD ~ 
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod8, exponentiate = TRUE)
drug_ipw8 <- augment_columns(logmod8, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCNSD/.fitted) + ((1-NACCNSD)/(1-.fitted)))
summary(drug_ipw8$ipw) ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~ NACCNSD + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, #trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw8$ipw) 
summary(drugSMw)



#9 PDMD
logmod9 <- glm(NACCPDMD ~ PARK +
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod9, exponentiate = TRUE)
drug_ipw9 <- augment_columns(logmod9, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCPDMD/.fitted) + ((1-NACCPDMD)/(1-.fitted)))
summary(drug_ipw9$ipw) ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~ NACCPDMD + PARK +
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw9$ipw) 
summary(drugSMw)

park <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$PARK,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(PARK, Diagnosis))))
pdmd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCPDMD,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCPDMD, Diagnosis))))
write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}
write.excel(park)
write.excel(pdmd)

#9 ADMD
logmod10 <- glm(NACCADMD ~ 
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod10, exponentiate = TRUE)
drug_ipw10 <- augment_columns(logmod10, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCADMD/.fitted) + ((1-NACCADMD)/(1-.fitted)))
summary(drug_ipw10$ipw) ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~ NACCADMD + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, #trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw10$ipw) 
summary(drugSMw)

admd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCADMD,
 lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCADMD, Diagnosis))))
write.excel(admd)

####################### Drug Subcategories ##################
#### Antihypertensives ##################################

rm(list = ls())
library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)
library(car)
library(coxphw)
library(survival)
library(ggfortify)
library(ggplot2)
library(survminer)
library(sjPlot)

File <- file.path("C:/Users/Daman Kaur/Desktop/Ulster/Research Projects/NACC drug classes analysis paper/Cox regression csv")
setwd(File)
nacc <- fread("MDhypercox.csv")
df <- nacc

df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1

xy <- df %>%
  mutate(date_visit = as.Date(Date, "%Y/%m/%d"))
xy <- xy %>% 
  group_by(NACCID) %>%
  summarise(first_visit = min(date_visit), 
            last_visit = max(date_visit), 
            exposure = last_visit - first_visit)
xy$exposure <- as.numeric(xy$exposure)
#xy$years <- xy$exposure / 365.25

df <- df[df$NEW_VISITNUM == 1, ] 

new <- left_join(df,xy)
df <- as.data.frame(new)

df <- subset(df, NACCAGE>=40,
             select=c(NACCID:exposure))
df <- as.data.frame(df)
df$years <- df$exposure / 365.25

#df <- df %>% filter(SEX == 1)

#df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
#df <- df[df$Diagnosis != 1, ]
#df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% 
 replace_with_na(replace = list(NACCHTNC = 9, NACCACEI = 9, NACCAAAS = 9,
                                 NACCBETA = 9, NACCCCBS = 9, NACCDIUR = 9, 
                                NACCVASD = 9, NACCANGI = 9))


df.baseDHH <- df[,c(1,8:10,14,26,107,110,114,116,117,119,120,124,127,130:137,159)]
df.baseDHH <- df.baseDHH[,c(1:6,7,15:24)] # select relevant drug class

df.baseDHH.compl <- df.baseDHH[complete.cases(df.baseDHH), ]
lg.baseDHH <- df.baseDHH.compl

#1 AAAS
logmod <- glm(NACCAAAS ~ 
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod, exponentiate = TRUE)
drug_ipw <- augment_columns(logmod, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCAAAS/.fitted) + ((1-NACCAAAS)/(1-.fitted)))
summary(drug_ipw$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  NACCAAAS +
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw$ipw) 
summary(drugSMw)


#2 ACEI
logmod2 <- glm(NACCACEI ~ 
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod2, exponentiate = TRUE)
drug_ipw2 <- augment_columns(logmod2, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCACEI/.fitted) + ((1-NACCACEI)/(1-.fitted)))
summary(drug_ipw2$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  NACCACEI +
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw2$ipw) 
summary(drugSMw)

#3 ANGI
logmod3 <- glm(NACCANGI  ~ 
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod3, exponentiate = TRUE)
drug_ipw3 <- augment_columns(logmod3, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCANGI/.fitted) + ((1-NACCANGI)/(1-.fitted)))
summary(drug_ipw3$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  NACCANGI +
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw3$ipw) 
summary(drugSMw)

#4 BETA
logmod4 <- glm(NACCBETA  ~ 
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod4, exponentiate = TRUE)
drug_ipw4 <- augment_columns(logmod4, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCBETA/.fitted) + ((1-NACCBETA)/(1-.fitted)))
summary(drug_ipw4$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  NACCBETA +
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, #trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw4$ipw) 
summary(drugSMw)


#5 CCBS
logmod5 <- glm(NACCCCBS  ~ 
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod5, exponentiate = TRUE)
drug_ipw5 <- augment_columns(logmod5, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCCCBS/.fitted) + ((1-NACCCCBS)/(1-.fitted)))
summary(drug_ipw5$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  NACCCCBS +
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD,# trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw5$ipw) 
summary(drugSMw)

#6 DIUR
logmod6 <- glm(NACCDIUR  ~ 
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod6, exponentiate = TRUE)
drug_ipw6 <- augment_columns(logmod6, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCDIUR/.fitted) + ((1-NACCDIUR)/(1-.fitted)))
summary(drug_ipw6$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  NACCDIUR +
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw6$ipw) 
summary(drugSMw)

#7 HTNC
logmod7 <- glm(NACCHTNC  ~
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod7, exponentiate = TRUE)
drug_ipw7 <- augment_columns(logmod7, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCHTNC/.fitted) + ((1-NACCHTNC)/(1-.fitted)))
summary(drug_ipw7$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  NACCHTNC +
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw7$ipw) 
summary(drugSMw)

#8 VASD
logmod8 <- glm(NACCVASD  ~ 
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod8, exponentiate = TRUE)
drug_ipw8 <- augment_columns(logmod8, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NACCVASD/.fitted) + ((1-NACCVASD)/(1-.fitted)))
summary(drug_ipw8$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  NACCVASD +
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw8$ipw) 
summary(drugSMw)

## Contingency tables
sex <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SEX,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SEX, Diagnosis))))

race <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$RACE,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(RACE, Diagnosis))))

tbi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TBI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TBI, Diagnosis))))

hear <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HEAR,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HEAR, Diagnosis))))

alcohol <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ALCOHOL,
                                                           lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ALCOHOL, Diagnosis))))

diab <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DIAB,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DIAB, Diagnosis))))

hyper <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HYPERT,
                                                         lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HYPERT, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

aaas <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCAAAS,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCAAAS, Diagnosis))))

acei <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCACEI,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCACEI, Diagnosis))))

angi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCANGI,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCANGI, Diagnosis))))

htnc <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCHTNC,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCHTNC, Diagnosis))))

beta <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCBETA,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCBETA, Diagnosis))))

ccbs <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCCCBS,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCCCBS, Diagnosis))))

diur <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCDIUR,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCDIUR, Diagnosis))))

vasd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCVASD,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCVASD, Diagnosis))))

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}

write.excel(sex)
write.excel(race)
write.excel(alcohol)
write.excel(hear)
write.excel(aaas)
write.excel(acei)
write.excel(angi)
write.excel(htnc)
write.excel(beta)
write.excel(ccbs)
write.excel(diur)
write.excel(vasd)
write.excel(diab)
write.excel(hyper)
write.excel(dep)
write.excel(bmi)
write.excel(tbi)
write.excel(cvd)


#### Anxiolytics ####

rm(list = ls())
library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)
library(car)
library(coxphw)
library(survival)
library(ggfortify)
library(ggplot2)
library(survminer)
library(sjPlot)
library(broom)

File <- file.path("C:/Users/Daman Kaur/Desktop/Ulster/Research Projects/NACC drug classes analysis paper/Cox regression csv")
setwd(File)
nacc <- fread("MDallmedscox.csv")
df <- nacc

df <- df %>%
  group_by(NACCID) %>%
  mutate(check1=is.na(ANX1) & is.na(ANX2) & is.na(ANX3))
df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))
df <- df[,-(148)] 
df <- as.data.frame(df)

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(ANX1, na.rm=TRUE))
df$ANX1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$ANX1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(ANX2, na.rm=TRUE))
df$ANX2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$ANX2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(ANX3, na.rm=TRUE))
df$ANX3[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$ANX3[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df[,-(148)] 
df <- as.data.frame(df)

df$ANX[(df$ANX1==1 | df$ANX2==1 | df$ANX3==1)] <- 1
df$ANX[(df$ANX1==0 & df$ANX2==0 & df$ANX3==0)] <- 0

df <- as.data.frame(df)
df <- df[,c(1,8:10,14,26,49:88,108,125:137,144,148)]

## Find drug terms based on category
df$BARB <- apply(df, 1, function(x)as.integer(any(grep("secobarbital|butalbital|amobarbital|
                    pentobarbital|mephobarbital|butabarbital|phenobarbital|thiopental",x, ignore.case = TRUE))))

df$BENZ <- apply(df, 1, function(x)as.integer(any(grep("oxazepam|diazepam|lorazepam|alprazolam|
                          chlordiazepoxide|clonazepam|flurazepam|temazepam|triazolam|
                                            halazepam|estazolam|quazepam|clobazam",x, ignore.case = TRUE))))

df$ATH <- apply(df, 1, function(x)as.integer(any(grep("diphenhydramine|pyrilamine|hydroxyzine|
                                                  doxylamine|promethazine",x, ignore.case = TRUE))))

df$HYP <- apply(df, 1, function(x)as.integer(any(grep("chloral hydrate|zolpidem|melatonin|
                          zaleplon|eszopiclone|ramelteon|sodium oxybate|
                                  acetylcarbromal",x, ignore.case = TRUE))))

df$CARB <- apply(df, 1, function(x)as.integer(any(grep("meprobamate",x, ignore.case = TRUE))))
df$DOX <- apply(df, 1, function(x)as.integer(any(grep("doxepin",x, ignore.case = TRUE))))
df$AZA <- apply(df, 1, function(x)as.integer(any(grep("buspirone",x, ignore.case = TRUE))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

### Adjusting anx values ###

#### BARB ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NBARB = sum(BARB, na.rm=TRUE))

df$BARB[(df$NBARB==0)]<- 0
df$BARB[(df$NBARB==1)]<- 9 
df$BARB[(df$NBARB>=2)]<- 1 
df <- as.data.frame(df)

#### BENZ ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NBENZ = sum(BENZ, na.rm=TRUE))

df$BENZ[(df$NBENZ==0)]<- 0
df$BENZ[(df$NBENZ==1)]<- 9 
df$BENZ[(df$NBENZ>=2)]<- 1 
df <- as.data.frame(df)

#### ATH ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NATH = sum(ATH, na.rm=TRUE))

df$ATH[(df$NATH==0)]<- 0
df$ATH[(df$NATH==1)]<- 9 
df$ATH[(df$NATH>=2)]<- 1 
df <- as.data.frame(df)

#### HYP ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NHYP = sum(HYP, na.rm=TRUE))

df$HYP[(df$NHYP==0)]<- 0
df$HYP[(df$NHYP==1)]<- 9 
df$HYP[(df$NHYP>=2)]<- 1 
df <- as.data.frame(df)

#### CARB ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCARB = sum(CARB, na.rm=TRUE))

df$CARB[(df$NCARB==0)]<- 0
df$CARB[(df$NCARB==1)]<- 9 
df$CARB[(df$NCARB>=2)]<- 1 
df <- as.data.frame(df)

#### DOX ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NDOX = sum(DOX, na.rm=TRUE))

df$DOX[(df$NDOX==0)]<- 0
df$DOX[(df$NDOX==1)]<- 9 
df$DOX[(df$NDOX>=2)]<- 1 
df <- as.data.frame(df)

#### AZA ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NAZA = sum(AZA, na.rm=TRUE))

df$AZA[(df$NAZA==0)]<- 0
df$AZA[(df$NAZA==1)]<- 9 
df$AZA[(df$NAZA>=2)]<- 1 
file <- as.data.frame(df)

df <- file
df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
xy <- df %>%
  mutate(date_visit = as.Date(Date, "%Y/%m/%d"))
xy <- xy %>% 
  group_by(NACCID) %>%
  summarise(first_visit = min(date_visit), 
            last_visit = max(date_visit), 
            exposure = last_visit - first_visit)
xy$exposure <- as.numeric(xy$exposure)

df <- df[df$NEW_VISITNUM == 1, ] 

new <- left_join(df,xy)
df <- as.data.frame(new)

df <- subset(df, NACCAGE>=40,
             select=c(NACCID:exposure))
df <- as.data.frame(df)
df$years <- df$exposure / 365.25
#df <- df %>% filter(SEX == 1)

#df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
#df <- df[df$Diagnosis != 1, ]
#df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% replace_with_na(replace = list(BARB = 9, BENZ = 9, ATH = 9, HYP = 9,
                                            CARB = 9, DOX = 9, AZA = 9))


df.baseDHH <- df[,c(1:6,10,13:20,22:29,40)] 
df.baseDHH <- df.baseDHH[,c(1:16,20,24)] 
df.baseDHH.compl <- df.baseDHH[complete.cases(df.baseDHH), ]
lg.baseDHH <- df.baseDHH.compl

#1 MED
logmod <- glm(HYP ~ ANX +
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod, exponentiate = TRUE)
drug_ipw <- augment_columns(logmod, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (HYP/.fitted) + ((1-HYP)/(1-.fitted)))
summary(drug_ipw$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  HYP + ANX +
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw$ipw) 
summary(drugSMw)

#### contingency tables ####
sex <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SEX,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SEX, Diagnosis))))

race <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$RACE,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(RACE, Diagnosis))))

tbi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TBI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TBI, Diagnosis))))

hear <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HEAR,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HEAR, Diagnosis))))

alcohol <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ALCOHOL,
                                                           lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ALCOHOL, Diagnosis))))

diab <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DIAB,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DIAB, Diagnosis))))

hyper <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HYPERT,
                                                         lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HYPERT, Diagnosis))))

anx <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ANX,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ANX, Diagnosis))))

bar <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$BARB,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(BARB, Diagnosis))))

ben <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$BENZ,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(BENZ, Diagnosis))))

ath <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ATH,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ATH, Diagnosis))))

hyp <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HYP,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HYP, Diagnosis))))

car <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CARB,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CARB, Diagnosis))))

dox <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DOX,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DOX, Diagnosis))))

aza <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$AZA,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(AZA, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

xy <- cbind(sex,alcohol,hear,diab,anx,hyp,dep,tbi,cvd)

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}

write.excel(aza)

write.excel(xy)
write.excel(race)
write.excel(bmi)

#### Anticoagulants ####
rm(list = ls())
library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)
library(car)
library(coxphw)
library(survival)
library(ggfortify)
library(ggplot2)
library(survminer)
library(sjPlot)
library(broom)

File <- file.path("C:/Users/Daman Kaur/Desktop/Ulster/Research Projects/NACC drug classes analysis paper/Cox regression csv")
setwd(File)
nacc <- fread("MDallmedscdrcox.csv")
df <- nacc

df <- as.data.frame(df)
df <- df[,c(1,8:10,14,26,49:88,109,125:137,141)]

#heparin, VitK, Direct Xa, Thrombin IIa, Phosphodiesterase, ADP, Aspirin

## Find drug terms based on category
df$HEPA <- apply(df, 1, function(x)as.integer(any(grep("heparin|enoxaparin|dalteparin|fondaparinux",x))))

df$VITK <- apply(df, 1, function(x)as.integer(any(grep("warfarin|dicumarol",x, ignore.case = TRUE))))

df$XAIN <- apply(df, 1, function(x)as.integer(any(grep("rivaroxaban|apixaban",x, ignore.case = TRUE))))

df$THII <- apply(df, 1, function(x)as.integer(any(grep("dabigatran",x, ignore.case = TRUE))))

df$PHOS <- apply(df, 1, function(x)as.integer(any(grep("dipyridamole|cilostazol",x, ignore.case = TRUE))))
df$ADP <- apply(df, 1, function(x)as.integer(any(grep("ticlopidine|clopidogrel|prasugrel|
                                                       ticagrelor",x, ignore.case = TRUE))))
df$ASPR <- apply(df, 1, function(x)as.integer(any(grep("aspirin|asa/|aspirin-",x, ignore.case = TRUE))))

df <- df[,-(7:46)] 
df[df==-4]<-NA


#### HEPA ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NHEPA = sum(HEPA, na.rm=TRUE))

df$HEPA[(df$NHEPA==0)]<- 0
df$HEPA[(df$NHEPA==1)]<- 9 
df$HEPA[(df$NHEPA>=2)]<- 1 
df <- as.data.frame(df)

#### VITK ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NVITK = sum(VITK, na.rm=TRUE))

df$VITK[(df$NVITK==0)]<- 0
df$VITK[(df$NVITK==1)]<- 9 
df$VITK[(df$NVITK>=2)]<- 1 
df <- as.data.frame(df)

#### XAIN ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NXAIN = sum(XAIN, na.rm=TRUE))

df$XAIN[(df$NXAIN==0)]<- 0
df$XAIN[(df$NXAIN==1)]<- 9 
df$XAIN[(df$NXAIN>=2)]<- 1 
df <- as.data.frame(df)

#### THII ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NTHII = sum(THII, na.rm=TRUE))

df$THII[(df$NTHII==0)]<- 0
df$THII[(df$NTHII==1)]<- 9 
df$THII[(df$NTHII>=2)]<- 1 
df <- as.data.frame(df)

#### PHOS ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NPHOS = sum(PHOS, na.rm=TRUE))

df$PHOS[(df$NPHOS==0)]<- 0
df$PHOS[(df$NPHOS==1)]<- 9 
df$PHOS[(df$NPHOS>=2)]<- 1 
df <- as.data.frame(df)

#### ADP ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NADP = sum(ADP, na.rm=TRUE))

df$ADP[(df$NADP==0)]<- 0
df$ADP[(df$NADP==1)]<- 9 
df$ADP[(df$NADP>=2)]<- 1 
df <- as.data.frame(df)

#### ASPR ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NASPR = sum(ASPR, na.rm=TRUE))

df$ASPR[(df$NASPR==0)]<- 0
df$ASPR[(df$NASPR==1)]<- 9 
df$ASPR[(df$NASPR>=2)]<- 1 
file <- as.data.frame(df)

df <- file
df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
xy <- df %>%
  mutate(date_visit = as.Date(Date, "%Y/%m/%d"))
xy <- xy %>% 
  group_by(NACCID) %>%
  summarise(first_visit = min(date_visit), 
            last_visit = max(date_visit), 
            exposure = last_visit - first_visit)
xy$exposure <- as.numeric(xy$exposure)

df <- df[df$NEW_VISITNUM == 1, ] 

new <- left_join(df,xy)
df <- as.data.frame(new)

df <- subset(df, NACCAGE>=40,
             select=c(NACCID:exposure))
df <- as.data.frame(df)
df$years <- df$exposure / 365.25
df <- as.data.frame(df)
#df <- df %>% filter(SEX == 1)

#df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
#df <- df[df$Diagnosis != 1, ]
#df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% replace_with_na(replace = list(HEPA = 9, VITK = 9, XAIN = 9, PHOS = 9,
                                            THII = 9, ADP = 9, ASPR = 9))

df.baseDHH <- df[,c(1:6,10,13:20,22:28,39)] 
df.baseDHH <- df.baseDHH[,c(1:15,22,23)] 
df.baseDHH.compl <- df.baseDHH[complete.cases(df.baseDHH), ]
lg.baseDHH <- df.baseDHH.compl

#1 MED
logmod <- glm(ASPR ~
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod, exponentiate = TRUE)
drug_ipw <- augment_columns(logmod, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (ASPR/.fitted) + ((1-ASPR)/(1-.fitted)))
summary(drug_ipw$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  ASPR + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD,# trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw$ipw) 
summary(drugSMw)


sex <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SEX,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SEX, Diagnosis))))

race <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$RACE,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(RACE, Diagnosis))))

tbi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TBI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TBI, Diagnosis))))

hear <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HEAR,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HEAR, Diagnosis))))

alcohol <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ALCOHOL,
                                                           lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ALCOHOL, Diagnosis))))

diab <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DIAB,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DIAB, Diagnosis))))

hyper <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HYPERT,
                                                         lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HYPERT, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

hepa <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HEPA,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HEPA, Diagnosis))))

vitk <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$VITK,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(VITK, Diagnosis))))

xain <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$XAIN,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(XAIN, Diagnosis))))

thii <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$THII,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(THII, Diagnosis))))

phos <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$PHOS,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(PHOS, Diagnosis))))

adp <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ADP,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ADP, Diagnosis))))

aspr <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ASPR,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ASPR, Diagnosis))))

xy <- cbind(sex,alcohol,hear,diab,aspr,hyper,dep,tbi,cvd)

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}

write.excel(aspr)

write.excel(xy)
write.excel(race)
write.excel(bmi)

#### Antidepressants ####

rm(list = ls())
library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)
library(car)
library(coxphw)
library(survival)
library(ggfortify)
library(ggplot2)
library(survminer)
library(sjPlot)
library(broom)

File <- file.path("C:/Users/Daman Kaur/Desktop/Ulster/Research Projects/NACC drug classes analysis paper/Cox regression csv")
setwd(File)
nacc <- fread("MDallmedscdrcox.csv")
df <- nacc

df <- as.data.frame(df)
df <- df[,c(1,8:10,14,26,49:88,111,125:137,142)]

## Find drug terms based on category
df$SSRI <- apply(df, 1, function(x)as.integer(any(grep("FLUOXETINE|SERTRALINE|PAROXETINE|FLUVOXAMINE|
                                                       CITALOPRAM|ESCITALOPRAM|5-HYDROXYTRYPTOPHAN",x))))

df$SARI <- apply(df, 1, function(x)as.integer(any(grep("TRAZODONE|NEFAZODONE",x))))

df$SNRI <- apply(df, 1, function(x)as.integer(any(grep("VENLAFAXINE|DULOXETINE|MILNACIPRAN|PRISTIQ|
                                                       DESVENLAFAXINE|LEVOMILNACIPRAN",x))))

df$TCA <- apply(df, 1, function(x)as.integer(any(grep("NORTRIPTYLINE|DESIPRAMINE|AMITRIPTYLINE|DOXEPIN|
                                                      IMIPRAMINE|PROTRIPTYLINE|CLOMIPRAMINE",x))))

df$TECA <- apply(df, 1, function(x)as.integer(any(grep("AMOXAPINE|MAPROTILINE|MIRTAZAPINE",x))))
df$NDRI <- apply(df, 1, function(x)as.integer(any(grep("BUPROPION",x))))
df$SMS <- apply(df, 1, function(x)as.integer(any(grep("VILAZODONE",x))))
df$MAOB <- apply(df, 1, function(x)as.integer(any(grep("SELEGILINE",x))))
df$COMBO <- apply(df, 1, function(x)as.integer(any(grep("fluoxetine-|bupropion-|amitriptyline-",x, ignore.case = TRUE))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

#### SSRI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSSRI = sum(SSRI, na.rm=TRUE))

df$SSRI[(df$NSSRI==0)]<- 0
df$SSRI[(df$NSSRI==1)]<- 9 
df$SSRI[(df$NSSRI>=2)]<- 1 
df <- as.data.frame(df)

#### SARI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSARI = sum(SARI, na.rm=TRUE))

df$SARI[(df$NSARI==0)]<- 0
df$SARI[(df$NSARI==1)]<- 9 
df$SARI[(df$NSARI>=2)]<- 1 
df <- as.data.frame(df)

#### SNRI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSNRI = sum(SNRI, na.rm=TRUE))

df$SNRI[(df$NSNRI==0)]<- 0
df$SNRI[(df$NSNRI==1)]<- 9 
df$SNRI[(df$NSNRI>=2)]<- 1 
df <- as.data.frame(df)

#### TCA ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NTCA = sum(TCA, na.rm=TRUE))

df$TCA[(df$NTCA==0)]<- 0
df$TCA[(df$NTCA==1)]<- 9 
df$TCA[(df$NTCA>=2)]<- 1 
df <- as.data.frame(df)

#### TECA ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NTECA = sum(TECA, na.rm=TRUE))

df$TECA[(df$NTECA==0)]<- 0
df$TECA[(df$NTECA==1)]<- 9 
df$TECA[(df$NTECA>=2)]<- 1 
df <- as.data.frame(df)

#### NDRI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NNDRI = sum(NDRI, na.rm=TRUE))

df$NDRI[(df$NNDRI==0)]<- 0
df$NDRI[(df$NNDRI==1)]<- 9 
df$NDRI[(df$NNDRI>=2)]<- 1 
df <- as.data.frame(df)

#### SMS ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSMS = sum(SMS, na.rm=TRUE))

df$SMS[(df$NSMS==0)]<- 0
df$SMS[(df$NSMS==1)]<- 9 
df$SMS[(df$NSMS>=2)]<- 1 
df <- as.data.frame(df)

#### MAOB ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NMAOB = sum(MAOB, na.rm=TRUE))

df$MAOB[(df$NMAOB==0)]<- 0
df$MAOB[(df$NMAOB==1)]<- 9 
df$MAOB[(df$NMAOB>=2)]<- 1 
df <- as.data.frame(df)

#### COMBO ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCOMBO = sum(COMBO, na.rm=TRUE))

df$COMBO[(df$NCOMBO==0)]<- 0
df$COMBO[(df$NCOMBO==1)]<- 9 
df$COMBO[(df$NCOMBO>=2)]<- 1 
file <- as.data.frame(df)

df <- file
df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
xy <- df %>%
  mutate(date_visit = as.Date(Date, "%Y/%m/%d"))
xy <- xy %>% 
  group_by(NACCID) %>%
  summarise(first_visit = min(date_visit), 
            last_visit = max(date_visit), 
            exposure = last_visit - first_visit)
xy$exposure <- as.numeric(xy$exposure)

df <- df[df$NEW_VISITNUM == 1, ] 

new <- left_join(df,xy)
df <- as.data.frame(new)

df <- subset(df, NACCAGE>=40,
             select=c(NACCID:exposure))
df <- as.data.frame(df)
df$years <- df$exposure / 365.25

df <- as.data.frame(df)
#df <- df %>% filter(SEX == 1)

#df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
#df <- df[df$Diagnosis != 1, ]
#df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% replace_with_na(replace = list(SSRI = 9, SARI = 9, SNRI = 9, TCA = 9,
                                            TECA = 9, NDRI = 9, SMS = 9, MAOB = 9, COMBO = 9))


df.baseDHH <- df[,c(1:6,10,13:20,22:30,43)] 
df.baseDHH <- df.baseDHH[,c(1:15,24,25)] 

df.baseDHH.compl <- df.baseDHH[complete.cases(df.baseDHH), ]
lg.baseDHH <- df.baseDHH.compl

#1 SSRI
logmod <- glm(SSRI ~
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod, exponentiate = TRUE)
drug_ipw <- augment_columns(logmod, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (SSRI/.fitted) + ((1-SSRI)/(1-.fitted)))
summary(drug_ipw$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  SSRI + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw$ipw) 
summary(drugSMw)

#2 SARI
logmod2 <- glm(SARI ~
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod2, exponentiate = TRUE)
drug_ipw2 <- augment_columns(logmod2, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (SARI/.fitted) + ((1-SARI)/(1-.fitted)))
summary(drug_ipw2$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  SARI + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw2$ipw) 
summary(drugSMw)

#3 SNRI
logmod3 <- glm(SNRI ~
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod3, exponentiate = TRUE)
drug_ipw3 <- augment_columns(logmod3, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (SNRI/.fitted) + ((1-SNRI)/(1-.fitted)))
summary(drug_ipw3$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  SNRI + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw3$ipw) 
summary(drugSMw)

#4 TCA
logmod4 <- glm(TCA ~
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod4, exponentiate = TRUE)
drug_ipw4 <- augment_columns(logmod4, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (TCA/.fitted) + ((1-TCA)/(1-.fitted)))
summary(drug_ipw4$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  TCA + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw4$ipw) 
summary(drugSMw)

#5 NDRI
logmod5 <- glm(NDRI ~
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod5, exponentiate = TRUE)
drug_ipw5 <- augment_columns(logmod5, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (NDRI/.fitted) + ((1-NDRI)/(1-.fitted)))
summary(drug_ipw5$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  NDRI + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, #trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw5$ipw) 
summary(drugSMw)

### Contingency tables
sex <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SEX,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SEX, Diagnosis))))

race <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$RACE,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(RACE, Diagnosis))))

tbi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TBI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TBI, Diagnosis))))

hear <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HEAR,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HEAR, Diagnosis))))

alcohol <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ALCOHOL,
                                                           lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ALCOHOL, Diagnosis))))

diab <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DIAB,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DIAB, Diagnosis))))

hyper <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HYPERT,
                                                         lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HYPERT, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

ssr <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SSRI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SSRI, Diagnosis))))

sar <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SARI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SARI, Diagnosis))))

snr <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SNRI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SNRI, Diagnosis))))

tca <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TCA,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TCA, Diagnosis))))

teca <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TECA,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TECA, Diagnosis))))

ndr <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NDRI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NDRI, Diagnosis))))

sms <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SMS,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SMS, Diagnosis))))

mao <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$MAOB,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(MAOB, Diagnosis))))

com <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$COMBO,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(COMBO, Diagnosis))))


xy <- cbind(sex,alcohol,hear,diab,ndr,hyper,dep,tbi,cvd)

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}

write.excel(com)

write.excel(sms)
write.excel(race)
write.excel(bmi)

#### Antidiabetic agents ####

rm(list = ls())
library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)
library(car)
library(coxphw)
library(survival)
library(ggfortify)
library(ggplot2)
library(survminer)
library(sjPlot)
library(broom)

File <- file.path("C:/Users/Daman Kaur/Desktop/Ulster/Research Projects/NACC drug classes analysis paper/Cox regression csv")
setwd(File)
nacc <- fread("MDallmedscdrcox.csv")
df <- nacc

df <- as.data.frame(df)
df <- df[,c(1,8:10,14,26,49:88,118,125:137,145)]

## Find drug terms based on category
df$SULF <- apply(df, 1, function(x)as.integer(any(grep("CHLORPROPAMIDE|ACETOHEXAMIDE|GLIPIZIDE|
                    GLYBURIDE|TOLAZAMIDE|TOLBUTAMIDE|GLIMEPIRIDE",x))))

df$THI <- apply(df, 1, function(x)as.integer(any(grep("TROGLITAZONE|ROSIGLITAZONE|PIOGLITAZONE",x))))

df$COM <- apply(df, 1, function(x)as.integer(any(grep("metformin-|-metformin|glimepiride-|
                       -sitagliptin|-linagliptin",x, ignore.case = TRUE))))

df$SGLT <- apply(df, 1, function(x)as.integer(any(grep("EMPAGLIFLOZIN|CANAGLIFLOZIN|DAPAGLIFLOZIN",x))))

df$DPP <- apply(df, 1, function(x)as.integer(any(grep("SITAGLIPTIN|SAXAGLIPTIN|LINAGLIPTIN|ALOGLIPTIN",
                                                      x))))

df$INS <- apply(df, 1, function(x)as.integer(any(grep("insulin",x, ignore.case = TRUE))))
df$MET <- apply(df, 1, function(x)as.integer(any(grep("METFORMIN",x))))
df$AGI <- apply(df, 1, function(x)as.integer(any(grep("ACARBOSE|MIGLITOL",x))))
df$IMM <- apply(df, 1, function(x)as.integer(any(grep("EXENATIDE|LIRAGLUTIDE|DULAGLUTIDE",x))))
df$MEG <- apply(df, 1, function(x)as.integer(any(grep("REPAGLINIDE|NATEGLINIDE",x))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

#### SULF ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSULF = sum(SULF, na.rm=TRUE))

df$SULF[(df$NSULF==0)]<- 0
df$SULF[(df$NSULF==1)]<- 9 
df$SULF[(df$NSULF>=2)]<- 1 
df <- as.data.frame(df)

#### THI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NTHI = sum(THI, na.rm=TRUE))

df$THI[(df$NTHI==0)]<- 0
df$THI[(df$NTHI==1)]<- 9 
df$THI[(df$NTHI>=2)]<- 1 
df <- as.data.frame(df)

#### COM ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCOM = sum(COM, na.rm=TRUE))

df$COM[(df$NCOM==0)]<- 0
df$COM[(df$NCOM==1)]<- 9 
df$COM[(df$NCOM>=2)]<- 1 
df <- as.data.frame(df)

#### SGLT ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSGLT = sum(SGLT, na.rm=TRUE))

df$SGLT[(df$NSGLT==0)]<- 0
df$SGLT[(df$NSGLT==1)]<- 9 
df$SGLT[(df$NSGLT>=2)]<- 1 
df <- as.data.frame(df)

#### DPP ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NDPP = sum(DPP, na.rm=TRUE))

df$DPP[(df$NDPP==0)]<- 0
df$DPP[(df$NDPP==1)]<- 9 
df$DPP[(df$NDPP>=2)]<- 1 
df <- as.data.frame(df)

#### INS ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NINS = sum(INS, na.rm=TRUE))

df$INS[(df$NINS==0)]<- 0
df$INS[(df$NINS==1)]<- 9 
df$INS[(df$NINS>=2)]<- 1 
df <- as.data.frame(df)

#### MET ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NMET = sum(MET, na.rm=TRUE))

df$MET[(df$NMET==0)]<- 0
df$MET[(df$NMET==1)]<- 9 
df$MET[(df$NMET>=2)]<- 1 
df <- as.data.frame(df)

#### AGI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NAGI = sum(AGI, na.rm=TRUE))

df$AGI[(df$NAGI==0)]<- 0
df$AGI[(df$NAGI==1)]<- 9 
df$AGI[(df$NAGI>=2)]<- 1 
df <- as.data.frame(df)

#### IMM ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NIMM = sum(IMM, na.rm=TRUE))

df$IMM[(df$NIMM==0)]<- 0
df$IMM[(df$NIMM==1)]<- 9 
df$IMM[(df$NIMM>=2)]<- 1 
df <- as.data.frame(df)

#### MEG ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NMEG = sum(MEG, na.rm=TRUE))

df$MEG[(df$NMEG==0)]<- 0
df$MEG[(df$NMEG==1)]<- 9 
df$MEG[(df$NMEG>=2)]<- 1 

file <- as.data.frame(df)

df <- file
df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
xy <- df %>%
  mutate(date_visit = as.Date(Date, "%Y/%m/%d"))
xy <- xy %>% 
  group_by(NACCID) %>%
  summarise(first_visit = min(date_visit), 
            last_visit = max(date_visit), 
            exposure = last_visit - first_visit)
xy$exposure <- as.numeric(xy$exposure)

df <- df[df$NEW_VISITNUM == 1, ] 

new <- left_join(df,xy)
df <- as.data.frame(new)

df <- subset(df, NACCAGE>=40,
             select=c(NACCID:exposure))
df <- as.data.frame(df)
df$years <- df$exposure / 365.25
df <- as.data.frame(df)
#df <- df %>% filter(SEX == 0)

#df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
#df <- df[df$Diagnosis != 1, ]
#df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% replace_with_na(replace = list(SULF = 9, THI = 9, COM = 9, SGLT = 9,
                                            DPP = 9, INS = 9, MET = 9, AGI = 9,
                                            IMM = 9, MEG = 9))


df.baseDHH <- df[,c(1:6,10,13:20,22:31,45)] 
df.baseDHH <- df.baseDHH[,c(1:15,25,26)] 

df.baseDHH.compl <- df.baseDHH[complete.cases(df.baseDHH), ]
lg.baseDHH <- df.baseDHH.compl

#1 SULF
logmod <- glm(SULF ~
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod, exponentiate = TRUE)
drug_ipw <- augment_columns(logmod, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (SULF/.fitted) + ((1-SULF)/(1-.fitted)))
summary(drug_ipw$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  SULF + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw$ipw) 
summary(drugSMw)

#2 MET
logmod2 <- glm(MET ~
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod2, exponentiate = TRUE)
drug_ipw2 <- augment_columns(logmod2, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (MET/.fitted) + ((1-MET)/(1-.fitted)))
summary(drug_ipw2$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  MET + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw2$ipw) 
summary(drugSMw)

### Contingency tables
sex <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SEX,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SEX, Diagnosis))))

race <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$RACE,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(RACE, Diagnosis))))

tbi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TBI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TBI, Diagnosis))))

hear <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HEAR,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HEAR, Diagnosis))))

alcohol <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ALCOHOL,
                                                           lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ALCOHOL, Diagnosis))))

diab <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DIAB,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DIAB, Diagnosis))))

hyper <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HYPERT,
                                                         lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HYPERT, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

sul <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SULF,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SULF, Diagnosis))))

thi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$THI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(THI, Diagnosis))))

com <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$COM,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(COM, Diagnosis))))

sglt <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SGLT,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SGLT, Diagnosis))))

dpp <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DPP,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DPP, Diagnosis))))

ins <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$INS,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(INS, Diagnosis))))

met <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$MET,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(MET, Diagnosis))))

agi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$AGI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(AGI, Diagnosis))))

imm <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$IMM,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(IMM, Diagnosis))))

meg <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$MEG,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(MEG, Diagnosis))))

write.excel(meg)

xy <- cbind(sex, alcohol, hear,diab,met,hyper,dep,tbi,cvd)

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}

write.excel(xy)
write.excel(race)
write.excel(bmi)

#### Hypolipidemic drugs ####

rm(list = ls())
library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)
library(car)
library(coxphw)
library(survival)
library(ggfortify)
library(ggplot2)
library(survminer)
library(sjPlot)
library(broom)

File <- file.path("C:/Users/Daman Kaur/Desktop/Ulster/Research Projects/NACC drug classes analysis paper/Cox regression csv")
setwd(File)
nacc <- fread("MDallmedscox.csv")
df <- nacc

df <- df %>%
  group_by(NACCID) %>%
  mutate(check1=is.na(HYPERC1) & is.na(HYPERC2))
df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))
df <- df[,-(148)] 
df <- as.data.frame(df)

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(HYPERC1, na.rm=TRUE))
df$HYPERC1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$HYPERC1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(HYPERC2, na.rm=TRUE))
df$HYPERC2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$HYPERC2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df[,-(148)] 
df <- as.data.frame(df)

df$HYPERC[(df$HYPERC1==1 | df$HYPERC2==1)] <- 1
df$HYPERC[(df$HYPERC1==0 & df$HYPERC2==0)] <- 0

df <- as.data.frame(df)
df <- df[,c(1,8:10,14,26,49:88,121,125:137,139,148)]

## Find drug terms based on category
df$STAT <- apply(df, 1, function(x)as.integer(any(grep("LOVASTATIN|PRAVASTATIN|SIMVASTATIN|
            FLUVASTATIN|ATORVASTATIN|CERIVASTATIN|ROSUVASTATIN|PITAVASTATIN|RED YEAST RICE",x))))

df$NIAC <- apply(df, 1, function(x)as.integer(any(grep("niacin|niacinamide",x, ignore.case = TRUE))))

df$FIBR <- apply(df, 1, function(x)as.integer(any(grep("gemfibrozil|fenofibrate|
                                                       fenofibric acid",x, ignore.case = TRUE))))

df$BILE <- apply(df, 1, function(x)as.integer(any(grep("colestipol|cholestyramine",x, ignore.case = TRUE))))

df$CAI <- apply(df, 1, function(x)as.integer(any(grep("ezetimibe",x, ignore.case = TRUE))))
df$COMB <- apply(df, 1, function(x)as.integer(any(grep("lovastatin-niacin|aspirin-pravastatin|ezetimibe-simvastatin|
                                                       amlodipine-atorvastatin|niacin-simvastatin",x, ignore.case = TRUE))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

#### STAT ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSTAT = sum(STAT, na.rm=TRUE))

df$STAT[(df$NSTAT==0)]<- 0
df$STAT[(df$NSTAT==1)]<- 9 
df$STAT[(df$NSTAT>=2)]<- 1 
df <- as.data.frame(df)

#### NIAC ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NNIAC = sum(NIAC, na.rm=TRUE))

df$NIAC[(df$NNIAC==0)]<- 0
df$NIAC[(df$NNIAC==1)]<- 9 
df$NIAC[(df$NNIAC>=2)]<- 1 
df <- as.data.frame(df)

#### FIBR ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NFIBR = sum(FIBR, na.rm=TRUE))

df$FIBR[(df$NFIBR==0)]<- 0
df$FIBR[(df$NFIBR==1)]<- 9 
df$FIBR[(df$NFIBR>=2)]<- 1 
df <- as.data.frame(df)

#### BILE ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NBILE = sum(BILE, na.rm=TRUE))

df$BILE[(df$NBILE==0)]<- 0
df$BILE[(df$NBILE==1)]<- 9 
df$BILE[(df$NBILE>=2)]<- 1 
df <- as.data.frame(df)

#### CAI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCAI = sum(CAI, na.rm=TRUE))

df$CAI[(df$NCAI==0)]<- 0
df$CAI[(df$NCAI==1)]<- 9 
df$CAI[(df$NCAI>=2)]<- 1 
df <- as.data.frame(df)

#### COMB ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCOMB = sum(COMB, na.rm=TRUE))

df$COMB[(df$NCOMB==0)]<- 0
df$COMB[(df$NCOMB==1)]<- 9 
df$COMB[(df$NCOMB>=2)]<- 1 
file <- as.data.frame(df)

df <- file
df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
xy <- df %>%
  mutate(date_visit = as.Date(Date, "%Y/%m/%d"))
xy <- xy %>% 
  group_by(NACCID) %>%
  summarise(first_visit = min(date_visit), 
            last_visit = max(date_visit), 
            exposure = last_visit - first_visit)
xy$exposure <- as.numeric(xy$exposure)

df <- df[df$NEW_VISITNUM == 1, ] 

new <- left_join(df,xy)
df <- as.data.frame(new)

df <- subset(df, NACCAGE>=40,
             select=c(NACCID:exposure))
df <- as.data.frame(df)
df$years <- df$exposure / 365.25

#df <- df %>% filter(SEX == 1)

#df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
#df <- df[df$Diagnosis != 1, ]
#df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% replace_with_na(replace = list(STAT = 9, NIAC = 9, FIBR = 9, BILE = 9,
                                            CAI = 9, COMB = 9))


df.baseDHH <- df[,c(1:6,10,13:20,22:28,38)] 
df.baseDHH <- df.baseDHH[,c(1:16,21,23)] 

df.baseDHH.compl <- df.baseDHH[complete.cases(df.baseDHH), ]
lg.baseDHH <- df.baseDHH.compl

#1 STAT
logmod <- glm(STAT ~ HYPERC +
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod, exponentiate = TRUE)
drug_ipw <- augment_columns(logmod, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (STAT/.fitted) + ((1-STAT)/(1-.fitted)))
summary(drug_ipw$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  STAT + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw$ipw) 
summary(drugSMw)

#2 CAI
logmod2 <- glm(CAI ~ HYPERC +
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod2, exponentiate = TRUE)
drug_ipw2 <- augment_columns(logmod2, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (CAI/.fitted) + ((1-CAI)/(1-.fitted)))
summary(drug_ipw2$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  CAI + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw2$ipw) 
summary(drugSMw)

### Contingency tables
sex <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SEX,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SEX, Diagnosis))))

race <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$RACE,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(RACE, Diagnosis))))

tbi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TBI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TBI, Diagnosis))))

hear <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HEAR,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HEAR, Diagnosis))))

alcohol <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ALCOHOL,
                                                           lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ALCOHOL, Diagnosis))))

diab <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DIAB,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DIAB, Diagnosis))))

hyper <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HYPERT,
                                                         lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HYPERT, Diagnosis))))

hyperc <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HYPERC,
                                                          lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HYPERC, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

stat <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$STAT,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(STAT, Diagnosis))))

niac <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NIAC,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NIAC, Diagnosis))))

fib <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$FIBR,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(FIBR, Diagnosis))))

bile <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$BILE,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(BILE, Diagnosis))))

cai <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CAI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CAI, Diagnosis))))

comb <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$COMB,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(COMB, Diagnosis))))

write.excel(comb)

xy <- cbind(sex,alcohol,hear,diab,hyperc,cai,hyper,dep,tbi,cvd)

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}

write.excel(xy)
write.excel(race)
write.excel(bmi)

#### NSAIDs ####

rm(list = ls())
library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)
library(car)
library(coxphw)
library(survival)
library(ggfortify)
library(ggplot2)
library(survminer)
library(sjPlot)
library(broom)

File <- file.path("C:/Users/Daman Kaur/Desktop/Ulster/Research Projects/NACC drug classes analysis paper/Cox regression csv")
setwd(File)
nacc <- fread("MDallmedscdrcox.csv")
df <- nacc

df <- as.data.frame(df)
df <- df[,c(1,8:10,14,26,49:88,122,125:137,140)]

df$SALC <- apply(df, 1, function(x)as.integer(any(grep("aspirin|asa/|diflunisal|choline sal|salsalate|
                                        sodium sal|magnesium sal|methyl sal|-salicylic|/salicylic|
                                        salicylic acid topical|salicylic acid-",x, ignore.case = TRUE))))

df$PROP <- apply(df, 1, function(x)as.integer(any(grep("IBUPROFEN|NAPROXEN|FENOPROFEN|KETOPROFEN|
                                                       FLURBIPROFEN|OXAPROZIN",x, ignore.case = TRUE))))

df$OXI <- apply(df, 1, function(x)as.integer(any(grep("oxicam",x, ignore.case = TRUE))))

df$COX <- apply(df, 1, function(x)as.integer(any(grep("coxib",x, ignore.case = TRUE))))

df$ACD <- apply(df, 1, function(x)as.integer(any(grep("SULINDAC|TOLMETIN|bromfenac|diclofenac|
                                indomethacin|etodolac|ketorolac|nabumetone",x, ignore.case = TRUE))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

#### SALC ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSALC = sum(SALC, na.rm=TRUE))

df$SALC[(df$NSALC==0)]<- 0
df$SALC[(df$NSALC==1)]<- 9 
df$SALC[(df$NSALC>=2)]<- 1 
df <- as.data.frame(df)

#### PROP ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NPROP = sum(PROP, na.rm=TRUE))

df$PROP[(df$NPROP==0)]<- 0
df$PROP[(df$NPROP==1)]<- 9 
df$PROP[(df$NPROP>=2)]<- 1 
df <- as.data.frame(df)

#### OXI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NOXI = sum(OXI, na.rm=TRUE))

df$OXI[(df$NOXI==0)]<- 0
df$OXI[(df$NOXI==1)]<- 9 
df$OXI[(df$NOXI>=2)]<- 1 
df <- as.data.frame(df)

#### COX ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCOX = sum(COX, na.rm=TRUE))

df$COX[(df$NCOX==0)]<- 0
df$COX[(df$NCOX==1)]<- 9 
df$COX[(df$NCOX>=2)]<- 1 
df <- as.data.frame(df)

#### ACD ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NACD = sum(ACD, na.rm=TRUE))

df$ACD[(df$NACD==0)]<- 0
df$ACD[(df$NACD==1)]<- 9 
df$ACD[(df$NACD>=2)]<- 1 
file <- as.data.frame(df)

df <- file
df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
xy <- df %>%
  mutate(date_visit = as.Date(Date, "%Y/%m/%d"))
xy <- xy %>% 
  group_by(NACCID) %>%
  summarise(first_visit = min(date_visit), 
            last_visit = max(date_visit), 
            exposure = last_visit - first_visit)
xy$exposure <- as.numeric(xy$exposure)

df <- df[df$NEW_VISITNUM == 1, ] 

new <- left_join(df,xy)
df <- as.data.frame(new)

df <- subset(df, NACCAGE>=40,
             select=c(NACCID:exposure))
df <- as.data.frame(df)
df$years <- df$exposure / 365.25
#df <- df %>% filter(SEX == 1)

#df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
#df <- df[df$Diagnosis != 1, ]
#df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% replace_with_na(replace = list(SALC = 9, PROP = 9, OXI = 9, COX = 9,
                                            ACD = 9))

df.baseDHH <- df[,c(1:6,10,13:20,22:26,35)] 
df.baseDHH <- df.baseDHH[,c(1:15,20,21)] 

df.baseDHH.compl <- df.baseDHH[complete.cases(df.baseDHH), ]
lg.baseDHH <- df.baseDHH.compl

#1 SALC
logmod <- glm(SALC ~
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod, exponentiate = TRUE)
drug_ipw <- augment_columns(logmod, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (SALC/.fitted) + ((1-SALC)/(1-.fitted)))
summary(drug_ipw$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  SALC + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD,# trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw$ipw) 
summary(drugSMw)

#2 PROP
logmod2 <- glm(PROP ~
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod2, exponentiate = TRUE)
drug_ipw2 <- augment_columns(logmod2, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (PROP/.fitted) + ((1-PROP)/(1-.fitted)))
summary(drug_ipw2$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  PROP + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, #trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw2$ipw) 
summary(drugSMw)

### Contingency tables
sex <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SEX,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SEX, Diagnosis))))

race <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$RACE,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(RACE, Diagnosis))))

tbi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TBI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TBI, Diagnosis))))

hear <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HEAR,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HEAR, Diagnosis))))

alcohol <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ALCOHOL,
                                                           lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ALCOHOL, Diagnosis))))

diab <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DIAB,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DIAB, Diagnosis))))

hyper <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HYPERT,
                                                         lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HYPERT, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

sal <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SALC,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SALC, Diagnosis))))

pro <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$PROP,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(PROP, Diagnosis))))

oxi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$OXI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(OXI, Diagnosis))))

cox <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$COX,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(COX, Diagnosis))))

acd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ACD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ACD, Diagnosis))))


xy <- cbind(sex,alcohol,hear,diab,pro,hyper,dep,tbi,cvd)

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}

write.excel(acd)

write.excel(xy)
write.excel(race)
write.excel(bmi)

#### Parkinson's drugs ####

rm(list = ls())
library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)
library(car)
library(coxphw)
library(survival)
library(ggfortify)
library(ggplot2)
library(survminer)
library(sjPlot)
library(broom)

File <- file.path("C:/Users/Daman Kaur/Desktop/Ulster/Research Projects/NACC drug classes analysis paper/Cox regression csv")
setwd(File)
nacc <- fread("MDallmedscdrcox.csv")
df <- nacc

df <- df %>%
  group_by(NACCID) %>%
  mutate(check1=is.na(PARK1) & is.na(PARK2))
df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))
df <- df[,-(148)] 
df <- as.data.frame(df)

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(PARK1, na.rm=TRUE))
df$PARK1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$PARK1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(PARK2, na.rm=TRUE))
df$PARK2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$PARK2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df[,-(148)] 
df <- as.data.frame(df)

df$PARK[(df$PARK2==1 | df$PARK2==1)] <- 1
df$PARK[(df$PARK2==0 & df$PARK2==0)] <- 0

df <- as.data.frame(df)
df <- df[,c(1,8:10,14,26,49:88,123,125:137,147,148)]

## Find drug terms based on category
df$LEV <- apply(df, 1, function(x)as.integer(any(grep("LEVODOPA",x, ignore.case = TRUE))))

df$DORA <- apply(df, 1, function(x)as.integer(any(grep("bromocriptine|pergolide|cabergoline|
                  apomorphine|pramipexole|ropinirole|rotigotine",x, ignore.case = TRUE))))

df$ACH <- apply(df, 1, function(x)as.integer(any(grep("benzotropine|procyclidine|biperiden|
                                           trihexylphenidyl|diphenhydramine",x, ignore.case = TRUE))))

df$COM <- apply(df, 1, function(x)as.integer(any(grep("TOLCAPONE|ENTACAPONE",x))))
df$MAO <- apply(df, 1, function(x)as.integer(any(grep("selegiline|rasagiline",x, ignore.case = TRUE))))
df$AMA <- apply(df, 1, function(x)as.integer(any(grep("amantidine",x, ignore.case = TRUE))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

#### LEV ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NLEV = sum(LEV, na.rm=TRUE))

df$LEV[(df$NLEV==0)]<- 0
df$LEV[(df$NLEV==1)]<- 9 
df$LEV[(df$NLEV>=2)]<- 1 
df <- as.data.frame(df)

#### DORA ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NDORA = sum(DORA, na.rm=TRUE))

df$DORA[(df$NDORA==0)]<- 0
df$DORA[(df$NDORA==1)]<- 9 
df$DORA[(df$NDORA>=2)]<- 1 
df <- as.data.frame(df)

#### ACH ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NACH = sum(ACH, na.rm=TRUE))

df$ACH[(df$NACH==0)]<- 0
df$ACH[(df$NACH==1)]<- 9 
df$ACH[(df$NACH>=2)]<- 1 
df <- as.data.frame(df)

#### MAO ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NMAO = sum(MAO, na.rm=TRUE))

df$MAO[(df$NMAO==0)]<- 0
df$MAO[(df$NMAO==1)]<- 9 
df$MAO[(df$NMAO>=2)]<- 1 
df <- as.data.frame(df)

#### AMA ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NAMA = sum(AMA, na.rm=TRUE))

df$AMA[(df$NAMA==0)]<- 0
df$AMA[(df$NAMA==1)]<- 9 
df$AMA[(df$NAMA>=2)]<- 1 
df <- as.data.frame(df)

#### COM ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCOM = sum(COM, na.rm=TRUE))

df$COM[(df$NCOM==0)]<- 0
df$COM[(df$NCOM==1)]<- 9 
df$COM[(df$NCOM>=2)]<- 1 
file <- as.data.frame(df)

df <- file
df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
xy <- df %>%
  mutate(date_visit = as.Date(Date, "%Y/%m/%d"))
xy <- xy %>% 
  group_by(NACCID) %>%
  summarise(first_visit = min(date_visit), 
            last_visit = max(date_visit), 
            exposure = last_visit - first_visit)
xy$exposure <- as.numeric(xy$exposure)

df <- df[df$NEW_VISITNUM == 1, ] 

new <- left_join(df,xy)
df <- as.data.frame(new)

df <- subset(df, NACCAGE>=40,
             select=c(NACCID:exposure))
df <- as.data.frame(df)
df$years <- df$exposure / 365.25

#df <- df %>% filter(SEX == 1)

#df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
#df <- df[df$Diagnosis != 1, ]
#df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% replace_with_na(replace = list(LEV = 9, DORA = 9, ACH = 9, COM = 9,
                                            MAO = 9, AMA = 9))

df.baseDHH <- df[,c(1:6,10,13:20,22:28,38)] 
df.baseDHH <- df.baseDHH[,c(1:16,17,23)] 

df.baseDHH.compl <- df.baseDHH[complete.cases(df.baseDHH), ]
lg.baseDHH <- df.baseDHH.compl

#1 LEV
logmod <- glm(LEV ~ PARK +
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod, exponentiate = TRUE)
drug_ipw <- augment_columns(logmod, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (LEV/.fitted) + ((1-LEV)/(1-.fitted)))
summary(drug_ipw$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  LEV + PARK +
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw$ipw) 
summary(drugSMw)

#2 DORA
logmod2 <- glm(DORA ~ PARK +
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod2, exponentiate = TRUE)
drug_ipw2 <- augment_columns(logmod2, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (DORA/.fitted) + ((1-DORA)/(1-.fitted)))
summary(drug_ipw2$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  DORA + PARK +
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw2$ipw) 
summary(drugSMw)

### Contingency tables
sex <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SEX,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SEX, Diagnosis))))

race <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$RACE,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(RACE, Diagnosis))))

tbi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TBI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TBI, Diagnosis))))

hear <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HEAR,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HEAR, Diagnosis))))

alcohol <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ALCOHOL,
                                                           lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ALCOHOL, Diagnosis))))

diab <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DIAB,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DIAB, Diagnosis))))

hyper <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HYPERT,
                                                         lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HYPERT, Diagnosis))))

park <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$PARK,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(PARK, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

lev <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$LEV,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(LEV, Diagnosis))))

dora <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DORA,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DORA, Diagnosis))))

ach <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ACH,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ACH, Diagnosis))))

com <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$COM,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(COM, Diagnosis))))

mao <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$MAO,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(MAO, Diagnosis))))

ama <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$AMA,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(AMA, Diagnosis))))

write.excel(ama)

xy <- cbind(sex,alcohol,hear,diab,park,lev,hyper,dep,tbi,cvd)

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}

write.excel(xy)
write.excel(race)
write.excel(bmi)

#### Alzheimer drugs ####

rm(list = ls())
library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)
library(car)
library(coxphw)
library(survival)
library(ggfortify)
library(ggplot2)
library(survminer)
library(sjPlot)
library(broom)

File <- file.path("C:/Users/Daman Kaur/Desktop/Ulster/Research Projects/NACC drug classes analysis paper/Cox regression csv")
setwd(File)
nacc <- fread("MDallmedscdrcox.csv")
df <- nacc

df <- as.data.frame(df)
df <- df[,c(1,8:10,14,26,49:88,112,125:137,146)]

## Find drug terms based on category
df$DON <- apply(df, 1, function(x)as.integer(any(grep("DONEPEZIL",x))))
df$MEM <- apply(df, 1, function(x)as.integer(any(grep("MEMANTINE",x))))
df$GAL <- apply(df, 1, function(x)as.integer(any(grep("GALANTAMINE",x, ignore.case = TRUE))))
df$RIV <- apply(df, 1, function(x)as.integer(any(grep("RIVASTIGMINE",x, ignore.case = TRUE))))
df$TAC <- apply(df, 1, function(x)as.integer(any(grep("TACRINE",x, ignore.case = TRUE))))
df$COM <- apply(df, 1, function(x)as.integer(any(grep("donepezil-memantine|donepezil and memantine",x, 
                                                      ignore.case = TRUE))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

#### DON ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NDON = sum(DON, na.rm=TRUE))

df$DON[(df$NDON==0)]<- 0
df$DON[(df$NDON==1)]<- 9 
df$DON[(df$NDON>=2)]<- 1 
df <- as.data.frame(df)

#### MEM ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NMEM = sum(MEM, na.rm=TRUE))

df$MEM[(df$NMEM==0)]<- 0
df$MEM[(df$NMEM==1)]<- 9 
df$MEM[(df$NMEM>=2)]<- 1 
df <- as.data.frame(df)

#### GAL ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NGAL = sum(GAL, na.rm=TRUE))

df$GAL[(df$NGAL==0)]<- 0
df$GAL[(df$NGAL==1)]<- 9 
df$GAL[(df$NGAL>=2)]<- 1 
df <- as.data.frame(df)

#### RIV ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NRIV = sum(RIV, na.rm=TRUE))

df$RIV[(df$NRIV==0)]<- 0
df$RIV[(df$NRIV==1)]<- 9 
df$RIV[(df$NRIV>=2)]<- 1 
df <- as.data.frame(df)

#### TAC ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NTAC = sum(TAC, na.rm=TRUE))

df$TAC[(df$NTAC==0)]<- 0
df$TAC[(df$NTAC==1)]<- 9 
df$TAC[(df$NTAC>=2)]<- 1 
df <- as.data.frame(df)

#### COM ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCOM = sum(COM, na.rm=TRUE))

df$COM[(df$NCOM==0)]<- 0
df$COM[(df$NCOM==1)]<- 9 
df$COM[(df$NCOM>=2)]<- 1 

file <- as.data.frame(df)

df <- as.data.frame(file) 
df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
xy <- df %>%
  mutate(date_visit = as.Date(Date, "%Y/%m/%d"))
xy <- xy %>% 
  group_by(NACCID) %>%
  summarise(first_visit = min(date_visit), 
            last_visit = max(date_visit), 
            exposure = last_visit - first_visit)
xy$exposure <- as.numeric(xy$exposure)

df <- df[df$NEW_VISITNUM == 1, ] 

new <- left_join(df,xy)
df <- as.data.frame(new)

df <- subset(df, NACCAGE>=40,
             select=c(NACCID:exposure))
df <- as.data.frame(df)
df$years <- df$exposure / 365.25

#df <- df %>% filter(SEX == 1)

#df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
#df <- df[df$Diagnosis != 1, ]
#df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% replace_with_na(replace = list(DON = 9, MEM = 9, COM = 9, TAC = 9,
                                            GAL = 9, RIV = 9))

df.baseDHH <- df[,c(1:6,10,13:20,22:27,37)] 
df.baseDHH <- df.baseDHH[,c(1:15,21,22)] 

df.baseDHH.compl <- df.baseDHH[complete.cases(df.baseDHH), ]
lg.baseDHH <- df.baseDHH.compl

#1 DON
logmod <- glm(DON ~
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod, exponentiate = TRUE)
drug_ipw <- augment_columns(logmod, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (DON/.fitted) + ((1-DON)/(1-.fitted)))
summary(drug_ipw$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  DON + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD,# trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw$ipw) 
summary(drugSMw)

#2 MEM
logmod2 <- glm(MEM ~
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod2, exponentiate = TRUE)
drug_ipw2 <- augment_columns(logmod2, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (MEM/.fitted) + ((1-MEM)/(1-.fitted)))
summary(drug_ipw2$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  MEM + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw2$ipw) 
summary(drugSMw)

#3 GAL
logmod3 <- glm(GAL ~
                SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                SMAX + EDMAX + ABMI + TBI + CVD, 
              data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod3, exponentiate = TRUE)
drug_ipw3 <- augment_columns(logmod3, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (GAL/.fitted) + ((1-GAL)/(1-.fitted)))
summary(drug_ipw3$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  GAL + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw3$ipw) 
summary(drugSMw)

#4 RIV
logmod4 <- glm(RIV ~
                 SEX + NACCAGE + RACE + ALCOHOL + HEAR + DIAB + HYPERT + DEPRSN + 
                 SMAX + EDMAX + ABMI + TBI + CVD, 
               data = lg.baseDHH, family =  binomial(link = "logit")) 
tidy(logmod4, exponentiate = TRUE)
drug_ipw4 <- augment_columns(logmod4, lg.baseDHH, type.predict = "response") %>%
  mutate(ipw = (RIV/.fitted) + ((1-RIV)/(1-.fitted)))
summary(drug_ipw4$ipw)   ### TRUNCATE ONLY IF NECESSARY

drugSMw <- coxphw(Surv(years, Diagnosis) ~  RIV + 
                    SEX + NACCAGE + RACE + ALCOHOL + HEAR + HYPERT + DIAB + DEPRSN + 
                    SMAX + EDMAX + ABMI + TBI + CVD, trunc.weights = 0.95, normalize=TRUE,
                  data = lg.baseDHH, template = "AHR", caseweights = drug_ipw4$ipw) 
summary(drugSMw)

### Contingency tables
sex <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SEX,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SEX, Diagnosis))))

race <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$RACE,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(RACE, Diagnosis))))

tbi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TBI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TBI, Diagnosis))))

hear <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HEAR,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HEAR, Diagnosis))))

alcohol <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ALCOHOL,
                                                           lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ALCOHOL, Diagnosis))))

diab <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DIAB,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DIAB, Diagnosis))))

hyper <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HYPERT,
                                                         lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HYPERT, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

don <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DON,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DON, Diagnosis))))

mem <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$MEM,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(MEM, Diagnosis))))

gal <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$GAL,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(GAL, Diagnosis))))

riv <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$RIV,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(RIV, Diagnosis))))

tac <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TAC,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TAC, Diagnosis))))

com <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$COM,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(COM, Diagnosis))))

write.excel(com)

xy <- cbind(sex,alcohol,hear,diab,riv,hyper,dep,tbi,cvd)

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}

write.excel(xy)
write.excel(race)
write.excel(bmi)

