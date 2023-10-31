rm(list = ls())
library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)
library(moments)

File <- file.path("C:/Users/Daman Kaur/Desktop/Ulster/Research Projects/NACC drug classes analysis paper/Cox regression csv")
setwd(File)
df <- fread("NMallmedscox.csv")

#### age at diagnosis ####
nacc <- fread("MtoDcdrcox.csv")
dh <- which(nacc$NACCID %in% df$NACCID)
data <- nacc[ c(dh), ]
data <- as.data.frame(data)

#write.csv(data, file = "NMc.csv",row.names=FALSE, na="")
#data <- nacc
data <- data[data$RACE < 50, ] 
file <- data %>% filter(Diagnosis == 1)
file <- file %>% filter(DIAG == 4) ## change according to progression group # 3 for MCI 4 for Dementia
mean(file$NACCAGE)
sd(file$NACCAGE)

###### duration of follow up / progression #####
file <- df %>% filter(Diagnosis == 1)

xy <- file %>%
  mutate(date_visit = as.Date(Date, "%Y/%m/%d"))
xy <- xy %>% 
  group_by(NACCID) %>%
  summarise(first_visit = min(date_visit), 
            last_visit = max(date_visit), 
            exposure = last_visit - first_visit)
xy$exposure <- as.numeric(xy$exposure)

xy$yrs <- (xy$exposure/365.2422)
mean(xy$yrs)
sd(xy$yrs)

############################################

df <- df[df$NEW_VISITNUM == 1, ] 
df <- df[df$RACE < 50, ] 
df <- df[df$NACCAGE >= 40, ] 

stbl <- df %>% filter(Diagnosis == 0)
prog <- df %>% filter(Diagnosis == 1)

#### numbers - ttest - normality ####

#normality
shapiro.test(stbl$NACCAGE)
shapiro.test(prog$NACCAGE)
agostino.test(stbl$NACCAGE)

#mean and SD
mean(stbl$NACCAGE)
sd(stbl$NACCAGE)
mean(prog$NACCAGE)
sd(prog$NACCAGE)
t.test(stbl$NACCAGE,prog$NACCAGE)
wilcox.test(stbl$NACCAGE,prog$NACCAGE) 

mean(stbl$ALVST)
sd(stbl$ALVST)
mean(prog$ALVST)
sd(prog$ALVST)

mean(stbl$EDUC)
sd(stbl$EDUC)
mean(prog$EDUC)
sd(prog$EDUC)
t.test(stbl$EDUC,prog$EDUC)
wilcox.test(stbl$EDUC,prog$EDUC) 

mean(stbl$SMOKE, na.rm = TRUE)
sd(stbl$SMOKE, na.rm = TRUE)
mean(prog$SMOKE, na.rm = TRUE)
sd(prog$SMOKE, na.rm = TRUE)
t.test(stbl$SMOKE,prog$SMOKE)
wilcox.test(stbl$SMOKE,prog$SMOKE) 


## marry / fam

stbl1 <- stbl %>% filter(MARI == 1)
stbl1 <- stbl %>% filter(FAM == 1)

stbl1 <- prog %>% filter(MARI == 1)
stbl1 <- prog %>% filter(FAM == 1)

#### other covariates
sum(stbl$MARI==1, na.rm=TRUE)    
sum(prog$MARI==1, na.rm=TRUE)

sum(stbl$FAM==1, na.rm=TRUE)
sum(prog$FAM==1, na.rm=TRUE)


sum(stbl$SEX==2, na.rm=TRUE)    ## sum((stbl$SEX==2)/9706*100) # 0 female
sum(prog$SEX==2, na.rm=TRUE)
sum(stbl$SEX==1, na.rm=TRUE)
sum(prog$SEX==1, na.rm=TRUE)

#sum(stbl$DIAB==0)
#sum(prog$DIAB==0)
sum(stbl$DIAB==1, na.rm=TRUE)
sum(prog$DIAB==1, na.rm=TRUE)

#sum(stbl$HYPERT==0)
#sum(prog$HYPERT==0)
sum(stbl$HYPERT==1, na.rm=TRUE)
sum(prog$HYPERT==1, na.rm=TRUE)

#sum(stbl$DEPRSN==0)
#sum(prog$DEPRSN==0)
sum(stbl$DEPRSN==1, na.rm=TRUE)
sum(prog$DEPRSN==1, na.rm=TRUE)

#sum(stbl$CVD==0)
#sum(prog$CVD==0)
sum(stbl$CVD==1, na.rm=TRUE)
sum(prog$CVD==1, na.rm=TRUE)

#sum(stbl$HEAR==0)
#sum(prog$HEAR==0)
sum(stbl$HEAR==1, na.rm=TRUE)
sum(prog$HEAR==1, na.rm=TRUE)

#sum(stbl$ALCOHOL==0)
#sum(prog$ALCOHOL==0)
sum(stbl$ALCOHOL==1, na.rm=TRUE)
sum(prog$ALCOHOL==1, na.rm=TRUE)

#sum(stbl$ABMI==0)
#sum(prog$ABMI==0)
sum(stbl$ABMI==1, na.rm=TRUE)
sum(prog$ABMI==1, na.rm=TRUE)
sum(stbl$ABMI==2, na.rm=TRUE)
sum(prog$ABMI==2, na.rm=TRUE)

sum(stbl$RACE==1, na.rm=TRUE)
sum(prog$RACE==1, na.rm=TRUE)
sum(stbl$RACE==2, na.rm=TRUE)
sum(prog$RACE==2, na.rm=TRUE)
sum(stbl$RACE==3, na.rm=TRUE)
sum(prog$RACE==3, na.rm=TRUE)
sum(stbl$RACE==4, na.rm=TRUE)
sum(prog$RACE==4, na.rm=TRUE)
sum(stbl$RACE==5, na.rm=TRUE)
sum(prog$RACE==5, na.rm=TRUE)

#### medication numbers ####

prog <- prog[,c(138:147)]
prog[prog==0]<-NA

stbl <- stbl[,c(138:147)]
stbl[stbl==0]<-NA

xy <- as.data.frame(colMeans(prog, na.rm = TRUE))
write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}

write.excel(xy)

sd(stbl$NAHTN, na.rm = TRUE)
sd(stbl$NLIPL, na.rm = TRUE)
sd(stbl$NNSD, na.rm = TRUE)
sd(stbl$NAC, na.rm = TRUE)
sd(stbl$NADEP, na.rm = TRUE)
sd(stbl$NAPSY, na.rm = TRUE)
sd(stbl$NANX, na.rm = TRUE)
sd(stbl$NDBM, na.rm = TRUE)
sd(stbl$NADM, na.rm = TRUE)
sd(stbl$NPDM, na.rm = TRUE)

sd(prog$NAHTN, na.rm = TRUE)
sd(prog$NLIPL, na.rm = TRUE)
sd(prog$NNSD, na.rm = TRUE)
sd(prog$NAC, na.rm = TRUE)
sd(prog$NADEP, na.rm = TRUE)
sd(prog$NAPSY, na.rm = TRUE)
sd(prog$NANX, na.rm = TRUE)
sd(prog$NDBM, na.rm = TRUE)
sd(prog$NADM, na.rm = TRUE)
sd(prog$NPDM, na.rm = TRUE)

###### chisq test #####

#data <- matrix(c(582,6492,359,3214), nrow = 2)
#chisq.test(data)

table <- table(df$Diagnosis, df$SEX)
chisq.test(table)

table <- table(df$Diagnosis, df$MARI)
chisq.test(table)

table <- table(df$Diagnosis, df$ABMI)
chisq.test(table)

table <- table(df$Diagnosis, df$FAM)
chisq.test(table)

table <- table(df$Diagnosis, df$DIAB)
chisq.test(table)

table <- table(df$Diagnosis, df$HYPERT)
chisq.test(table)

table <- table(df$Diagnosis, df$CVD)
chisq.test(table)

table <- table(df$Diagnosis, df$DEPRSN)
chisq.test(table)

table <- table(df$Diagnosis, df$HEAR)
chisq.test(table)

table <- table(df$Diagnosis, df$RACE)
chisq.test(table)

table <- table(df$Diagnosis, df$ALCOHOL)
chisq.test(table)
