### Source code for deriving UK Biobank input phenotypes for latent class analysis (step 1)
### as carried out in Thijssen et al., 2023, Molecular Psychiatry
### Prepared by Jeanne Savage (j.e.savage@vu.nl)
### 19 Jan 2024

### This script assumes that a data file "phenotypes.txt" exists in the working directory, 
# in which the user has extracted the UKB data fields listed in the file "ThijssenMP2023_fieldcodes.txt".
# The file "ThijssenMP2023_fieldcodes.txt" is also assumed to be in the working directory.


# ----------------------------------------------------
# Setup
# ----------------------------------------------------

rm(list=ls())
setwd("")
require(data.table)
require(psych)
require(stringr)
require(Hmisc)
require(car)
require(ggplot2)
require(MASS)
require(polycor)
require(GPArotation)


# Read phenotype files into R & recode varnames
alc <- fread("phenotypes.txt",na.strings = c("NA","-3","-1","-6","-818",'-819',"-121"))

codes <- fread("ThijssenMP2023_fieldcodes.txt",fill=T)
for(i in 1:length(codes$V1)) {
  names(alc)[grep(codes$V1[i],names(alc))] <- str_replace(names(alc)[grep(codes$V1[i],names(alc))], paste0("f.",codes$V1[i]), codes$V2[i])
}

dim(alc)
names(alc)
summary(alc)


# ----------------------------------------------------
# Alcohol Phenotype preparation
# ----------------------------------------------------

# ----------------------------------------------------
### Create total quantity variable (gmEth/month; non-drinkers missing)

# Grams=sum of monthly*amount ethanol  
# first convert drink units to standard drinks
#8gm per standard drink in uk
#beer measured as 1 pint = 2.3 standard drinks (assuming 4% abv)
#wine as a glass (1/6 of a bottle); 9 standard drinks per bottle (assuming 12% abv) = 1.5 drinks per glass
#spirits as a standard measure (1/25 of a bottle) = 1 standard drink
#fortified wine as a glass (1/12 of a bottle) ~= 1.25 standard drinks (assuming 20% abv)
#"other" can't be determined. alcopops are 1.1-1.5. leave as 1 standard drink without further information.


### Visit 1
## Drinks per month
head(alc[,grep("quant.*mo.0",names(alc),val=T)])
# Rescale to standard drinks
alc$quant_rwine_mo.0.0 <- alc$quant_rwine_mo.0.0*1.5
alc$quant_wwine_mo.0.0 <- alc$quant_wwine_mo.0.0*1.5
alc$quant_beer_mo.0.0 <- alc$quant_beer_mo.0.0*2.3
alc$quant_fwine_mo.0.0 <- alc$quant_fwine_mo.0.0*1.25
# Sum drinks per month (exclude missing alls)
alc$quant_allmiss_mo <- apply(apply(alc[,grep("quant.*mo.0",names(alc),val=T)],1,is.na),2,sum)
alc$quant_mo.0.0 <- NA
alc$quant_mo.0.0[alc$quant_allmiss_mo != 6] <- apply(alc[alc$quant_allmiss_mo != 6,grep("quant.*mo.0",names(alc),val=T)],1,sum,na.rm=T)
summary(alc[,grep("quant.*mo.0",names(alc),val=T)])
hist(alc$quant_mo.0.0)
table(alc$quant_mo.0.0 > 600)
alc$quant_mo.0.0[which(alc$quant_mo.0.0 > 600)] <- 600 #Cap max drinks/mo at 600, or 20 drinks per day

## Drinks per week
head(alc[,grep("quant.*wk.0",names(alc),val=T)])
# Rescale to standard drinks
alc$quant_rwine_wk.0.0 <- alc$quant_rwine_wk.0.0*1.5
alc$quant_wwine_wk.0.0 <- alc$quant_wwine_wk.0.0*1.5
alc$quant_beer_wk.0.0 <- alc$quant_beer_wk.0.0*2.3
alc$quant_fwine_wk.0.0 <- alc$quant_fwine_wk.0.0*1.25
# Sum drinks per week (exclude missing alls)
alc$quant_allmiss_wk <- apply(apply(alc[,grep("quant.*wk.0",names(alc),val=T)],1,is.na),2,sum)
alc$quant_wk.0.0 <- NA
alc$quant_wk.0.0[alc$quant_allmiss_wk != 6] <- apply(alc[alc$quant_allmiss_wk != 6,grep("quant.*wk.0",names(alc),val=T)],1,sum,na.rm=T)
# Rescale drinks per week -> per month
for (i in grep("quant.*wk.0",names(alc),val=T)) {
  alc[,i] <- alc[,i]*4
}
summary(alc[,grep("quant.*wk.0",names(alc),val=T)])
hist(alc$quant_wk.0.0)
table(alc$quant_wk.0.0 > 600)
alc$quant_wk.0.0[which(alc$quant_wk.0.0 > 600)] <- 600 #Cap max drinks/mo at 600, or 20 drinks per day
table(is.na(alc$quant_mo.0.0), is.na(alc$quant_wk.0.0))

#calculate gm/eth per day
#8gm per standard drink in uk
alc$gm_eth0.0 <- NA
alc$gm_eth0.0[!(is.na(alc$quant_mo.0.0) & is.na(alc$quant_wk.0.0))] <- apply(alc[!(is.na(alc$quant_mo.0.0) & is.na(alc$quant_wk.0.0)),c("quant_mo.0.0","quant_wk.0.0")],1,sum,na.rm=T)*8/30
table(is.na(alc$gm_eth0.0))
summary(alc$gm_eth0.0)
hist(alc$gm_eth0.0)



### Repeat for visit 2
head(alc[,grep("quant.*mo.1",names(alc),val=T)])
alc$quant_rwine_mo.1.0 <- alc$quant_rwine_mo.1.0*1.5
alc$quant_wwine_mo.1.0 <- alc$quant_wwine_mo.1.0*1.5
alc$quant_beer_mo.1.0 <- alc$quant_beer_mo.1.0*2.3
alc$quant_fwine_mo.1.0 <- alc$quant_fwine_mo.1.0*1.25
alc$quant_allmiss_mo <- apply(apply(alc[,grep("quant.*mo.1",names(alc),val=T)],1,is.na),2,sum)
alc$quant_mo.1.0 <- NA
alc$quant_mo.1.0[alc$quant_allmiss_mo != 6] <- apply(alc[alc$quant_allmiss_mo != 6,grep("quant.*mo.1",names(alc),val=T)],1,sum,na.rm=T)
summary(alc[,grep("quant.*mo.1",names(alc),val=T)])
hist(alc$quant_mo.1.0)
table(alc$quant_mo.1.0 > 600)
alc$quant_mo.1.0[which(alc$quant_mo.1.0 > 600)] <- 600 #Cap max drinks/mo at 900, or 20 drinks per day

head(alc[,grep("quant.*wk.1",names(alc),val=T)])
alc$quant_rwine_wk.1.0 <- alc$quant_rwine_wk.1.0*1.5
alc$quant_wwine_wk.1.0 <- alc$quant_wwine_wk.1.0*1.5
alc$quant_beer_wk.1.0 <- alc$quant_beer_wk.1.0*2.3
alc$quant_fwine_wk.1.0 <- alc$quant_fwine_wk.1.0*1.25
alc$quant_allmiss_wk <- apply(apply(alc[,grep("quant.*wk.1",names(alc),val=T)],1,is.na),2,sum)
alc$quant_wk.1.0 <- NA
alc$quant_wk.1.0[alc$quant_allmiss_wk != 6] <- apply(alc[alc$quant_allmiss_wk != 6,grep("quant.*wk.1",names(alc),val=T)],1,sum,na.rm=T)
for (i in grep("quant.*wk.1",names(alc),val=T)) {
  alc[,i] <- alc[,i]*4
}
summary(alc[,grep("quant.*wk.1",names(alc),val=T)])
hist(alc$quant_wk.1.0)
table(alc$quant_wk.1.0 > 600)
alc$quant_wk.1.0[which(alc$quant_wk.1.0 > 600)] <- 600 #Cap max drinks/mo at 900, or 20 drinks per day
table(is.na(alc$quant_mo.1.0), is.na(alc$quant_wk.1.0))

alc$gm_eth1.0 <- NA
alc$gm_eth1.0[!(is.na(alc$quant_mo.1.0) & is.na(alc$quant_wk.1.0))] <- apply(alc[!(is.na(alc$quant_mo.1.0) & is.na(alc$quant_wk.1.0)),c("quant_mo.1.0","quant_wk.1.0")],1,sum,na.rm=T)*8/30
table(is.na(alc$gm_eth1.0))
summary(alc$gm_eth1.0)
hist(alc$gm_eth1.0)


### Repeat for visit 3
head(alc[,grep("quant.*mo.2",names(alc),val=T)])
alc$quant_rwine_mo.2.0 <- alc$quant_rwine_mo.2.0*1.5
alc$quant_wwine_mo.2.0 <- alc$quant_wwine_mo.2.0*1.5
alc$quant_beer_mo.2.0 <- alc$quant_beer_mo.2.0*2.3
alc$quant_fwine_mo.2.0 <- alc$quant_fwine_mo.2.0*1.25
alc$quant_allmiss_mo <- apply(apply(alc[,grep("quant.*mo.2",names(alc),val=T)],1,is.na),2,sum)
alc$quant_mo.2.0 <- NA
alc$quant_mo.2.0[alc$quant_allmiss_mo != 6] <- apply(alc[alc$quant_allmiss_mo != 6,grep("quant.*mo.2",names(alc),val=T)],1,sum,na.rm=T)
summary(alc[,grep("quant.*mo.2",names(alc),val=T)])
hist(alc$quant_mo.2.0)
table(alc$quant_mo.2.0 > 600)
alc$quant_mo.2.0[which(alc$quant_mo.2.0 > 600)] <- 600 #Cap max drinks/mo at 900, or 20 drinks per day

head(alc[,grep("quant.*wk.2",names(alc),val=T)])
alc$quant_rwine_wk.2.0 <- alc$quant_rwine_wk.2.0*1.5
alc$quant_wwine_wk.2.0 <- alc$quant_wwine_wk.2.0*1.5
alc$quant_beer_wk.2.0 <- alc$quant_beer_wk.2.0*2.3
alc$quant_fwine_wk.2.0 <- alc$quant_fwine_wk.2.0*1.25
alc$quant_allmiss_wk <- apply(apply(alc[,grep("quant.*wk.2",names(alc),val=T)],1,is.na),2,sum)
alc$quant_wk.2.0 <- NA
alc$quant_wk.2.0[alc$quant_allmiss_wk != 6] <- apply(alc[alc$quant_allmiss_wk != 6,grep("quant.*wk.2",names(alc),val=T)],1,sum,na.rm=T)
for (i in grep("quant.*wk.2",names(alc),val=T)) {
  alc[,i] <- alc[,i]*4
}
summary(alc[,grep("quant.*wk.2",names(alc),val=T)])
hist(alc$quant_wk.2.0)
table(alc$quant_wk.2.0 > 600)
alc$quant_wk.2.0[which(alc$quant_wk.2.0 > 600)] <- 600 #Cap max drinks/mo at 900, or 20 drinks per day
table(is.na(alc$quant_mo.2.0), is.na(alc$quant_wk.2.0))

alc$gm_eth2.0 <- NA
alc$gm_eth2.0[!(is.na(alc$quant_mo.2.0) & is.na(alc$quant_wk.2.0))] <- apply(alc[!(is.na(alc$quant_mo.2.0) & is.na(alc$quant_wk.2.0)),c("quant_mo.2.0","quant_wk.2.0")],1,sum,na.rm=T)*8/30
table(is.na(alc$gm_eth2.0))
summary(alc$gm_eth2.0)
hist(alc$gm_eth2.0)


# fill in missing data with valid data from later assessments, if available
table(is.na(alc$gm_eth0.0),is.na(alc$gm_eth1.0),is.na(alc$gm_eth2.0))
alc$gm_eth <- alc$gm_eth0.0
alc$gm_eth[is.na(alc$gm_eth)] <- alc$gm_eth1.0[is.na(alc$gm_eth)]  
alc$gm_eth[is.na(alc$gm_eth)] <- alc$gm_eth2.0[is.na(alc$gm_eth)]  
table(is.na(alc$gm_eth))

# log transform skewed data
alc$gm_eth_ln <- log(alc$gm_eth+1)
hist(alc$gm_eth_ln)



# ----------------------------------------------------
### Recode typical frequency vars to days/month (non-drinkers missing)

alc$drinkfreq.0.0 <- recode(alc$drinkfreq.0.0, "1=20;2=14;3=6;4=2;5=.5;6=NA")
alc$drinkfreq.1.0 <- recode(alc$drinkfreq.1.0, "1=20;2=14;3=6;4=2;5=.5;6=NA")
alc$drinkfreq.2.0 <- recode(alc$drinkfreq.2.0, "1=20;2=14;3=6;4=2;5=.5;6=NA")

# fill in missing data with valid data from later assessments, if available
alc$drinkfreq <- alc$drinkfreq.0.0
alc$drinkfreq[is.na(alc$drinkfreq)] <- alc$drinkfreq.1.0[is.na(alc$drinkfreq)]  
alc$drinkfreq[is.na(alc$drinkfreq)] <- alc$drinkfreq.2.0[is.na(alc$drinkfreq)]  
table(is.na(alc$drinkfreq))
hist(alc$drinkfreq)


# ----------------------------------------------------
### Recode binge frequency vars to days/month

table(alc$mh_bingefreq.0.0)
#recode freq as 0 in those who answered mental health frequency question as "never" and exclude abstainers
alc$mh_bingefreq.0.0[alc$mh_freq.0.0==0 & !(alc$drinkerstatus.0.0 %in% c(0,1))] <- 1
alc$mh_bingefreq.0.0[alc$mh_bingefreq.0.0==1 & (alc$drinkerstatus.0.0 %in% c(0,1))] <- NA
alc$mh_binge_ln <- log(recode(alc$mh_bingefreq.0.0, "1=0;2=.5;3=1;4=4;5=20")+1)
hist(alc$mh_binge_ln)
# log binge frequency is still very skewed - create dichotomized version
alc$mh_binge_D <- recode(alc$mh_bingefreq.0.0, "1=1;2:5=2")
table(alc$mh_binge_D)


# ----------------------------------------------------
### Individual AUD symptoms from AUDIT (mental health questionnaire)

table(alc$aud_blackouts.0.0)
table(alc$aud_responsibilities.0.0)
table(alc$aud_cantstop.0.0)
table(alc$aud_guilt.0.0)
table(alc$aud_withdrawal.0.0)
table(alc$concerned_friends.0.0)
table(alc$injured_drinking.0.0)
alc$blackouts <- recode(alc$aud_blackouts.0.0,"1=1;2:5=2;else=NA")
alc$responsibilities <- recode(alc$aud_responsibilities.0.0,"1=1;2:5=2;else=NA")
alc$cantstop <- recode(alc$aud_cantstop.0.0,"1=1;2:5=2;else=NA")
alc$guilt <- recode(alc$aud_guilt.0.0,"1=1;2:5=2;else=NA")
alc$withdrawal <- recode(alc$aud_withdrawal.0.0,"1=1;2:5=2;else=NA")
alc$concern <- recode(alc$concerned_friends.0.0,"0=1;1=2;2=2;else=NA")
alc$injury <- recode(alc$injured_drinking.0.0,"0=1;1=2;2=2;else=NA")
table(alc$blackouts)
table(alc$responsibilities)
table(alc$cantstop)
table(alc$guilt)
table(alc$withdrawal)
table(alc$concern)
table(alc$injury)



# ----------------------------------------------------
# Dummy variables for abstainer vs. current and former vs. current drinker

table(alc$drinkerstatus.0.0)
table(alc$drinkerstatus.0.0,alc$former_drinker.0.0)
alc$abstainer.0.0 <- NA
alc$abstainer.0.0[alc$drinkerstatus.0.0==0] <- 1
alc$abstainer.0.0[alc$drinkerstatus.0.0==2] <- 0
alc$former.0.0 <- NA
alc$former.0.0[alc$drinkerstatus.0.0==1] <- 1
alc$former.0.0[alc$drinkerstatus.0.0==2] <- 0
alc$abstainer.1.0 <- NA
alc$abstainer.1.0[alc$drinkerstatus.1.0==0] <- 1
alc$abstainer.1.0[alc$drinkerstatus.1.0==2] <- 0
alc$former.1.0 <- NA
alc$former.1.0[alc$drinkerstatus.1.0==1] <- 1
alc$former.1.0[alc$drinkerstatus.1.0==2] <- 0
alc$abstainer.2.0 <- NA
alc$abstainer.2.0[alc$drinkerstatus.2.0==0] <- 1
alc$abstainer.2.0[alc$drinkerstatus.2.0==2] <- 0
alc$former.2.0 <- NA
alc$former.2.0[alc$drinkerstatus.2.0==1] <- 1
alc$former.2.0[alc$drinkerstatus.2.0==2] <- 0
table(alc$abstainer.0.0)
table(alc$former.0.0)
table(alc$abstainer.0.0,alc$abstainer.1.0,alc$abstainer.2.0)
table(alc$former.0.0,alc$former.1.0,alc$former.2.0)

alc$abstainer <- alc$abstainer.0.0
alc$abstainer[is.na(alc$abstainer)] <- alc$abstainer.1.0[is.na(alc$abstainer)]
alc$abstainer[is.na(alc$abstainer)] <- alc$abstainer.2.0[is.na(alc$abstainer)]

alc$former <- alc$former.0.0
alc$former[is.na(alc$former)] <- alc$former.1.0[is.na(alc$former)]
alc$former[is.na(alc$former)] <- alc$former.2.0[is.na(alc$former)]

table(alc$abstainer)
table(alc$former)

# ----------------------------------------------------
# Dummy variables for increase vs. stable and decrease vs. stable quantity (past 10 years)
table(alc$drink_diff_10yrs_ago.0.0)
alc$increasedrink.0.0 <- NA
alc$increasedrink.0.0[alc$drink_diff_10yrs_ago.0.0==1] <- 1
alc$increasedrink.0.0[alc$drink_diff_10yrs_ago.0.0==2] <- 0
alc$decreasedrink.0.0 <- NA
alc$decreasedrink.0.0[alc$drink_diff_10yrs_ago.0.0==3] <- 1
alc$decreasedrink.0.0[alc$drink_diff_10yrs_ago.0.0==2] <- 0
alc$increasedrink.1.0 <- NA
alc$increasedrink.1.0[alc$drink_diff_10yrs_ago.1.0==1] <- 1
alc$increasedrink.1.0[alc$drink_diff_10yrs_ago.1.0==2] <- 0
alc$decreasedrink.1.0 <- NA
alc$decreasedrink.1.0[alc$drink_diff_10yrs_ago.1.0==3] <- 1
alc$decreasedrink.1.0[alc$drink_diff_10yrs_ago.1.0==2] <- 0
alc$increasedrink.2.0 <- NA
alc$increasedrink.2.0[alc$drink_diff_10yrs_ago.2.0==1] <- 1
alc$increasedrink.2.0[alc$drink_diff_10yrs_ago.2.0==2] <- 0
alc$decreasedrink.2.0 <- NA
alc$decreasedrink.2.0[alc$drink_diff_10yrs_ago.2.0==3] <- 1
alc$decreasedrink.2.0[alc$drink_diff_10yrs_ago.2.0==2] <- 0

table(alc$increasedrink.0.0,alc$increasedrink.1.0,alc$increasedrink.2.0)
table(alc$decreasedrink.0.0,alc$decreasedrink.1.0,alc$decreasedrink.2.0)

alc$increasedrink <- alc$increasedrink.0.0
alc$increasedrink[is.na(alc$increasedrink)] <- alc$increasedrink.1.0[is.na(alc$increasedrink)]
alc$increasedrink[is.na(alc$increasedrink)] <- alc$increasedrink.2.0[is.na(alc$increasedrink)]

alc$decreasedrink <- alc$decreasedrink.0.0
alc$decreasedrink[is.na(alc$decreasedrink)] <- alc$decreasedrink.1.0[is.na(alc$decreasedrink)]
alc$decreasedrink[is.na(alc$decreasedrink)] <- alc$decreasedrink.2.0[is.na(alc$decreasedrink)]

table(alc$increasedrink)
table(alc$decreasedrink)

# ----------------------------------------------------

alcvars <- c("mh_binge_D", "gm_eth_ln","drinkfreq",
             "blackouts","responsibilities",
             "cantstop","guilt","withdrawal","concern","injury",
             "former","abstainer","increasedrink","decreasedrink")




# ----------------------------------------------------
# Internalizing Phenotype preparation
# ----------------------------------------------------

# Select variables
anxvars2w <- grep("a2_",names(alc),val=T)
depvars2w <- grep("d2_",names(alc),val=T)
neuvars <- grep("neu",names(alc),val=T)

summary(alc[,anxvars2w])
summary(alc[,depvars2w])

# Create sum scores for symptom counts (excluding those missing data on all questions)
alc$allmiss <- apply(apply(alc[,anxvars2w],1,is.na),2,sum)
summary(alc$allmiss)
alc$anx_sxs_2w <- NA
alc$anx_sxs_2w[alc$allmiss != max(alc$allmiss)] <- apply(alc[alc$allmiss != max(alc$allmiss),anxvars2w],1,sum,na.rm=T)
summary(alc$anx_sxs_2w)
hist(alc$anx_sxs_2w)
alc$allmiss <- apply(apply(alc[,depvars2w],1,is.na),2,sum)
summary(alc$allmiss)
alc$dep_sxs_2w <- NA
alc$dep_sxs_2w[alc$allmiss != max(alc$allmiss)] <- apply(alc[alc$allmiss != max(alc$allmiss),depvars2w],1,sum,na.rm=T)
summary(alc$dep_sxs_2w)
hist(alc$dep_sxs_2w)

hist(alc$neu_sum.0.0)

intvars <- c("anx_sxs_2w","dep_sxs_2w","neu_sum.0.0")




# ----------------------------------------------------
# Externalizing Phenotype preparation
# ----------------------------------------------------

table(alc$ever_cannabis.0.0)
alc$ever_cann_D <- recode(alc$ever_cannabis.0.0, "0=1;1:4=2")
table(alc$ever_cann_D)

table(alc$ever_addicted.0.0)



# ----------------------------------------------------
# Psychiatric disorder phenotype preparation
# ----------------------------------------------------

#update derived first occurrence diagnosis variables
alc$alcohol_use_dx[alc$alcohol_use_dx > 0 | alc$alc_cirrhosis > 0] <- 1
alc$alcohol_use_dx[is.na(alc$alcohol_use_dx)] <- 0
table(alc$alcohol_use_dx)
alc$alcohol_use_dx[alc$drinkerstatus.0.0==0 & !alc$alcohol_use_dx==1] <- NA         # Exclude non-initiators as AUD controls

alc$tobacco_use_dx[alc$tobacco_use_dx > 0] <- 1
alc$tobacco_use_dx[is.na(alc$tobacco_use_dx)] <- 0
table(alc$tobacco_use_dx)

alc$illicit_substance_use_dx <- 0
alc$illicit_substance_use_dx[alc$opioid_use_dx>0 | alc$cannabis_use_dx>0 | alc$sedative_use_dx>0 | 
                              alc$cocaine_use_dx>0 | alc$stim_use_dx>0 | alc$halluc_use_dx>0 | 
                              alc$solvent_use_dx>0 | alc$multdrug_use_dx>0] <- 1
table(alc$illicit_substance_use_dx)

#coding scheme for diagnoses from mental health questionnaire
vars.mh <- grep("mh_diagnosis",names(alc),val=T)
codes <- list(mdd=11,
              panic=6,
              generalized_anxiety=15,
              phobia=c(1,5,17))


# Extract MHQ diagnostic codes for each phenotype
isin <- function(data, target) {
  if (!FALSE %in% is.na(data)) {return(0)}
  if (!FALSE %in% (data=="")) {return(0)}
  return(as.numeric(TRUE %in% (data %in% target)))
}

alc$mdd <- apply(alc[,vars.mh],1, isin, target=codes$mdd)
alc$panic <- apply(alc[,vars.mh],1, isin, target=codes$panic)
alc$generalized_anxiety <- apply(alc[,vars.mh],1, isin, target=codes$generalized_anxiety)
alc$phobia <- apply(alc[,vars.mh],1, isin, target=codes$phobia)

#combine FO and MHQ diagnoses
alc$mdd[alc$mdd.0.0>0] <- 1
alc$panic[alc$panic.0.0>0] <- 1
alc$phobia[alc$phobia.0.0>0] <- 1

summary(alc)

psyvars <- c("mdd","panic","generalized_anxiety","phobia","alcohol_use_dx","tobacco_use_dx","illicit_substance_use_dx")


################
# Output file for Mplus mixture model analyses

#select variables
mpvars <- c("f.eid",psyvars,alcvars,intvars,"ever_addicted.0.0","ever_cann_D")

#exclude relatives
mplus <- alc[alc$relative_exclusion==F,mpvars]

#save mplus input file
fwrite(mplus,"ukb_alc.dat",sep="\t",quote=F,row.names=F, col.names=F, na=9999)





