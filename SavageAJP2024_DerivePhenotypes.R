### Source code for deriving UK Biobank input phenotypes for GWAS/genomic SEM
### as carried out in Savage et al., American Journal of Psychiatry, 2024
### Prepared by Jeanne Savage (j.e.savage@vu.nl)
### 8 Feb 2025

### This script assumes that two data files exist in the working directory:
# (1) "all_lkps_maps.xlsx" - UKB provided spreadsheet listing all clinical codes and mapping schemes for electronic health data
#       File can be downloaded from https://biobank.ctsu.ox.ac.uk/showcase/label.cgi?id=3000
# (2) "SavageAJP2024_clinical_alcohol_maps.xlsx" - subset of (1) with selected alcohol-related clinical codes
#       File is included with this data return



# ----------------------------------------------------
# Setup
# ----------------------------------------------------

rm(list=ls())
setwd("")
library(readxl)
library(data.table)
library(car)
library(stringr)


# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------
# Define clinical event/drug codes for use in extracting data from UKB GP database
### As multiple coding systems have been used, we first check correspondence of different codes
### and identify all possible codes for the same event/diagnosis/prescription drug
### We create list of codes in both READ2 and READ3 formats for extraction (where applicable)
# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------

#icd to read2 mapping
ri9 <- read_excel("all_lkps_maps.xlsx",sheet=11)
ri10 <- read_excel("all_lkps_maps.xlsx",sheet=12)
                                                                                                                                       
#read2 to read3/ctv                                                                                                                                        
r2r3 <- read_excel("all_lkps_maps.xlsx",sheet=14)
                                                                                                                                        
#preselected icd/read2 codes to extract                                                                                                                                        
i9 <- read_excel("SavageAJP2024_clinical_alcohol_maps.xlsx",sheet=5)
i10 <- read_excel("SavageAJP2024_clinical_alcohol_maps.xlsx",sheet=6)
r2 <- read_excel("SavageAJP2024_clinical_alcohol_maps.xlsx",sheet=8)

#get read2 codes for icd + convert to read3
table(i9$ICD9 %in% ri9$icd9_code)
mri9 <- merge(ri9,i9,by.x="icd9_code",by.y="ICD9")
table(i10$ALT_CODE%in% ri10$icd10_code)
mri10 <- merge(ri10,i10[,c(1,2,4,5)],by.x="icd10_code",by.y="ALT_CODE")

table(mri9$read_code %in% r2r3$READV2_CODE)
mr2r3i9 <- merge(r2r3[,c(2,7,8)],mri9,by.y="read_code",by.x="READV2_CODE")
table(mri10$read_code %in% r2r3$READV2_CODE)
mr2r3i10 <- merge(r2r3[,c(2,7,8)],mri10,by.y="read_code",by.x="READV2_CODE")

#collate single list of clinical event read codes
table(r2$read_code %in% c(mri10$read_code,mri9$read_code) )
allr2 <- unique(c(r2$read_code, mri10$read_code,mri9$read_code) )
mr2r3 <- merge(r2r3[,c(2,7,8)],r2,by.y="read_code",by.x="READV2_CODE")
allr3 <- unique(c(mr2r3$READV3_CODE, mr2r3i10$READV3_CODE,mr2r3i9$READV3_CODE) )


# read in drug codes
bnf <- read_excel("SavageAJP2024_clinical_alcohol_maps.xlsx",sheet=3)
dmd <- read_excel("SavageAJP2024_clinical_alcohol_maps.xlsx",sheet=4)
r2d <- read_excel("SavageAJP2024_clinical_alcohol_maps.xlsx",sheet=9)


# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------
# Extract relevant codes from UKB databases
### This step requires access to individual level record data (https://biobank.ctsu.ox.ac.uk/showcase/label.cgi?id=3001) via an approved Project
### from the databases "GP clinical event records" (field 42040, gp_clinical) and "GP prescription records" (field 42039, gp_scripts)
### The R code below generates two SQL commands that can be pasted into the database to generate downloadable files
### Save these files as "gp_clinical_alcohol.txt" and "gp_scripts_alcohol.txt", respectively
### Note that these databases are continuously updated so results will not exactly match the original publication
# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------

# Generate sql codes to extract data
paste0("SELECT * FROM gp_scripts WHERE read_2 = '", paste0(r2d$read_code, collapse="' OR read_2 = '"), "' OR bnf_code = '", paste0(bnf$BNF_Presentation_Code, collapse="' OR bnf_code = '"), "' OR dmd_code = '", paste0(dmd$concept_id, collapse="' OR dmd_code = '"), "'")
paste0("SELECT * FROM gp_clinical WHERE read_2 = '", paste0(allr2, collapse="' OR read_2 = '"), "' OR read_3 = '", paste0(allr3, collapse="' OR read_3 = '"),  "'")




# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------
# Check and clean extracted files
# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------

# annotate codes from sql output dataset
###check duplicated codes and remove
as.data.frame(r2[r2$read_code %in% r2$read_code[duplicated(r2$read_code)],])
ar2 <- r2[!duplicated(r2$read_code),]

r3 <- r2r3[r2r3$READV2_CODE %in% allr2,c(7,8,10)]
r3 <- r3[!is.na(r3$TERMV3_DESC),]
as.data.frame(r3[r3$READV3_CODE %in% r3$READV3_CODE[duplicated(r3$READV3_CODE)],])
ar3 <- r3[!duplicated(r3$READV3_CODE),]

###merge clinical
clin <- fread("gp_clinical_alcohol.txt")
clin2 <- merge(clin,ar2,all.x=T,by.x="read_2",by.y="read_code")
clin3 <- merge(clin2,ar3,all.x=T,by.x="read_3",by.y="READV3_CODE")
clin3$desc <- apply(clin3[,c('term_description','TERMV3_DESC')],1,max,na.rm=T)

###fix errors (!!! check for errors manually as these are likely to differ across iterations of the database)
clin3$desc[clin3$read_2=="F1440"] <- "Cerebral ataxia due to alcoholism"
clin3$desc[clin3$read_2=="J61.."] <- "Cirrhosis and chronic liver disease"
clin3$desc[clin3$read_2=="SM0.."] <- "Alcohol causing toxic effect"
clin3$desc[clin3$read_2=="ZV4KC"] <- "Alcohol use"
clin3 <- clin3[!is.na(clin3$desc),]
clin3 <- clin3[!substr(clin3$read_3,1,2)=="A9",]
sort(table(clin3$desc))

fwrite(clin3, "gp_clinical_alcohol_annotated.txt", sep="\t",quote=F,na=NA)


### merge prescriptions
rx <- fread("gp_scripts_alcohol.txt")
rx <- merge(rx,r2d,all.x=T,by.x="read_2",by.y="read_code")
rx$desc <- apply(rx[,c('term_description','drug_name')],1,max,na.rm=T)
sort(table(rx$desc))
rx$drug <- character()
rx$drug[rx$desc %in% c("CAMPRAL EC tab 333mg","Campral EC 333mg tablets (Merck Serono Ltd)","CAMPRAL EC tabs 333mg","CAMPRAL EC 333mg e/c tablets","Acamprosate 333mg gastro-resistant tablets","ACAMPROSATE CALCIUM 333mg e/c tablets")] <- "acamprosate"
rx$drug[rx$desc %in% c("ANTABUSE TAB 200mg","Antabuse 200mg tablets (Actavis UK Ltd)","ANTABUSE 200mg tablets","Disulfiram 200mg tablets","DISULFIRAM 200mg tablets")] <- "disulfiram"
rx$drug[rx$desc %in% c("Nalmefene 18mg tablets")] <- "nalmefene"
table(rx$drug)
fwrite(rx, "gp_scripts_alcohol_annotated.txt", sep="\t",quote=F,na=NA)



###group phenotype categories

clin3$value1[clin3$value1=="Y"] <- NA

liver <- clin3[clin3$desc %in% c("Alcoholic hepatic failure","Portal cirrhosis unspecified","Chronic alcoholic hepatitis","Alcoholic fibrosis and sclerosis of liver","Acute alcoholic hepatitis","Alcoholic hepatitis","Cirrhosis and chronic liver disease","Alcoholic cirrhosis of liver","Alcoholic fatty liver","Alcoholic liver damage unspecified"),]
korsakoff <- clin3[clin3$desc %in% c("Cerebral degeneration due to alcoholism","Chronic alcoholic brain syndrome","Cerebral ataxia due to alcoholism","Alcoholic encephalopathy","Alcoholic dementia NOS","Cerebellar ataxia due to alcoholism","Korsakov psychosis","Korsakov's alcoholic psychosis with peripheral neuritis","Other alcoholic psychosis","[X]Mental and behavioural disorders due to use of alcohol: residual and late-onset psychotic disorder","Wernicke-Korsakov syndrome","Alcoholic psychoses","Korsakov's alcoholic psychosis","Alcohol amnestic syndrome NOS","[X]Mental and behavioural disorders due to use of alcohol: amnesic syndrome","[X]Mental and behavioural disorders due to use of alcohol: psychotic disorder","Alcoholic psychosis NOS","Other alcoholic dementia"),]
rehab <- clin3[clin3$desc %in% c("Delivery of rehabilitation for alcohol addiction","In-house alcohol detoxification","Alcoholics anonymous","Alcohol counselling by other agencies","[V]Alcohol rehabilitation","Aversion therapy - alcoholism","Referral to alcohol brief intervention service","Admitted to alcohol detoxification centre","Under care of community alcohol team","Referral to specialist alcohol treatment service","Referral to community alcohol team","Alcohol detoxification","Community detoxification registered","Referral to community alcohol team declined","Declined referral to specialist alcohol treatment service"),]
withdrawal <- clin3[clin3$desc %in% c("[X]Mental and behavioural disorders due to use of alcohol: withdrawal state","Alcohol withdrawal hallucinosis","[X]Mental and behavioural disorders due to use of alcohol: withdrawal state with delirium","Alcohol withdrawal delirium","[X]Alcohol withdrawal-induced seizure","Alcohol withdrawal syndrome"),]
abuse <- clin3[clin3$desc %in% c("[X]Mental and behavioural disorders due to use of alcohol: harmful use","Nondependent alcohol abuse, continuous","Nondependent alcohol abuse, episodic","Nondependent alcohol abuse NOS","Nondependent alcohol abuse, unspecified","Harmful alcohol use","Hazardous alcohol use","Nondependent alcohol abuse","Alcohol misuse","Alcohol misuse - enhanced services administration","Alcohol misuse enhanced service completed","Alcohol misuse enhanced services administration","Alcohol misuse - enhanced service completed","Increasing risk drinking","Higher risk drinking","Binge drinker"),]
excess <- clin3[clin3$desc %in% c("Alcohol intake above recommended sensible limits","Alcohol intake within recommended sensible limits","Feels should cut down drinking","Brief intervention for excessive alcohol consumption completed","Extended intervention for excessive alcohol consumption completed","Brief intervention for excessive alcohol consumption declined","Extended intervention for excessive alcohol consumption declined"),]
aud <- clin3[clin3$desc %in% c("[V]Personal history of alcoholism","[X]Mental and behavioural disorders due to use of alcohol","Alcohol dependence syndrome NOS","Alcohol dependence syndrome","[X]Mental and behavioural disorders due to use of alcohol: dependence syndrome","Unspecified chronic alcoholism","Continuous chronic alcoholism","Chronic alcoholism NOS","Chronic alcoholism","Episodic chronic alcoholism","Acute alcoholic intoxication in alcoholism NOS","Acute alcoholic intoxication in alcoholism","Continuous acute alcoholic intoxication in alcoholism","Acute alcoholic intoxication, unspecified, in alcoholism"),]
chronicaud <- clin3[clin3$desc %in% c("Unspecified chronic alcoholism","Continuous chronic alcoholism","Chronic alcoholism NOS","Chronic alcoholism","Episodic chronic alcoholism"),]
recovery <- clin3[clin3$desc %in% c("Alcohol dependence resolved","Nondependent alcohol abuse in remission","Chronic alcoholism in remission","Ex-very heavy drinker-(>9u/d)","Ex-heavy drinker - (7-9u/day)","Ex-moderate drinker - (3-6u/d)","Ex-light drinker - (1-2u/day)","Ex-trivial drinker (<1u/day)"),]
quant_gp <- clin3[clin3$desc %in% c("Alcohol consumption screen","Alcohol consumption NOS","Alcohol consumption","Alcohol units per week"),]
drinktype <- clin3[clin3$desc %in% c("Drinks beer and spirits","Spirit drinker","Beer drinker","Drinks wine"),]
drinkerstatus <- clin3[clin3$desc %in% c("Current non-drinker","Current non drinker","Light drinker","Moderate drinker","Social drinker","Heavy drinker","Very heavy drinker","Binge drinker"),]
advice <- clin3[clin3$desc %in% c("Lifestyle advice regarding alcohol","Patient advised about alcohol","Advised to abstain from alcohol consumption","Advised to contact primary care alcohol worker","Alcohol leaflet given","Alcohol consumption counselling","Alcohol disorder monitoring","Alcohol abuse monitoring","[V]Alcohol abuse counselling and surveillance"),]
audit <- clin3[clin3$desc %in% c("Alcohol use disorders identification test","Alcohol screen - alcohol use disorder identification test completed"),]
auditc <- clin3[clin3$desc %in% c("Alcohol use disorder identification test consumption questionnaire","Alcohol screen - alcohol use disorder identification test consumption questions completed","Alcohol screen - alcohol use disorder identification test completed"),]
auditpic <- clin3[clin3$desc %in% c("Alcohol use disorder identification test Piccinelli consumption questionnaire","Alcohol screen - alcohol use disorder identification test Piccinelli consumption questions completed"),]
fast <- clin3[clin3$desc %in% c("Fast alcohol screening test","FAST - Fast Alcohol Screening Test","Alcohol screen - fast alcohol screening test completed"),]
acute <- clin3[clin3$desc %in% c("Inebriety NOS","[X]Mental and behavioural disorders due to use of alcohol: acute intoxication","Drunkenness NOS","O/E - breath smell NOS","Hangover from alcohol","Pathological alcohol intoxication","Substance level in breath","[X]Evidence of alcohol involvement determined by blood alcohol level","[X]Evidence of alcohol involvement determined by blood alcohol level of 240 mg/100 ml or more","[X]Evidence of alcohol involvement determined by level of intoxication","[X]Evidence of alcohol involvement determined by blood alcohol level of 20-39 mg/100 ml","Breath ethanol level","[D]Alcohol blood level excessive"),]
otherclin <- clin3[clin3$desc %in% c("Ethanol causing toxic effect","Ethyl alcohol causing toxic effect","Alcohol causing toxic effect","HO/RTS - police: venesection alcohol","HO/RTS - police: venesection alcohol","Disqualified from driving due to excess alcohol","Hospital alcohol liaison team report received","Emergency department attendance related to personal alcohol consumption","Alcoholic paranoia","[X]Mental and behavioural disorders due to use of alcohol: unspecified mental and behavioural disorder","[X]Mental and behavioural disorders due to use of alcohol: other mental and behavioural disorders"),]

#define categorical phenotypes
liver$liver <- korsakoff$korsakoff <- rehab$rehab <- withdrawal$withdrawal <- abuse$abuse <- aud$aud <- chronicaud$chronicaud <- advice$advice <- rx$alc_rx <- acute$acute <- otherclin$otherclin <- 1
excess$excess <- recode(excess$desc, "'Alcohol intake within recommended sensible limits'=0;else=1")
drinkerstatus$drinkerstatus <- recode(drinkerstatus$desc, "'Current non drinker'=0;'Current non-drinker'=0;'Light drinker'=1;'Social drinker'=2;'Moderate drinker'=3;'Binge drinker'=4;'Heavy drinker'=5;'Very heavy drinker'=6")

#define and clean quantitative phenotypes
quant_gp$value1[quant_gp$value1 %in% c("ADV001", "^", "136..>] priority=1","136..>] priority=2","136Z.>] priority=1", "136Z.>] priority=2", "D", "N", "Y          80", "Y0100805", "YBO005", "YZZZZZZZ", "<", ">")] <- NA
quant_gp$value2[quant_gp$value2 %in% c("ADV001", "^", "136..>] priority=1","136..>] priority=2","136Z.>] priority=1", "136Z.>] priority=2", "D", "N", "Y          80", "Y0100805", "YBO005", "YZZZZZZZ", "<", ">")] <- NA
quant_gp$quant_gp <- as.numeric(apply(quant_gp[,c("value1","value2")],1,max,na.rm=T))
quant_gp$quant_gp[quant_gp$quant_gp > 600] <- NA #remove extreme/implausible outliers
quant_gp$quant_gp_ln <- log(quant_gp$quant_gp+1)
hist(quant_gp$quant_gp_ln)

audit$audit <- as.numeric(audit$value1); audit$audit[audit$audit > 40] <- NA
audit$audit_ln <- log(audit$audit+1)
hist(audit$audit_ln)

auditc$auditc <- round(as.numeric(auditc$value1)); auditc$auditc[auditc$auditc > 12] <- NA
auditc$auditc_ln <- log(auditc$auditc+1)
hist(auditc$auditc_ln)

auditpic$auditpic <- as.numeric(auditpic$value1); auditpic$auditpic[auditpic$auditpic > 16] <- NA
auditpic$auditpic_ln <- log(auditpic$auditpic+1)
hist(auditpic$auditpic_ln)

fast$fast <- as.numeric(fast$value1); fast$fast[fast$fast > 16] <- NA
fast$fast_ln <- log(fast$fast+1)
hist(fast$fast_ln)


#remove duplicates
quant_gp <- quant_gp[!duplicated(eid),]
audit <- audit[!duplicated(eid)]
auditc <- auditc[!duplicated(eid)]
auditpic <- auditpic[!duplicated(eid)]
fast <- fast[!duplicated(eid)]
liver <- liver[!duplicated(eid)]
korsakoff <- korsakoff[!duplicated(eid)]
rehab <- rehab[!duplicated(eid)]
withdrawal <- withdrawal[!duplicated(eid)]
abuse <- abuse[!duplicated(eid)]
aud <- aud[!duplicated(eid)]
chronicaud <- chronicaud[!duplicated(eid)]
advice <- advice[!duplicated(eid)]
excess <- excess[!duplicated(eid)]
drinkerstatus <- drinkerstatus[!duplicated(eid)]
acute <- acute[!duplicated(eid)]
otherclin <- otherclin[!duplicated(eid)]
rx <- rx[!duplicated(eid)]



# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------
# Combine GP data with clinical info from other sources (hospital records; interviews)
### Assumes the following files exist in the working directory:
### (3) rx_phenos - individual level UKB data for field 20003 (Treatment/medication code)
### (4) icd_phenos - individual level UKB data for fields 41202	(Diagnoses - main ICD10), 41204	(Diagnoses - secondary ICD10), 40001 (Primary cause of death - ICD10), 40002 (Secondary cause of death - ICD10)
### (5) iv_dxs_phenos -  individual level UKB data for fields 2473 (Medical conditions - self report) and 20002 (Non-cancer illness code)
### (6) icd_coding.tsv - UKB data coding 19 (https://biobank.ctsu.ox.ac.uk/showcase/coding.cgi?id=19)
### (7) iv_coding.tsv - UKB data coding 6 (https://biobank.ctsu.ox.ac.uk/showcase/coding.cgi?id=6)
# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------

# function to identify whether a certain diagnostic code is present across a set of fields/array value
isin <- function(data, target) {
  if (!FALSE %in% is.na(data)) {return(0)}
  if (!FALSE %in% (data=="")) {return(0)}
  return(as.numeric(TRUE %in% (data %in% target)))
}

# Read in prescription drug information reported in medical interview
ivrx <- fread("rx_phenos")
vars.ivrx <- grep("20003",names(ivrx),val=T)

# Identify prescriptions for alcohol-related conditions (see UKB data coding scheme 4 - https://biobank.ctsu.ox.ac.uk/showcase/coding.cgi?id=4&nl=1)
codes <- list(alc_rx=list(ivrx=c(1140926990,1140926994,1141157338,1140872480,1140872484)))
for (i in 1:length(codes)) {
  ivrx[,names(codes)[i]:=apply(ivrx[,..vars.ivrx],1, isin, target=codes[[i]]$ivrx)]
}
ivrx <- ivrx[alc_rx==1,]
table(ivrx$f.eid %in% rx$eid)


# Read in ICD codes reported in hospital records and medical interview
icd <- fread("icd_phenos")
icd <- icd[!(f.40001.0.0=="" & f.41202.0.0==""),]     # remove individuals with missing records for this source
dim(icd)
iv <- fread("iv_dxs_phenos")
iv <- iv[!is.na(iv$f.2473.0.0),]     # remove individuals with missing records for this source
dim(iv)

vars.icd <- grep(paste("41202","41204","40001","40002",sep="|"),names(icd),val=T)
vars.iv <- grep("20002",names(iv),val=T)
codes.icd <- fread("icd_coding.tsv")
codes.iv <- fread("iv_coding.tsv")

# Define list of alcohol-related ICD codes to search for (matching GP clinical code definitions where relevant)
codes <- list(acute=list(icd="F100"),harmful=list(icd="F101"),dependence=list(icd="F102",iv=c(1408)),withdrawal=list(icd=c("F103","F104")),korsakoff=list(icd=c("F105","F106","F107","E510","E512","E518","E519")),mentbehav=list(icd=c("F108","F109")),otheralc=list(icd="F10"),liver=list(icd=c("K70","K700","K701","K702","K703","K704","K709"),iv=c(1604)),rehab=list(icd="Z502"),pershistory=list(icd="Z864"),famhistory=list(icd="Z811"),toxicity=list(icd=c("T510","T519","X450","X451","X452","X453","X454","X455","X458","X459")),bac=list(icd="R780"))

for (i in 1:length(codes)) {
  icd[,names(codes)[i]:=apply(icd[,..vars.icd],1, isin, target=codes[[i]]$icd)]
  iv[,names(codes)[i]:=apply(iv[,..vars.iv],1, isin, target=codes[[i]]$iv)]
}

for(j in names(codes)){print(j);print(table(icd[,..j]))}
icd <- icd[,c(1,299:311)]
iv <- iv[,c(1,144,149)]




# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------
# Combine clinical data with self-report alcohol measures
### Assumes the following files exist in the working directory:
### (8) survey_phenos - individual level UKB data for fields listed in (9)
### (9) SavageAJP2024_fieldcodes.txt - field codes and name mappings for (8)
# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------


# Read phenotype files into R & recode varnames
alc <- fread("survey_phenos",na.strings = c("NA","-3","-1","-6","-818",'-819',"-121"))

codes <- fread("SavageAJP2024_fieldcodes.txt",fill=T)
for(i in 1:length(codes$V1)) {
  names(alc)[grep(codes$V1[i],names(alc))] <- str_replace(names(alc)[grep(codes$V1[i],names(alc))], paste0("f.",codes$V1[i]), codes$V2[i])
}


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
# Rescale to standard drinks
alc$quant_rwine_mo.0.0 <- alc$quant_rwine_mo.0.0*1.5
alc$quant_wwine_mo.0.0 <- alc$quant_wwine_mo.0.0*1.5
alc$quant_beer_mo.0.0 <- alc$quant_beer_mo.0.0*2.3
alc$quant_fwine_mo.0.0 <- alc$quant_fwine_mo.0.0*1.25
# Sum drinks per month (exclude missing alls)
alc$quant_allmiss_mo <- apply(apply(alc[,grep("quant.*mo.0",names(alc),val=T)],1,is.na),2,sum)
alc$quant_mo.0.0 <- NA
alc$quant_mo.0.0[alc$quant_allmiss_mo != 6] <- apply(alc[alc$quant_allmiss_mo != 6,grep("quant.*mo.0",names(alc),val=T)],1,sum,na.rm=T)
alc$quant_mo.0.0[which(alc$quant_mo.0.0 > 600)] <- 600 #Cap max drinks/mo at 600, or 20 drinks per day

## Drinks per week
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
alc$quant_wk.0.0[which(alc$quant_wk.0.0 > 600)] <- 600 #Cap max drinks/mo at 600, or 20 drinks per day

#calculate gm/eth per day
#8gm per standard drink in uk
alc$quant_ts0.0 <- NA
alc$quant_ts0.0[!(is.na(alc$quant_mo.0.0) & is.na(alc$quant_wk.0.0))] <- apply(alc[!(is.na(alc$quant_mo.0.0) & is.na(alc$quant_wk.0.0)),c("quant_mo.0.0","quant_wk.0.0")],1,sum,na.rm=T)*8/30



### Repeat for visit 2
alc$quant_rwine_mo.1.0 <- alc$quant_rwine_mo.1.0*1.5
alc$quant_wwine_mo.1.0 <- alc$quant_wwine_mo.1.0*1.5
alc$quant_beer_mo.1.0 <- alc$quant_beer_mo.1.0*2.3
alc$quant_fwine_mo.1.0 <- alc$quant_fwine_mo.1.0*1.25
alc$quant_allmiss_mo <- apply(apply(alc[,grep("quant.*mo.1",names(alc),val=T)],1,is.na),2,sum)
alc$quant_mo.1.0 <- NA
alc$quant_mo.1.0[alc$quant_allmiss_mo != 6] <- apply(alc[alc$quant_allmiss_mo != 6,grep("quant.*mo.1",names(alc),val=T)],1,sum,na.rm=T)
alc$quant_mo.1.0[which(alc$quant_mo.1.0 > 600)] <- 600 #Cap max drinks/mo at 900, or 20 drinks per day

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
alc$quant_wk.1.0[which(alc$quant_wk.1.0 > 600)] <- 600 #Cap max drinks/mo at 900, or 20 drinks per day

alc$quant_ts1.0 <- NA
alc$quant_ts1.0[!(is.na(alc$quant_mo.1.0) & is.na(alc$quant_wk.1.0))] <- apply(alc[!(is.na(alc$quant_mo.1.0) & is.na(alc$quant_wk.1.0)),c("quant_mo.1.0","quant_wk.1.0")],1,sum,na.rm=T)*8/30


### Repeat for visit 3
alc$quant_rwine_mo.2.0 <- alc$quant_rwine_mo.2.0*1.5
alc$quant_wwine_mo.2.0 <- alc$quant_wwine_mo.2.0*1.5
alc$quant_beer_mo.2.0 <- alc$quant_beer_mo.2.0*2.3
alc$quant_fwine_mo.2.0 <- alc$quant_fwine_mo.2.0*1.25
alc$quant_allmiss_mo <- apply(apply(alc[,grep("quant.*mo.2",names(alc),val=T)],1,is.na),2,sum)
alc$quant_mo.2.0 <- NA
alc$quant_mo.2.0[alc$quant_allmiss_mo != 6] <- apply(alc[alc$quant_allmiss_mo != 6,grep("quant.*mo.2",names(alc),val=T)],1,sum,na.rm=T)
alc$quant_mo.2.0[which(alc$quant_mo.2.0 > 600)] <- 600 #Cap max drinks/mo at 900, or 20 drinks per day

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
alc$quant_wk.2.0[which(alc$quant_wk.2.0 > 600)] <- 600 #Cap max drinks/mo at 900, or 20 drinks per day

alc$quant_ts2.0 <- NA
alc$quant_ts2.0[!(is.na(alc$quant_mo.2.0) & is.na(alc$quant_wk.2.0))] <- apply(alc[!(is.na(alc$quant_mo.2.0) & is.na(alc$quant_wk.2.0)),c("quant_mo.2.0","quant_wk.2.0")],1,sum,na.rm=T)*8/30


# fill in missing data with valid data from later assessments, if available
table(is.na(alc$quant_ts0.0),is.na(alc$quant_ts1.0),is.na(alc$quant_ts2.0))
alc$quant_ts <- alc$quant_ts0.0
alc$quant_ts[is.na(alc$quant_ts)] <- alc$quant_ts1.0[is.na(alc$quant_ts)]  
alc$quant_ts[is.na(alc$quant_ts)] <- alc$quant_ts2.0[is.na(alc$quant_ts)]  

# log transform skewed data
alc$quant_ts_ln <- log(alc$quant_ts+1)
hist(alc$quant_ts_ln)



# ----------------------------------------------------
### Beverage-specific quantity variables

alc$quantbeer <- log(alc$quant_beer_mo.0.0+1)
alc$quantbeer[is.na(alc$quantbeer)] <- log((alc$quant_beer_wk.0.0[is.na(alc$quantbeer)]*4)+1)
alc$quantrwine <- log(alc$quant_rwine_mo.0.0+1)
alc$quantrwine[is.na(alc$quantrwine)] <- log((alc$quant_rwine_wk.0.0[is.na(alc$quantrwine)]*4)+1)
alc$quantwwine <- log(alc$quant_wwine_mo.0.0+1) 
alc$quantwwine[is.na(alc$quantwwine)] <- log((alc$quant_wwine_wk.0.0[is.na(alc$quantwwine)]*4)+1)
alc$quantfwine <- log(alc$quant_fwine_mo.0.0+1)
alc$quantfwine[is.na(alc$quantfwine)] <- log((alc$quant_fwine_wk.0.0[is.na(alc$quantfwine)]*4)+1)
alc$quantspirit <- log(alc$quant_spirit_mo.0.0+1) 
alc$quantspirit[is.na(alc$quantspirit)] <- log((alc$quant_spirit_wk.0.0[is.na(alc$quantspirit)]*4)+1)
alc$quantother <- log(alc$quant_other_mo.0.0+1)
alc$quantsother[is.na(alc$quantother)] <- log((alc$quant_other_wk.0.0[is.na(alc$quantother)]*4)+1)



# ----------------------------------------------------
### Recode typical frequency vars to days/month (non-drinkers missing)

alc$drinkfreq.0.0 <- recode(alc$drinkfreq.0.0, "1=20;2=14;3=6;4=2;5=.5;6=NA")
alc$drinkfreq.1.0 <- recode(alc$drinkfreq.1.0, "1=20;2=14;3=6;4=2;5=.5;6=NA")
alc$drinkfreq.2.0 <- recode(alc$drinkfreq.2.0, "1=20;2=14;3=6;4=2;5=.5;6=NA")

# fill in missing data with valid data from later assessments, if available
alc$drinkfreq <- alc$drinkfreq.0.0
alc$drinkfreq[is.na(alc$drinkfreq)] <- alc$drinkfreq.1.0[is.na(alc$drinkfreq)]  
alc$drinkfreq[is.na(alc$drinkfreq)] <- alc$drinkfreq.2.0[is.na(alc$drinkfreq)]  



# ----------------------------------------------------
### AUDIT-C and AUDIT-P items/scores from mental health questionnaire

### Create total grams/month quantity variable for mh questionnaire measure
alc$mh_quantR <- recode(alc$mh_quant.0.0, "1=1.5;2=3.5;3=5.5;4=8;5=10")
alc$mh_quantR[alc$mh_freq.0.0==0 & !(alc$drinkerstatus.0.0 %in% c(0,1))] <- 0
hist(alc$mh_quantR)
alc$mh_freqR <- recode(alc$mh_freq.0.0, "0=0;1=1;2=3;3=10;4=16")
hist(alc$mh_freqR,breaks = 30)
alc$mh_quant_ln <- log((alc$mh_freqR*alc$mh_quantR*8)+1)
hist(alc$mh_quant_ln)

### Recode binge frequency vars to days/month
#recode freq as 0 in those who answered mental health frequency question as "never" and exclude abstainers
alc$mh_bingefreq.0.0[alc$mh_freq.0.0==0 & !(alc$drinkerstatus.0.0 %in% c(0,1))] <- 1
alc$mh_bingefreq.0.0[alc$mh_bingefreq.0.0==1 & (alc$drinkerstatus.0.0 %in% c(0,1))] <- NA
alc$mh_binge_ln <- log(recode(alc$mh_bingefreq.0.0, "1=0;2=.5;3=1;4=4;5=20")+1)

alc$mh_auditc <- NA
alc$mh_auditp <- NA
#rescale to match audit scoring rules
alc$mh_quant.0.0 <- alc$mh_quant.0.0-1
alc$mh_bingefreq.0.0 <- alc$mh_bingefreq.0.0-1
alc$aud_responsibilities.0.0 <- alc$aud_responsibilities.0.0-1
alc$aud_blackouts.0.0 <- alc$aud_blackouts.0.0-1
alc$aud_guilt.0.0 <- alc$aud_guilt.0.0-1
alc$aud_withdrawal.0.0 <- alc$aud_withdrawal.0.0-1
alc$aud_cantstop.0.0 <- alc$aud_cantstop.0.0-1
alc$concerned_friends.0.0 <- alc$concerned_friends.0.0*2
alc$injured_drinking.0.0 <- alc$injured_drinking.0.0*2
audit.c <- c("mh_quant.0.0","mh_freq.0.0","mh_bingefreq.0.0")
audit.p <- c("aud_responsibilities.0.0", "aud_blackouts.0.0", "aud_guilt.0.0", "concerned_friends.0.0", "injured_drinking.0.0", "aud_withdrawal.0.0", "aud_cantstop.0.0")
alc <- as.data.frame(alc)
# Recode individuals who took the MH survey but were skipped out of AUDIT-P questions 
for (i in auditp) {
  alc[which(is.na(alc[,i]) & !is.na(alc$mh_freq.0.0) ),i] <- 0
}
# Sum 
alc$mh_auditp[which(!is.na(alc$mh_freq.0.0))] <- apply(alc[which(!is.na(alc$mh_freq.0.0)), audit.p],1,sum,na.rm=T)
alc$mh_auditc[which(!is.na(alc$mh_freq.0.0))] <- apply(alc[which(!is.na(alc$mh_freq.0.0)), audit.c],1,sum,na.rm=T)





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


# ----------------------------------------------------
# Self-report alcohol addiction phenotypes
# Recode  for individuals who took the MH survey but were skipped out of AUDIT-P questions 
alc[which(is.na(alc$phys_alcohol.0.0) & !is.na(alc$mh_freq.0.0) ),'phys_alcohol.0.0'] <- 0
alc[which(is.na(alc$ever_alcohol.0.0) & !is.na(alc$mh_freq.0.0) ),'ever_alcohol.0.0'] <- 0
alc[which(is.na(alc$current_alcohol.0.0) & !is.na(alc$mh_freq.0.0) ),'current_alcohol.0.0'] <- 0



# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------
### create phenotype file
# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------

phen <- alc[,c("f.eid","quant_ts_ln","quantbeer", "quantrwine","quantwwine","quantfwine","quantspirit","quantother",
               "mh_auditc","mh_auditp","mh_binge_ln","mh_quant_ln",
               "drinkfreq","drink_w_meals.0.0","increasedrink", "decreasedrink", 
               "phys_alcohol.0.0", "ever_alcohol.0.0", "current_alcohol.0.0")]
phen <- merge(phen, quant_gp[,c("eid","quant_gp_ln")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, audit[,c("eid","audit_ln")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, auditc[,c("eid","auditc_ln")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, auditpic[,c("eid","auditpic_ln")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, fast[,c("eid","fast_ln")], by.x="f.eid",by.y="eid",all=T)


phen <- merge(phen, liver[,c("eid","liver")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, korsakoff[,c("eid","korsakoff")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, rehab[,c("eid","rehab")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, withdrawal[,c("eid","withdrawal")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, abuse[,c("eid","abuse")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, aud[,c("eid","aud")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, chronicaud[,c("eid","chronicaud")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, advice[,c("eid","advice")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, audit[,c("eid","audit")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, excess[,c("eid","excess")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, drinkerstatus[,c("eid","drinkerstatus")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, acute[,c("eid","acute")], by.x="f.eid",by.y="eid",all=T)
phen <- merge(phen, otherclin[,c("eid","otherclin")], by.x="f.eid",by.y="eid",all=T)

phen <- merge(phen, icd[,c("f.eid","pershistory")], by="f.eid",all=T)
phen <- merge(phen, icd[,c("f.eid","toxicity")], by="f.eid",all=T)

phen <- merge(phen, rx[,c("eid","alc_rx")], by.x="f.eid",by.y="eid",all=T)

# Code cases as 1 where needed/missing
phen$acute[phen$f.eid %in% icd$f.eid[which(icd$acute==1)]] <- 1
phen$abuse[phen$f.eid %in% icd$f.eid[which(icd$harmful==1)]] <- 1
phen$aud[phen$f.eid %in% icd$f.eid[which(icd$dependence==1)]] <- 1
phen$withdrawal[phen$f.eid %in% icd$f.eid[which(icd$withdrawal==1)]] <- 1
phen$korsakoff[phen$f.eid %in% icd$f.eid[which(icd$korsakoff==1)]] <- 1
phen$otherclin[phen$f.eid %in% icd$f.eid[which(icd$mentbehav==1)]] <- 1
phen$otherclin[phen$f.eid %in% icd$f.eid[which(icd$otheralc==1)]] <- 1
phen$liver[phen$f.eid %in% icd$f.eid[which(icd$liver==1)]] <- 1
phen$rehab[phen$f.eid %in% icd$f.eid[which(icd$rehab==1)]] <- 1

phen$toxicity[phen$f.eid %in% clin3$eid[which(clin3$desc %in% c("Ethanol causing toxic effect","Ethyl alcohol causing toxic effect","Alcohol causing toxic effect"))]] <- 1

phen$aud[phen$f.eid %in% iv$f.eid[which(iv$dependence==1)]] <- 1
phen$liver[phen$f.eid %in% iv$f.eid[which(iv$liver==1)]] <- 1

phen$alc_rx[phen$f.eid %in% ivrx$f.eid[which(ivrx$alc_rx==1)]] <- 1


###broad aud case definition to also include prescriptions and hospital diagnoses and self-report (and optionally personal history of alcoholism ICD code)
phen$broad_aud <- phen$aud
phen$broad_aud[phen$liver==1|phen$korsakoff==1|phen$rehab==1|phen$withdrawal==1|phen$alc_rx==1|phen$phys_alcohol.0.0==1|phen$ever_alcohol.0.0==1] <- 1
phen$broad_aud_ph <- phen$broad_aud
phen$broad_aud_ph[phen$pershistory==1] <- 1

###broad clinical consequences of alc phenotype
phen$anyclin <- phen$broad_aud
phen$anyclin[phen$abuse==1|phen$advice==1|phen$excess==1|phen$acute==1|phen$alc_rx==1|phen$otherclin==1|phen$toxicity==1] <- 1


###cobmine self-report and gp-report of weekly quantity to maximize sample size
phen$quant_ts_gp <- apply(phen[,c("quant_ts_ln","quant_gp_ln")],1,max,na.rm=T)

