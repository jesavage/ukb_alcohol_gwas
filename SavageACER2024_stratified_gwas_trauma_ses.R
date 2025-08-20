require(data.table)
require(stringr)
require(car)

strat <- fread("strat_phenos.gz",na.strings = c("NA","-3","-2","-1","-6","-7","-818","-121"))
codes <- fread("strat_codes")
for (i in 1:length(codes$V1)) {
  names(strat)[grep(codes$V1[i],names(strat))] <- str_replace(names(strat)[grep(codes$V1[i],names(strat))], paste0("f.",codes$V1[i]), codes$V3[i])
}
summary(strat)

#recode reverse-coded trauma items
strat$tr_loved.0.0_R <- recode(strat$tr_loved.0.0, "1=4;2=3;3=2;4=1")
strat$tr_doctor.0.0_R <- recode(strat$tr_doctor.0.0, "1=4;2=3;3=2;4=1")
strat$tr_confide.0.0_R <- recode(strat$tr_confide.0.0, "1=4;2=3;3=2;4=1")
strat$tr_rent.0.0_R <- recode(strat$tr_rent.0.0, "1=4;2=3;3=2;4=1")


#sum & dichotomize trauma exposure (among MH questionnaire respondents)
trauma <- c("tr_child_physab.0.0","tr_child_sexab.0.0","tr_part_phys.0.0","tr_part_sex.0.0","tr_assault.0.0","tr_crime.0.0","tr_accident.0.0","tr_witness.0.0","tr_illness.0.0","tr_combat.0.0")
for (i in trauma) { print(table(strat[,..i]))}
strat$trauma_sum <- integer()
strat$trauma_sum[which(!is.na(strat$tr_loved.0.0))] <- apply(strat[which(!is.na(strat$tr_loved.0.0)),..trauma],1,sum,na.rm=T)
table(strat$trauma_sum)
hist(strat$trauma_sum)
hist(log(strat$trauma_sum+1))

strat$trauma_dich <- as.integer(strat$trauma_sum>0)
table(strat$trauma_dich)

#add alcohol vars
alc <- fread("alc_phenos") #see manuscript and prior publications for details on alc phenotype derivations
phen <- merge(alc[,c("f.eid","quant_ts_gp","broad_aud","broad_aud_ph","mh_binge_ln","drinkfreq","drink_w_meals.0.0")],strat[,c("f.eid","towndepr.0.0","trauma_sum","trauma_child","trauma_dich","tdi_dich")],by="f.eid",all=T)

#normalize interaction terms
phen$Zquant_ts_gp <- scale(phen$quant_ts_gp)
phen$Zmh_binge_ln <- scale(phen$mh_binge_ln)
phen$Zdrinkfreq <- scale(phen$drinkfreq)
phen$Ztowndepr <- scale(phen$towndepr.0.0)
phen$Ztrauma_sum_ln <- scale(log(phen$trauma_sum+1))

fwrite(phen[,c("IID","quant_ts_gp","broad_aud","broad_aud_ph","mh_binge_ln","drinkfreq","drink_w_meals.0.0","Zquant_ts_gp","Zmh_binge_ln","Zdrinkfreq")],"allstrat.pheno",sep="\t",quote=F,row.names=F,na=NA)
covs <- fread("/gwa_covs.txt.gz") #see manuscript for details on other covariates
covs <- covs[,1:36]
phen$FID <- phen$IID
covs <- merge(covs,phen[,c("FID","IID","towndepr.0.0","trauma_sum","trauma_child","trauma_dich","tdi_dich","Ztowndepr","Ztrauma_sum_ln")],by=c("FID","IID"))
fwrite(covs,"allstrat.covs",sep="\t",quote=F,row.names=F,na=NA)

#Run GWAS for each GxE pair
#example GWAS syntax:
# plink2 --bfile <chunk-prefix> \
#  --remove defaultSubjectExclusions.txt \
#  --exclude defaultVariantExclusions.txt \
#  --covar allstrat.covs \
#  --covar-name Ztowndepr,sex,pop_pc1-pop_pc5 \
#  --maf 0.001 \
#  --geno 0.05 \
#  --freq \
#  --ci 0.95 \
#  --pheno-name Zdrinkfreq \
#  --pheno allstrat.pheno \
#  --glm no-x-sex interaction \
#  --memory 2000

  
  
  
