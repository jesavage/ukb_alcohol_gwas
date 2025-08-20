# Genomic SEM script carried out after item-level GWAS has been run for all candidate phenotypes

### FORMAT SUMSTAT FILES

#for binary variables, effective N needs to be calculated from c/c allele count columns (=4/((1/(N_case_alleles/2)) + (1/(N_control_alleles/2)). then sample prevalence for h2/rg should be set at 0.5 when using the effective N

for i in audit_ln auditc_ln consumption_ln drinkerstatus drinkfreq fast_ln gm_eth_ln mh_auditc mh_auditp mh_binge_ln quantbeer quantrwine quantwwine quantfwine quantspirit quant_ts_gp; do
awk '{print $16"\t"$2"\t"$3"\t"$4"\t"$6"\t"$8"\t"$9"\t"$10"\t"$14"\t"$15}' ${i}/${i}_EUR_sumstats.txt > ${i}/${i}_EUR_sumstats_short.txt
sed -i '1d' ${i}/${i}_EUR_sumstats_short.txt
sed -i '1iSNP\tCHR\tPOS\tA2\tA1\tN\tBETA\tSE\tP\tMAF' ${i}/${i}_EUR_sumstats_short.txt
done

for i in abuse advice anyclin aud broad_aud broad_aud_ph decreasedrink drink_w_meals.0.0 ever_alcohol.0.0 excess increasedrink pershistory; do
cp ${i}/${i}_EUR_sumstats.txt ${i}/tmp
sed -i '1d' ${i}/tmp
awk '{print $18"\t"$2"\t"$3"\t"$4"\t"$6"\t"$11"\t"$12"\t"$16"\t"$17"\t"($7/2)"\t"($8/2)"\t"4/((1/($7/2))+(1/($8/2)))}' ${i}/tmp > ${i}/${i}_EUR_sumstats_short.txt
sed -i '1iSNP\tCHR\tPOS\tA2\tA1\tOR\tSE\tP\tMAF\tNcase\tNcontrol\tNeff' ${i}/${i}_EUR_sumstats_short.txt
rm ${i}/tmp
done


#munge sumstats 

gzip */*sumstats*.txt

LDSC_PATH=/programs/ldsc
LDSC_REF_PATH=/programs/LDSC_reffiles
SSDIR=${LDSC_REF_PATH}/SUMSTATS_OTHERTRAITS
POP="EUR"
LDSCOREFILE="${LDSC_REF_PATH}/eur_w_ld_chr/"

for i in abuse aud advice anyclin broad_aud broad_aud_ph decreasedrink drink_w_meals.0.0 excess ever_alcohol.0.0 increasedrink pershistory audit_ln auditc_ln consumption_ln drinkerstatus drinkfreq fast_ln gm_eth_ln mh_auditc mh_auditp mh_binge_ln quantbeer quantrwine quantwwine quantfwine quantspirit quant_ts_gp; do 
PHEN=${i}
${LDSC_PATH}/munge_sumstats.py \
--sumstats ${PHEN}/${PHEN}_EUR_sumstats_short.txt.gz \
--merge-alleles ${LDSC_REF_PATH}/w_hm3.snplist \
--out ${PHEN}/${PHEN}_EUR_munge 
done

#check h2 to identify low-performing items
#quantitative phenotypes
for i in audit_ln auditc_ln consumption_ln drinkerstatus drinkfreq fast_ln gm_eth_ln mh_auditc mh_auditp mh_binge_ln quantbeer quantrwine quantwwine quantfwine quantspirit quant_ts_gp; do
PHEN=${i}
${LDSC_PATH}/ldsc.py \
--h2 ${PHEN}/${PHEN}_EUR_munge.sumstats.gz \
--ref-ld-chr ${LDSCOREFILE} \
--w-ld-chr ${LDSCOREFILE} \
--out ${PHEN}/${PHEN}_EUR_h2 
done

#case control phenotypes
# assumes text file 'popprev' with the sample and population prevalence estimates per phenotype; see supplementary tables
while IFS=$'\t' read -r PHEN SAMPPREV POPPREV SOURCE; do
echo "PHEN: $PHEN"
echo "SAMPPREV: $SAMPPREV"
echo "POPPREV: $POPPREV"
PHEN=${i}

#liability scale est
${LDSC_PATH}/ldsc.py \
--h2 ${PHEN}/${PHEN}_EUR_munge.sumstats.gz \
--ref-ld-chr ${LDSCOREFILE} \
--w-ld-chr ${LDSCOREFILE} \
--samp-prev 0.5 \
--pop-prev ${POPPREV} \
--out ${PHEN}/${PHEN}_EUR_h2

#observed scale est
${LDSC_PATH}/ldsc.py \
--h2 ${PHEN}/${PHEN}_EUR_munge.sumstats.gz \
--ref-ld-chr ${LDSCOREFILE} \
--w-ld-chr ${LDSCOREFILE} \
--out ${PHEN}/${PHEN}_EUR_h2_obsscale

done < popprev




### CALCULATE GENETIC COVARIANCE

Rscript - << 'END'
#require(devtools)
#install_github("GenomicSEM/GenomicSEM")
require(GenomicSEM)
print(packageVersion("GenomicSEM"))

#using already munged sumstats from ldsc
#otherwise go through munging steps within gSEM or in ldsc first

#calculate multivariable genetic covariance matrix
traits <- c("abuse/abuse_EUR_munge.sumstats.gz",
"advice/advice_EUR_munge.sumstats.gz",
"anyclin/anyclin_EUR_munge.sumstats.gz",
"aud/aud_EUR_munge.sumstats.gz",
"audit_ln/audit_ln_EUR_munge.sumstats.gz",
"auditc_ln/auditc_ln_EUR_munge.sumstats.gz",
"broad_aud/broad_aud_EUR_munge.sumstats.gz",
"broad_aud_ph/broad_aud_ph_EUR_munge.sumstats.gz",
"consumption_ln/consumption_ln_EUR_munge.sumstats.gz",
"decreasedrink/decreasedrink_EUR_munge.sumstats.gz",
"drink_w_meals.0.0/drink_w_meals.0.0_EUR_munge.sumstats.gz",
"drinkerstatus/drinkerstatus_EUR_munge.sumstats.gz",
"drinkfreq/drinkfreq_EUR_munge.sumstats.gz",
"ever_alcohol.0.0/ever_alcohol.0.0_EUR_munge.sumstats.gz",
"excess/excess_EUR_munge.sumstats.gz",
"fast_ln/fast_ln_EUR_munge.sumstats.gz",
"gm_eth_ln/gm_eth_ln_EUR_munge.sumstats.gz",
"increasedrink/increasedrink_EUR_munge.sumstats.gz",
"mh_auditc/mh_auditc_EUR_munge.sumstats.gz",
"mh_auditp/mh_auditp_EUR_munge.sumstats.gz",
"mh_binge_ln/mh_binge_ln_EUR_munge.sumstats.gz",
"pershistory/pershistory_EUR_munge.sumstats.gz",
"quantbeer/quantbeer_EUR_munge.sumstats.gz",
"quantfwine/quantfwine_EUR_munge.sumstats.gz",
"quantrwine/quantrwine_EUR_munge.sumstats.gz",
"quantwwine/quantwwine_EUR_munge.sumstats.gz",
"quantspirit/quantspirit_EUR_munge.sumstats.gz")

pop.prev <- c(.13,.145,.3,.038,NA,NA,.047,.142,NA,.545,.675,NA,NA,.038,.24,NA,NA,.301,NA,NA,NA,.138,NA,NA,NA,NA,NA)
samp.prev <- c(.5,.5,.5,.5,NA,NA,.5,.5,NA,.5,.5,NA,NA,.5,.5,NA,NA,.5,NA,NA,NA,.5,NA,NA,NA,NA,NA)
trait.names <- c("abuse","advice","anyclin","aud","audit_ln","auditc_ln","broad_aud","broad_aud_ph","consumption_ln","decreasedrink","drink_w_meals.0.0","drinkerstatus","drinkfreq","ever_alcohol.0.0","excess","fast_ln","gm_eth_ln","increasedrink","mh_auditc","mh_auditp","mh_binge_ln","pershistory","quantbeer","quantfwine","quantrwine","quantwwine","quantspirit")
ld <- "/programs/LDSC_reffiles/eur_w_ld_chr/"
wld <- "/programs/LDSC_reffiles/eur_w_ld_chr/"
 
LDSCoutput <- ldsc(traits, samp.prev, pop.prev, ld, wld, trait.names)
save(LDSCoutput, file="gsem/alc_items_gencov_matrix_allitems.RData")




### Based on LDSC results, remove audit_ln, drinkerstatus, ever_alcohol.0.0, excess, and fast_ln (n.s. h2), as well as consumption_ln+gm_eth_ln (rg=.98, use combined quant_ts_gp instead), aud (rg=.96, use broad_aud instead), broad_aud_ph (rg=.99, use pershistory instead), mh_auditc (rg=.93 with gm_eth_ln and .94 with mh_binge_ln; quant, freq, and binge are all already captured by individual items)

#calculate multivariable genetic covariance matrix
traits <- c("abuse/abuse_EUR_munge.sumstats.gz",
"advice/advice_EUR_munge.sumstats.gz",
"anyclin/anyclin_EUR_munge.sumstats.gz",
"auditc_ln/auditc_ln_EUR_munge.sumstats.gz",
"broad_aud/broad_aud_EUR_munge.sumstats.gz",
"decreasedrink/decreasedrink_EUR_munge.sumstats.gz",
"drink_w_meals.0.0/drink_w_meals.0.0_EUR_munge.sumstats.gz",
"drinkfreq/drinkfreq_EUR_munge.sumstats.gz",
"increasedrink/increasedrink_EUR_munge.sumstats.gz",
"mh_auditp/mh_auditp_EUR_munge.sumstats.gz",
"mh_binge_ln/mh_binge_ln_EUR_munge.sumstats.gz",
"pershistory/pershistory_EUR_munge.sumstats.gz",
"quantbeer/quantbeer_EUR_munge.sumstats.gz",
"quantfwine/quantfwine_EUR_munge.sumstats.gz",
"quantrwine/quantrwine_EUR_munge.sumstats.gz",
"quant_ts_gp/quant_ts_gp_EUR_munge.sumstats.gz",
"quantwwine/quantwwine_EUR_munge.sumstats.gz",
"quantspirit/quantspirit_EUR_munge.sumstats.gz")

pop.prev <- c(.13,.145,.3,NA,.047,.545,.675,NA,.301,NA,NA,.138,NA,NA,NA,NA,NA,NA)
samp.prev <- c(.5,.5,.5,NA,.5,.5,.5,NA,.5,NA,NA,.5,NA,NA,NA,NA,NA,NA)
trait.names <- c("abuse","advice","anyclin","auditc_ln","broad_aud","decreasedrink","drink_w_meals.0.0","drinkfreq","increasedrink","mh_auditp","mh_binge_ln","pershistory","quantbeer","quantfwine","quantrwine","quant_ts_gp","quantwwine","quantspirit")
ld <- "/programs/LDSC_reffiles/eur_w_ld_chr/"
wld <- "/programs/LDSC_reffiles/eur_w_ld_chr/"
 
LDSCoutput <- ldsc(traits, samp.prev, pop.prev, ld, wld, trait.names)
save(LDSCoutput, file="/gsem/alc_items_gencov_matrix.RData")




### MODELS

#smooth cov matrix to make positive definite
require(Matrix)
Ssmooth<-as.matrix((nearPD(LDSCoutput$S, corr = FALSE))$mat)

#eigenvalues
eigen(cov2cor(Ssmooth))$values
#[1] 7.401667e+00 4.956589e+00 1.330556e+00 1.074839e+00 7.028927e-01 5.636834e-01 4.706933e-01 3.919352e-01
# [9] 3.465654e-01 2.690050e-01 1.705004e-01 1.392628e-01 1.093177e-01 5.704051e-02 1.545259e-02 2.730841e-07
#[17] 1.901485e-07 1.213296e-07
plot(1:18,eigen(cov2cor(Ssmooth))$values,type="b",ylab="Eigenvalues",xlab="")
abline(h=1,col="red")

#efa
require(stats)
EFA1<-factanal(covmat = Ssmooth, factors = 1)
EFA2<-factanal(covmat = Ssmooth, factors = 2, rotation = "promax")
EFA3<-factanal(covmat = Ssmooth, factors = 3, rotation = "promax")
EFA4<-factanal(covmat = Ssmooth, factors = 4, rotation = "promax")
EFA5<-factanal(covmat = Ssmooth, factors = 5, rotation = "promax")
EFA2orth<-factanal(covmat = Ssmooth, factors = 2, rotation = "varimax")
EFA3orth<-factanal(covmat = Ssmooth, factors = 3, rotation = "varimax")
EFA4orth<-factanal(covmat = Ssmooth, factors = 4, rotation = "varimax")
EFA5orth<-factanal(covmat = Ssmooth, factors = 5, rotation = "varimax")

# cfa - final model with unit factor loadings and highest h2 item listed first per factor - other specifications artificially inflate Neff estimates
CFA4unit_arg <- 
'F1 =~ mh_auditp + abuse + anyclin + broad_aud + increasedrink + mh_binge_ln + quant_ts_gp
F2 =~ drink_w_meals.0.0 + advice + anyclin + decreasedrink + pershistory + quantbeer + quantfwine + quantrwine + quantwwine
F3 =~ drinkfreq + auditc_ln + mh_auditp + mh_binge_ln + quantbeer + quantrwine + quant_ts_gp + quantwwine
F4 =~ quantspirit + auditc_ln + quantfwine'
CFA4unit <-usermodel(LDSCoutput, estimation = "DWLS", model = CFA4unit_arg, CFIcalc = TRUE, std.lv = FALSE, imp_cov = FALSE)
CFA4unit



### GWAS of factor models

# 1. Prep sumstats 

setwd("../")
files <- c("abuse/abuse_EUR_sumstats_short.txt.gz",
"advice/advice_EUR_sumstats_short.txt.gz",
"anyclin/anyclin_EUR_sumstats_short.txt.gz",
"auditc_ln/auditc_ln_EUR_sumstats_short.txt.gz",
"broad_aud/broad_aud_EUR_sumstats_short.txt.gz",
"decreasedrink/decreasedrink_EUR_sumstats_short.txt.gz",
"drink_w_meals.0.0/drink_w_meals.0.0_EUR_sumstats_short.txt.gz",
"drinkfreq/drinkfreq_EUR_sumstats_short.txt.gz",
"increasedrink/increasedrink_EUR_sumstats_short.txt.gz",
"mh_auditp/mh_auditp_EUR_sumstats_short.txt.gz",
"mh_binge_ln/mh_binge_ln_EUR_sumstats_short.txt.gz",
"pershistory/pershistory_EUR_sumstats_short.txt.gz",
"quantbeer/quantbeer_EUR_sumstats_short.txt.gz",
"quantfwine/quantfwine_EUR_sumstats_short.txt.gz",
"quantrwine/quantrwine_EUR_sumstats_short.txt.gz",
"quant_ts_gp/quant_ts_gp_EUR_sumstats_short.txt.gz",
"quantwwine/quantwwine_EUR_sumstats_short.txt.gz",
"quantspirit/quantspirit_EUR_sumstats_short.txt.gz")

ref="/ref/reference.1000G.maf.0.005.txt" #download from 1kG website

se.logit <- c(T,T,T,F,T,T,T,F,T,F,F,T,F,F,F,F,F,F)
ols <- c(F,F,F,T,F,F,F,T,F,T,T,F,T,T,T,T,T,T)

#give effective sample size for dichotomous traits, NA to use the column in sumstats file for continuous traits
#all sumstats have Neff calculated per-snp within file, so give NA for all phenos
ns <- rep(NA,length(files))

p_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS=ols,linprob=NULL,N=ns,maf.filter=.01,parallel=F)
save.image(file="gsem/alc_4facGWAS.RData")



## 2. Run GWAS 

model <-
'F1 =~ mh_binge_ln + abuse + anyclin + broad_aud + increasedrink + mh_auditp + quant_ts_gp
F2 =~ drink_w_meals.0.0 + advice + anyclin + decreasedrink + pershistory + quantbeer + quantfwine + quantrwine + quantwwine
F3 =~ drinkfreq + auditc_ln + mh_auditp + mh_binge_ln + quantbeer + quantrwine + quant_ts_gp + quantwwine
F4 =~ quantspirit + quantfwine + auditc_ln
F1 ~ SNP
F2~SNP
F3~SNP
F4~SNP
F1~~F2
F1~~F3
F1~~F4
F2~~F3
F2~~F4
F3~~F4'
fourFactorGWAS <- userGWAS(covstruc = LDSCoutput, SNPs = p_sumstats, estimation="DWLS", model=model, printwarn=TRUE, SNPSE = FALSE, smooth_check=T, parallel=F, cores=16)

write.table(fourFactorGWAS,"alc_4facGWAS.txt",na="NA",row.names=F)

END

# Extract and format sumstats for each latent factor
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' alc_4facGWASout.txt > alc_4facGWAS_f1.out
awk '{print $23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44}' alc_4facGWASout.txt > alc_4facGWAS_f2.out 
awk '{print $45,$46,$47,$48,$49,$50,$51,$52,$53,$54,$55,$56,$57,$58,$59,$60,$61,$62,$63,$64,$65,$66}' alc_4facGWASout.txt > alc_4facGWAS_f3.out
awk '{print $67,$68,$69,$70,$71,$72,$73,$74,$75,$76,$77,$78,$79,$80,$81,$82,$83,$84,$85,$86,$87,$88}' alc_4facGWASout.txt > alc_4facGWAS_f4.out 



