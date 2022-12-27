library(data.table)
library(TwoSampleMR)


Hypo_exp <-fread('20002_1226_logistic.EUR.sumstats.MACfilt.txt',header=T,sep='\t')
exp_dat <- format_data(
    Hypo_exp,
    type = 'exposure',
    snp_col = "SNPID_UKB",
    beta_col = "OR",
    se_col = "SE",
    effect_allele_col = "A1_UKB",
    other_allele_col = "A2_UKB",
    eaf_col = "MAF_UKB",
    pval_col = "P",
    chr = "CHR",
    pos = "BP",
    ncase = "18740",
    ncontrol = "270567",
)

exp_dat$exposure <- "Hypothyroidism"

exp_dat_ob <- exp_dat[exp_dat$pval.exposure < 5e-8,]
exp_dat_ob <- na.omit(exp_dat_ob)


Hypo_exp_dat <- clump_data(exp_dat_ob, clump_r2=0.01, clump_kb = 5000, pop = "EUR")
Hypo_exp_dat_new <- Hypo_exp_dat[Hypo_exp_dat$beta.exposure > -0.07,]
dat <- harmonise_data(exposure_dat = Hypo_exp_dat_new, outcome_dat = outcome_dat, action=2)



#N <- 270567+18740
#r_2 <- 1/(1+(N*Hypo_exp_dat$se.exposure^2/Hypo_exp_dat$beta.exposure^2))
#F_score <- r_2*(N-2)/(1-r_2)


outcome_dat <- read_outcome_data(
    filename = "34594039-GCST90018690-EFO_0000685.h.tsv",
    snps = Hypo_exp_dat$SNP,
    sep = "\t",
    snp_col = "variant_id",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col ="effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_freq",
    pval_col = "p-value",
    chr = "chromosome",
    pos = "base_pair_location",
    ncase = "5348",
    ncontrol = "173268",
    phenotype = "Rheumatoid arthritis"
    
)
outcome_dat$outcome <- "Rheumatoid arthritis"

dat <- harmonise_data(exposure_dat = Hypo_exp_dat, outcome_dat = outcome_dat, action=2)
res <- mr(dat)
generate_odds_ratios(res)



