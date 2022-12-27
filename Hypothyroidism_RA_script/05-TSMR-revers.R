library(TwoSampleMR)
library(ggplot2)
exp_dat <- extract_instruments(outcomes = 'finn-b-M13_RHEUMA',clump = TRUE,r2 = 0.01,kb = 5000,access_token = NULL)
outcome_dat <- extract_outcome_data(exp_dat$SNP, outcomes = 'finn-b-HYPOTHY_REIMB', proxies = 1, rsq = 0.8, align_alleles = 1,  palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exp_dat, outcome_dat, action = 2)
dat$outcome <- "Hypothyroidism, drug reimbursement"
dat$exposure <- "Rheumatoid arthritis"
mr_result <- mr(dat)
odds_res_finn <- generate_odds_ratios(mr_result)
plt_finn <- mr_pleiotropy_test(dat)



exp_dat <- extract_instruments(outcomes = 'ebi-a-GCST000679',clump = TRUE,r2 = 0.01,kb = 5000,access_token = NULL)
outcome_dat <- extract_outcome_data(exp_dat$SNP, outcomes = 'ukb-a-77', proxies = 1, rsq = 0.8, align_alleles = 1,  palindromes = 1, maf_threshold = 0.3)
correction_factor <- abs(mean(outcome_dat$beta.outcome)/mean(exp_dat$beta.exposure))
exp_dat$beta.exposure <- exp_dat$beta.exposure * correction_factor
dat <- harmonise_data(exp_dat, outcome_dat, action = 2)
dat$outcome <- "Non-cancer illness code self-reported: hypothyroidism/myxoedema"
dat$exposure <- "Rheumatoid arthritis"
mr_result <- mr(dat)
odds_res_UKB77 <- generate_odds_ratios(mr_result)
plt_UKB77 <- mr_pleiotropy_test(dat)


exp_dat <- extract_instruments(outcomes = 'ebi-a-GCST000679',clump = TRUE,r2 = 0.01,kb = 5000,access_token = NULL)
outcome_dat <- extract_outcome_data(exp_dat$SNP, outcomes = 'ukb-a-190', proxies = 1, rsq = 0.8, align_alleles = 1,  palindromes = 1, maf_threshold = 0.3)
correction_factor <- abs(mean(outcome_dat$beta.outcome)/mean(exp_dat$beta.exposure))
exp_dat$beta.exposure <- exp_dat$beta.exposure * correction_factor
dat <- harmonise_data(exp_dat, outcome_dat, action = 2)
dat$outcome <- "Treatment/medication code: levothyroxine sodium"
dat$exposure <- "Rheumatoid arthritis"
mr_result <- mr(dat)
odds_res_UKB190 <- generate_odds_ratios(mr_result)
plt_UKB190 <- mr_pleiotropy_test(dat)

plt_result <- rbind(plt_finn,plt_UKB77,plt_UKB190)
odds_result <- rbind(odds_res_finn,odds_res_UKB77,odds_res_UKB190)
write.table(odds_result, "dat_MR_reverse_results.txt", row.names=FALSE, sep='\t')
write.table(plt_result, "dat_MR_reverse_plt_results.txt", row.names=FALSE, sep='\t')
