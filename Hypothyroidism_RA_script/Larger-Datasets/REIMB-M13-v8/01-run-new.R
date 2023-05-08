library(data.table)
library(TwoSampleMR)
library(ggplot2)

Hypo_exp <-fread('finngen_R8_HYPOTHY_REIMB_use.tsv',header=T,sep='\t')
##Read Data
exp_dat <- format_data(
    Hypo_exp,
    type = 'exposure',
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    eaf_col = "af_alt",
    pval_col = "pval",
    chr = "#chrom",
    pos = "pos",
    ncase = "10943",
    ncontrol = "81899",
)

exp_dat$exposure <- "Hypothyroidism, drug reimbursement"

#exp_dat_ob <- exp_dat[exp_dat$pval.exposure < 5e-8,]
Hypo_exp_dat <- clump_data(exp_dat, clump_r2=0.01, clump_kb = 5000, pop = "EUR")

N <- 92842
r_2 <- 1/(1+(N*Hypo_exp_dat$se.exposure^2/Hypo_exp_dat$beta.exposure^2))
F_score <- r_2*(N-2)/(1-r_2)


outcome_dat <- read_outcome_data(
    filename = "finngen_R8_M13_RHEUMA_use.tsv",
    snps = Hypo_exp_dat$SNP,
    sep = "\t",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col ="alt",
    other_allele_col = "ref",
    eaf_col = "af_alt",
    pval_col = "pval",
    chr = "chrom",
    pos = "pos",
    ncase = "11178",
    ncontrol = "221323",
    phenotype = "Rheumatoid arthritis"
    
)
outcome_dat$outcome <- "Rheumatoid arthritis"

dat <- harmonise_data(exposure_dat = Hypo_exp_dat, outcome_dat = outcome_dat, action=2)

#heterogeneity
het <- mr_heterogeneity(dat)
#pleiotropy
plt <- mr_pleiotropy_test(dat)

#leave-one-out analysis
res_loo <- mr_leaveoneout(dat)

if (plt$pval < 0.05){
        print("This data has pleiotropy!")
        library(MRPRESSO)

        mr_presso_result <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000, SignifThreshold = 0.05)

        out_ind <- mr_presso_result$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
        dat_new <- dat[-out_ind,]
        dat <- dat_new

}


mr_result <- mr(dat)

#Analysis

#heterogeneity
het <- mr_heterogeneity(dat)
#pleiotropy
plt <- mr_pleiotropy_test(dat)
#leave-one-out analysis
res_loo <- mr_leaveoneout(dat)



#scatte
p1 <- mr_scatter_plot(mr_result,dat)

#forest
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)

#Leave-one-out plot
p3 <- mr_leaveoneout_plot(res_loo)

#funnel
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)

ggsave(p1[[1]], file="mr_scatter_plot.pdf", width = 8.4, height=8.4,unit='cm')
ggsave(p2[[1]], file="mr_forest_plot.pdf", width = 8.4, height=8.4,unit='cm')
ggsave(p3[[1]], file="mr_leaveoneout_plot.pdf", width = 8.4, height=8.4,unit='cm')
ggsave(p4[[1]], file="mr_funnel_plot.pdf", width = 8.4, height=8.4,unit='cm')

#Output all result
write.table(dat$SNP, "dat_SNP.txt", row.names=FALSE, col.names=FALSE)
odds_res <- generate_odds_ratios(mr_result)
write.table(odds_res, "dat_MR_results.txt", row.names=FALSE, sep='\t')
write.table(res_single, "dat_res_single.txt", row.names=FALSE, sep='\t')
write.table(res_loo,"dat_res_leaveoneout.txt",row.names=FALSE, sep='\t')
write.table(het,"dat_instrument_het.txt",row.names=FALSE, sep='\t')
write.table(plt,"dat_instrument_plt.txt",row.names=FALSE, sep='\t')



