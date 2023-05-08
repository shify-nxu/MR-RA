library(TwoSampleMR)
library(ggplot2)
library(MRPRESSO)

exp_dat <- extract_instruments(outcomes = 'ukb-a-77',clump = TRUE,r2 = 0.01,kb = 5000,access_token = NULL)

#Out_RH <-fread('RA_GWASmeta_European_v2_usetorepare.tsv',header=T,sep='\t')
#colnames(Out_RH) <- c('SNPID','Chr','Position','A1','A2','OR','OR_low','OR_up','Pval','Numbers')
#se <- abs(log(Out_RH$OR)/qnorm(Out_RH$Pval/2))
#beta <- log(Out_RH$OR)
#write.table(Out_RH,"RA_GWASmeta_European_v2_usetorepare_bs.tsv",sep='\t',row.names=FALSE,quote=FALSE)


outcome_dat <- read_outcome_data(
    filename = "../RA_GWASmeta_European_v2_usetorepare_bs.tsv",
    snps = exp_dat$SNP,
    sep = "\t",
    snp_col = "SNPID",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col ="A1",
    other_allele_col = "A2",
    pval_col = "Pval",
    chr = "Chr",
    pos = "Position",
    ncase = "14361",
    ncontrol = "42923",
    phenotype = "Rheumatoid arthritis"

)
outcome_dat$outcome <- "Rheumatoid arthritis"

correction_factor <- abs(mean(outcome_dat$beta.outcome)/mean(exp_dat$beta.exposure))
exp_dat$beta.exposure <- exp_dat$beta.exposure * correction_factor

dat <- harmonise_data(exp_dat, outcome_dat, action = 2)

dat$exposure <- "Non-cancer illness code self-reported: hypothyroidism/myxoedema"
dat$outcome <- "Rheumatoid arthritis"


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

#heterogeneity
het <- mr_heterogeneity(dat)
#pleiotropy
plt <- mr_pleiotropy_test(dat)
#leave-one-out analysis
res_loo <- mr_leaveoneout(dat)



#scatter
p1 <- mr_scatter_plot(mr_result,dat)

#forest
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)

#Leave-one-out plot
p3 <- mr_leaveoneout_plot(res_loo)

#funnel
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)

ggsave(p1[[1]], file="mr_scatter_plot.png", width = 7, height=7, dpi=900)
ggsave(p2[[1]], file="mr_forest_plot.png", width = 7, height=7, dpi=900)
ggsave(p3[[1]], file="mr_leaveoneout_plot.png", width = 7, height=7, dpi=900)
ggsave(p4[[1]], file="mr_funnel_plot.png", width = 7, height=7, dpi=900)

#Output all result
write.table(dat$SNP, "dat_SNP.txt", row.names=FALSE, col.names=FALSE)
odds_res <- generate_odds_ratios(mr_result)
write.table(odds_res, "dat_MR_results.txt", row.names=FALSE, sep='\t')
write.table(res_single, "dat_res_single.txt", row.names=FALSE, sep='\t')
write.table(res_loo,"dat_res_leaveoneout.txt",row.names=FALSE, sep='\t')
write.table(het,"dat_instrument_het.txt",row.names=FALSE, sep='\t')
write.table(plt,"dat_instrument_plt.txt",row.names=FALSE, sep='\t')


