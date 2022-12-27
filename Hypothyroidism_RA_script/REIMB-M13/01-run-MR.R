library(TwoSampleMR)
library(ggplot2)
exp_dat <- extract_instruments(outcomes = 'finn-b-HYPOTHY_REIMB',clump = TRUE,r2 = 0.01,kb = 5000,access_token = NULL)

outcome_dat <- extract_outcome_data(exp_dat$SNP, outcomes = 'finn-b-M13_RHEUMA', proxies = 1, rsq = 0.8, align_alleles = 1,  palindromes = 1, maf_threshold = 0.3)

dat <- harmonise_data(exp_dat, outcome_dat, action = 2)
dat$exposure <- "Hypothyroidism, drug reimbursement"
dat$outcome <- "Rheumatoid arthritis"


#异质性检验
het <- mr_heterogeneity(dat)
#水平多效性检验
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


#可视化
#散点图
p1 <- mr_scatter_plot(mr_result,dat)

#森林图
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)

#Leave-one-out plot
p3 <- mr_leaveoneout_plot(res_loo)

#漏斗图
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)

#ggsave(p1[[1]], file="mr_scatter_plot.png", width = 7, height=7, dpi=900)
#ggsave(p2[[1]], file="mr_forest_plot.png", width = 7, height=7, dpi=900)
#ggsave(p3[[1]], file="mr_leaveoneout_plot.png", width = 7, height=7, dpi=900)
#ggsave(p4[[1]], file="mr_funnel_plot.png", width = 7, height=7, dpi=900)


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


