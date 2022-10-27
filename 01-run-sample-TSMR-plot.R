library(TwoSampleMR)
library(ggplot2)
bmi_exp <- extract_instruments(outcomes = 'finn-b-E4_HYTHY_AI_STRICT_PURCH',clump = TRUE,r2 = 0.01,kb = 5000,access_token = NULL)

outcome_dat <- extract_outcome_data(bmi_exp$SNP, outcomes = 'ieu-a-833', proxies = 1, rsq = 0.8, align_alleles = 1,  palindromes = 1, maf_threshold = 0.3)

dat <- harmonise_data(bmi_exp, outcome_dat, action = 2)

mr_result <- mr(dat)

#异质性检验
mr_heterogeneity(dat)
#水平多效性检验
mr_pleiotropy_test(dat)

library(MRPRESSO)

#devtools::install_github("rondolab/MR-PRESSO",force = TRUE)

mr_presso_result <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000, SignifThreshold = 0.05)

out_ind <- mr_presso_result$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`

dat_new <- dat[-out_ind,]
mr_result <- mr(dat_new)


#可视化
#散点图
p1 <- mr_scatter_plot(mr_result,dat_new)

#森林图
res_single <- mr_singlesnp(dat_new)
p2 <- mr_forest_plot(res_single)

#漏斗图
res_single <- mr_singlesnp(dat_new)
p4 <- mr_funnel_plot(res_single)

ggsave(p1[[1]], file="mr_scatter_plot.png", width = 7, height=7, dpi=900)
ggsave(p2[[1]], file="mr_forest_plot.png", width = 7, height=7, dpi=900)
ggsave(p4[[1]], file="mr_funnel_plot.png", width = 7, height=7, dpi=900)

write.table(dat_new$SNP,dat_new_SNP.txt")
