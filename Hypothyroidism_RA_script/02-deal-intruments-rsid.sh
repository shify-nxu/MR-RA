# Merge intrument from three datasets
cat ./REIMB-M13/dat_SNP.txt ./UKB-190-IEU-RA/dat_SNP.txt ./UKB-77-IEU-RA/dat_SNP.txt | cut -d ' ' -f 2 | grep -v "x" |sort | uniq -c | sort -k1rn | awk '{print $2}' | sed 's/\"//g' > MR-intrument-SNPs-list

# Merge single SNP results from MR of three datasets
cat ./REIMB-M13/dat_res_single.txt ./UKB-190-IEU-RA/dat_res_single.txt ./UKB-77-IEU-RA/dat_res_single.txt | cut -f 3,4,6,7 | grep -v 'id.exposure' | grep -v 'All' | sort -k4rn | awk '{print $1"-"$2, $3,$4}' | tr ' ' '\t' > Merge-SNP-singleanalysis-result.txt


