import sys
import pandas as pd
import scipy.stats as stats

# Load MR singleSNP b score
single_SNP_MR_result = pd.read_table("Merge-SNP-singleanalysis-result.txt", header=None, sep='\t')
single_SNP_MR_result.columns = ['exposure_result','SNP_id','b_score']

# Merge SNP b score, for the same rsid has mean b_score
b_score_mean = single_SNP_MR_result.groupby(['SNP_id'])['b_score'].mean()
# Merge SNP exposure and result ID together
ID_exposure = single_SNP_MR_result.groupby('SNP_id')['exposure_result'].apply(lambda x:x.str.cat(sep=';')).reset_index()
ID_exposure.index = ID_exposure.SNP_id
ID_exposure['b_score'] = b_score_mean

# Load Sei score
sei_result = pd.read_table("/home/sdc1/Shify/Cooperation_project/MR/01-TwoSanmpleMR-ID-data/result/Sei-predict-variants/predict2/9f39337b-6788-4522-a77b-a3d1328be632_MR_intrument_SNPs_list.reformat.GRCh37_sequence-class-scores.tsv",sep='\t')

sei_result_extract = pd.concat([sei_result.iloc[:,5], sei_result.iloc[:,9:]], axis=1)
sei_result_extract['id_max'] = sei_result_extract.iloc[:,1:].abs().idxmax(axis=1)

def Get_max_id_value(each_line):
    sei_value = each_line[each_line["id_max"]]
    snp = each_line['id']
    return sei_value

sei_result_extract_maxscore = sei_result_extract.apply(Get_max_id_value,axis=1)
sei_result_extract['value_max'] = sei_result_extract_maxscore
sei_result_use = sei_result_extract.loc[:,['id','value_max']]
sei_result_use_merge = sei_result_use.groupby('id')['value_max'].mean()

ID_exposure['sei_score'] = sei_result_use_merge[ID_exposure['SNP_id']]
ID_exposure_use = ID_exposure.dropna(axis=0,how='any')

ID_exposure_use.to_csv("Merge-SNP-singleanalysis-result-Sei-score.txt",sep='\t')

# Make file with pick variants
sei_pick_variants = pd.read_table("/home/sdc1/Shify/Cooperation_project/MR/01-TwoSanmpleMR-ID-data/result/Sei-predict-variants/predict2/pick_snp.txt",sep='\t')
sei_pick_variants_result_extract = pd.concat([sei_pick_variants.iloc[:,2], sei_pick_variants.iloc[:,7:]], axis=1)
# Merge the differnces of same variant
sei_pick_variants_result_extract_merge = sei_pick_variants_result_extract.groupby('id').mean()

# Add the b score of each variant to the merge result
sei_pick_variants_result_extract_merge['b_score'] = ID_exposure_use.loc[sei_pick_variants_result_extract_merge.index,'b_score']

sei_pick_variants_result_extract_merge.to_csv("Merge-SNP-singleanalysis-Sei_pick_vairiants.txt",sep='\t')
