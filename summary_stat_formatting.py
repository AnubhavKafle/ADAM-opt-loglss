#To clean and format summary data, with already columns having L2
import numpy as np
import pandas as pd



def summ_data_format(sum_data_file, hapmap_loc):

    hapmap_snps_file=hapmap_loc
    hapmap_snps = np.loadtxt(hapmap_snps_file, dtype="str")
    
    sum_data = pd.read_csv(sum_data_file,sep="\t", header=0)
    sum_data.set_index('SNP', inplace=True)

    SNPs_summary_data = set(sum_data.index)
    hapmap_snps = set(hapmap_snps)

    common_hapmap_summary_data_SNP = hapmap_snps.intersection(SNPs_summary_data)

    sum_data_hapmap = sum_data.loc[common_hapmap_summary_data_SNP,:]


    remove_outliers = sum_data_hapmap.Stat <= sum_data_hapmap.n/99
    sum_data_clean = sum_data_hapmap.loc[np.array(remove_outliers, dtype=bool),:]

    print("Summary data formatting completed")

    return sum_data_clean

def format_annotation(annot_file):
    annot_categor = pd.read_csv(annot_file, sep="\t",header=0,index_col=[1])
    print(annot_categor.columns)
    #annot_categor = annot_categor['SNP'].astype(str)
    #annot_categor.set_index('SNP', inplace=True)
    print("Annotation formatting completed")
    return annot_categor
