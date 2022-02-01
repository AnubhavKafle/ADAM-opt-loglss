import math
import numpy as np
import pandas as pd
import gc

###################################################################

# Compute logl_ss for GWAS summary data based on ADAM optimization

####################################################################



def calc_logl (sum_data_batch, t_p, annot_ldscores):

    #Maybe apply multiprocessing here
    ldscores_category = annot_ldscores.iloc[:,2:]
    ldscores_category = np.transpose(ldscores_category)

    taus = t_p.reshape((ldscores_category.shape[0],1))

    #taus = t_p[0:-1].reshape((ldscores_category.shape[0],1))
    #snps = np.array(sum_data_batch.index.values) ----------> Should not be this cuz there are missing SNPs due to filtering
    snps_ldscores  = set(np.array(ldscores_category.columns))
    snps_sumdata   = set(np.array(sum_data_batch.index.values))
    snps = np.array(list(snps_ldscores.intersection(snps_sumdata)))

    matrix_category = np.array(ldscores_category.loc[:,snps]).reshape((ldscores_category.shape[0],snps.shape[0]))
    matrix_category[np.where(matrix_category<=0)] = 1e-6
    matrix_category = np.log(matrix_category)

    exp_log_tau_log_ld = np.array(np.e**(taus+matrix_category))

    sigma_taus_sigma_ldsc = np.sum(exp_log_tau_log_ld, axis=0)
    sum_data_batch = sum_data_batch.loc[snps,:]

    no_of_samples = np.array(sum_data_batch.n)

    expected_sj = 1 + no_of_samples * sigma_taus_sigma_ldsc
    expected_sj[np.where(expected_sj<=0)]=1e-6

    sum_stat = np.array(sum_data_batch.Stat)
    full_ld_scores = np.array(sum_data_batch.L2)

    logl =  (- sum_stat / (2*expected_sj) - 0.5*np.log(sum_stat) - 0.5*np.log(2*expected_sj) - 0.5*np.log(math.pi))/full_ld_scores
    logl = np.sum(logl)

    gc.collect()

    return logl
