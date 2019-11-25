import math
import numpy as np
import pandas as pd
import gc


def neg_gradient (sum_data_batch, t_p, annot_ldscores):

    #Maybe apply multiprocessing here

    ldscores_category = annot_ldscores.iloc[:,2:]

    
  
    ldscores_category = np.transpose(ldscores_category)

 
    taus = t_p.reshape((ldscores_category.shape[0],1))


    #snps = np.array(sum_data_batch.index.values) ----------> Should not be this cuz there are missing SNPs due to filtering
    snps_ldscores  = set(np.char.strip(np.array(ldscores_category.columns, dtype='str')))
    snps_sumdata   = set(np.char.strip(np.array(sum_data_batch.index.values,dtype='str')))

 

    snps = np.array(list(snps_ldscores.intersection(snps_sumdata)))

    print(snps)

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
    first_part  =  (0.5 * sum_stat / expected_sj**2 - 0.5  / expected_sj) / full_ld_scores

    first_part = first_part * no_of_samples

    first_part = first_part.reshape((1, snps.shape[0]))

    exp_tau_ldscores_category = np.multiply(np.transpose(np.e**taus), np.transpose(ldscores_category.loc[:,snps]))

    all_tau_deri = np.dot(first_part, exp_tau_ldscores_category)

    logl =  (- sum_stat / (2*expected_sj) - 0.5*np.log(sum_stat) - 0.5*np.log(2*expected_sj) - 0.5*np.log(math.pi))/full_ld_scores
    logl = np.sum(logl)

    gradient_vec = -1 * all_tau_deri

    gc.collect()

    return (gradient_vec, logl)

def stochastic_opt (sum_data_batch, alpha_learn, tol, dim_annot, no_iteration, annot_prefix):
    dim_annot = dim_annot + 1

    theta_p_curr = np.linspace(-15,1,dim_annot, dtype='float64') #np.random.normal(0,1,dim_annot)

    epoch=1

# ADAM part
#initialize the values of the parameters
    beta_1 = 0.9
    beta_2 = 0.999
    m_t = np.zeros(dim_annot)
    v_t = np.zeros(dim_annot)
    epoch=0
    loglikeli=list()

    #Start iteration here
    while epoch <= no_iteration :
        #print("current epoch {0}".format(epoch))
        epoch+=1
        theta_p_previous = theta_p_curr

        (grad, logl) = neg_gradient(sum_data_batch,theta_p_curr, annot_prefix)

        loglikeli.append(logl)

        grad = grad #/ sum_data_batch.shape[0]
        print("Gradient is",grad)                                   # should it be divided by batch size?? check it out
        m_t = beta_1 * m_t + (1-beta_1)*grad                        #updates the moving averages of the gradient
        v_t = beta_2 * v_t + (1-beta_2)*(np.power(grad,2))          #updates the moving averages of the squared gradient

        #calculate the bias-corrected estimates
        m_corrected = m_t/(1-(beta_1**epoch))   + (1 - beta_1) * grad / (1 - np.power(beta_1, epoch))
        v_corrected = v_t/(1-(beta_2**epoch))

        test=(alpha_learn*m_corrected)/(np.sqrt(v_corrected)+tol)

        theta_p_curr = theta_p_curr - (alpha_learn*m_corrected)/(np.sqrt(v_corrected)+tol)
       

        print("loglikelihood is",logl)
        #check for convergence -- better way of doing this?

        #if sum(abs(theta_p_curr - theta_p_previous)) == 0.0 : #checks if it is converged or not
            #break

        if epoch > 3000:
            tols=0.5e-2
            diff_likelihood = np.abs(np.max(loglikeli[::-1][0:10]) - np.min(loglikeli[::-1][0:10]))
            #left = -1*np.log(1+tols)
            #right = 1*np.log(1+tols)
            #log_ratio = np.log(np.abs(np.max(loglikeli[::-1][0:10]))) - np.log(np.abs(np.min(loglikeli[::-1][0:10])))
            #if left < log_ratio and log_ratio < right :
            if diff_likelihood < tols:
                print("convergence achieved at the tolerance level")
                print(" Achieved at epoch number {0}".format(epoch))
                break

    if epoch >= no_iteration:
        print("Could not be optimized inside the iteration. Current value is saved ..")
        print(" Epoch number {0}".format(epoch))
        return theta_p_curr,loglikeli
    #else:
    return theta_p_curr, loglikeli

# In[4]:
