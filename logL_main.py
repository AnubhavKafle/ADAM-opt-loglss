import numpy as np
import pandas as pd
import argparse
import os

#importing custom modules here
from optimization import *
from summary_stat_formatting import *
from likelihood_estimation import *

def parse_args():
    parser = argparse.ArgumentParser(description='Optimization of the new loglikelihood')

    parser.add_argument('--annot_file', default=None, type=str,
        help='Provide annotation file ')

    parser.add_argument('--sum_data', default=None, type=str,
        help='Summary data file')

    parser.add_argument('--alpha', default=None, type=float,
        help='learning rate')

    parser.add_argument('--no_annot', default=61, type=str,
        help='no of annotations')

    parser.add_argument('--tol', default=None, type=float,
        help='tolerance')
    parser.add_argument('--no_iteration', default=None, type=int,
        help='no_of_iterations')
    parser.add_argument('--which_maf', default=None, type=str,
        help='what maf threshold used for LDscores scaling')
    parser.add_argument('--out', default=None, type=str,
        help='input output directory')
    parser.add_argument('--hapmap', default=None, type=str,
        help='input hapmap snps file')

    opts = parser.parse_args()
    return opts


###########################
#      MPI parallelizing  #
###########################
#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#ncore = comm.Get_size()

if __name__ == '__main__':
    args = parse_args()
    sum_data_file = args.sum_data
    no_annot = int(args.no_annot)
    annot_file = args.annot_file
    alpha=args.alpha
    tol = args.tol
    no_iteration = args.no_iteration
    which_maf   = args.which_maf
    hapmap = args.hapmap

    summ_data_clean = summ_data_format(sum_data_file,hapmap)
    annot = format_annotation(annot_file)

    #perform optimizaiton here

    tau, loglikeli = stochastic_opt(summ_data_clean, float(alpha), float(tol), no_annot, int(no_iteration), annot )

    savefolder=args.out

    outfile_tau = savefolder + "best_tau_maf_{0}".format(which_maf)
    outfile_logl = savefolder + "loglikeli_maf_{0}".format(which_maf)
    np.save(outfile_tau,tau)

    np.save(outfile_logl,loglikeli)

    print("The program has ended")
