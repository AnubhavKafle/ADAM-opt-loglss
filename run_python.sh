#!/bin/bash


SCRIPT="/home/student.unimelb.edu.au/akaphle/Desktop/work/Research/eQTLs_heritability/Neals_lab_UKBB/new_likelihood/Profile_Likelihood_based_optimization/scripts/new/"
SUMDATA="/home/student.unimelb.edu.au/akaphle/Desktop/work/Research/eQTLs_heritability/Neals_lab_UKBB/new_likelihood/Profile_Likelihood_based_optimization/scripts/UKBB_BMI_with_L2.txt"
ANNOT="/home/student.unimelb.edu.au/akaphle/Desktop/work/Research/eQTLs_heritability/Neals_lab_UKBB/new_likelihood/Profile_Likelihood_based_optimization/results/6_LD_annotations_with_base/All_chr_ldscores_0.00_header.txt"

#ANNOT="/home/student.unimelb.edu.au/akaphle/Desktop/work/Research/eQTLs_heritability/Neals_lab_UKBB/new_likelihood/Profile_Likelihood_based_optimization/results/24_funct_ldscores_maf_6_LD_annotations/All_chr_ldscores_2.00_header.txt"
MAF=1.00
OUTFILENAME="/home/student.unimelb.edu.au/akaphle/Desktop/work/Research/eQTLs_heritability/Neals_lab_UKBB/new_likelihood/Profile_Likelihood_based_optimization/scripts/"
HAPMAP="/home/student.unimelb.edu.au/akaphle/Desktop/work/Research/eQTLs_heritability/Neals_lab_UKBB/new_likelihood/data/hapmap3_snps/hm.all_chr.snp"
DIMENSION=6

python ${SCRIPT}/logL_main.py --sum_data ${SUMDATA} --annot_file ${ANNOT} --alpha 0.05 --no_annot ${DIMENSION} --tol 1e-8 --no_iteration 6500 --which_maf ${MAF}_ldsc --hapmap ${HAPMAP} --out ${OUTFILENAME}
