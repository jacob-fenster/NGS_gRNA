import os, sys, glob, pdb, time

counts = "/Users/jacobfenster/Documents/NGS_gRNA/output/20231126_exact_counts/20231126_Carbon_exact_MAGeCK.txt"
NT = "/Users/jacobfenster/Documents/NGS_gRNA/databases/Lib3_NT.txt"
matrix_dir = "/Users/jacobfenster/Documents/NGS_gRNA/design_matricies"
output_dir = "/Users/jacobfenster/Documents/NGS_gRNA/output/20231126_exact_counts/20231206_MAGeCK_test"
threads = 11
log_file = f"{output_dir}/20231206_test.log"
os.system(f"> {log_file}")
#testing different parameters initially
#start with one sample LB_t2_R1 vs start LB_t0_R1
matrix = "Large_design_matrix_test20231206.txt"
#MLE algorithm on one exp vs one control with and without a control group specified
time_i=time.time()
os.system("docker run -v ${HOME}:${HOME} -w /work jacobfenster/mageck:v1.0 /bin/bash -c " f"\"mageck mle -k {counts} -d {matrix_dir}/{matrix} --control-sgrna {NT} -n {output_dir}/{matrix.split('.')[0]}_mle_perm10 --permutation-round 10 --threads {threads} >> {log_file}\"")
time_f=time.time()
os.system(f'\"echo this script took {time_f-time_i} seconds with {threads} threads\" >> design_matricies/Large_design_matrix_test20231206.txt')
