import os, sys, glob, pdb, time

counts = "/Users/jacobfenster/Documents/NGS_gRNA/output/20231126_exact_counts/20231126_Carbon_exact_MAGeCK.txt"
NT = "/Users/jacobfenster/Documents/NGS_gRNA/databases/Lib3_NT.txt"
matrix_dir = "/Users/jacobfenster/Documents/NGS_gRNA/design_matricies"
output_dir = "/Users/jacobfenster/Documents/NGS_gRNA/output/20231126_exact_counts/MAGeCK_test"
threads = 11
log_file = f"{output_dir}/20231127_LB_tests.log"
os.system(f"> {log_file}")
#testing different parameters initially
#start with one sample LB_t2_R1 vs start LB_t0_R1
matrix = "LB_t0_t2_Rep1.txt"
#MLE algorithm on one exp vs one control with and without a control group specified
os.system("docker run -v ${HOME}:${HOME} -w /work jacobfenster/mageck:v1.0 /bin/bash -c " f"\"mageck mle -k {counts} -d {matrix_dir}/{matrix} --control-sgrna {NT} -n {output_dir}/{matrix.split('.')[0]}_mle-noday0_perm2 --permutation-round 2 --threads {threads} >> {log_file}\"")
os.system("docker run -v ${HOME}:${HOME} -w /work jacobfenster/mageck:v1.0 /bin/bash -c " f"\"mageck mle -k {counts} --day0-label LB_t0_R1 -d {matrix_dir}/{matrix} --control-sgrna {NT} -n {output_dir}/{matrix.split('.')[0]}_mle-withday0_perm2 --permutation-round 2 --threads {threads} >> {log_file}\"")
#test algoritm
os.system("docker run -v ${HOME}:${HOME} -w /work jacobfenster/mageck:v1.0 /bin/bash -c " f"\"mageck test -k {counts} -t LB_t2_R1 -c LB_t0_R1 -n {output_dir}/{matrix.split('.')[0]}_test_noday0 --control-sgrna {NT} --pdf-report >> {log_file}\"")
os.system("docker run -v ${HOME}:${HOME} -w /work jacobfenster/mageck:v1.0 /bin/bash -c " f"\"mageck test -k {counts} -t LB_t2_R1 -c LB_t0_R1 --day0-label LB_t0_R1 -n {output_dir}/{matrix.split('.')[0]}_test_withday0 --control-sgrna {NT} --pdf-report >> {log_file}\"")
#working with replicates
#start with paired test
matrix = 'LB_t2_vs_t0.txt'
#this outputs a separate value for each replicate _r0 and _r1
os.system("docker run -v ${HOME}:${HOME} -w /work jacobfenster/mageck:v1.0 /bin/bash -c " f"\"mageck test -k {counts} -t LB_t2_R1,LB_t2_R2 -c LB_t0_R1,LB_t0_R2 -n {output_dir}/{matrix.split('.')[0]}_test_2rep_noday0 --control-sgrna {NT} --pdf-report\"")
os.system("docker run -v ${HOME}:${HOME} -w /work jacobfenster/mageck:v1.0 /bin/bash -c " f"\"mageck test -k {counts} -t LB_t2_R1,LB_t2_R2 -c LB_t0_R1,LB_t0_R2 --paired -n {output_dir}/{matrix.split('.')[0]}_test_2rep_noday0_paired --control-sgrna {NT} --pdf-report\"")
#now MLE
os.system("docker run -v ${HOME}:${HOME} -w /work jacobfenster/mageck:v1.0 /bin/bash -c " f"\"mageck mle -k {counts} -d {matrix_dir}/{matrix} --control-sgrna {NT} -n {output_dir}/{matrix.split('.')[0]}_2rep_mle-noday0_perm2 --permutation-round 2 --threads {threads}\"")
os.system("docker run -v ${HOME}:${HOME} -w /work jacobfenster/mageck:v1.0 /bin/bash -c " f"\"mageck mle -k {counts} --day0-label LB_t0_R1 -d {matrix_dir}/{matrix} --control-sgrna {NT} -n {output_dir}/{matrix.split('.')[0]}_2rep_mle-withday0_perm2 --permutation-round 2 --threads {threads}\"")

if False: #going to loop once I know what setup I want to compute
    matrix_names = ['LB_t1_vs_t0.txt', 'LB_t0_t1_t2.txt', 'LB_t0_t1_Rep1.txt', 'LB_t2_vs_t1.txt', 'LB_t2_vs_t0.txt']
    paired_names = ['LB_t1_vs_t0.txt', 'LB_t0_t1_t2.txt', 'LB_t2_vs_t1.txt', 'LB_t2_vs_t0.txt']
    for matrix in matrix_names:
        breakpoint()
        #this is the command to run the mle module with control sgRNA input and 10 permutations
        #started Tue, 28 Nov 2023 02:06:21
        os.system("docker run -v ${HOME}:${HOME} -w /work jacobfenster/mageck:v1.0 /bin/bash -c " f"\"mageck mle -k {counts} -d {matrix_dir}/{matrix} --control-sgrna {NT} -n {output_dir}/{matrix.split('.')[0]}_mle --permutation-round 10 --threads {threads}\"")
    for matrix in paired_names:
        os.system(f"")
