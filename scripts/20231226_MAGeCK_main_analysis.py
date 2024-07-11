import os, sys, glob, pdb, time

counts = "/Users/jacobfenster/Documents/NGS_gRNA/output/20231126_exact_counts/20231126_Carbon_exact_MAGeCK.txt"
NT = "/Users/jacobfenster/Documents/NGS_gRNA/databases/Lib3_NT.txt"
output_dir = '/Users/jacobfenster/Documents/NGS_gRNA/output/20231226_main_analysis'
threads = 11
log_file = f"{output_dir}/20231226_main_analysis_script.log"
os.system(f"> {log_file}")
# LB_t0_R1	LB_t0_R2	LB_t1_R1	LB_t1_R2	LB_t2_R1	LB_t2_R2	gluc_t1_R1	gluc_t1_R2	gluc_t2_R1	gluc_t2_R2	ace_t1_R1	ace_t1_R2	ace_t2_R1	ace_t2_R2	10pCA_t1_R1	10pCA_t1_R2	10pCA_t2_R1	10pCA_t2_R2	50pCA_t1_R1	50pCA_t1_R2	50pCA_t2_R1	50pCA_t2_R2
# each tuple has the format exp[0-1] are -t (test) and exp[2-3] are -c (control) [end - initial]
experiments = [
               ('LB_t2_R1', 'LB_t2_R2', 'LB_t1_R1', 'LB_t1_R2'),
               ('gluc_t2_R1', 'gluc_t2_R2', 'gluc_t1_R1', 'gluc_t1_R2'),
               ('ace_t2_R1', 'ace_t2_R2', 'ace_t1_R1', 'ace_t1_R2'),
               ('10pCA_t2_R1', '10pCA_t2_R2', '10pCA_t1_R1', '10pCA_t1_R2'),
               ('50pCA_t2_R1', '50pCA_t2_R2', '50pCA_t1_R1', '50pCA_t1_R2'),
               ('LB_t1_R1', 'LB_t1_R2', 'LB_t0_R1', 'LB_t0_R2'),
               ('gluc_t1_R1', 'gluc_t1_R2', 'LB_t0_R1', 'LB_t0_R2'),
               ('ace_t1_R1', 'ace_t1_R2', 'LB_t0_R1', 'LB_t0_R2'),
               ('10pCA_t1_R1', '10pCA_t1_R2', 'LB_t0_R1', 'LB_t0_R2'),
               ('50pCA_t1_R1', '50pCA_t1_R2', 'LB_t0_R1', 'LB_t0_R2')
               ]
# there is some error where it is not permutating the NT gRNAs. troubleshoot tomorrow. 
for exp in experiments:
    timei = time.time()
    os.system(f"echo \"experiment {exp[0].split('_R')[0]}-vs-{exp[2].split('_R')[0]} has begun\" >> {log_file}")
    with open(log_file, 'a') as f:
        f.write("The command is: docker run -v ${HOME}:${HOME} -w /work jacobfenster/mageck:v1.0 /bin/bash -c " f"\"mageck test -k {counts} -t {exp[0]},{exp[1]} -c {exp[2]},{exp[3]} -n {output_dir}/{exp[0].split('_R')[0]}-vs-{exp[2].split('_R')[0]} --control-sgrna {NT} --pdf-report\"")
    os.system("docker run -v ${HOME}:${HOME} -w /work jacobfenster/mageck:v1.0 /bin/bash -c " f"\"mageck test -k {counts} -t {exp[0]},{exp[1]} -c {exp[2]},{exp[3]} -n {output_dir}/{exp[0].split('_R')[0]}-vs-{exp[2].split('_R')[0]} --control-sgrna {NT} --pdf-report\"")
    timef = time.time()
    os.system(f"echo \"experiment {exp[0].split('_R')[0]}-vs-{exp[2].split('_R')[0]} has finished in {(timef-timei)/60:.1f} minutes\" >> {log_file}")
    os.system(f"echo \"\" >> {log_file}")
os.system(f"echo \"the script has completed\" >> {log_file}")