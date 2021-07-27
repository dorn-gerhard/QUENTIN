# Gerhard Dorn, 21st July, 2021

We work with two directories:
[1] program directory: bm_program - contains the program
[2] job directory: run/$ID - collects all data for a certain system, which has a name $ID specified in the setup_file, e.g. run/benzene

Changing and compiling the program is done within the bm_program folder

Setup of the system and calling of the program with various parameters is done from the run/$ID directory


ad 1)
a) run "compile.m" after changes in the program (e.g. enabling CPT), check that paths are correct (e.g. GROUNDPATH = "~/QUENTIN/bm_program/")
b) the bash file "start_matlab.sh" calls the automatically generated "run_exec_bm.sh" file which starts the compiled program "exec_bm" with the parameters in the correct order
-> check that the path in "start_matlab.sh" links to "bm_program/run_exec_bm.sh"


ad 2)
c) Choose a short name for the system you want to study (e.g. TEST) and create a job directory for the system you want to study in /run/ (e.g. /run/TEST) (copy the TEST directory and change all TEST names to your name)
d) Setup for the open quantum system is done within the setup_file_$ID
running it will create the file "setup_file_$JOB_ID" in the ./data/ subfolder. The program takes all  of its input parameters except the job_id from this file. The created $JOB_ID contains the day of its creation and an index. 

e) You have to write the $JOB_ID in the batch_config file (e.g. Conder_sub.q) and define how many runs shall be started

f) call your batch job (using e.g. HTCondor or SLURM) from the job directory [2] and get results saved in ./data/ and ./log_files.

g) use "post_fun_fuse_files($RUN_ID)" (in bm_program/functions/) to create a combined_data_file that contains all return paramters from all runs of a job for further data evaluation.



