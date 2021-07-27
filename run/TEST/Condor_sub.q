
Universe   = vanilla
#Universe   = scheduler


# =======================Change here========================

# change for each run:

IDENTIFIER=TEST
SETUP_ID=191120_1
#Queues the job x times with process id ranging from 0 to x-1
Queue 101

# change once for your setup:
Executable = /QUENTIN/bm_program/start_matlab.sh 
# path to your Matlab installation:
PATH=/opt/matlab/R2018b/

notify_user = @gmail.com
#===========================================================

Output     = ./log_files/Console_$(SETUP_ID)_$(Process).out
Error      = ./log_files/Error_$(SETUP_ID)_$(Process).err
Log        = ./log_files/Condor_$(SETUP_ID)_$(Process).log
run_as_owner = true


notification = Error
nice_user = false

#requirements = (TARGET.Machine =!= "129.27.161.28") && (TARGET.Machine =!= "129.27.161.29") && (TARGET.Machine =!= "129.27.161.64") && (TARGET.Machine =!= "129.27.161.132")
#requirements = (TARGET.Machine =!= "faepop68.tu-graz.ac.at") && (TARGET.Machine =!= "faepop23.tu-graz.ac.at") && (TARGET.Machine =!= "faepop30.tu-graz.ac.at") && (TARGET.Machine =!= "faepop69.tu-graz.ac.at") && (TARGET.Machine =!= "faepop71.tu-graz.ac.at")
requirements = (TARGET.Machine =!= "faepop04.tu-graz.ac.at") && (TARGET.Machine =!= "faepop03.tu-graz.ac.at") && (TARGET.Machine =!= "faepop02.tu-graz.ac.at") && (TARGET.Machine =!= "faepop01.tu-graz.ac.at") && (TARGET.Machine =!= "faepop67.tu-graz.ac.at")

request_memory = max({12, Target.TotalSlotMemory})
rank = Memory
request_cpus = 4



SETUP_FILE=./data/setup_$(IDENTIFIER)_$(SETUP_ID).mat


#1st ARG: Process number; 2nd ARG: number of input arguments line in CMD_ARG file; 3rd ARG: Outputfolder;
#if 4th (optional): if neccessary we can also give the path to the Parameter file
#intialdir = this would set an intial directory
Arguments  = $(Process) $(SETUP_FILE) $(PATH)


