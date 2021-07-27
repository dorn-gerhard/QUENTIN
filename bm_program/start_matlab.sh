#!/bin/bash
n_run=$1;
setup_file=$2;
path=$3;

#Omega=$(cat CMD_line_ARG.txt | head -n $((n_run+2)) | tail -n 1);
echo "+++++++++++SCRIPT OUTPUT: Started run for $n_run +++++++++++++++";
#nohup ~/bm_program/run_exec_bm.sh $path $setup_file $n_run
nohup ~/QUENTIN/bm_program/run_exec_bm.sh $path $setup_file $n_run
echo "+++++++++++SCRIPT OUTPUT: Ended run for $n_run +++++++++++++++++";
