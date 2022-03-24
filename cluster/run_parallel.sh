#!/bin/bash
#Submission options:
#$-S /bin/bash 
#$-N seeds 
#$-cwd 
#$-t 1-120
#$-cwd 

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/lib/
#/share/apps/bin/python2.7 2017_07_21_fit_v1.py $SGE_TASK_ID $niters $outdir 
/share/apps/bin/python2.7 2018_10_24_rerun_2018_05_31_DpbyGi_multseed_v16.py $SGE_TASK_ID 
