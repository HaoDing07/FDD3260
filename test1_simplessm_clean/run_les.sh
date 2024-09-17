#!/bin/bash
#SBATCH -A NAISS2024-1-3
#SBATCH -n 16
#SBATCH -t 120:00:00
#SBATCH -J test1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hao.ding@misu.su.se
#

mpprun ./mimicav5_* 2>&1 | tee myrun01.log
# End of script
