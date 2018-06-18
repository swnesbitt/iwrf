#!/bin/bash

#SBATCH --job-name=plot
#SBATCH --array=0-48
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --time=1:00:00
#SBATCH --partition=g
#SBATCH --mem-per-cpu=10000
export LD_LIBRARY_PATH=

#sleep $[ ( $RANDOM % 10 )  + 1 ]s

run=${1}

cd /data/keeling/a/snesbitt/wrf_realtime/scripts/

#./get_smn.sh ${run}

cd /data/keeling/a/snesbitt/wrf_realtime/py-postplot

export fhr=`printf "%02d\n" ${SLURM_ARRAY_TASK_ID}`

python wrf_plot_FULL_domain.py $run $fhr full
python wrf_plot_FULL_domain.py $run $fhr cordoba_zoom
python wrf_plot_FULL_domain.py $run $fhr mendoza_zoom


#python iwrf_postproc.py $run $fhr cba "Servicio Meteorologico Nacional GFS-WRF" "SMN_WRF"
#python iwrf_postproc.py $run $fhr full
#python iwrf_postproc.py $run $fhr mdz

