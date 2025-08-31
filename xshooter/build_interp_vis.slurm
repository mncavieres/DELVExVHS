
#!/bin/bash
#SBATCH --job=vacuum_rvspecfit
#SBATCH --output=%x_%j.out
#SBATCH --mail-user="cavierescarrera@strw.leiden.nl"
#SBATCH --mail-type="ALL"
#SBATCH --mem=16000M               # 64 GB of memory for all tasks
#SBATCH --time=100:00:00           # Maximum runtime
#SBATCH --account=gpu_strw
#SBATCH --partition=gpu_strw       # Specify the partition
#SBATCH --ntasks=1                 # 24 MPI tasks
#SBATCH --cpus-per-task=20          # 1 CPU per task (for pure MPI jobs)
#SBATCH --gpus=1

# Load modules

# Activate the Python virtual environment
source /home/cavierescarreramc/data1/gemini_3/environment/rvspecfit_env_torch/bin/activate

# Job information
echo "[$SHELL] #### Starting Python test"
echo "[$SHELL] ## This is $SLURM_JOB_USER on $HOSTNAME and this job has the ID $SLURM_JOB_ID"
export CWD=$(pwd)
echo "[$SHELL] ## Current working directory: $CWD"


PREFIX=/home/cavierescarreramc/data1/XSHOOTER/
TEMPLPREF=/home/cavierescarreramc/data1/PHOENIX/phoenix.astro.physik.uni-goettingen.de/data/v2.0/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/
WAVEFILE=/home/cavierescarreramc/data1/PHOENIX/phoenix.astro.physik.uni-goettingen.de/data/v2.0/HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits
RLAM0=5300
RLAM1=10300
VSINIS=0.01,400
RESOL=11333
STEP=0.02
SETUP=VIS_XSHOOTER

rvs_read_grid --prefix $TEMPLPREF --templdb files.db

# changed this to vacuum just to check if this makes my measurement match Kreuzer and brown. 

rvs_make_interpol --setup $SETUP --lambda0 $RLAM0 --lambda1 $RLAM1 --resol $RESOL --step $STEP --templdb ${PREFIX}/files.db --oprefix ${PREFIX}/ --templprefix $TEMPLPREF --wavefile $WAVEFILE --nthreads 19

rvs_make_nd --prefix ${PREFIX}/ --setup $SETUP # 

rvs_train_nn_interpolator  --dir ${PREFIX}/ --setup $SETUP

rvs_make_ccf --setup $SETUP --lambda0 $RLAM0 --lambda1 $RLAM1 --every 30 --vsinis $VSINIS --prefix ${PREFIX}/ --oprefix=${PREFIX} --step $STEP


echo "[$SHELL] ## Script finished"

echo "[$SHELL] #### Finished Python test. Have a nice day"
