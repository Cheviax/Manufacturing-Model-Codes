#!/bin/sh
#SBATCH --job-name=ml_job
#SBATCH --partition=midmem                    # submission queue
#SBATCH --time=7-0:00:00                        # 1-1 means one day and one hour -> here 20min
#SBATCH --mail-type=END
##SBATCH --mail-user=alain.ngandjong@u-picardie.fr       #e-mail notification
#SBATCH --output=job_ml-%j.out         # if --error is absent, includes alsothe errors
#SBATCH --nodes=1      # 2 cpus
#SBATCH --ntasks-per-node=40
##SBATCH --cpus=10
#SBATCH --mem=5G                                # T-tera, G-giga, M-mega
#SBATCH --reservation=midmem
#SBATCH -w midnode001

##==================================================================================
# this job run on 2 nodes with 28 cpus on each
##==================================================================================

echo "------------------------------------------------------------------------------"
echo "hostname                     =   $(hostname)"
echo "SLURM_JOB_NAME               =   $SLURM_JOB_NAME"
echo "SLURM_SUBMIT_DIR             =   $SLURM_SUBMIT_DIR"
echo "SLURM_JOBID                  =   $SLURM_JOBID"
echo "SLURM_JOB_ID                 =   $SLURM_JOB_ID"
echo "SLURM_NODELIST               =   $SLURM_NODELIST"
echo "SLURM_JOB_NODELIST           =   $SLURM_JOB_NODELIST"
echo "SLURM_TASKS_PER_NODE         =   $SLURM_TASKS_PER_NODE"
echo "SLURM_JOB_CPUS_PER_NODE      =   $SLURM_JOB_CPUS_PER_NODE"
echo "SLURM_TOPOLOGY_ADDR_PATTERN  =   $SLURM_TOPOLOGY_ADDR_PATTERN"
echo "SLURM_TOPOLOGY_ADDR          =   $SLURM_TOPOLOGY_ADDR"
echo "SLURM_CPUS_ON_NODE           =   $SLURM_CPUS_ON_NODE"
echo "SLURM_NNODES                 =   $SLURM_NNODES"
echo "SLURM_JOB_NUM_NODES          =   $SLURM_JOB_NUM_NODES"
echo "SLURMD_NODENAME              =   $SLURMD_NODENAME"
echo "SLURM_NTASKS                 =   $SLURM_NTASKS"
echo "SLURM_NPROCS                 =   $SLURM_NPROCS"
echo "SLURM_MEM_PER_NODE           =   $SLURM_MEM_PER_NODE"
echo "SLURM_PRIO_PROCESS           =   $SLURM_PRIO_PROCESS"
echo "------------------------------------------------------------------------------"

#python3.6 test_viscosity_Keras.py
python3.6 test_viscosity.py
touch mlfinish
