To run a job on the cluster with <n> cpus execute the following commands:

export OMP_NUM_THREADS=<n>
export OMP_PROC_BIND=FALSE
bsub -n <n> R --vanilla --slave -f runCluster.R --args <n>

Wait for the job to start before changing the variables OMP_NUM_THREADS and OMP_PROC_BIND.

