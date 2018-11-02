export OMP_NUM_THREADS=1
export OMP_PROC_BIND=FALSE
R --vanilla --slave -f runCluster.R --args 1 1
R --vanilla --slave -f runCluster.R --args 1 0

export OMP_NUM_THREADS=2
export OMP_PROC_BIND=FALSE
R --vanilla --slave -f runCluster.R --args 2 1
R --vanilla --slave -f runCluster.R --args 2 0

for n in `seq 4 4 24`
do
export OMP_NUM_THREADS="$n"
export OMP_PROC_BIND=FALSE
R --vanilla --slave -f runCluster.R --args "$n" 1
R --vanilla --slave -f runCluster.R --args "$n" 0
done


