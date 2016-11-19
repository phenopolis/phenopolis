#$ -S /bin/sh
#$ -l h_vmem=5G
#$ -l tmem=5G
#$ -l h_rt=50:00:0
#$ -V
#$ -R y
#$ -pe smp 1
#$ -cwd
#$ -o /cluster/project8/vyp/cian/scripts/cluster/output/
#$ -e /cluster/project8/vyp/cian/scripts/cluster/error/
/cluster/project8/vyp/vincent/Software/R-3.3.0/bin/R CMD BATCH --no-save makeGoodPhenoLists.R
