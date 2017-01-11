#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -cwd
#$ -V
#$ -l tmem=6G,h_vmem=6G
#$ -l h_rt=240:0:0
#$ -t 1-25
set -u
set -x
scriptname=macular_dystrophy
mkdir -p ${scriptname}.qsub.out ${scriptname}.qsub.err
exec >${scriptname}.qsub.out/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.out 2>${scriptname}.qsub.err/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.err
args=( header 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
chrom=${args[$SGE_TASK_ID]}

Rscript="/cluster/project8/vyp/vincent/Software/R-3.3.0/bin/Rscript"
SKAT="/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/Gene_Based_Tests/SKAT/SKAT_function.R"
oDir="/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Macular_Dystrophy/macular_dystrophy"
genes="/SAN/vyplab/UCLex/mainset_July2016/cian/gene.names"
control_list="/SAN/vyplab/UCLex/support/IoO/Opthal.control.list"
case_list="/SAN/vyplab/UCLex/support/IoO/Macular_Dystrophy/macular_dystrophy_cases"
$Rscript $SKAT --case.list $case_list --TargetGenes $genes --control.list $control_list \
 --oDir ${oDir}_chr_${chrom} --SavePrep TRUE --chrom $chrom --compoundHets TRUE --qcPREP TRUE 
