oDir=CaseControlScripts/
mkdir $oDir

/cluster/project8/vyp/vincent/Software/R-3.3.0/bin/R CMD BATCH prep.R
first=/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/Gene_Based_Tests/SKAT/HPO/first.txt
second=/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/Gene_Based_Tests/SKAT/HPO/second.txt

phenoDir=/SAN/vyplab/UCLex/support/HPO/Phenos/

while read p; do
  	case_list=${phenoDir}${p}_cases
  	ctrl_list=${phenoDir}${p}_ctrls
	cleanP=$(echo "$p" | tr : _)
	oScript=${oDir}${cleanP}.sh

  	SKATout=/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Gene_Based_Tests/${p}

  	cat $first > $oScript
  	echo 'case_list='"$case_list" >> $oScript
  	echo 'control_list='"$ctrl_list" >> $oScript
  	echo 'oDir='"$SKATout" >> $oScript
  	cat $second >> $oScript
  	qsub $oScript
done < GoodSets