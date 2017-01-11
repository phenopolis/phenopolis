python=/share/apps/python/bin/python

$python /SAN/vyplab/UCLex/scripts/all_hpo_terms.py > HPO.terms

while read p; do
  echo $p
	$python /SAN/vyplab/UCLex/scripts/case_control_hpo_db.py $p > ${p}.pheno
done < HPO.terms
