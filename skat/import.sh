

#for hpo in /slms/gee/research/vyplab/UCLex/mainset_July2016/cian/HPO/Gene_Based_Tests/HP*
#for hpo in `head -n1 hpo.txt`
#for hpo in `tail -n1 hpo.txt`
for hpo in HP:0000556
do
    python import.py --hpo $hpo --basedir /slms/gee/research/vyplab/UCLex/mainset_July2016/cian/HPO/Gene_Based_Tests/
done
