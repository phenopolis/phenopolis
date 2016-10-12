# need to set primary key on table

for f in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
    echo chr $f
    if [ ! -f VEP_${f}.json ]
    then
        gzip -f -d VEP_${f}.json.gz > VEP_${f}.json
    fi
    mongoimport --host localhost --port 27017 --collection variants --db uclex --file VEP_${f}.json
    # -username user --password pass
done

