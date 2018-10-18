#! /bin/bash

for (( i=1; i<=$1; i++))
do
        ./run.sh $i $2 $3 $4
done 
mv result_* ../results/
./clean.sh
