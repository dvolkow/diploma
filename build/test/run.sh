#! /bin/bash
for (( i=1; i<=$1; i++))
do
        mkdir ${i}
        ../diploma --file ../apogee-ucac.txt --solution f --low-bound 24 --err 25 --upper-bound 26 --profile --mksize 10 --ord ${i}
        ../../sources/scripts/plottest.py r .
        ../../sources/scripts/plottest.py th ./R0Theta0.txt
        ../../sources/scripts/plottest.py mk ./mk_results.txt
        ../../sources/scripts/plottest.py xy .
        gnuplot -c ../../sources/scripts/main_cmd.gnu "cur.eps" 'objs.txt' 'ERROR_LIMITED_R_THETA'
        mv *.txt ${i}
        mv *.png ${i}
done 

