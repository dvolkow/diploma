#! /bin/bash
for (( i=$1; i<=$2; i++))
do
        mkdir ${i}
        ../diploma --file ../apogee-ucac.txt --solution f --low-bound 10 --err 25 --upper-bound 40 --profile --bolter 1 --mksize $3 --ord ${i}
        ../../sources/scripts/plottest.py r .
        ../../sources/scripts/plottest.py th ./R0Theta0.txt
        ../../sources/scripts/plottest.py mk ./mk_results.txt
        ../../sources/scripts/plottest.py xy . ./missing_xyz.txt ./err_xyz.txt
        ../../sources/scripts/cur_plot.py 'objs.txt' 'rotc.txt' 'sun.txt'
        mv *.txt ${i}
        mv *.png ${i}
        mv *.eps ${i}
done 

#        gnuplot -c ../../sources/scripts/main_cmd.gnu "cur.eps" 'objs.txt' 'ERROR_LIMITED_R_THETA'
