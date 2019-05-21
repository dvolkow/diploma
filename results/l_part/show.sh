#! /bin/bash
for (( i=$1; i<=$2; i++))
do
        cat ${i}/r0_bounds.txt
done 

#        gnuplot -c ../../sources/scripts/main_cmd.gnu "cur.eps" 'objs.txt' 'ERROR_LIMITED_R_THETA'
#        mv *.eps ${i}
