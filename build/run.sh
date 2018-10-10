#! /bin/bash

ORD=$1
IN_FILE=$2
L=0
H=$3
./diploma -f ${IN_FILE} --ord ${ORD} --filter L ${L} ${H}
gnuplot ../sources/scripts/main.gnu
NEW_DIR=$(echo "result_${ORD}_${L}_${H}_${IN_FILE}")
mkdir ${NEW_DIR}

mv ./sun.txt ${NEW_DIR}
mv ./averages.txt ${NEW_DIR}
mv ./rotc.eps ${NEW_DIR}/${IN_FILE}_${ORD}_${L}_${H}.eps
mv ./rotc.txt ${NEW_DIR}
mv ./result.txt ${NEW_DIR}
mv ./objs.txt ${NEW_DIR}
