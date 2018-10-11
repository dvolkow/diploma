#! /bin/bash

ORD=$1
IN_FILE=$2
L=0
H=$3
./diploma -f ${IN_FILE} --ord ${ORD} --filter L ${L} ${H}
gnuplot ../sources/scripts/main.gnu
gnuplot ../sources/scripts/background.gnu
NEW_DIR=$(echo "result_${ORD}_${L}_${H}_${IN_FILE}")
mkdir ${NEW_DIR}

mv ./sun.txt ${NEW_DIR}
mv ./averages.txt ${NEW_DIR}
mv ./rotc.eps ${NEW_DIR}/${IN_FILE}_${ORD}_${L}_${H}.eps
mv ./rotc.txt ${NEW_DIR}
mv ./result.txt ${NEW_DIR}
mv ./objs.txt ${NEW_DIR}
mv ./bk_sd.txt ${NEW_DIR}
mv ./back.eps ${NEW_DIR}/${IN_FILE}_${ORD}_${L}_${H}_back.eps

cat ${NEW_DIR}/result.txt
cat ${NEW_DIR}/bk_sd.txt
