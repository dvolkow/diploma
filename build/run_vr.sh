#! /bin/bash

ORD=$1
IN_FILE=$2
H=$3

if [[ -n "$4" ]]; then
        ./diploma -f ${IN_FILE} --ord ${ORD} --mksize ${H} -b $4
        NEW_DIR=$(echo "result_${ORD}_${H}_$(echo ${IN_FILE}_$4 | sed 's/.txt//')")
elif [[ -n "$3" ]]; then
        ./diploma -f ${IN_FILE} --ord ${ORD} --mksize ${H} -b 1
        NEW_DIR=$(echo "result_${ORD}_${H}_$(echo ${IN_FILE} | sed 's/.txt//')")
else 
        ./diploma -f ${IN_FILE} --ord ${ORD} -b 1 
        NEW_DIR=$(echo "result_${ORD}_$(echo ${IN_FILE} | sed 's/.txt//')")
fi

gnuplot -c ../sources/scripts/main_cmd.gnu 'rotc.eps' 'vr_objs.txt' 'vr_objs_err.txt'
gnuplot ../sources/scripts/AR0.gnu
mkdir ${NEW_DIR}

PREFIX="${IN_FILE}_${ORD}"

mv ./ERROR_LIMITED ${NEW_DIR}

mv ./R0Theta0.txt ${NEW_DIR}
mv ./R0Theta0_main.txt ${NEW_DIR}
mv ./AR0.eps ${NEW_DIR}/${PREFIX}_AR0.eps

mv ./result.txt ${NEW_DIR}
mv ./sun.txt ${NEW_DIR}

mv ./rotc.eps ${NEW_DIR}/${PREFIX}.eps
mv ./rotc.txt ${NEW_DIR}

mv ./vr_objs.txt ${NEW_DIR}
mv ./vr_objs_err.txt ${NEW_DIR}

cat ${NEW_DIR}/result.txt

