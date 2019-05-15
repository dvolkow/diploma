#! /bin/bash

ORD=$1
IN_FILE=$2
H=$3

./diploma -f ${IN_FILE} --ord ${ORD} --mksize ${H} -s $4 -p 1
NEW_DIR=$(echo "result_${ORD}_${H}_$(echo ${IN_FILE}_$4 | sed 's/.txt//')")

gnuplot -c ../sources/scripts/main_cmd.gnu 'rotc.eps' 'vr_objs.txt' 'vr_objs_err.txt'
gnuplot ../sources/scripts/AR0.gnu
gnuplot -c ../sources/scripts/sigma.gnu "profile.eps" 'vr_profile.txt' 'red'
mkdir ${NEW_DIR}

PREFIX="${IN_FILE}_${ORD}"

mv ./ERROR_LIMITED ${NEW_DIR}

mv ./R0Theta0.txt ${NEW_DIR}
mv ./R0Theta0_main.txt ${NEW_DIR}
mv ./AR0.eps ${NEW_DIR}/${PREFIX}_AR0.eps

mv ./result.txt ${NEW_DIR}
mv ./unfresult.txt ${NEW_DIR}
mv ./sun.txt ${NEW_DIR}

mv ./vr_profile.txt ${NEW_DIR}
mv ./profile.eps ${NEW_DIR}/${PREFIX}_vr_profile.eps

mv ./rotc.eps ${NEW_DIR}/${PREFIX}.eps
mv ./rotc.txt ${NEW_DIR}

mv ./vr_objs.txt ${NEW_DIR}
mv ./vr_objs_err.txt ${NEW_DIR}

cat ${NEW_DIR}/result.txt

