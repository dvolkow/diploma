#! /bin/bash

ORD=$1
IN_FILE=$2
H=$3

if [[ -n "$4" ]]; then
        ./diploma -f ${IN_FILE} --ord ${ORD} --mksize ${H} -s $4 -b 0 -p 1
        NEW_DIR=$(echo "result_${ORD}_${H}_$(echo ${IN_FILE}_$4 | sed 's/.txt//')")
fi

#mv get_untited_solution_178_0 get_solution_178_0

./plot.sh

mkdir ${NEW_DIR}

PREFIX="${IN_FILE}_${ORD}"

mv ./b_cur.txt ${NEW_DIR}
mv ./B_PART_OBJ.txt ${NEW_DIR}
mv ./b_rotc.eps ${NEW_DIR}/${PREFIX}_b_rotc.eps

mv ./ERROR_LIMITED ${NEW_DIR}
mv ./ERROR_LIMITED_R_THETA ${NEW_DIR}

mv ./get_solution_178_0 ${NEW_DIR}

mv ./l_cur.txt ${NEW_DIR}
mv ./L_PART_OBJ.txt ${NEW_DIR}
mv ./l_rotc.eps ${NEW_DIR}/${PREFIX}_l_rotc.eps

mv ./objs.txt ${NEW_DIR}

mv ./R0Theta0.txt ${NEW_DIR}
mv ./R0Theta0_main.txt ${NEW_DIR}
mv ./R0Theta0.eps ${NEW_DIR}/${PREFIX}_R0Theta0.eps

mv ./result.txt ${NEW_DIR}
mv ./sun.txt ${NEW_DIR}

mv ./rotc.eps ${NEW_DIR}/${PREFIX}.eps
mv ./profile.eps ${NEW_DIR}/${PREFIX}_profile.eps
mv ./rotc.txt ${NEW_DIR}

mv ./unfresult.txt ${NEW_DIR}

mv ./vr_cur.txt ${NEW_DIR}
mv ./VR_PART_OBJ.txt ${NEW_DIR}
mv ./vr_rotc.eps ${NEW_DIR}/${PREFIX}_vr_rotc.eps

mv ./xy.eps ${NEW_DIR}/${PREFIX}_xy.eps
mv ./xz.eps ${NEW_DIR}/${PREFIX}_xz.eps
mv ./yz.eps ${NEW_DIR}/${PREFIX}_yz.eps

cat ${NEW_DIR}/result.txt
