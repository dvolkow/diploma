#! /bin/bash

ORD=$1
IN_FILE=$2
H=$3

if [[ -n "$4" ]]; then
        ./diploma -f ${IN_FILE} --ord ${ORD} --mksize ${H} -b $4
        NEW_DIR=$(echo "result_${ORD}__${H}_$(echo ${IN_FILE}_$4 | sed 's/.txt//')")
elif [[ -n "$3" ]]; then
        ./diploma -f ${IN_FILE} --ord ${ORD} --mksize ${H} -b 1
        NEW_DIR=$(echo "result_${ORD}_${H}_$(echo ${IN_FILE} | sed 's/.txt//')")
else 
        ./diploma -f ${IN_FILE} --ord ${ORD} 
        NEW_DIR=$(echo "result_${ORD}_$(echo ${IN_FILE} | sed 's/.txt//')")
fi
gnuplot ../sources/scripts/main.gnu
#gnuplot ../sources/scripts/background.gnu
gnuplot ../sources/scripts/xy.gnu
gnuplot ../sources/scripts/yz.gnu
gnuplot ../sources/scripts/xz.gnu
#gnuplot ../sources/scripts/lb.gnu

mkdir ${NEW_DIR}

mv ./sun.txt ${NEW_DIR}
#mv ./averages.txt ${NEW_DIR}
#mv ./dump_table.txt ${NEW_DIR}
mv ./rotc.eps ${NEW_DIR}/${IN_FILE}_${ORD}.eps
mv ./rotc.txt ${NEW_DIR}
#mv ./background.txt ${NEW_DIR}
mv ./result.txt ${NEW_DIR}
mv ./unfresult.txt ${NEW_DIR}
mv ./objs.txt ${NEW_DIR}
#mv ./bk_sd.txt ${NEW_DIR}
#mv ./xyz_obj.txt ${NEW_DIR}
mv ./get_solution_178_0 ${NEW_DIR}
mv ./ERROR_LIMITED ${NEW_DIR}
mv ./ERROR_LIMITED_R_THETA ${NEW_DIR}
#mv ./back.eps ${NEW_DIR}/${IN_FILE}_${ORD}_${H}$4_back.eps
mv ./xy.eps ${NEW_DIR}/${IN_FILE}_${ORD}_xy.eps
mv ./xz.eps ${NEW_DIR}/${IN_FILE}_${ORD}_xz.eps
mv ./yz.eps ${NEW_DIR}/${IN_FILE}_${ORD}_yz.eps
#mv ./lb.eps ${NEW_DIR}/${IN_FILE}_${ORD}_${H}$4_lb.eps

cat ${NEW_DIR}/result.txt
#cat ${NEW_DIR}/bk_sd.txt
