#! /bin/bash

ORD=$1
IN_FILE=$2
L=0
H=$3

if [[ -n "$4" ]]; then
        ./diploma -f ${IN_FILE} --ord ${ORD} --filter L ${L} ${H} --mode $4
elif [[ -n "$3" ]]; then
        ./diploma -f ${IN_FILE} --ord ${ORD} --filter L ${L} ${H}
else 
        ./diploma -f ${IN_FILE} --ord ${ORD} 
fi
gnuplot ../sources/scripts/main.gnu
gnuplot ../sources/scripts/background.gnu
gnuplot ../sources/scripts/xy.gnu
gnuplot ../sources/scripts/yz.gnu
gnuplot ../sources/scripts/xz.gnu
gnuplot ../sources/scripts/lb.gnu

NEW_DIR=$(echo "result_${ORD}_${L}_${H}_$(echo ${IN_FILE} | sed 's/.txt//')")
mkdir ${NEW_DIR}

mv ./sun.txt ${NEW_DIR}
mv ./averages.txt ${NEW_DIR}
cp ./dump_table.txt ${NEW_DIR}
mv ./rotc.eps ${NEW_DIR}/${IN_FILE}_${ORD}_${L}_${H}$4.eps
mv ./rotc.txt ${NEW_DIR}
mv ./result.txt ${NEW_DIR}
mv ./objs.txt ${NEW_DIR}
mv ./bk_sd.txt ${NEW_DIR}
mv ./xyz_obj.txt ${NEW_DIR}
mv ./back.eps ${NEW_DIR}/${IN_FILE}_${ORD}_${L}_${H}$4_back.eps
mv ./xy.eps ${NEW_DIR}/${IN_FILE}_${ORD}_${L}_${H}$4_xy.eps
mv ./xz.eps ${NEW_DIR}/${IN_FILE}_${ORD}_${L}_${H}$4_xz.eps
mv ./yz.eps ${NEW_DIR}/${IN_FILE}_${ORD}_${L}_${H}$4_yz.eps
mv ./lb.eps ${NEW_DIR}/${IN_FILE}_${ORD}_${L}_${H}$4_lb.eps

cat ${NEW_DIR}/result.txt
cat ${NEW_DIR}/bk_sd.txt
