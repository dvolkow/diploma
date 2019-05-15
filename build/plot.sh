#! /bin/bash
gnuplot ../sources/scripts/main.gnu
gnuplot ../sources/scripts/xy.gnu
gnuplot ../sources/scripts/xz.gnu
gnuplot ../sources/scripts/yz.gnu
gnuplot ../sources/scripts/vr.gnu
gnuplot ../sources/scripts/b.gnu
gnuplot ../sources/scripts/l.gnu
gnuplot ../sources/scripts/R0Theta0.gnu
gnuplot -c ../sources/scripts/sigma.gnu "profile.eps" 'uni_profile.txt' 'red'
