set terminal postscript eps enhanced color 20
set output "r_0.eps"
unset key
set linestyle 1 lt 1 lw 3
set xlabel "R, kpc"
set ylabel "R_0"
set yrange [5:10]
set mxtics 4
set mytics 5
plot 'r_0.txt' using 1:2:($2-$3):($2+$4) with yerrorbars ls 10 ps 0.5 pt 7 lt rgb 'red'
reset

