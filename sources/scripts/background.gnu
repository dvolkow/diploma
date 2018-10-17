set terminal postscript eps enhanced color 20
set output "back.eps"
unset key
set linestyle 1 lt 1 lw 3
set xlabel "R, kpc"
set ylabel "{/Symbol Q}, km/s"
set yrange [0:100]
set mxtics 4
set mytics 5
plot 'background.txt' using 1:4:2:3:($4-$5):($4+$5) with xyerrorbars ls 10 ps 0.5 pt 7 lt rgb 'blue'
reset

