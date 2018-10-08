set terminal postscript eps enhanced color 20
set output "rotc.eps"
unset key
set linestyle 1 lt 1 lw 3
set xlabel "R, kpc"
set ylabel "{/Symbol Q}, km/s"
set mxtics 4
set mytics 5
set yrange [-400:700]
set xrange [2:15]
plot 'objs.txt' using 1:2 with p ps 0.1 lt rgb 'black', 'rotc.txt' using 1:2 with l ls 1 lt rgb 'blue', 'sun.txt' using 1:2 with p ps 0.7 pt 6 lt rgb 'yellow', 'averages.txt' using 3:6:($3-$4):($3+$5):($6-$7):($6+$7) with xyerrorbars ls 10 ps 0.5 pt 7 lt rgb 'red'
reset
   










