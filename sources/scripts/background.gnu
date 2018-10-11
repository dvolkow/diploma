set terminal postscript eps enhanced color 20
set output "back.eps"
unset key
set linestyle 1 lt 1 lw 3
set xlabel "R, kpc"
set ylabel "{/Symbol Q}, km/s"
set yrange [-100:100]
set mxtics 4
set mytics 5
plot 'background.txt' using 1:4 with l ls 1 lt rgb 'blue', 'background.txt' using 1:3 with l ls 1 lt rgb 'red'
reset

