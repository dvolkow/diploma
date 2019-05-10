set terminal postscript eps enhanced color 20
set output "b_rotc.eps"
unset key
set linestyle 1 lt 1 lw 3
set xlabel "R, kpc"
set ylabel "{/Symbol Q}, km/s"
set mxtics 4
set mytics 5
set yrange [-400:700]
set xrange [2:15]
plot 'B_PART_OBJ.txt' using 1:2 with p ps 0.1 lt rgb 'gray', 'b_cur.txt' using 1:2 with l ls 1 lt rgb 'blue' #,'sun.txt' using 1:2 with p ps 1.4 pt 7 lt rgb 'yellow', 'rotc.txt' using 1:3 with l ls 2 lt rgb 'red', 'rotc.txt' using 1:4 with l ls 2 lt rgb 'red', 
reset
   












