set terminal postscript eps enhanced color 20
set output "lb.eps"
unset key
set linestyle 1 lt 1 lw 3
set size ratio -1
set yrange [-90:90]
set xrange [0:360]
set mxtics 10
set mytics 10
plot 'dump_table.txt' using 1:2 with p ps 0.1 lt rgb 'black'
reset

