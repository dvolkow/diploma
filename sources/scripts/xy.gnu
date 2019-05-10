set terminal postscript eps enhanced color 20
set output "xy.eps"
set size ratio -1
unset key
set yrange [-14:14]
set xrange [-14:14]
set mxtics 5
set mytics 5
set xlabel "X, kpc" 
set ylabel "Y, kpc"
plot 'get_solution_178_0' using 1:2 with p ps 0.1 lt rgb 'black', 'ERROR_LIMITED' using 1:2 with p ps 0.2 lt rgb 'red'
reset

