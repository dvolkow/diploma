set terminal postscript eps enhanced color 20
set output "yz.eps"
set size ratio -1
unset key
set yrange [-14:14]
set xrange [-14:14]
set mxtics 5
set mytics 5
set xlabel "Y, kpc" 
set ylabel "Z, kpc"
plot 'xyz_obj.txt' using 2:3 with p ps 0.1 lt rgb 'black'
reset

