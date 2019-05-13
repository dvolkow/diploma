set terminal postscript eps enhanced color 20
set output "AR0.eps"
unset key
set yrange [85:115]
set xrange [6.5:8.5]
set mxtics 5
set mytics 5
set xlabel "R_0, kpc" 
set ylabel "{/Symbol T}, kmps"
plot 'R0Theta0.txt' using 1:2 with p ps 0.3 lt rgb 'black', 'R0Theta0_main.txt' using 1:2 with p ps 0.5 lt rgb 'blue'
reset



