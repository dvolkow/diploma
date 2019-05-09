set terminal postscript eps enhanced color 20
set output "sigma.eps"
unset key
set linestyle 1 lt 1 lw 3
set xlabel "N"
set ylabel "{/Symbol s}"
set mxtics 4
set mytics 5
plot 'sigma.txt' using 1:2 with l ls 10 lt rgb 'red'
reset

