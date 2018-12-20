set terminal postscript eps enhanced color 20
set output "cds.eps"
unset key
set linestyle 1 lt 1 lw 3
set xlabel "r_2, kpc"
set ylabel "r_1, kpc"
set size ratio -1
set mxtics 5
set mytics 5
set yrange [0:8]
set xrange [0:8]
f(x) = x
g(x) = 1.007775 * x
g_up(x) = (1.007775 + 0.021401) * x
g_b(x) = (1.007775 - 0.021401) * x
plot 'ejected.dat' using 2:1 with p ps 0.1 lt rgb 'red', 'delta.dat' using 2:1 with p ps 0.1 lt rgb 'black', f(x) lt rgb 'yellow', g(x) lt rgb 'green', g_up(x) lt rgb 'blue', g_b(x) lt rgb 'blue'
reset
   










