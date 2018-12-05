set terminal postscript eps enhanced color 20
set output "cds_without_ej.eps"
unset key
set linestyle 1 lt 1 lw 3
set xlabel "r_1, kpc"
set ylabel "r_2, kpc"
set mxtics 5
set mytics 5
set yrange [0:15]
set xrange [0:15]
f(x) = x
g(x) = 1.007775 * x
g_up(x) = (1.007775 + 0.021401) * x
g_b(x) = (1.007775 - 0.021401) * x
plot 'delta.dat' using 1:2 with p ps 0.2 lt rgb 'black', f(x) lt rgb 'black', g(x) lt rgb 'green', g_up(x) lt rgb 'blue', g_b(x) lt rgb 'blue'
reset
   


