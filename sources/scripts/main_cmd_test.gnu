set terminal postscript eps enhanced color 20
set output ARG1
unset key

set multiplot title "Rotation"

set lmargin at screen 0.15
set rmargin at screen 0.85
set bmargin at screen 0.15
set tmargin at screen 0.90

set linestyle 1 lt 1 lw 3
set xlabel "R, kpc"
set ylabel "{/Symbol Q}, km/s"
set mxtics 4
set mytics 5
set yrange [-400:700]
set xrange [2:15]
plot ARG2 using 1:2 with p ps 0.1 lt rgb 'gray', ARG3 using 1:2 with p ps 0.2 lt rgb 'green', 'rotc.txt' using 1:2 with l ls 1 lt rgb 'blue','sun.txt' using 1:2 with p ps 1.4 pt 7 lt rgb 'yellow', 'rotc.txt' using 1:3 with l ls 2 lt rgb 'red', 'rotc.txt' using 1:4 with l ls 2 lt rgb 'red', 
reset
   

# second plot  (tall and narrow; at right of main plot)
#
set lmargin at screen 0.85
set rmargin at screen 1.00
set bmargin at screen 0.15
set tmargin at screen 0.90

set parametric
set dummy u,v
set view map

f(h) = sin(sqrt(h**2))/sqrt(h**2)

set surface

splot f(u), u, 0 with lines lc rgb "green"

unset parametric
unset multiplot











