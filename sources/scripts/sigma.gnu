set terminal postscript eps enhanced color 20
set output ARG1
unset key
set linestyle 1 lt 1 lw 3
set xlabel "N"
set ylabel "{/Symbol s}"
set mxtics 4
set mytics 5
plot ARG2 using 1:2 with l ls 10 lt rgb ARG3
reset

