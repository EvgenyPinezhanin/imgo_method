#! /usr/bin/gnuplot

reset

datafile="peano_operational_characteristics.txt"

set xlabel "K"
set ylabel "P_s(k)"
set grid

set title "Operational characteristics for imgo with peano" font "Helvetica Bold, 20"
plot datafile using 1:2 with lines ls 5 lc rgb "red" title "P_s(k)", \

bind all "alt-End" "exit gnuplot"
pause mouse close
