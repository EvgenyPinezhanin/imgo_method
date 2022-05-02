#! /usr/bin/gnuplot

reset

datafile="peano_operational_characteristics.txt"

array Name[4]
array R[3*4]
load "peano_operational_characteristics_opt.txt"

ind=3*ARG1

set xlabel "K"
set ylabel "P_s(k)"
set grid

set title "Operational characteristics for imgo on a family of tasks ".Name[ARG1+1] font "Helvetica Bold, 20"
plot datafile index ind using 1:2 with lines ls 5 lc rgb "red" title "r=".R[ind+1], \
     datafile index ind+1 using 1:2 with lines ls 5 lc rgb "green" title "r=".R[ind+2], \
     datafile index ind+2 using 1:2 with lines ls 5 lc rgb "blue" title "r=".R[ind+3], \

bind all "alt-End" "exit gnuplot"
pause mouse close
