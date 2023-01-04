#! /usr/bin/gnuplot

datafile="output_data/mggsa_operational_characteristics_r_test.txt"

array Name[4]
load "output_data/mggsa_operational_characteristics_r_test_opt.txt"

ind=1*ARG1

set xlabel "r" font ", 15"
set ylabel "P_s^{max}(r)" font ", 15"
set grid

set tics font ", 11"
set key font ", 15"

set title "Operational characteristics for mggsa on a family of tasks ".Name[ARG1+1] font "Helvetica Bold, 20"
plot datafile index ind using 1:2 with lines ls 5 lc rgb "red" title "P_s^{max}(r)"

bind all "alt-End" "exit gnuplot"
pause mouse close
