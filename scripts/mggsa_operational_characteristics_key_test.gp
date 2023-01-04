#! /usr/bin/gnuplot

datafile="output_data/mggsa_operational_characteristics_key_test.txt"

load "output_data/mggsa_operational_characteristics_key_test_opt.txt"

set linetype 1 lc rgb "red"
set linetype 2 lc rgb "green"
set linetype 3 lc rgb "blue"
set linetype 4 lc rgb "orange"

set linetype cycle 4

ind = count_key * ARG1

set xlabel "K" font ", 15"
set ylabel "P_s(k)" font ", 15"
set grid

set tics font ", 11"
set key font ", 15"

set title "Operational characteristics for mggsa on a family of tasks ".Name[ARG1+1] font "Helvetica Bold, 20"
plot for [i = 1:count_key] datafile index ind + i - 1 using 1:2 with lines lt i title "r=".R[ind+i].", key = ".Key[i]

bind all "alt-End" "exit gnuplot"
pause mouse close
