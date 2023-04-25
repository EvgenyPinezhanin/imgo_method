#! /usr/bin/gnuplot

datafile="output_data/mggsa_operational_characteristics_test_incr.txt"

load "output_data/mggsa_operational_characteristics_test_incr_opt.txt"

set linetype 1 lc rgb "red"
set linetype 2 lc rgb "green"
set linetype 3 lc rgb "blue"
set linetype 4 lc rgb "orange"

set linetype cycle 4

ind = count_key * ARG2

set xlabel "increment" font ", 15"
set grid

set tics font ", 11"
set key font ", 15"

set title "Operational characteristics for mggsa on a family of tasks ".Name[ARG2 + 1]."(incr test)" font "Helvetica Bold, 20"
if (ARG1 == 0) {
    set ylabel "trials" font ", 15"
    plot for [i = 1:count_key] datafile index ind + i - 1 using 1:2 with lines lt i title "r = ".R[ind+i]
}
if (ARG1 == 1) {
    set ylabel "points" font ", 15"
    plot for [i = 1:count_key] datafile index ind + i - 1 using 1:3 with lines lt i title "r = ".R[ind+i]
}
if (ARG1 == 2) {
    set ylabel "points / trials" font ", 15"
    plot for [i = 1:count_key] datafile index ind + i - 1 using 1:($3 / $2) with lines lt i title "r = ".R[ind+i]
}

bind all "alt-End" "exit gnuplot"
pause mouse close
