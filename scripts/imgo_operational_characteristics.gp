#! /usr/bin/gnuplot

datafile="output_data/imgo_operational_characteristics.txt"

load "output_data/imgo_operational_characteristics_opt.txt"

set linetype 1 lc rgb "red"
set linetype 2 lc rgb "green"
set linetype 3 lc rgb "blue"
set linetype 4 lc rgb "orange"
set linetype 5 lc rgb "black"
set linetype 6 lc rgb "violet"

set linetype cycle 6

set grid

set xlabel "K" font ", 13"
set ylabel "P_s(k)" font ", 13"

set tics font ", 13"
set key font ", 15"

if (ARG1 == 0 || ARG1 == 1) {
    set title "Operational characteristics on a family of tasks ".Name[ARG1 + 1] font "Helvetica Bold, 25"
    plot for [i = 1:3] datafile index ARG1 * 3 + i - 1 using 1:2 with lines lt i title "r = ".R[ARG1 * 3 + i]
}
if (ARG1 == 2) {
    set title "Comparison of operational characteristics for the Hill and Shekel families" font "Helvetica Bold, 20"
    plot for [i = 1:6] datafile index i - 1 using 1:2 with lines lt i title "r = ".R[i]."(".Name[(i - 1) / 3 + 1].")"
}

bind all "alt-End" "exit gnuplot"
pause mouse close
