#! /usr/bin/gnuplot

datafile = "output_data/direct_operational_characteristics.txt"

load "output_data/direct_operational_characteristics_opt.txt"

set linetype 1 lc rgb "red"
set linetype 2 lc rgb "green"
set linetype 3 lc rgb "blue"
set linetype 4 lc rgb "orange"
set linetype 5 lc rgb "black"
set linetype 6 lc rgb "violet"

set linetype cycle 6

set grid

set xlabel "K" font ", 15"
set ylabel "P_s(k)" font ", 15"

set tics font ", 11"
set key font ", 15"

if (ARG1 < 4) {
     set title "Operational characteristics for DIRECT on a family of tasks ".Name[ARG1 + 1] font "Helvetica Bold, 20"
     plot for [i = 1:3] datafile index 3 * ARG1 + i - 1 using 1:2 with lines lt i title "r = ".R[3 * ARG1 + i]
}
if (ARG1 == 4) {
     set title "Comparison of operational characteristics Grishagin and GKLS for DIRECT" font "Helvetica Bold, 20"
     plot for [i = 1:6] datafile index i - 1 using 1:2 with lines lt i title "r = ".R[i]."(".Name[(i - 1) / 3 + 1].")"
}
if (ARG1 == 5) {
     set title "Comparison of operational characteristics Grishagin and GKLS(constrained) for DIRECT" font "Helvetica Bold, 20"
     plot for [i = 7:12] datafile index i - 1 using 1:2 with lines lt i - 6 title "r = ".R[i]."(".Name[(i - 1) / 3 + 1].")"
}

bind all "alt-End" "exit gnuplot"
pause mouse close
