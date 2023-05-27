#! /usr/bin/gnuplot

datafile = "output_data/direct_operational_characteristics.txt"

load "output_data/direct_operational_characteristics_opt.txt"

set linetype 1 lc rgb "red"
set linetype 2 lc rgb "green"
set linetype 3 lc rgb "blue"
set linetype 4 lc rgb "orange"

set linetype cycle 4

set grid

set xlabel "K" font ", 15"
set ylabel "P_s(k)" font ", 15"

set tics font ", 11"
set key font ", 15"

array alghorithm[2]
alghorithm[1] = "DIRECT ORIGINAL"
alghorithm[2] = "DIRECT GABLONSKY"

if (ARG1 < 4) {
     set title "Operational characteristics for DIRECT on a family of tasks ".familyNames[ARG1 + 1] font "Helvetica Bold, 20"
     plot for [i = 1:2] datafile index 2 * ARG1 + i - 1 using 1:2 with lines lt i title alghorithm[i]
}
if (ARG1 == 4) {
     set title "Comparison of operational characteristics Grishagin and GKLS for DIRECT" font "Helvetica Bold, 20"
     plot for [i = 1:4] datafile index i - 1 using 1:2 with lines lt i title alghorithm[i % 2 + 1]."(".familyNames[(i - 1) / 2 + 1].")"
}
if (ARG1 == 5) {
     set title "Comparison of operational characteristics Grishagin and GKLS(constrained) for DIRECT" font "Helvetica Bold, 20"
     plot for [i = 5:8] datafile index i - 1 using 1:2 with lines lt i - 4 title alghorithm[i % 2 + 1]."(".familyNames[(i - 1) / 2 + 1].")"
}

bind all "alt-End" "exit gnuplot"
pause mouse close
