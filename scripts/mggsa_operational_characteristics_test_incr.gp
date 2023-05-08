#! /usr/bin/gnuplot

datafile = "output_data/mggsa_operational_characteristics_test_incr.txt"

load "output_data/mggsa_operational_characteristics_test_incr_opt.txt"

set linetype 1 lc rgb "red"
set linetype 2 lc rgb "green"
set linetype 3 lc rgb "blue"
set linetype 4 lc rgb "orange"

set linetype cycle 4

set xlabel "increment" font ", 15"
set grid

set tics font ", 11"
set key font ", 15"

ind = numberKey * ARG2

set title "Operational characteristics for mggsa on a family of tasks ".familyNames[ARG2 + 1]."(incr test)" font "Helvetica Bold, 20"
if (ARG1 == 0) {
    set ylabel "trials" font ", 15"
    plot for [i = 1:numberKey] datafile index ind + i - 1 using 1:2 with lines lt i title "r = ".r[ind + i]
}
if (ARG1 == 1) {
    set ylabel "trial points" font ", 15"
    plot for [i = 1:numberKey] datafile index ind + i - 1 using 1:3 with lines lt i title "r = ".r[ind + i]
}
if (ARG1 == 2) {
    set ylabel "trial points / trials" font ", 15"
    plot for [i = 1:numberKey] datafile index ind + i - 1 using 1:($3 / $2) with lines lt i title "r = ".r[ind + i]
}

bind all "alt-End" "exit gnuplot"
pause mouse close
