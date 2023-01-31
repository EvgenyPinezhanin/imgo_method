#! /usr/bin/gnuplot

datafile = "output_data/mggsa_operational_characteristics_r_test.txt"

array Name[4]
load "output_data/mggsa_operational_characteristics_r_test_opt.txt"

set linetype 1 lc rgb "red"
set linetype 2 lc rgb "green"
set linetype 3 lc rgb "blue"
set linetype 4 lc rgb "orange"

set linetype cycle 4

ind = 1 * ARG2

set xlabel "r" font ", 15"
set tics font ", 11"
set key font ", 15"
set grid

set title "Operational characteristics for mggsa on a family of tasks ".Name[ARG2 + 1]."(r test)" font "Helvetica Bold, 20"
if (ARG1 == 0) {
    set ylabel "P_s^{max}" font ", 15"
    plot for [i = 1:1] datafile index ind + i - 1 using 1:2 with lines lt i title "incr = 30"
}
if (ARG1 == 1) {
    set y2tics
    set ylabel "P_s^{max}" font ", 15"
    set y2label "trials" font ", 15"
    plot for [i = 1:1] datafile index ind + i - 1 using 1:2 with lines lt i title "P_s^{max}(incr = 30)" axes x1y1, \
         for [i = 1:1] datafile index ind + i - 1 using 1:3 with lines lt i + 1 title "trials(incr = 30)" axes x1y2
}
if (ARG1 == 2) {
    set y2tics
    set ylabel "P_s^{max}" font ", 15"
    set y2label "points" font ", 15"
    plot for [i = 1:1] datafile index ind + i - 1 using 1:2 with lines lt i title "P_s^{max}(incr = 30)" axes x1y1, \
         for [i = 1:1] datafile index ind + i - 1 using 1:4 with lines lt i + 1 title "points(incr = 30)" axes x1y2 
}
if (ARG1 == 3) {
    set y2tics
    set ylabel "P_s^{max}" font ", 15"
    set y2label "points / trials" font ", 15"
    plot for [i = 1:1] datafile index ind + i - 1 using 1:2 with lines lt i title "P_s^{max}(incr = 30)" axes x1y1, \
         for [i = 1:1] datafile index ind + i - 1 using 1:($4 / $3) with lines lt i + 1 title "points / trials(incr = 30)" axes x1y2
}

bind all "alt-End" "exit gnuplot"
pause mouse close
