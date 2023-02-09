#! /usr/bin/gnuplot

datafile = "output_data/mggsa_operational_characteristics_r_test.txt"

load "output_data/mggsa_operational_characteristics_r_test_opt.txt"

set linetype 1 lc rgb "red" lw 3
set linetype 2 lc rgb "green" lw 3
set linetype 3 lc rgb "orange" lw 2
set linetype 4 lc rgb "black" lw 2
set linetype 5 lc rgb "blue" lw 1

array Parameters[2]
Parameters[1] = "key = 1"
Parameters[2] = "key = 3, incr = 20"

ind = 2 * ARG2

set grid

set xlabel "r" font ", 15"

set tics font ", 11"
set key font ", 15"

set title "Operational characteristics for mggsa on a family of tasks ".Name[ARG2 + 1]."(r test)" font "Helvetica Bold, 20"
if (ARG1 == 0) {
    set ylabel "P_s^{max}" font ", 15"
    plot for [i = 1:2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title Parameters[i]
}
if (ARG1 == 1) {
    set y2tics
    set ylabel "P_s^{max}" font ", 15"
    set y2label "trials" font ", 15"
    plot for [i = 1:2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title "P_s^{max}(".Parameters[i].")" axes x1y1, \
         for [i = 1:2] datafile index ind + i - 1 using 1:3 with lines lt i * 2 title "trials(".Parameters[i].")" axes x1y2, \
         datafile index ind using 1:(K_max[ARG2 + 1]) with lines lt 5 title "K^{max}" axes x1y2
}
if (ARG1 == 2) {
    set y2tics
    set ylabel "P_s^{max}" font ", 15"
    set y2label "points" font ", 15"
    plot for [i = 1:2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title "P_s^{max}(".Parameters[i].")" axes x1y1, \
         for [i = 1:2] datafile index ind + i - 1 using 1:4 with lines lt i * 2 title "points(".Parameters[i].")" axes x1y2
}
if (ARG1 == 3) {
    set y2tics
    set ylabel "P_s^{max}" font ", 15"
    set y2label "points / trials" font ", 15"
    plot for [i = 1:2] datafile index ind + i - 1 using 1:2 with lines lt i * 2 - 1 title "P_s^{max}(".Parameters[i].")" axes x1y1, \
         for [i = 1:2] datafile index ind + i - 1 using 1:($4 / $3) with lines lt i * 2 title "points / trials(".Parameters[i].")" axes x1y2
}

bind all "alt-End" "exit gnuplot"
pause mouse close
