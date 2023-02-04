#! /usr/bin/gnuplot

datafile = "output_data/incr_test.txt"

load "output_data/incr_test_opt.txt"

set linetype 1 lc rgb "#0000FF" lw 5 pt 1 dt 1
set linetype 2 lc rgb "#0020FF" lw 4 pt 1 dt 2
set linetype 3 lc rgb "#0040FF" lw 3 pt 1 dt 3
set linetype 4 lc rgb "#0060FF" lw 2 pt 1 dt 4
set linetype 5 lc rgb "#0080FF" lw 1 pt 1 dt 5

set linetype cycle 5

set xlabel "increment"
set xrange [I_MIN[ARG2 - ARG3 + 1]:I_MAX[ARG2 - ARG3 + 1]]

set title 'Increment test for Rastrigin function(N = '.ARG2.')' font "Helvetica Bold, 20"

if (ARG1 == 0) {
    set ylabel "trials"
    plot for [i = ARG5:ARG6] datafile index (ARG2 - ARG3) * (ARG6 - ARG5 + 1) + i - ARG5 using 1:2 with lines lt i - ARG5 + 1 title sprintf("m = %d", i)
}
if (ARG1 == 1) {
    set ylabel "points"
    plot for [i = ARG5:ARG6] datafile index (ARG2 - ARG3) * (ARG6 - ARG5 + 1) + i - ARG5 using 1:3 with lines lt i - ARG5 + 1 title sprintf("m = %d", i)
}
if (ARG1 == 2) {
    set ylabel "accuracy"
    plot for [i = ARG5:ARG6] datafile index (ARG2 - ARG3) * (ARG6 - ARG5 + 1) + i - ARG5 using 1:4 with lines lt i - ARG5 + 1 title sprintf("m = %d", i)
}
if (ARG1 == 3) {
    set ylabel "points / trials"
    plot for [i = ARG5:ARG6] datafile index (ARG2 - ARG3) * (ARG6 - ARG5 + 1) + i - ARG5 using 1:($3 / $2) with lines lt i - ARG5 + 1 title sprintf("m = %d", i)
}

bind all "alt-End" "exit gnuplot"
pause mouse close
