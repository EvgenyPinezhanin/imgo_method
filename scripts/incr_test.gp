#! /usr/bin/gnuplot

reset

datafile = "incr_test.txt"

set linetype 1 lc rgb "#000FFF" lw 2 pt 1
set linetype 2 lc rgb "#002FFF" lw 2 pt 1
set linetype 3 lc rgb "#004FFF" lw 2 pt 1
set linetype 4 lc rgb "#006FFF" lw 2 pt 1
set linetype 5 lc rgb "#008FFF" lw 2 pt 1
set linetype 6 lc rgb "#00AFFF" lw 2 pt 1
set linetype 7 lc rgb "#00CFFF" lw 2 pt 1
set linetype 8 lc rgb "#00EFFF" lw 2 pt 1

set linetype cycle 8

set xlabel "incr"
set xrange [ARG2:ARG3]

set title 'increment test'

if (ARG1 == 0) {
    set ylabel "trials"
    plot for [i = ARG4:ARG5] datafile index i - ARG4 using 1:2 with lines lt i - ARG4 + 1 title sprintf("m = %d", i), \
         for [i = ARG4:ARG5] datafile index i - ARG4 using 1:2:(sprintf("%d", i)) with labels offset 0.0, 0.5 notitle
} else {
    set ylabel "accuracy"
    plot for [i = ARG4:ARG5] datafile index i - ARG4 using 1:3 with lines lt i - ARG4 + 1 title sprintf("m = %d", i), \
         for [i = ARG4:ARG5] datafile index i - ARG4 using 1:3:(sprintf("%d", i)) with labels offset 0.0, 0.5 notitle
}

bind all "alt-End" "exit gnuplot"
pause mouse close
