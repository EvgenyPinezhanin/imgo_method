#! /usr/bin/gnuplot

datafile = "output_data/mggsa_test_lipschitz.txt"

load "output_data/mggsa_test_lipschitz_opt.txt"

set linetype 1 lc rgb "#0000FF" lw 5 pt 1 dt 1
set linetype 2 lc rgb "#0020FF" lw 5 pt 1 dt 2
set linetype 3 lc rgb "#0040FF" lw 5 pt 1 dt 3
set linetype 4 lc rgb "#0060FF" lw 4 pt 1 dt 4
set linetype 5 lc rgb "#0080FF" lw 4 pt 1 dt 5
set linetype 6 lc rgb "#00A0FF" lw 4 pt 1 dt 4
set linetype 7 lc rgb "#00C0FF" lw 3 pt 1 dt 3
set linetype 8 lc rgb "#00E0FF" lw 3 pt 1 dt 2
set linetype 9 lc rgb "#00FFFF" lw 3 pt 1 dt 1

set linetype 10 lc rgb "#00FF00" lw 4 pt 1 dt 1
set linetype 11 lc rgb "#20FF00" lw 4 pt 1 dt 2
set linetype 12 lc rgb "#40FF00" lw 4 pt 1 dt 3
set linetype 13 lc rgb "#60FF00" lw 3 pt 1 dt 4
set linetype 14 lc rgb "#80FF00" lw 3 pt 1 dt 5
set linetype 15 lc rgb "#A0FF00" lw 3 pt 1 dt 4
set linetype 16 lc rgb "#C0FF00" lw 2 pt 1 dt 3
set linetype 17 lc rgb "#E0FF00" lw 2 pt 1 dt 2
set linetype 18 lc rgb "#FFFF00" lw 2 pt 1 dt 1

set linetype 19 lc rgb "#FF0000" lw 3 pt 1 dt 1
set linetype 20 lc rgb "#FF0020" lw 3 pt 1 dt 2
set linetype 21 lc rgb "#FF0040" lw 3 pt 1 dt 3
set linetype 22 lc rgb "#FF0060" lw 2 pt 1 dt 4
set linetype 23 lc rgb "#FF0080" lw 2 pt 1 dt 5
set linetype 24 lc rgb "#FF00A0" lw 2 pt 1 dt 4
set linetype 25 lc rgb "#FF00C0" lw 1 pt 1 dt 3
set linetype 26 lc rgb "#FF00E0" lw 1 pt 1 dt 2
set linetype 27 lc rgb "#FF00FF" lw 1 pt 1 dt 1

set linetype cycle 27

set xlabel "increment"
set ylabel "Lipschitz constant"
set xrange [ARG6:ARG7]

set title 'Lipschitz test for the function '.familyNames[ARG1 + 1] font "Helvetica Bold, 20"

plot for [i = ARG2:ARG3] for [j = ARG4:ARG5] datafile index (ARG5 - ARG4 + 1) * (i - 1) + j - ARG4 \
     using 1:2 with lines lt 9 * (i - 1) + j - ARG4 + 1 title sprintf("key = %d, m = %d", i, j)

bind all "alt-End" "exit gnuplot"
pause mouse close
