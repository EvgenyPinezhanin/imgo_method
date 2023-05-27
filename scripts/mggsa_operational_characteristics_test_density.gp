#! /usr/bin/gnuplot

datafile = "output_data/mggsa_operational_characteristics_test_density.txt"

load "output_data/mggsa_operational_characteristics_test_density_opt.txt"

set style line 1 lt 1 lc rgb "red" lw 4
set style line 2 lt 2 lc rgb "blue" lw 4
set style line 3 lt 3 lc rgb "orange" lw 4
set style line 4 lt 4 lc rgb "green" lw 4
set style line 5 lt 5 lc rgb "dark-violet" lw 4

set style line 6 lt 1 lc rgb "brown" lw 2
set style line 7 lt 2 lc rgb "dark-yellow" lw 2
set style line 8 lt 3 lc rgb "dark-violet" lw 2
set style line 9 lt 4 lc rgb "black" lw 2
set style line 10 lt 5 lc rgb "cyan" lw 2

set grid
set xlabel "K" font ", 15"
set ylabel "P_s(k)" font ", 15"
set tics font ", 11"
set key font ", 15"

set title "Operational characteristics for mggsa on a family of tasks ".familyNames[ARG1 + 1]." (density test)" font "Helvetica Bold, 20"
plot for [i = 1:sizeDen] datafile index sizeDen * ARG1 + i - 1 using 1:2 with linespoints ls i title "den = ".den[i]

bind all "alt-End" "exit gnuplot"
pause mouse close
